#include "eamxx_gocart_surface_emission.hpp"

#include "share/grid/remap/identity_remapper.hpp"
#include "share/grid/remap/refining_remapper_p2p.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"
#include "share/util/eamxx_timing.hpp"

#include "share/field/field.hpp"

#include "share/grid/remap/coarsening_remapper.hpp"
#include "share/grid/remap/horiz_interp_remapper_base.hpp"


namespace scream::gocart
{

    // ------- Static variable intial values --------
    std::shared_ptr<Kokkos::View<Real*>> SurfaceEmission::s_io_buffer        = nullptr;
    bool                                 SurfaceEmission::s_io_buffer_locked = false;


    /*---------------------------------------------------------------------------------------------
     * Constructor
     *-------------------------------------------------------------------------------------------*/ 
    SurfaceEmission::SurfaceEmission(std::string                                species_name,
                                     std::string                                data_file_name,
                                     std::vector<std::string>                   sectors,
                                     const std::shared_ptr<const AbstractGrid>& model_grid,
                                     const std::string&                         map_file)
        : m_species_name(species_name), m_data_file(data_file_name), m_sectors(sectors)            
    {
        using Kokkos::View;

        //std::cout<< "in SurfaceEmission (Constructor)..." << std::endl;
        // Spatial remapping may store shallow copy of the grid
        m_horiz_remapper = create_horiz_remapper(model_grid,map_file);

        // Create the data reader
        const int n_fields = m_horiz_remapper->get_num_fields();
        std::vector<Field> field_emiss_sectors;
        field_emiss_sectors.reserve(n_fields);
        for(int i = 0; i < n_fields; ++i) 
            field_emiss_sectors.push_back(m_horiz_remapper->get_src_field(i));

        m_data_reader = std::make_shared<AtmosphereInput>(m_data_file, 
                                                          m_horiz_remapper->get_src_grid(),
                                                          field_emiss_sectors, 
                                                          true);

        // Initialize views for I/O data
        int ncols    = model_grid->get_num_local_dofs();
        int nsectors = m_sectors.size();
        m_data_start = view_2d("AllSectors", nsectors, ncols);
        m_data_end   = view_2d("AllSectors", nsectors, ncols);
        m_data_out   = view_2d("AllSectors", 1, ncols);




        // To reduce the memory requirement, minimal buffers are used instead of a scorpio 
        // AtmosphereInput, which requires each input/variable to have a preallocated Field
        // (as far as I can tell). Input buffers are 1D Views with size "ncols" according
        // to the input file, while remap buffers are 1D Views with size "ndof" according
        // to the grid. Two buffers of each type are used to support time interpolation.
        // Each emission file has 12 snapshots (1 per month).

        m_ncols_gbl = model_grid->get_num_global_dofs();
        m_ncols_loc = model_grid->get_num_local_dofs();
        m_ncols_in  = scorpio::get_dimlen(m_data_file, "ncol");
        m_ntime     = scorpio::get_time_len(m_data_file);
        
        // Allocate buffers
        m_input_buffer1 = View<Real*>("input_buffer1",m_ncols_in);
        m_input_buffer2 = View<Real*>("input_buffer2",m_ncols_in);

        // Only allocate remap buffers if grids are different
        // (Possibly remappers take care of mapping global data to rank-local grid, 
        // in which case even for same grid we should have a buffer sized to fit the local rank.)
        //if (m_ncols_gbl != m_ncols_in)
        //{
            m_remap_buffer1 = View<Real*>("remap_buffer1",m_ncols_loc);
            m_remap_buffer2 = View<Real*>("remap_buffer2",m_ncols_loc);
        //}

        if (s_io_buffer == nullptr) 
            s_io_buffer = std::make_shared<Kokkos::View<Real*>>("SurfaceEmission_io_buffer",m_ncols_in);
    } // Constructor    
    
    

    SurfaceEmission::~SurfaceEmission() {}
    //{scorpio::release_file(m_data_file);} // Destructor (not necessary??)






    /*---------------------------------------------------------------------------------------------
     * std::shared_ptr<AbstractRemapper> create_horiz_remapper(...)
     *
     * Internal method to create the spatial interpolation remapper
     *-------------------------------------------------------------------------------------------*/ 
    std::shared_ptr<AbstractRemapper> 
    SurfaceEmission::create_horiz_remapper(const std::shared_ptr<const AbstractGrid> &model_grid,
                                           const std::string                         &map_file)
    {
        using grid_type = scream::GridsManager::grid_type;
    
        // Declare remapper (output)   
        std::shared_ptr<AbstractRemapper> remapper;

        // Read number of columns in data file
        int ncols_data;
        try
        { 
            scorpio::register_file(m_data_file, scorpio::Read);
            //const 
            ncols_data = scorpio::get_dimlen(m_data_file, "ncol");
            std::cout<< "    ncols_data = " << ncols_data << std::endl;
            scorpio::release_file(m_data_file);
        }
        catch (...)
        {
            std::cerr << "Error in SurfaceEmission::create_horiz_remapper:" 
                      << "Could not read parameter 'ncol' from data file: "
                      << m_data_file << ". Data file may be in an older EAM format. "
                      << "Consider remapping data prior to the run."
                      << std::endl;          
        }    

        // number of columns in the model grid    
        const int ncols_model = model_grid->get_num_global_dofs();

        // Create a shallow copy of the grid to interpolate to
        // Shallow copies do not copy the entire grid, they just share
        // a view, so they are cheap (just don't modify the grid)
        std::shared_ptr<grid_type> horiz_interp_target_grid =
            model_grid->clone("surf_emission_horiz_interp_target_grid", true); 
        //auto horiz_interp_target_grid =
        //    model_grid->clone("surf_emission_horiz_interp_target_grid", true);     
        
        // If the data file and model use the same grid, no interpolation is necessary
        if(ncols_data == ncols_model) 
        {
            // Interpolation target grid is an alias of the model grid
            remapper = std::make_shared<IdentityRemapper>(horiz_interp_target_grid, 
                                        IdentityRemapper::SrcAliasTgt);
        } 
        else 
        {
            std::cout << "    WARNING: runtime remapping process may be bugged. "
                      << "Conisder remapping data to the data to the desired grid prior to the run." 
                      << std::endl;
            // Only refining remappers allowed (I don't know why)
            EKAT_REQUIRE_MSG(ncols_data <= ncols_model,
                            "Error! We do not allow to coarsen srfEmiss data to fit "
                            "the model. We only allow\n"
                            "srfEmiss data to be at the same or coarser resolution as "
                            "the model.\n");
            // We must have a valid map file
            EKAT_REQUIRE_MSG(
                map_file != "",
                "ERROR: srfEmiss data is on a different grid than the model one,\n"
                "but srfEmiss_remap_file is missing from srfEmiss parameter "
                "list.");
      
            remapper = std::make_shared<RefiningRemapperP2P>(horiz_interp_target_grid, map_file);
        }
        
        return remapper;
    } // create_horiz_remapper



    /*---------------------------------------------------------------------------------------------
     * void update_data_from_file(...)
     *
     * Adapted from 'srfEmissFunctions::update_srfEmiss_data_from_file' function from the 'mam'
     * process. For improved clarity, this has been changed to a method of SurfaceEmission.
     * 
     * Parameters:
     *   time_step
     *   time_index
     *-------------------------------------------------------------------------------------------*/
    void SurfaceEmission::update_data_from_file(Field&                 tracer,
                                                const util::TimeStamp& time_step, 
                                                const int              time_index)
    {
        //std::cout << "in SurfaceEmission::update_data_from_file..." << std::endl;
        start_timer("EAMxx::SurfaceEmission::update_data_from_file");

        // For now, data from the emissions file at one time level is simply summed into
        // one of the input buffers. Once development has progressed in the
        // GOCART bridge, this should be changed to include time interpolation
        read_all_sectors(m_input_buffer1, m_data_file, m_sectors, time_index);

        // Add emissions to tracer at lowest elevation cell
        Kokkos::View<Real**> data = tracer.get_view<Real**>();
        for (int i=0; i<m_ncols_gbl; ++i) data(i,0) += m_input_buffer1(i);

        // Somehow the remapper will need to interpolate these to the grid. For
        // now, we are using single-processor and the same grid, so no remapping
        // is necessary

        //scorpio::release_file(m_filename); << not sure if this is done here or just at finalize


        stop_timer("EAMxx::SurfaceEmission::update_data_from_file");


        #if 0
        std::cout << "in SurfaceEmission::update_data_from_file..." << std::endl;
        start_timer("EAMxx::SurfaceEmission::update_data_from_file");

        // 1. Read from file
        start_timer("EAMxx::SurfaceEmission::update_data_from_file::read_data");
        m_data_reader->read_variables(time_index);
        stop_timer("EAMxx::SurfaceEmission::update_data_from_file::read_data");

        // 2. Run the horiz remapper (it is a do-nothing op if srfEmiss data is on
        // same grid as model)
        start_timer("EAMxx::SurfaceEmission::update_data_from_file::horiz_remap");
        m_horiz_remapper->remap_fwd();
        stop_timer("EAMxx::SurfaceEmission::update_data_from_file::horiz_remap");

        // 3. Copy from the tgt field of the remapper into the srfEmiss_data, padding
        // data if necessary
        start_timer("EAMxx::SurfaceEmission::update_data_from_file::copy_and_pad");
        // Recall, the fields are registered in the order: ps, ccn3, g_sw, ssa_sw,
        // tau_sw, tau_lw

        // Read fields from the file
        const int num_fields = m_horiz_remapper->get_num_fields();
        for(int i = 0; i < num_fields; ++i) 
        {
            /*
            template<typename DT, HostOrDevice HD = Device>
            get_view_type<DT,HD>
            get_view () const;
            */
            
            // Copy each field into its own 'subview' stored in the 'm_data_end' view

            // 'sector' and 'emiss' are view types. These are left as auto as the Kokkos
            // types are similarly unclear due to heavy template metaprogramming
            auto sector      = m_horiz_remapper->get_tgt_field(i).get_view<const Real *>();
            int n1 = sector.extent(0);
            std::cout << "sector.extent(0) = " << n1 << std::endl;
            int n2 = sector.extent(1);
            std::cout << "sector.extent(1) = " << n2 << std::endl;

            float sum = 0.0;
            for (int ival=0; i<n1; ++i) sum+=sector(ival);
            std::cout << "sum over sector = " << sum << std::endl;


            const auto emiss = Kokkos::subview(m_data_end, i, Kokkos::ALL());
            Kokkos::deep_copy(emiss, sector);
            // deep_copy(dest, src)
            // I think the argument has to be 'const' for some Kokkos metaprogramming witchcraft
            // but the data itself must have a 'mutable' tag
        }

        Kokkos::fence();
        stop_timer("EAMxx::SurfaceEmission::update_data_from_file::copy_and_pad");
        stop_timer("EAMxx::SurfaceEmission::update_data_from_file");
        #endif
    }



    void my_update_from_file()
    {
        //void read_var (const std::string &filename, const std::string &varname, T* buf, const int time_index)
    }




    /*---------------------------------------------------------------------------------------------
     * void read_all_sectors(...)
     *
     * Reads all emissions sectors from a netcdf file into the data buffer.
     * 
     * Parameters:
     *   buffer.........array to store sum total of all emissions sectors
     *   file...........name of the data file
     *   sectors........names of emissions sectors
     *   time_index.....time index (typically month) to read from
     *-------------------------------------------------------------------------------------------*/
    void SurfaceEmission::read_all_sectors(Kokkos::View<Real*>&            buffer,
                                           const std::string&              file, 
                                           const std::vector<std::string>& sectors,
                                           const int                       time_index)
    {
        // This should match "ncols" from the file due to the Constructor
        const int N = m_ncols_in;

        // Reset the buffer
        for (int i=0; i<N; ++i) buffer(i)=0.0;

        // Lock the IO buffer
        if (s_io_buffer_locked)
        {
            std::cerr << "Error in SurfaceEmission::read_all_sectors: IO buffer is in use." 
                      << std::endl;        
        }
        else {s_io_buffer_locked = true;}

        // Resizing the buffer may be necessary
        if (s_io_buffer->extent(0) < N) resize(*s_io_buffer, N);

        // Read each sector, summing into the input buffer
        for (const std::string& sector : sectors)
        {
            scorpio::read_var(file, sector, s_io_buffer->data(), time_index);
            for (int i=0; i<N; ++i) buffer(i) += (*s_io_buffer)(i);
        }

        // Unlock the IO buffer to finish
        s_io_buffer_locked = false;
    }



}