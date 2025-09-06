/*
 * Add any #include statements pointing to headers needed by this process.
 * Typical headers would be
 *     eamxx_$PROCESS_process_interface.hpp  (the eamxx interface header file for this process)
 *     physics/share/physics_constants.hpp   (a header storing all physical constants in eamxx)
 *     $PROCESS.hpp                          (a header file for the process itself)
 */
//#include "physics/gocart/eamxx_gocart_process_interface.hpp"

#include <functional> // std::ref
#include <algorithm>
#include <iostream>

#include "physics/gocart/eamxx_gocart_process_interface.hpp"
#include "physics/gocart/sslt-interface/sslt_sections.hpp"

#include "eam-interface/eam_bridge.hpp"
#include "sslt-interface/sslt_sections.hpp"
#include "setsox-interface/mo_setsox.hpp"
#include "gocart-interface/mo_gocart.hpp"

#include "ekat/ekat.hpp"
#include "share/util/eamxx_data_interpolation.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"
#include "physics/share/physics_constants.hpp"

//#include "externals/haero/haero/constants.hpp"


namespace scream
{

    void gocart_test_io(gocart::SurfaceEmission& surf_emiss)
    {
        std::string species = surf_emiss.species_name();
        std::string file    = surf_emiss.data_file();
        int time_index      = -1;

        int time_snapshots = scorpio::get_time_len(file);
        std::cout << "number of time snapshots: " << time_snapshots << std::endl;
        int buffer_size = scorpio::get_dimlen(file, "ncol");
        std::cout << "has_dim(file, ncol):" << scorpio::has_dim(file,"ncol") << std::endl;

        auto iotype = scorpio::str2iotype("default");
        scorpio::register_file(file,scorpio::Read,iotype);
        std::cout << "file registered..." << std::endl;
        Kokkos::View<Real*> buffer("emission buffer",buffer_size);
        std::cout << "buffer allocated..." << std::endl;

        for (const std::string& sector : surf_emiss.sectors())
        {
            int N = scorpio::get_dimlen(file, "ncol");
            std::cout << "dimension size: " << N << std::endl;;
            scorpio::read_var(file, sector, buffer.data(), time_index);
            
            float sum = 0.0;
            for (int i=0; i<N; ++i) sum += buffer(i);

            std::cout << "sector: " << sector << ", total: " << sum << std::endl;
        }

        scorpio::release_file(file);
        
    }

    /*---------------------------------------------------------------------------------------------
     * The Constructor for this interface
     *
     * Inputs:
     *     comm - an EKAT communication group
     *     params - a parameter list of options for the process.
     *-------------------------------------------------------------------------------------------*/
    GOCART::GOCART (const ekat::Comm& comm, const ekat::ParameterList& params)
        : AtmosphereProcess(comm,params)
    {
        // The 'params' variable will hold all runtime options.  
        // Note that `params.get<X>("Y") is the syntax to get the parameter labeled "Y" from the list
        // which is of type X.  
        //
        // ex:
        // m_Y = params.get<X>("Y");

        // -- Elevated Emissions Files --    
        EKAT_REQUIRE_MSG(m_params.isParameter("gocart_so2_elevated_emiss_file_name"),
            "ERROR: gocart_so2_elevated_emiss_file_name is missing from GOCART parameter list.");

        // -- Surface Emissions Files --    
        EKAT_REQUIRE_MSG(m_params.isParameter("gocart_so2_surf_emiss_file_name"),
            "ERROR: gocart_so2_surf_emiss_file_name is missing from GOCART parameter list.");

        // -- Surface Emission Remap File --        
        EKAT_REQUIRE_MSG(m_params.isParameter("srf_remap_file"),
            "ERROR: srf_remap_file is missing from GOCART parameter list.");

    }

    /*-----------------------------------------------------------------------------------------------
     * set_grids(grids_manager)
     *
     * set_grids establishes which grid this process will be run on.  
     * It is also where all fields (variables) that the process will need access to.
     *
     * Inputs:
     *     grids_manager - a grids manager object which stores all grids used in this simulation
     *-------------------------------------------------------------------------------------------*/
    void GOCART::
    set_grids (const std::shared_ptr<const GridsManager> grids_manager)
    {
        std::cout << "In GOCART::set_grids" << std::endl;
        using namespace gocart;
        using ekat::units::Units;
        
        constexpr Units kg     = ekat::units::kg;
        constexpr Units K      = ekat::units::K;
        constexpr Units m      = ekat::units::m;
        constexpr Units m3     = pow(m,3); 
        constexpr Units s      = ekat::units::s;
        constexpr Units Pa     = ekat::units::Pa;
        constexpr Units nondim = Units::nondimensional();


        // Some typical namespaces used in set_grids,
        //   using namespace ShortFieldTagsNames;
        //   using PC = scream::physics::Constants<Real>;

        // Specify which grid this process will act upon, typical options are "Dynamics" or "Physics".
        //auto m_grid = grids_manager->get_grid("Physics");
        m_grid = grids_manager->get_grid("physics");
        ncol_ = m_grid->get_num_local_dofs();       // Number of columns on this rank
        nlev_ = m_grid->get_num_vertical_levels();  // Number of levels per column


        // EAM infrastructure muse be initialized as it contains tracer info
        // EAM user input namelist for setup of wet/dry deposition
        std::string eam_namelist = m_params.get<std::string>("gocart_eam_namelist_file","");
        int nlon=-1; int nlat=-1; // Hopefully old grid structure can be pruned
        eam_bridge::eam_bridge_init(ncol_, nlev_, nlon, nlat, eam_namelist);


        // Initialize m_surface_emissions vector
        init_surface_emissions();

        
        // Add gocart tracers
        // Not sure which units to use (guessing kg/kg for now)
        // NOTE: These will likely need to change to "Updated", however in the standalone test,
        //       the executable seg faults when reading the IC file if it does not contain these tracers
        //       (which none do as far as I am aware)
        register_tracers(m_grid);
        

        // -- Layouts --
        const FieldLayout scalar2d = m_grid->get_2d_scalar_layout();
        const FieldLayout scalar3d = m_grid->get_3d_scalar_layout(true);

        // -- Surface Fields --
        add_field<Required>("sst", scalar2d, K, m_grid->name()); // Sea surface temperature [K]
        add_field<Required>("ocnfrac", scalar2d, nondim, m_grid->name()); // Ocean fraction in cell
        add_field<Computed>("u10_cubed", scalar2d, m/s, m_grid->name()); // 10m windspeed "cubed" 

        // -- 3D Fields --
        add_field<Required>("p_mid",            scalar3d, Pa,    m_grid->name()); // Midpoint pressure
        add_field<Required>("T_mid",            scalar3d, K,     m_grid->name()); // Midpoint temperature
        add_field<Required>("cldfrac_tot",      scalar3d, nondim,m_grid->name()); // Cloud fraction
        add_field<Required>("AtmosphereDensity",scalar3d, kg/m3, m_grid->name()); // Atmospheric density <-- Need to verify units
        add_field<Required>("qv",               scalar3d, kg/kg, m_grid->name()); // Specific humidity
        add_field<Required>("nc",               scalar3d, 1/kg,  m_grid->name()); // number concentration of cloud liquid particles [#/kg]
        add_field<Required>("pseudo_density",   scalar3d, Pa,    m_grid->name()); // Layer thickness (pdel) [Pa] at midpoints
        add_tracer<Required>("qc", m_grid, kg/kg);                                // cloud water mixing ratio

        // GOCART must be able to retrieve the windspeed at the midpoint of the grid cell
        constexpr int pack_size = SPack<Real>::n; // i.e. SCREAM_SMALL_PACK_SIZE
        FieldLayout vector3d_mid = m_grid->get_3d_vector_layout(true,2); // Packed as (u,v) components
        add_field<Updated>("horiz_winds", vector3d_mid, m/s, m_grid->name(), pack_size);
    }

    /*-----------------------------------------------------------------------------------------------
     * intialize_impl(run_type)
     *
     * called once for each process at initialization of EAMxx.  This impl can be defined with any
     * actions or functions that are needed at initialization. 
     *
     * Inputs:
     *     run_type - an enum which describes the run type.  Initial or Restart
     *
     * can also be empty
     *-------------------------------------------------------------------------------------------*/
    void GOCART::initialize_impl (const RunType /* run_type */)
    {
        using namespace gocart;

        // NOTE: run_type tells us if this is an initial or restarted run,
        std::cout << "in GOCART::initialize_impl..." << std::endl;


        // Build map between SurfaceEmissions and tracer species for more efficient looping
        // Must be initialized here, as Fields cannot be accessed during 'set_grids'
        for (SurfaceEmission& emission : m_surface_emissions)
        {
            Field& species = get_field_out(emission.species_name());
            m_surf_emiss_tracer_pairs.emplace_back(std::make_tuple(std::ref(emission),
                                                                   std::ref(species)));
        }
        std::cout << "map built... " << std::endl;


        // -- Allocate the giant tracer arrays required by bridged EAM routines --
        gas_pcnst_ = eam_bridge::gas_pcnst();
        m_eamf90_tracers        = Kokkos::View<Real***>("eamf90_tracers",        gas_pcnst_, nlev_, ncol_);
        m_eamf90_cloud_aerosols = Kokkos::View<Real***>("eamf90_cloud_aerosols", gas_pcnst_, nlev_, ncol_); 

        // -- Allocate other EAM-specific fields --
        m_eamf90_mbar   = Kokkos::View<Real**>("eamf90_mbar",   nlev_, ncol_);
        m_eamf90_T_mid  = Kokkos::View<Real**>("eamf90_T_mid",  nlev_, ncol_);
        m_eamf90_p_mid  = Kokkos::View<Real**>("eamf90_p_mid",  nlev_, ncol_);
        m_eamf90_cldfrc = Kokkos::View<Real**>("eamf90_cldfrc", nlev_, ncol_);
        m_eamf90_cldnum = Kokkos::View<Real**>("eamf90_cldnum", nlev_, ncol_);
        m_eamf90_xhnm   = Kokkos::View<Real**>("eamf90_xhnm",   nlev_, ncol_);
        m_eamf90_pdel   = Kokkos::View<Real**>("eamf90_pdel",   nlev_, ncol_);
        m_eamf90_lwc    = Kokkos::View<Real**>("eamf90_lwc",    nlev_, ncol_);

        // -- Initialize sslt_sections (sea salt) module --
        sslt_sections::sslt_sections_init();

        // -- Initialize mo_setsox (cloud chemistry) module --
        mo_setsox::sox_inti_c2f(true,false,false); // Not sure what the 2nd two arguments should actually be

        // Allocate an array to hold seasalt emissions
        const int sslt_nsec = sslt_sections::nsections();
        m_sslt_flux = view_2d<Real>("sslt_flux", sslt_nsec, ncol_);

        // -- Update surface emissions from file --
        const TimeStamp& time_stamp = start_of_step_ts();
        const int        curr_month = time_stamp.get_month() - 1; // 0-indexed
        std::cout << "m_surface_emissions.size() = " << m_surface_emissions.size() << std::endl;
        
        for (std::tuple<SurfaceEmission&, Field&>& pair : m_surf_emiss_tracer_pairs)
        {
            SurfaceEmission& emission = std::get<0>(pair);
            Field&           species  = std::get<1>(pair);
            //emission.update_data_from_file(species,time_stamp,curr_month);
        }

        const int gas_pcnst = eam_bridge::gas_pcnst();
        for (int i=0; i<gas_pcnst; ++i)
        {
            std::string spc_name  = eam_bridge::gas_species_name(i);
            int         spc_index = eam_bridge::get_spc_ndx(spc_name);
            std::cout << "name: " << spc_name << "  index:" << spc_index << "\n";
        }


    } //GOCART::initialize_impl

    /*---------------------------------------------------------------------------------------------
     * run_impl(dt)
     *
     * The main run call for the process.  This is where most of the interface to the underlying 
     * process takes place.  This impl is called every timestep. 
     *
     * Inputs:
     *     dt - the timestep for this run step.
     *-------------------------------------------------------------------------------------------*/
    void GOCART::run_impl (const double dt)
    {
        using namespace gocart;
        std::cout << "In GOCART::run_impl..." << std::endl;

        // -- Update surface emissions from file --
        const TimeStamp& time_stamp = start_of_step_ts();
        const int        curr_month = time_stamp.get_month() - 1; // 0-indexed
        
        for (std::tuple<SurfaceEmission&, Field&>& pair : m_surf_emiss_tracer_pairs)
        {
            SurfaceEmission& emission = std::get<0>(pair);
            Field&           species  = std::get<1>(pair);
            //emission.update_data_from_file(species,time_stamp,curr_month);
        }


        // Transfer advected EAMxx tracer data to EAM-compatible array
        const int ncol=ncol_;           // <-- Compilers are probably smart enough to deduce that the bounds won't change
        const int nlev=nlev_;           //     during the loop, but just in case, the 'const' local variables make it 
        const int gas_pcnst=gas_pcnst_; //     clear that more aggressive optimizations are allowed.
        for (int ispc=0; ispc<gas_pcnst; ++ispc)
            for (int icol=0; icol<ncol; ++icol)
                for (int ilev=0; ilev<nlev; ++ilev)
                {
                    std::string spc_name               = eam_bridge::gas_species_name(ispc);
                    Kokkos::View<const Real**> species = get_field_in(spc_name).get_view<const Real**>();

                    // Interpreted in Fortran as q(icol,ilev,ispc)
                    m_eamf90_tracers(ispc,ilev,icol) = species(icol,ilev);
                }

        // Populate other EAM-compatible arrays and recover molecular mass of atmosphere (amu)
        // Note: The conversion from MMR -> VMR could be done on the EAMxx side
        //       using `calculate_vmr_from_mmr` from `PhysicsFunctions`.
        //       This could be slightly more efficient, though we would need to 
        //       retrieve the molecular weight of each tracer species (`adv_mass`)
        //       from the bridged EAM. For now, the conversion is done on the EAM
        //       side for consistency.
        const Kokkos::View<const Real**> T_mid          = get_field_in("T_mid").get_view<const Real**>();
        const Kokkos::View<const Real**> p_mid          = get_field_in("p_mid").get_view<const Real**>();
        const Kokkos::View<const Real**> cldfrc         = get_field_in("cldfrac_tot").get_view<const Real**>(); // <-- could be 'cldfrac_liq'
        const Kokkos::View<const Real**> cldnum         = get_field_in("nc").get_view<const Real**>();
        const Kokkos::View<const Real**> qc             = get_field_in("qc").get_view<const Real**>();
        const Kokkos::View<const Real**> qv             = get_field_in("qv").get_view<const Real**>();
        const Kokkos::View<const Real**> atm_density    = get_field_in("AtmosphereDensity").get_view<const Real**>();
        const Kokkos::View<const Real**> pseudo_density = get_field_in("pseudo_density").get_view<const Real**>(); 
        constexpr Real air_mol_weight = scream::physics::Constants<Real>::MWdry;
        constexpr Real avogadro       = 6.022214076e23; // from haero::Constants
        for (int icol=0; icol<ncol; ++icol)
            for (int ilev=0; ilev<nlev; ++ilev)
            {
                m_eamf90_T_mid(ilev,icol)  = T_mid(icol,ilev);
                m_eamf90_p_mid(ilev,icol)  = p_mid(icol,ilev);
                m_eamf90_cldfrc(ilev,icol) = cldfrc(icol,ilev);
                m_eamf90_cldnum(ilev,icol) = cldnum(icol,ilev);
                m_eamf90_lwc(ilev,icol)    = qc(icol,ilev);
                m_eamf90_pdel(ilev,icol)   = pseudo_density(icol,ilev);
                m_eamf90_mbar(ilev,icol) = air_mol_weight / (1.0-qv(icol,ilev));
                m_eamf90_xhnm(ilev,icol) = atm_density(icol,ilev) * avogadro / (air_mol_weight * 1.0e06); // density in (#/cm^3)
            }

 
        
        // Convert from MMR to VMR (in place) for chemical processes
        eam_bridge::mmr2vmr(m_eamf90_tracers.data(), // mmr
                            m_eamf90_tracers.data(), // -> vmr
                            m_eamf90_mbar.data(),    // mean wet atmospheric mass (amu)
                            ncol_);
               

        //view_3d<const SPack<Real>> horiz_winds = get_field_in("horiz_winds").get_view<SPack<Real>***>();
        view_3d<SPack<Real>> horiz_winds   = get_field_out("horiz_winds").get_view<SPack<Real>***>();
        view_1d<const Real>  sea_surf_temp = get_field_in("sst").get_view<const Real*>(); 
        view_1d<Real>        u10_cubed     = get_field_out("u10_cubed").get_view<Real*>();

        update_10m_cubed_windspeed(u10_cubed,horiz_winds);
        std::cout << "returned from update_10m_cubed_windspeed...\n";

        sslt_sections::sslt_fluxes(m_sslt_flux, sea_surf_temp, u10_cubed);
        std::cout << "returned from sslt_fluxes...\n";

        /*---- Cloud Aqueous Chemistry Section ----*/
        /*
        mo_setsox::setsox_c2f(ncol_, dt, 
                              m_eamf90_p_mid.data(),  
                              m_eamf90_pdel.data(),
                              m_eamf90_T_mid.data(),   
                              m_eamf90_mbar.data(),   
                              m_eamf90_lwc.data(),
                              m_eamf90_cldfrc.data(), 
                              m_eamf90_cldnum.data(), 
                              m_eamf90_xhnm.data(), 
                              m_eamf90_cloud_aerosols.data(), 
                              m_eamf90_tracers.data());
        std::cout << "returned from setsox_c2f...\n";
        */
        mo_setsox::setsox(dt, 
                          m_eamf90_p_mid,  
                          m_eamf90_pdel,
                          m_eamf90_T_mid,   
                          m_eamf90_mbar,   
                          m_eamf90_lwc,
                          m_eamf90_cldfrc, 
                          m_eamf90_cldnum, 
                          m_eamf90_xhnm, 
                          m_eamf90_cloud_aerosols, 
                          m_eamf90_tracers);


        // Convert (in place) back to MMR from VMR for advection
        eam_bridge::mmr2vmr(m_eamf90_tracers.data(), // vmr
                            m_eamf90_tracers.data(), // -> mmr
                            m_eamf90_mbar.data(),    // mean wet atmospheric mass (amu)
                            ncol_);

        // Transfer EAM tracer data back to EAMxx Fields
        for (int ispc=0; ispc<gas_pcnst; ++ispc)
            for (int icol=0; icol<ncol; ++icol)
                for (int ilev=0; ilev<nlev; ++ilev)
                {
                    std::string spc_name         = eam_bridge::gas_species_name(ispc);
                    Kokkos::View<Real**> species = get_field_out(spc_name).get_view<Real**>();

                    species(icol,ilev) = m_eamf90_tracers(ispc,ilev,icol);
                }

    } // run_impl

/*-----------------------------------------------------------------------------------------------
 * finalize_impl()
 *
 * Called at the end of the simulation, handles all finalization of the process.
 *
 * In most cases this is left blank, as EAMxx takes care of most finalization steps.
 * But just in case a process has specific needs the option is available.
 */
void GOCART::finalize_impl ()
{
    std::cout <<  "In GOCART::finalize_impl..." << std::endl;
    SurfaceEmission::finalize();
  // Usually blank
}
/*-----------------------------------------------------------------------------------------------*/


    /*---------------------------------------------------------------------------------------------
     * void register_tracers(...)
     *
     * Adds required tracers to the physics grid. This method is called only during the 
     * 'set_grids' and is separated to provide greater clarity during initialization of the 
     * GOCART process.
     * 
     * Parameters:
     *   grid.........physics grid
     *-------------------------------------------------------------------------------------------*/
    void GOCART::register_tracers(std::shared_ptr<const AbstractGrid> grid)
    {
        using ekat::units::Units;
        using namespace gocart;
        constexpr Units kg = ekat::units::kg;

        std::cout << "in register_tracers...\n";

        const int gas_pcnst = eam_bridge::gas_pcnst();
        for (int i=0; i< gas_pcnst; ++i)
        {
            std::string spc_name = eam_bridge::gas_species_name(i);
            add_tracer<Updated>(spc_name, m_grid, kg/kg);
        }

        /*
        // Should already be taken care of by the loop
        add_tracer<Updated>("C10H16", m_grid, emiss_units);
        add_tracer<Updated>("C2H4", m_grid, emiss_units);
        add_tracer<Updated>("C2H6", m_grid, emiss_units);
        add_tracer<Updated>("C3H8", m_grid, emiss_units);
        add_tracer<Updated>("CH2O", m_grid, emiss_units);
        add_tracer<Updated>("CH3CHO", m_grid, emiss_units);
        add_tracer<Updated>("CH3COCH3", m_grid, emiss_units);
        add_tracer<Updated>("CO", m_grid, emiss_units);
        add_tracer<Updated>("dms", m_grid, emiss_units);
        add_tracer<Updated>("e90", m_grid, emiss_units);
        add_tracer<Updated>("isop", m_grid, emiss_units);
        add_tracer<Updated>("isop_vbs", m_grid, emiss_units);
        add_tracer<Updated>("no", m_grid, emiss_units);
        add_tracer<Updated>("SO2", m_grid, emiss_units);
        add_tracer<Updated>("soag0", m_grid, emiss_units);
        add_tracer<Updated>("bc_a4", m_grid, emiss_units);
        add_tracer<Updated>("num_a1", m_grid, emiss_units);
        add_tracer<Updated>("num_a2", m_grid, emiss_units);
        add_tracer<Updated>("num_a4", m_grid, emiss_units);
        add_tracer<Updated>("pom_a4", m_grid, emiss_units);
        add_tracer<Updated>("so4_a1", m_grid, emiss_units);
        add_tracer<Updated>("so4_a2", m_grid, emiss_units);
        */

        // Additional tracers for gocart
        add_tracer<Updated>("gcso4_a", m_grid,kg/kg);
        add_tracer<Updated>("gcsea_a1",m_grid,kg/kg);
        add_tracer<Updated>("gcsea_a2",m_grid,kg/kg);
        add_tracer<Updated>("gcsea_a3",m_grid,kg/kg);
        add_tracer<Updated>("gcsea_a4",m_grid,kg/kg);
        add_tracer<Updated>("gcsea_a5",m_grid,kg/kg);
        add_tracer<Updated>("gcdst_a1",m_grid,kg/kg);
        add_tracer<Updated>("gcdst_a2",m_grid,kg/kg);
        add_tracer<Updated>("gcdst_a3",m_grid,kg/kg);
        add_tracer<Updated>("gcdst_a4",m_grid,kg/kg);
        add_tracer<Updated>("gcdst_a5",m_grid,kg/kg);
        add_tracer<Updated>("gcbco_a", m_grid,kg/kg);
        add_tracer<Updated>("gcbci_a", m_grid,kg/kg);
        add_tracer<Updated>("gcoco_a", m_grid,kg/kg);
        add_tracer<Updated>("gcoci_a", m_grid,kg/kg);
        add_tracer<Updated>("gcnh3_a", m_grid,kg/kg);
        add_tracer<Updated>("gcnh4_a", m_grid,kg/kg);
        add_tracer<Updated>("gcnit_a", m_grid,kg/kg);
        add_tracer<Updated>("gcsoa_a", m_grid,kg/kg);
    } // register_tracers


    /*---------------------------------------------------------------------------------------------
     * void init_surface_emissions()
     *
     * Adds initializes surface emissions vector. This method is called only during the 
     * 'set_grids' and is separated to provide greater clarity during initialization of the 
     * GOCART process.
     * 
     * Parameters:
     *-------------------------------------------------------------------------------------------*/
    void GOCART::init_surface_emissions()
    {
        std::cout << "In GOCART::init_surface_emissions..." << std::endl;
        using ekat::units::Units;
        constexpr Units kg = ekat::units::kg;
        auto emiss_units = kg/kg;

        // Surface emissions remapping file
        std::string surface_remap_file = m_params.get<std::string>("srf_remap_file", "");

        // -- Initialize vector of surface emissions --
        std::string              species_name;
        std::string              data_file_name;
        std::vector<std::string> sectors;

        // Check that data has not already been initialized already
        if (m_surface_emissions.size() != 0)
        {
            std::cerr << "ERROR! GOCART::init_surface_emissions: surface emission "
                      << "vector has already been initialized!" << std::endl;
        }

        // C10H16 Surface emissions
        species_name   = "C10H16";
        sectors        = {"emiss_bb", "emiss_bio"};
        data_file_name = m_params.get<std::string>("gocart_c10h16_surf_emiss_file_name","");
        m_surface_emissions.emplace_back(SurfaceEmission(species_name,
                                                         data_file_name,
                                                         sectors,
                                                         m_grid,
                                                         surface_remap_file));

        // C2H4 Surface emissions
        species_name   = "C2H4";
        sectors        = {"emiss_anthro", "emiss_bb", "emiss_bio", "emiss_oceans"};
        data_file_name = m_params.get<std::string>("gocart_c2h4_surf_emiss_file_name","");
        m_surface_emissions.emplace_back(SurfaceEmission(species_name,
                                                         data_file_name,
                                                         sectors,
                                                         m_grid,
                                                         surface_remap_file));

        // C2H6 Surface emissions
        species_name   = "C2H6";
        sectors        = {"emiss_anthro", "emiss_bb", "emiss_bio", "emiss_oceans"};
        data_file_name = m_params.get<std::string>("gocart_c2h6_surf_emiss_file_name","");
        m_surface_emissions.emplace_back(SurfaceEmission(species_name,
                                                         data_file_name,
                                                         sectors,
                                                         m_grid,
                                                         surface_remap_file));

        // C3H8 Surface emissions
        species_name   = "C3H8";
        sectors        = {"emiss_anthro", "emiss_bb", "emiss_bio", "emiss_oceans"};
        data_file_name = m_params.get<std::string>("gocart_c3h8_surf_emiss_file_name","");
        m_surface_emissions.emplace_back(SurfaceEmission(species_name,
                                                         data_file_name,
                                                         sectors,
                                                         m_grid,
                                                         surface_remap_file));

        // CH2O Surface emissions
        species_name   = "CH2O";
        sectors        = {"emiss_anthro", "emiss_bb", "emiss_bio"};
        data_file_name = m_params.get<std::string>("gocart_ch2o_surf_emiss_file_name","");
        m_surface_emissions.emplace_back(SurfaceEmission(species_name,
                                                         data_file_name,
                                                         sectors,
                                                         m_grid,
                                                         surface_remap_file));

        // CH3CHO Surface emissions
        species_name   = "CH3CHO";
        sectors        = {"emiss_anthro", "emiss_bb", "emiss_bio"};
        data_file_name = m_params.get<std::string>("gocart_ch3cho_surf_emiss_file_name","");
        m_surface_emissions.emplace_back(SurfaceEmission(species_name,
                                                         data_file_name,
                                                         sectors,
                                                         m_grid,
                                                         surface_remap_file));

        // CH3COCH3 Surface emissions
        species_name   = "CH3COCH3";
        sectors        = {"emiss_anthro", "emiss_bb", "emiss_bio"};
        data_file_name = m_params.get<std::string>("gocart_ch3coch3_surf_emiss_file_name","");
        m_surface_emissions.emplace_back(SurfaceEmission(species_name,
                                                         data_file_name,
                                                         sectors,
                                                         m_grid,
                                                         surface_remap_file));

        // CO Surface emissions
        species_name   = "CO";
        sectors        = {"emiss_anthro", "emiss_bb", "emiss_bio", "emiss_oceans"};
        data_file_name = m_params.get<std::string>("gocart_co_surf_emiss_file_name","");
        m_surface_emissions.emplace_back(SurfaceEmission(species_name,
                                                         data_file_name,
                                                         sectors,
                                                         m_grid,
                                                         surface_remap_file));
        // DMS Surface emissions
        species_name   = "DMS";
        sectors        = {"DMS"};
        data_file_name = m_params.get<std::string>("gocart_dms_surf_emiss_file_name",""); // <-- ncdump can't read dms file for some reason!?!
        m_surface_emissions.emplace_back(SurfaceEmission(species_name,
                                                         data_file_name,
                                                         sectors,
                                                         m_grid,
                                                         surface_remap_file));

        // E90 Surface emissions
        species_name   = "E90";
        sectors        = {"emiss_anthro"};
        data_file_name = m_params.get<std::string>("gocart_e90_surf_emiss_file_name","");
        m_surface_emissions.emplace_back(SurfaceEmission(species_name,
                                                         data_file_name,
                                                         sectors,
                                                         m_grid,
                                                         surface_remap_file));

        // ISOP Surface emissions
        species_name   = "ISOP";
        sectors        = {"emiss_bb", "emiss_bio"};
        data_file_name = m_params.get<std::string>("gocart_isop_surf_emiss_file_name","");
        m_surface_emissions.emplace_back(SurfaceEmission(species_name,
                                                         data_file_name,
                                                         sectors,
                                                         m_grid,
                                                         surface_remap_file));

        // ISOP VBS Surface emissions 
        species_name   = "ISOP_VBS";
        sectors        = {"emiss_bb", "emiss_bio"};
        data_file_name = m_params.get<std::string>("gocart_isop_vbs_surf_emiss_file_name","");
        m_surface_emissions.emplace_back(SurfaceEmission(species_name,
                                                         data_file_name,
                                                         sectors,
                                                         m_grid,
                                                         surface_remap_file));

        // NO Surface emissions
        species_name   = "NO";
        sectors        = {"emiss_anthro", "emiss_bb", "emiss_soils"};
        data_file_name = m_params.get<std::string>("gocart_no_surf_emiss_file_name","");
        m_surface_emissions.emplace_back(SurfaceEmission(species_name,
                                                         data_file_name,
                                                         sectors,
                                                         m_grid,
                                                         surface_remap_file));

        // SO2 Surface emissions
        species_name   = "SO2";
        sectors        = {"AGR", "RCO", "SHP", "SLV", "TRA", "WST"};
        data_file_name = m_params.get<std::string>("gocart_so2_surf_emiss_file_name","");
        m_surface_emissions.emplace_back(SurfaceEmission(species_name,
                                                         data_file_name,
                                                         sectors,
                                                         m_grid,
                                                         surface_remap_file));
                                                        
        // SOAG0 Surface emissions
        species_name   = "SOAG0";
        sectors        = {"AGR", "ENE", "IND", "RCO", "SHP", "SLV", "TRA", "WST"};
        data_file_name = m_params.get<std::string>("gocart_soag0_surf_emiss_file_name","");
        m_surface_emissions.emplace_back(SurfaceEmission(species_name,
                                                         data_file_name,
                                                         sectors,
                                                         m_grid,
                                                         surface_remap_file));

        // BC A4 Surface emissions
        species_name   = "bc_a4";
        sectors        = {"AGR", "ENE", "IND", "RCO", "SHP", "SLV", "TRA", "WST"};
        data_file_name = m_params.get<std::string>("gocart_bc_a4_surf_emiss_file_name","");
        m_surface_emissions.emplace_back(SurfaceEmission(species_name,
                                                         data_file_name,
                                                         sectors,
                                                         m_grid,
                                                         surface_remap_file));

        // NUM A1 Surface emissions
        species_name   = "num_a1";
        sectors        = {"num_a1_SO4_AGR", "num_a1_SO4_SHP", "num_a1_SO4_SLV", "num_a1_SO4_WST"};
        data_file_name = m_params.get<std::string>("gocart_num_a1_surf_emiss_file_name","");
        m_surface_emissions.emplace_back(SurfaceEmission(species_name,
                                                         data_file_name,
                                                         sectors,
                                                         m_grid,
                                                         surface_remap_file));

        // NUM A2 Surface emissions
        species_name   = "num_a2";
        sectors        = {"num_a2_SO4_RCO", "num_a2_SO4_TRA"};
        data_file_name = m_params.get<std::string>("gocart_num_a2_surf_emiss_file_name","");
        m_surface_emissions.emplace_back(SurfaceEmission(species_name,
                                                         data_file_name,
                                                         sectors,
                                                         m_grid,
                                                         surface_remap_file));

        // NUM A4 Surface emissions
        species_name   = "num_a4";
        sectors        = {"num_a1_BC_AGR", "num_a1_BC_ENE", "num_a1_BC_IND", "num_a1_BC_RCO", 
                          "num_a1_BC_SHP", "num_a1_BC_SLV", "num_a1_BC_TRA", "num_a1_BC_WST",
                          "num_a1_POM_AGR", "num_a1_POM_ENE", "num_a1_POM_IND", "num_a1_POM_RCO",
                          "num_a1_POM_SHP", "num_a1_POM_SLV", "num_a1_POM_TRA", "num_a1_POM_WST"};
        data_file_name = m_params.get<std::string>("gocart_num_a4_surf_emiss_file_name","");
        m_surface_emissions.emplace_back(SurfaceEmission(species_name,
                                                         data_file_name,
                                                         sectors,
                                                         m_grid,
                                                         surface_remap_file));

        // POM A4 Surface emissions
        species_name   = "pom_a4";
        sectors        = {"AGR", "ENE", "IND", "RCO", "SHP", "SLV", "TRA", "WST"};
        data_file_name = m_params.get<std::string>("gocart_pom_a4_surf_emiss_file_name","");
        m_surface_emissions.emplace_back(SurfaceEmission(species_name,
                                                         data_file_name,
                                                         sectors,
                                                         m_grid,
                                                         surface_remap_file));

        // SO4 A1 Surface emissions
        species_name   = "so4_a1";
        sectors        = {"AGR", "SHP", "SLV", "WST"};
        data_file_name = m_params.get<std::string>("gocart_so4_a1_surf_emiss_file_name","");
        m_surface_emissions.emplace_back(SurfaceEmission(species_name,
                                                         data_file_name,
                                                         sectors,
                                                         m_grid,
                                                         surface_remap_file));

        // SO4 A2 Surface emissions
        species_name   = "so4_a2";
        sectors        = {"RCO", "TRA"};
        data_file_name = m_params.get<std::string>("gocart_so4_a2_surf_emiss_file_name","");
        m_surface_emissions.emplace_back(SurfaceEmission(species_name,
                                                         data_file_name,
                                                         sectors,
                                                         m_grid,
                                                         surface_remap_file));

    } // init_surface_emissions




    /*---------------------------------------------------------------------------------------------
     * void update_10m_cubed_windspeed(...)
     *
     * Computes the 10m windspeed (accoding to curve fit by Tie and Seinfeld and Pandis, p.859)
     * "cubed" (to the power of 3.41 according to according to Gong et al., 1997) needed for
     * seasalt emissions model
     * 
     * Arguments:
     *   u10_cubed.....2D Field of 10m horizontal windspeed "cubed" 
     *   horiz_winds...3D Field of (u,v) velocity vector
     *-------------------------------------------------------------------------------------------*/
    void GOCART::update_10m_cubed_windspeed(view_1d<Real>&                    u10_cubed,
                                            const_view_3d<SPack<Real>> const& horiz_winds)
    {
        using std::pow;
        using std::sqrt;
        using std::log;

        const Real z0 = 0.0001;  // (m) roughness length over oceans--from ocean model
        
        auto elevation_data = m_grid->get_geometry_data("lev").get_view<const Real*>();
        const Real z        = elevation_data(0);      // First grid cell midpoint elevation
        const Real coef_10m = log(10.0/z0)/log(z/z0); // coefficient to estimate 10m windspeed from first gridcell

        const int ncol = ncol_;
        for (int icol=0; icol<ncol; ++icol)
        {
            // Get horizontal windspeed from lowest cell center
            Real u = horiz_winds(icol,0,0)[0];
            Real v = horiz_winds(icol,1,0)[0];

            // Estimate at 10m (Tie and Seinfeld and Pandis)
            Real u10 = sqrt(u*u + v*v) * coef_10m;

            // exponent from Gong et al., 1997
            u10_cubed(icol) = pow(u10,3.41);
        }
    }




} // namespace scream
