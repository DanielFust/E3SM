/*-------------------------------------------------------------------------------------------------
 * gocart_data_structs.hpp
 *
 * This header file contains definitions of lightweight data structures used internally by the
 * GOCART AtmosphereProcess
 *-----------------------------------------------------------------------------------------------*/
#ifndef GOCART_SURFACE_EMISSION_HPP
#define GOCART_SURFACE_EMISSION_HPP

#include <string>
#include <vector>
#include <memory> // shared_ptr

#include "share/grid/remap/abstract_remapper.hpp"
#include "share/io/scorpio_input.hpp"

#include "physics/gocart/eamxx_gocart_data_structs.hpp" // TimeState

namespace scream::gocart
{
    /*---------------------------------------------------------------------------------------------
     * class SurfaceEmission
     *
     * The SurfaceEmission struct aggregates data required for I/O of surface emission data.
     * It is based on the 'surf_emiss_ struct' of the MAMSrfOnlineEmiss AtmosphereProcess
     * which requires similar surface emission data
     *-------------------------------------------------------------------------------------------*/
    class SurfaceEmission
    {
        //using KT            = ekat::KokkosTypes<DefaultDevice>;
        //using view_2d = typename KT::template view_2d<Real>;
        using view_2d = ekat::KokkosTypes<DefaultDevice>::view_2d<Real>;
      //protected:  
      public: // For debugging
        // species name
        std::string m_species_name;

        // Data file name
        std::string m_data_file;

        // Sector names in file
        std::vector<std::string> m_sectors;

        // Spatial interpolation
        std::shared_ptr<AbstractRemapper> m_horiz_remapper;
        
        // Reads netcdf input data
        std::shared_ptr<AtmosphereInput>  m_data_reader;

        // Simple struct for tracking date/time
        TimeState m_time_state;

        // Inputs and outputs are both 2D Views with the same number of columns
        // however, when working with inputs, one must consider the number of
        // sectors.

        // Inputs
        view_2d m_data_start;
        view_2d m_data_end;
        // Output
        view_2d m_data_out;


        // Okay, I don't know what the MAM crew is doing, but it doesn't seem to work so I'm doing my own thing now
        // Data buffers for 2 months for time interpolation
        Kokkos::View<Real*> m_input_buffer1;
        Kokkos::View<Real*> m_input_buffer2;

        Kokkos::View<Real*> m_remap_buffer1;
        Kokkos::View<Real*> m_remap_buffer2;

        // Only one SurfaceEmission object per rank reads or writes at a given instant,
        // so they may share an input buffer for reading the files
        static std::shared_ptr<Kokkos::View<Real*>> s_io_buffer;
        static bool                                 s_io_buffer_locked;

        int m_ncols_in;   // columns in input file
        int m_ncols_gbl;  // columns in global physics grid
        int m_ncols_loc;  // columns in physics grid (on this rank)
        int m_ntime;      // number of time snapshots
      
      // ------------ Public methods ---------------  
      public:  

        SurfaceEmission(std::string                                species_name,
                        std::string                                data_file_name,
                        std::vector<std::string>                   sectors,
                        const std::shared_ptr<const AbstractGrid>& model_grid,
                        const std::string&                         map_file);

        ~SurfaceEmission();

        // Deallocate shared buffers
        static void finalize() {s_io_buffer = nullptr;}

        void update_data_from_file(Field&, const util::TimeStamp &, const int);

        // Simple Getters
        std::string                     species_name() const {return m_species_name;}
        std::string                     data_file()    const {return m_data_file;}
        const std::vector<std::string>& sectors()      const {return m_sectors;}

      // ----------- Protected methods -------------
      protected:

        // Method to create the spatial interpolation remapper
        std::shared_ptr<AbstractRemapper> create_horiz_remapper(
            const std::shared_ptr<const AbstractGrid>& model_grid,
            const std::string& map_file);

      // Reads/sums all sectors from emission file into a 1D buffer
      void read_all_sectors(Kokkos::View<Real*>&            buffer,
                            const std::string&              file, 
                            const std::vector<std::string>& sectors,
                            const int                       time_index);      
        
      

    };


                  
}
#endif // GOCART_SURFACE_EMISSION_HPP




