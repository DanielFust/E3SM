#ifndef EAMXX_GOCART_PROCESS_INTERFACE_HPP
#define EAMXX_GOCART_PROCESS_INTERFACE_HPP



#include <vector>
#include "share/atm_process/atmosphere_process.hpp"
#include "share/util/eamxx_data_interpolation.hpp"

#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_workspace.hpp"

#include "surface-emission/eamxx_gocart_surface_emission.hpp"

namespace scream
{

    /*---------------------------------------------------------------------------------------------
     * class GOCART
     *
     * This class bridges the GOCART-2G model implemented in EAM as an EAMxx AtmosphereProcess.
     * This process is encompasses:
     *   1. Surface emissions of constituent species from file I/O
     *   2. Elevated emissions of constituent speices from file I/O
     *   3. Seasalt online emissions
     *   4. Aqueous cloud chemistry of SOx
     *   5. Dry deposition
     *   6. Wet deposition
     * 
     * Transport is provided by the EAMxx (HOMME), while radiation chemistry is presumably 
     * handled by RRTMGP or an equivalent process.
     *-------------------------------------------------------------------------------------------*/
    class GOCART : public AtmosphereProcess
    {
        /*
        using KT            = ekat::KokkosTypes<DefaultDevice>;
        using view_1d       = typename KT::template view_1d<Real>;
        using view_2d       = typename KT::template view_2d<Real>;
        using const_view_1d = typename KT::template view_1d<const Real>;
        using const_view_2d = typename KT::template view_2d<const Real>;
        */

        using KT            = ekat::KokkosTypes<DefaultDevice>;
        template<typename S>
        using SPack = ekat::Pack<S,SCREAM_SMALL_PACK_SIZE>;

        template <typename S>
        using view_1d       = typename KT::template view_1d<S>;
        template <typename S>
        using const_view_1d = typename KT::template view_1d<const S>;
        template <typename S>
        using view_2d       = typename KT::template view_2d<S>;
        template <typename S>
        using const_view_2d = typename KT::template view_2d<const S>;
        template <typename S>
        using view_3d       = typename KT::template view_3d<S>;
        template <typename S>
        using const_view_3d = typename KT::template view_3d<const S>;
        
        // For some reason the C++ council doesn't allow 'using' statements at class scope
        // using gocart::SurfaceEmission;
        typedef gocart::SurfaceEmission SurfaceEmission;
      
      // --------------------- Public methods ------------------------
      public:

        GOCART (const ekat::Comm& comm, const ekat::ParameterList& params);

        AtmosphereProcessType type () const override { return AtmosphereProcessType::Physics; }

        std::string name () const override { return "GOCART"; }

        void set_grids (const std::shared_ptr<const GridsManager> grids_manager) override;
      
      
      // -------------------- Protected Methods ------------------------
      protected:

        // AtmosphereProcess overrides
        void initialize_impl (const RunType run_type) override;
        void run_impl        (const double dt)        override;
        void finalize_impl   ()                       override;


      // -------------------- Private Methods ------------------------  
      private: 

        // Add tracers
        void register_tracers(std::shared_ptr<const AbstractGrid>); 

        // Initilize surface emission sources
        void init_surface_emissions();

        // Compute 10m windspeed to the 3.41 power
        void update_10m_cubed_windspeed(view_1d<Real>&                    u10_cubed, 
                                        const_view_3d<SPack<Real>> const& horiz_wind);
 

      //-------- Member variables for GOCART process ----------

      protected:
        // Grid pointer
        std::shared_ptr<const AbstractGrid> m_grid;

        // Surface emissions
        std::vector<SurfaceEmission> m_surface_emissions;

        // Seasalt emissions fluxes (nsec,ncol) <-- reversed indexing to match fortran storage order
        Kokkos::View<Real**> m_sslt_flux;

        // Map between each SurfaceEmission and species
        std::vector<std::tuple<SurfaceEmission&, Field&>> m_surf_emiss_tracer_pairs;

        // Arrays with EAM-compatible storage order
        Kokkos::View<Real***> m_eamf90_tracers;         // Contiguous array of EAM tracers
        Kokkos::View<Real***> m_eamf90_cloud_aerosols;  // Cloud-borne aerosols
        Kokkos::View<Real**>  m_eamf90_mbar;            // Molecular weight of atmosphere (amu)
        Kokkos::View<Real**>  m_eamf90_T_mid;           // Temperature (K)
        Kokkos::View<Real**>  m_eamf90_p_mid;           // Pressure (Pa)
        Kokkos::View<Real**>  m_eamf90_cldfrc;          // Cloud fraction
        Kokkos::View<Real**>  m_eamf90_cldnum;          // Cloud droplet concentration (#/kg)
        Kokkos::View<Real**>  m_eamf90_xhnm;            // Atmospheric Density in (#/cm^3)
        Kokkos::View<Real**>  m_eamf90_pdel;            // Layer pressure thickness at the midpoint (Pa)
        Kokkos::View<Real**>  m_eamf90_lwc;             // cloud liquid water content (kg/kg)

        std::shared_ptr<DataInterpolation>    m_data_interpolation;
      
      private:
        int ncol_;      // Number of columns on this rank
        int nlev_;      // Number of levels per column
        int gas_pcnst_; // Number of species expected by bridged EAM routines
        
        

      //--------------- utility functions  --------------------
        //void update_tracer_timestate

    };
    void gocart_test_io(gocart::SurfaceEmission& surf_emiss);


} // namespace scream

#endif // EAMXX_GOCART_PROCESS_INTERFACE_HPP
