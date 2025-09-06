/*-------------------------------------------------------------------------------------------------
 * sslt_sections.hpp
 *
 * C++ function declarations for external 'sslt_sections' Fortran module procedues for seasalt 
 * aerosol emissions
 *-------------------------------------------------------------------------------------------------*/
#ifndef GOCART_SSLT_SECTIONS_HPP
#define GOCART_SSLT_SECTIONS_HPP


namespace scream::gocart::sslt_sections
{
    /*=================  External Procedures =================*/

    // Module variable getters 
    extern "C" int  sslt_sections_nsections(); 

    // Initialization of module variables
    extern "C" void sslt_sections_init();

    // (???) fluxes due to seasalt
    extern "C" void sslt_fluxes_c2f(const int     ncol,
                                    const int     nsec,
                                    double*       fi, 
                                    const double* sst, 
                                    const double* u10cubed);

    /*---------------------------------------------------------------------------------------------
     * The 'sslt_sections' namespace can now be used to assign names that shadow the more natural
     * EAM equivalents while using cleaner call statements on the C++ side
     *-------------------------------------------------------------------------------------------*/
    int nsections() {return sslt_sections_nsections();}



    /*---------------------------------------------------------------------------------------------
     * void sslt_fluxes(...)
     *
     * This procedure provides an additional layer of safety between EAMxx types and the 
     * pointer-based operations of the C++ Fortran bridge of seasalt emissions.
     * 
     * Kokkos Views are assumed to be in row-major storage order, whereas Fortran arrays are
     * always in column-major order.
     * 
     * Arguments:
     *   flux............(ncols,nsections) Flux due to seasalt in each "section"
     *   sea_surf_temp...(ncols) Sea surface temperature in (K)
     *   u10_cubed.......(ncols) 10m windspeed with (3.41) exponent according to empircal fit by
     *                           Gong et al., 1997
     *-------------------------------------------------------------------------------------------*/
    void sslt_fluxes(Kokkos::View<Real**>&            flux,
                     Kokkos::View<const Real*> const& sea_surf_temp,
                     Kokkos::View<const Real*> const& u10_cubed)
    {
        const int ncols = sea_surf_temp.extent(0);
        const int nsecs = sslt_sections_nsections();
        
        // Dimension check for debugging
        #ifdef GOCART_DEBUG
        if (flux.extent(0) != ncols) || (flux.extent(1) != nsecs) || (u10_cubed.extent(0) != ncols)
        {
            std::cerr << "Error in gocart::sslt_fluxes:  Incompatible arrays sizes\n"
                      << "ncols,nsecs: " << ncols << "," << nsecs << "\n"
                      << "flux: " << flux.extent(0) << "x" << flux.extent(1) << "\n"
                      << "sea_surf_temp: " << sea_surf_temp.extent(0) << "\n"
                      << "u10_cubed: " << u10_cubed.extent(0) 
                      << std::endl;
        }
        #endif

        std::cout << "In sslt_fluxes..." << std::endl;
        //sslt_fluxes_c2f(flux.data(), sea_surf_temp.data(), u10_cubed.data(), ncols);
        sslt_fluxes_c2f(ncols, nsecs, flux.data(), sea_surf_temp.data(), u10_cubed.data());
    }



}
#endif // GOCART_SSLT_SECTIONS_HPP