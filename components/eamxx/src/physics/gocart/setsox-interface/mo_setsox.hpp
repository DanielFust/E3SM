#include <iostream>
namespace scream::gocart::mo_setsox
{

    // Initialization subroutine
    extern "C" void sox_inti_c2f(bool use_modal_aerosols, 
                                 bool use_mmf, 
                                 bool use_ecpp);

    // Primary `mo_setsox` aqueous chemistry subroutine
    extern "C" void setsox_c2f(int           ncol,   // num of columns
                               int           nlev,   // num of elevation levels
                               double        dtime,  // time step (sec)
                               const double* press,  // midpoint pressure ( Pa )
                               const double* pdel,   // pressure thickness of levels (Pa)
                               const double* tfld,   // temperature
                               const double* mbar,   // mean wet atmospheric mass ( amu )
                               const double* lwc,    // cloud liquid water content (kg/kg)
                               const double* cldfrc, // cloud fraction
                               const double* cldnum, // droplet number concentration (#/kg)
                               const double* xhnm,   // total atms density ( /cm**3)
                               double*       qcw,    // cloud-borne aerosol (vmr)
                               double*       qin);   // transported species ( vmr )

    /*---------------------------------------------------------------------------------------------
     * void setsox(...)
     *
     * EAMxx side interface for `mo_setsox` `SETSOX` subroutine for cloud aqueous chemistry
     * problem
     *-------------------------------------------------------------------------------------------*/
    void setsox(Real                         dtime,  // time step (sec)
                Kokkos::View<Real**> const&  press,  // midpoint pressure ( Pa )
                Kokkos::View<Real**> const&  pdel,   // pressure thickness of levels (Pa)
                Kokkos::View<Real**> const&  tfld,   // temperature
                Kokkos::View<Real**> const&  mbar,   // mean wet atmospheric mass ( amu )
                Kokkos::View<Real**> const&  lwc,    // cloud liquid water content (kg/kg)
                Kokkos::View<Real**> const&  cldfrc, // cloud fraction
                Kokkos::View<Real**> const&  cldnum, // droplet number concentration (#/kg)
                Kokkos::View<Real**> const&  xhnm,   // total atms density ( /cm**3)
                Kokkos::View<Real***>&       qcw,    // cloud-borne aerosol (vmr)
                Kokkos::View<Real***>&       qin)    // transported species ( vmr )
    {
        std::cout << "In mo_setsox::setsox...\n";
        int ncol = press.extent(1);
        int nlev = press.extent(0);
        setsox_c2f(ncol,nlev,dtime,
                   press.data(), 
                   pdel.data(),
                   tfld.data(), 
                   mbar.data(),
                   lwc.data(), 
                   cldfrc.data(),
                   cldnum.data(), 
                   xhnm.data(), 
                   qcw.data(), 
                   qin.data()); 
    }
}