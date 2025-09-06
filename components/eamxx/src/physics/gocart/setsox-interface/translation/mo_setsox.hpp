#ifndef SCREAM_GOCART_MO_SETSOX_HPP
#define SCREAM_GOCART_MO_SETSOX_HPP


#include "ekat/ekat.hpp"
#include "share/eamxx_types.hpp"

namespace scream::gocart::mo_setsox
{

    /*---------------------------------------------------------------------------------------------
     * At present, constants and hyperparameters used in 'setsox' are set at compile time
     * in the 'setsox_config' namespace. While this is unlikely to change in the near future,
     * it may be desireable o refactor this as a 'struct' or simply as static, non-const 
     * variables if user input is required.
     *-------------------------------------------------------------------------------------------*/
    namespace setsox_config
    {
        // Maximum number of iterations for root-finding algorithm
        constexpr int ITER_MAX = 20;
        
        // Molar gas constant  [(J/K*mol) * (atm/Pa)]? 
        constexpr Real Ra = 8314.0 / 101325.0;
        
        // water acidity
        constexpr Real xkw = 1.0e-14;

        // FIXME: AWFUL CONSTANT
        // [cm3/L / Avogadro constant]
        constexpr Real const0 = 1.0e3 / 6.023e23;
        
        // 330 ppm = 330.0e-6 atm
        constexpr Real co2g = 330.0e-6;

        // Liquid Water Content (LWC) threshold for aqueous chemistry problem
        constexpr Real SMALL_VALUE_LWC = 1.0e-08;

        // Small/minimum value for mixing ratio concentrations
        constexpr Real SMALL_VALUE_CONC = 1.0e-08;

        // ??? Reference pressure ??? [Pa]
        constexpr Real p0 = 101300.0;

        // Reaction rate coefficients??
        constexpr Real kh0 = 9.0e3;   // HO2(g)          -> Ho2(a)
        constexpr Real kh1 = 2.05e-5; // HO2(a)          -> H+ + O2-
        constexpr Real kh2 = 8.6e5;   // HO2(a) + ho2(a) -> h2o2(a) + o2
        constexpr Real kh3 = 1.0e8;   // HO2(a) + o2-    -> h2o2(a) + o2

        // FIX ME: It't not clear what this value should typically be
        constexpr bool cloud_borne = true;

        // FIX ME: It't not clear what this value should typically be
        constexpr bool modal_aerosols = false;

        // Tolerance charge balance problem
        constexpr Real TOLERANCE = 0.005;

        // Number of gas species
        // This is a Fortran parameter in chem_mods however, there are multiple competing modules/parameters 
        // with this name (presumably for running mutually exclusive chemistry models)
        constexpr int gas_pcnst = 1; // FIX ME: Not actually 1


        // FIX ME: Not sure how to choose this correctly
        /*
        #if ( defined MODAL_AERO_7MODE )
            integer, parameter :: ntot_amode = 7
        #elif ( defined MODAL_AERO_5MODE)
            integer, parameter :: ntot_amode = 5
        #elif ( defined MODAL_AERO_9MODE )
            integer, parameter :: ntot_amode = 9
        #elif (( defined MODAL_AERO_4MODE ) || ( defined MODAL_AERO_4MODE_MOM ))
            integer, parameter :: ntot_amode = 4
        #elif ( defined MODAL_AERO_3MODE )
            integer, parameter :: ntot_amode = 3
        #endif
        */
        constexpr int ntot_amode = 7;

    } 





    //=============================================================================================
    //                            Function Prototypes
    //=============================================================================================
    



    // ---- Set SOx for an individual atmospheric cell ----
    void setsox_for_cell(const Real dtime, 
                         const Real press,
                         const Real pdel, 
                         const Real tfld, 
                         const Real mbar,
                         const Real lwc, 
                         const Real cldfrc, 
                         const Real cldnum,
                         const Real xhnm,
                         Real& hno3, 
                         Real& nh3, 
                         Real& o3, 
                         Real& ho2, 
                         Real& h2so4, // if (cloud_borne)
                         // out (it appears that no other species are modified)
                         Real& so2, 
                         Real& so4,
                         Real& h2o2);

    // ---- Calculate pH of cloud droplet aqueous solution ----
    void calc_ph_values(const Real tfld, 
                        const Real patm, 
                        const Real xlwc,
                        const Real t_factor, 
                        const Real xhno3, 
                        const Real xnh3, 
                        const Real xnh4,
                        const Real xso2, 
                        const Real xso4,
                        const Real xhnm, 
                        const Real fact_so4,
                        // out
                        bool &converged, 
                        Real &xph); 

    // ---- Calculate net charge of aqueous solution ----
    Real compute_aqueous_charge(Real const& fact1_hno3, // HNO3 factor 1
                                Real const& fact2_hno3, //             2
                                Real const& fact3_hno3, //             3
                                Real const& fact1_so2,  // SO2 factor 1
                                Real const& fact2_so2,  //            2
                                Real const& fact3_so2,  //            3
                                Real const& fact4_so2,  //            4
                                Real const& fact1_nh3,  // NH3 factor 1
                                Real const& fact2_nh3,  //            2
                                Real const& fact3_nh3,  //            3
                                Real const& fact_so4,   // SO4 factor
                                Real const& Eh2o,       // H2O effect
                                Real const& Eco2,       // CO2 effect
                                Real const& Eso4,       // SO4 effect
                                Real  pH);        // pH of solution

    // ---- Bisection root finder ----
    template<typename FUNC>
    void find_root_bisection_method(FUNC  func, 
                                    Real  x_low, 
                                    Real  x_high, 
                                    Real& x_root, 
                                    bool& converged);









    //=============================================================================================
    //                            Implementations
    //=============================================================================================




    /*---------------------------------------------------------------------------------------------
     * void setsox(...)
     *
     * This function is intended as an EAMxx compatible refactoring of the subroutine 'SETSOX'
     * located in 'mo_setsox.F90'
     * 
     * An attempt was made to keep variable names and operations faithful to the source with 
     * the exception of a few notable optimizations:
     *   - Aqueous chemistry is performed cell-by-cell to support parallelism.
     *   - Unnecessary dynamic array allocations have been removed.
     *   - The 'cldaero' type is no longer used. It's required member variables have been
     *     refactored as local variables. Unused member variables have been removed.
     *   - The solution of the charge balance problem has been modularized to allow future 
     *     replacememt of the bisection method with a more efficient/accurate root-finding 
     *     algorithm.
     *-------------------------------------------------------------------------------------------*/
    void setsox(const Real                  dtime,   // time step (sec)
                const Kokkos::View<Real**>& press,   // midpoint pressure (Pa)
                const Kokkos::View<Real**>& pdel,    // pressure thickness levels (Pa)
                const Kokkos::View<Real**>& tfld,    // temperature
                const Kokkos::View<Real**>& mbar,    // mean wet atmospheric mass (amu)
                const Kokkos::View<Real**>& lwc,     // cloud liquid water content (kg/kg)  
                const Kokkos::View<Real**>& cldfrc,  // cloud fraction
                const Kokkos::View<Real**>& cldnum,  // droplet number concentration
                const Kokkos::View<Real**>& xhnm,    // total atms density ( /cm^3)
                // Tracked species
                Kokkos::View<Real**>&       hno3, 
                Kokkos::View<Real**>&       nh3, 
                Kokkos::View<Real**>&       o3, 
                Kokkos::View<Real**>&       ho2, 
                Kokkos::View<Real**>&       h2so4, // if (cloud_borne)
                Kokkos::View<Real**>&       so2, 
                Kokkos::View<Real**>&       so4,
                Kokkos::View<Real**>&       h2o2)
    {
        // FIX ME: This should be made into a Kokkos parallel_for once
        //         it has been tested
        const int ncol = press.extent(0);
        const int nlev = press.extent(1);
        
        for (int icol=0; icol<ncol; ++icol)
        {
            for (int ilev=0; ilev<nlev; ++ilev)
            {
                setsox_for_cell(dtime,  
                                press(icol,ilev),  pdel(icol,ilev), tfld(icol,ilev), 
                                mbar(icol,ilev),   lwc(icol,ilev),  cldfrc(icol,ilev), 
                                cldnum(icol,ilev), xhnm(icol,ilev),
                                hno3(icol,ilev),   nh3(icol,ilev),  o3(icol,ilev), 
                                ho2(icol,ilev),    h2so4(icol,ilev), 
                                // out
                                so2(icol,ilev), so4(icol,ilev), h2o2(icol,ilev));
            }
        } 
    }




    /*---------------------------------------------------------------------------------------------
     *-----------------------------------------------------------------------
     *      ... Compute heterogeneous reactions of SOX
     *
     *   (0) using initial PH to calculate PH
     *       (a) HENRYs law constants
     *       (b) PARTIONING
     *       (c) PH values
     *
     *   (1) using new PH to repeat
     *       (a) HENRYs law constants
     *       (b) PARTIONING
     *       (c) REACTION rates
     *       (d) PREDICTION
     *-----------------------------------------------------------------------
     *-----------------------------------------------------------------------
     *  ... Local variables
     *
     *  FORTRAN refactoring: the units are a little messy here
     *  my understanding (may not be right) is that, the PH value xph, shown
     *  in [H+] concentration, is (mol H+)/(L water), which can be
     *  transfered to kg/L or kg/kg the variables xso2, xso4, xo3 etc have
     *  units of [mol/mol] (maybe corresponding to kg/kg above?) the
     *  variable xhnm has unit of [#/cm3]. Some units may changes to
     *  different formats across modules
     *  Shuaiqi Tang  4/18/2023
     *-----------------------------------------------------------------------
     *-------------------------------------------------------------------------------------------*/
    KOKKOS_INLINE_FUNCTION
    void setsox_for_cell(const Real dtime, 
                         const Real press,
                         const Real pdel, 
                         const Real tfld, 
                         const Real mbar,
                         const Real lwc, 
                         const Real cldfrc, 
                         const Real cldnum,
                         const Real xhnm,
                         Real& hno3, 
                         Real& nh3, 
                         Real& o3, 
                         Real& ho2, 
                         Real& h2so4, // if (cloud_borne)
                         // out (it appears that no other species are modified)
                         Real& so2, 
                         Real& so4,
                         Real& h2o2)
    {
        using std::pow; // alternative to haero::pow
        using std::max;

        //xnhm := total atms density [#/cm3]
        constexpr Real one = 1.0;
        constexpr Real zero = 0.0;
        constexpr Real t298K = 298.0;

        constexpr Real ph0 = 5.0; // Initial PH values
        // initial PH value, in H+ concentration
        Real xph0 = pow(10, -ph0);
        // cfact := total atms density [kg/L]
        // FIXME: BAD CONSTANTS
        //           cm-3 * m-3    * kg/m3            * kg/L;
        Real cfact = xhnm * 1.0e6 * 1.38e-23 / 287.0 * 1.0e-3;


        /*---------------------------------------------
         * There is no need for the Couldconc struct.
         * The original Fortran version has several
         * uninitialized arrays in addition to trivial
         * computable quantities. In this refactoring
         * of the mam4xx implementation, necessary 
         * data are local variables in this function
         * 
         * Daniel Fust 07/31/2025
         *--------------------------------------------*/
        Real so4c = 0.0;
        Real xlwc = 0.0;

        // This step is the only non-allocation action performed by 'sox_cldaero_create_obj'
        // xlwc is in-cloud LWC with the unit of [kg/L]
        if (cldfrc > 0.0) {
            // cloud water L(water)/L(air)
            xlwc = lwc * cfact;
            // liquid water in the cloudy fraction of cell
            xlwc = xlwc / cldfrc;
        } else {
            xlwc = 0.0;
        }

        /*--------------------------------------------------
         * Note: The EAM version of SETSOX does not appear
         *       to distinguish between variants of SO4
         *       such as SO4_1a, SO4_2a, SO4_3a which
         *       are found in the mam4xx refactor, found
         *       as
         * so4c = qcw[id_so4_1a] + qcw[id_so4_2a] + qcw[id_so4_3a];
         * so4c = so4_1a + so4_2a + so4_4a;
         *-------------------------------------------------*/
        
        // Set as 2.0 in cldaero_allocate() in cldaero_mod.F90
        Real fact_so4 = 2.0;

        // species molar mixing ratios(?) [mol/mol]
        Real xso4 = zero;
        // initial PH value
        Real xph    = xph0;
        // I/O Data is copied then updated later to be more consistent with the old EAM version
        Real xhno3  = hno3;  //qin[setsox_config_.id_hno3];
        Real xso2   = so2;   //qin[setsox_config_.id_so2];
        Real xh2o2  = h2o2;  //qin[setsox_config_.id_h2o2];
        Real xo3    = o3;    //qin[setsox_config_.id_o3];
        Real xh2so4 = h2so4; //qin[setsox_config_.id_h2so4];
        Real xnh3   = nh3;   //qin[setsox_config_.id_nh3];
        Real xho2   = ho2;   //qin[setsox_config_.id_ho2];

        // Local species
        Real xnh4 = 0.0;
        Real xno3 = 0.0;

        // there doesn't appear to be any reason for doing this
        //Real xso4c = cldconc.so4c;

        //Real xlwc = cldconc.xlwc;

        Real t_factor = (one / tfld) - (one / t298K);
        // calculate press in atm
        Real patm = press / setsox_config::p0;

        // Compute pH if liquid water content is high enough
        if (xlwc >= setsox_config::SMALL_VALUE_LWC) 
        {
            if (setsox_config::cloud_borne && (cldfrc > zero)) 
            {xso4 = so4 / cldfrc;}

            bool converged = false;
            calc_ph_values(tfld, patm, xlwc, t_factor, xhno3,  xnh3, 
                           xnh4, xso2, xso4, xhnm, fact_so4,
                           // out
                           converged, xph); 
    
        }// end if (xlwc >= small_value_xlwc)

        //==============================================================
        //          ... Now use the actual pH
        //==============================================================
        Real xk,xe,x2;
        Real xam = press / (tfld * 1.38e-23); // air density /M3

        //  HNO3
        xk = 2.1e5 * exp( 8700.0*t_factor);
        xe = 15.4;
        Real hehno3 = xk*(1.0 + xe/xph);

        // H2O2
        xk = 7.4e4 * exp(6621.0*t_factor);
        xe = 2.2e-12 * exp(-3730.0*t_factor);
        Real heh2o2 = xk*(1.0 + xe/xph);

        // SO2
        xk = 1.23   * exp(3120.0 * t_factor);
        xe = 1.7e-2 * exp(2090.0 * t_factor);
        x2 = 6.0e-8 * exp(1120.0 * t_factor);

        Real heso2 = xk*(1.0 + (xe/xph)*(1.0 + x2/xph));

        // NH3
        constexpr Real xkw = setsox_config::xkw; 
        xk = 58.0   * exp( 4085.0 * t_factor);
        xe = 1.7e-5 * exp(-4325.0 * t_factor);
        Real henh3 = xk*(1.0 + xe*xph/xkw);

        // O3
        xk   = 1.15e-2 * exp(2560.0*t_factor);
        Real heo3 = xk;


        //------------------------------------------------------------------------
        //       ... for Ho2(g) -> H2o2(a) formation 
        //           schwartz JGR, 1984, 11589
        //------------------------------------------------------------------------
        constexpr Real kh0 = setsox_config::kh0;
        constexpr Real kh1 = setsox_config::kh1;
        constexpr Real kh2 = setsox_config::kh2;
        constexpr Real kh3 = setsox_config::kh3;
        constexpr Real Ra  = setsox_config::Ra;
        Real r2h2o2;
        Real xl = xlwc; // Relic of old version in which 'xlwc' was an array and 'xl = xlwc(icol,ilev)'
        Real kh4    = (kh2 + kh3*kh1/xph) / pow(1.0 + kh1/xph, 2);
        Real ho2s   = kh0*xho2*patm*(1.0 + kh1/xph);  // ho2s = ho2(a)+o2-
        Real r1h2o2 = kh4*ho2s*ho2s;                  // prod(h2o2) in mole/L(w)/s

        if ( setsox_config::cloud_borne )
        {
             r2h2o2 = r1h2o2*xl                  // mole/L(w)/s   * L(w)/fm3(a) = mole/fm3(a)/s
                  / setsox_config::const0*1.0e6  // correct a bug here ????
                  / xam;
        }          
        else
        {
            r2h2o2 = r1h2o2*xl           //  mole/L(w)/s   * L(w)/fm3(a) = mole/fm3(a)/s
                * setsox_config::const0  //  mole/fm3(a)/s * 1.e-3       = mole/cm3(a)/s
                / xam;                   //  /cm3(a)/s    / air-den      = mix-ratio/s
        } // if ( setsox_config::cloud_borne )


        if (!setsox_config::modal_aerosols) {xh2o2 += r2h2o2*dtime;} // updated h2o2 by het production
        
        //------------ Partitioning --------------
        Real px;

        // HNO3
        px    = hehno3 * Ra * tfld * xl;
        Real hno3g = (xhno3+xno3)/(1.0 + px);

        // H2O2
        px    = heh2o2 * Ra * tfld * xl;
        Real h2o2g = xh2o2/(1.0 + px);

        // SO2
        px   = heso2 * Ra * tfld * xl;
        Real so2g = xso2/(1.0 + px);

        // O3
        px  = heo3 * Ra * tfld * xl;
        Real o3g = xo3/(1.0 + px);

        // NH3
        /*
        OLD CODE:

        px = henh3(i,k) * Ra * tz * xl
        if (id_nh3>0) then
            nh3g(i,k) = (xnh3(i,k)+xnh4(i,k))/(1._r8+ px)
        else
            nh3g(i,k) = 0._r8
        endif

        It appears that NH3 may not always have been present and the author
        was using the "ID" or index within the tracer array to check for
        its presence. For this refactor, NH3 will be assumed present in the
        model

        Daniel Fust (08/7/2025)
        */
        px        = henh3 * Ra * tfld * xl;
        Real nh3g = (xnh3 + xnh4)/(1.0 + px);
        
        /*-----------------------------------------------
         *       ... Aqueous phase reaction rates
         *           SO2 + H2O2 -> SO4
         *           SO2 + O3   -> SO4
         *---------------------------------------------*/

        /*------------------------------------------------------------------------
         *       ... S(IV) (HSO3) + H2O2
         *----------------------------------------------------------------------*/
        Real rah2o2 = 8.e4 * exp( -3650.0*t_factor ) / (0.1 + xph);

        /*------------------------------------------------------------------------
         *        ... S(IV)+ O3
         *----------------------------------------------------------------------*/
        Real rao3   = 4.39e11 * exp(-4131.0 / tfld) + 2.56e3 * exp(-996.0 / tfld) / xph;


        /*-----------------------------------------------------------------
         *       ... Prediction after aqueous phase
         *       so4
         *       When Cloud is present 
         *   
         *       S(IV) + H2O2 = S(VI)
         *       S(IV) + O3   = S(VI)
         *
         *       reference:
         *           (1) Seinfeld
         *           (2) Benkovitz
         *---------------------------------------------------------------*/
        
        /*............................
         *       S(IV) + H2O2 = S(VI)
         *..........................*/

        if (xl >= setsox_config::SMALL_VALUE_LWC) // when cloud is present
        {
            Real pso4, patm_x, ccc, xdelso4hp;

            if (setsox_config::cloud_borne) {patm_x = patm;}
            else {patm_x = 1.0;}

            if (setsox_config::modal_aerosols)
            {
                pso4 = rah2o2 * 7.4e4 * exp(6621.0 * t_factor) * h2o2g * patm_x
                      * 1.23 * exp(3120.0 * t_factor) * so2g * patm_x;
            }
            else
            {
                pso4 = rah2o2 * heh2o2 * h2o2g * patm_x 
                       * heso2 * so2g  * patm_x;        // [M/s]
            }

            // [M/s] = [mole/L(w)/s]
            pso4 *= xl                      // ! [mole/L(a)/s]
                    / setsox_config::const0 // [/L(a)/s]
                    / xhnm;

            ccc = pso4*dtime;
            ccc = max(ccc, 1.0e-30);        

            Real xso4_init = xso4;

            if (xh2o2 > xso2)
            {
                if (ccc > xso2)
                { 
                    xso4 += xso2;
                    if (setsox_config::cloud_borne)
                    {
                        xh2o2 -= xso2;
                        xso2  =  1.0e-20;
                    }   
                    else //        ???? bug ????
                    {
                        xso2  =  1.0e-20;
                        xh2o2 -= xso2;
                    }  // if (setsox_config::cloud_borne)
                }   
                else
                {
                    xso4  += ccc;
                    xh2o2 -= ccc;
                    xso2  -= ccc;
                } // if (ccc > xso2)
            }
            else
            {
                if (ccc > xh2o2) 
                {
                   xso4  += xh2o2;
                   xso2  -= xh2o2;
                   xh2o2 =  1.0e-20;
                }   
                else
                {
                   xso4  += ccc;
                   xh2o2 -= ccc;
                   xso2  -= ccc;
                } // if (ccc > xh2o2)
            } // if (xh2o2 > xso2)

            
            if (setsox_config::modal_aerosols) {xdelso4hp = xso4 - xso4_init;}

            /*...........................
             *       S(IV) + O3 = S(VI)
             *.........................*/

            pso4 = rao3 * heo3 * o3g * patm_x * heso2 * so2g * patm_x;  // [M/s]

            // [M/s] =  [mole/L(w)/s]
            pso4 *= xl                    // [mole/L(a)/s]
                  / setsox_config::const0 // [/L(a)/s]
                  / xhnm;                 // [mixing ratio/s]
             
            ccc = max(pso4*dtime, 1.0e-30);

            xso4_init=xso4;

            if (ccc > xso2)
            {
                xso4 += xso2;
                xso2 =  1.0e-20;
            }    
            else
            {
                xso4 += ccc;
                xso2 -= ccc;
            } // if (ccc > xso2)
        } // if (xl >= setsox_config::SMALL_VALUE_LWC) // when cloud is present



        // FIX ME:  sox_cldaero_update from chemistry/modal_aero/sox_cldaero_mod.F90 is much more complicated
        //          but it's not clear which should be used.
        // equivalent to: call sox_cldaero_update (outside loop) from chemistry/bulk_aero/sox_cldaero_mod.F90
        so2  = max(xso2, setsox_config::SMALL_VALUE_CONC);
        h2o2 = max(xh2o2,setsox_config::SMALL_VALUE_CONC);
        so4  = max(xso4, setsox_config::SMALL_VALUE_CONC);


        // The following is only required if some replacement to the
        // 'outfld' subroutine call is needed.
        /*
        Real xphlwc
        if ((cldfrc >= 1.0e-5) && (lwc >= 1.0e-8))
        {xphlwc = -1.0 * log10(xph) * lwc;}
        else {xphlwc = 0.0;}
        */

        //call outfld( 'XPH_LWC', xphlwc(:ncol,:), ncol , lchnk )
        //call sox_cldaero_destroy_obj(cldconc)


        // FIX ME
        // Not sure if this matters in EAMxx. Check what 'outfld' routine does...
        /* 
        xphlwc(:,:) = 0._r8
        do k = 1, pver
            do i = 1, ncol
                if (cldfrc(i,k)>=1.e-5_r8 .and. lwc(i,k)>=1.e-8_r8) then
                    xphlwc(i,k) = -1._r8*log10(xph(i,k)) * lwc(i,k)
                endif
            end do
        end do
        call outfld( 'XPH_LWC', xphlwc(:ncol,:), ncol , lchnk )
        */
    }




    KOKKOS_INLINE_FUNCTION
    void calc_ph_values(const Real tfld, 
                        const Real patm, 
                        const Real xlwc,
                        const Real t_factor, 
                        const Real xhno3, 
                        const Real xnh3, 
                        const Real xnh4,
                        const Real xso2, 
                        const Real xso4,
                        const Real xhnm, 
                        const Real fact_so4,
                        // out
                        bool &converged, 
                        Real &xph) 
    {
        using std::exp;

        // Temperature dependent Henry constants (reused per species)
        Real xe, xk, x2;
        constexpr Real Ra   = setsox_config::Ra;
        constexpr Real xkw  = setsox_config::xkw;
        constexpr Real co2g = setsox_config::co2g;
        constexpr Real const0 = setsox_config::const0;

        // Factors for HNO3
        xk         = 2.1e5 * exp(8700.0*t_factor);
        xe         = 15.4;
        Real fact1_hno3 = xk*xe*patm*xhno3;
        Real fact2_hno3 = xk*Ra*tfld*xlwc;
        Real fact3_hno3 = xe;

        // Temperature dependent Henry's law factors for SO2
        xk = 1.23e-3 * exp(3120.0*t_factor);
        xe = 1.7e-2  * exp(2090.0*t_factor);
        x2 = 6.0e-8  * exp(1120.0*t_factor);
        Real fact1_so2 = xk*xe*patm*xso2;
        Real fact2_so2 = xk*Ra*tfld*xlwc;
        Real fact3_so2 = xe;
        Real fact4_so2 = x2;

        // Temperature dependent Henry's law factors for NH3
        xk = 58.0   * exp(4085.0*t_factor);
        xe = 1.7e-5 * exp(-4325.0*t_factor);
        Real fact1_nh3 = (xk*xe*patm/xkw)*(xnh3 + xnh4);
        Real fact2_nh3 = xk*Ra*tfld*xlwc;
        Real fact3_nh3 = xe/xkw;

        // H2O effects
        Real Eh2o = xkw;

        // CO2 effects
        //const Real co2g = 330.0e-6;            //330 ppm = 330.e-6 atm
        xk = 3.1e-2 * exp( 2423.0*t_factor);
        xe = 4.3e-7 * exp(-913.0 *t_factor);
        Real Eco2 = xk*xe*co2g  * patm;

        // SO4 effects
        
        Real Eso4 = xso4*xhnm   //  /cm3(a)
                    *const0/xlwc;

        // Lambda is used to wrap net charge computation to fit with more generic
        // root-finding algorithms (Not sure if it's safe to capture by refernce here or not)
        auto compute_ynetpos = [&fact1_hno3, &fact2_hno3, &fact3_hno3, &fact1_so2, 
                                &fact2_so2,  &fact3_so2,  &fact4_so2,  &fact1_nh3,  
                                &fact2_nh3,  &fact3_nh3,  &fact_so4,   
                                &Eh2o,       &Eco2,       &Eso4]
            (Real pH) -> Real
        {
            return compute_aqueous_charge(fact1_hno3,fact2_hno3,fact3_hno3,fact1_so2, 
                                          fact2_so2, fact3_so2, fact4_so2, fact1_nh3,  
                                          fact2_nh3, fact3_nh3, fact_so4,   
                                          Eh2o,      Eco2,      Eso4,      pH); 
        }; // compute_ynetpos

        // Solve electro-neutrality problem with root finding algorithm
        find_root_bisection_method(compute_ynetpos, 2.0, 7.0, xph, converged);
    }

        
        




    /*
    ---------------------------------------------------------------------------
    calculate PH value and H+ concentration

    21-mar-2011 changes by rce
    now uses bisection method to solve the electro-neutrality equation
    3-mode aerosols (where so4 is assumed to be nh4hso4)
            old code set xnh4c = so4c
            new code sets xnh4c = 0, then uses a -1 charge (instead of -2)
        for so4 when solving the electro-neutrality equation
    ---------------------------------------------------------------------------

    ----------------------------------------
    effect of chemical species
    ----------------------------------------
        
    */






    /*---------------------------------------------------------------------------------------------
     * Real compute_aqueous_charge(...)
     *
     * FIXME: CHEMISTRY GOES HERE
     * 
     * Returns:
     *   Net positive charge of aqueous chemistry problem
     *-------------------------------------------------------------------------------------------*/
    Real compute_aqueous_charge(Real const& fact1_hno3, // HNO3 factor 1
                                Real const& fact2_hno3, //             2
                                Real const& fact3_hno3, //             3
                                Real const& fact1_so2,  // SO2 factor 1
                                Real const& fact2_so2,  //            2
                                Real const& fact3_so2,  //            3
                                Real const& fact4_so2,  //            4
                                Real const& fact1_nh3,  // NH3 factor 1
                                Real const& fact2_nh3,  //            2
                                Real const& fact3_nh3,  //            3
                                Real const& fact_so4,   // SO4 factor
                                Real const& Eh2o,       // H2O effect
                                Real const& Eco2,       // CO2 effect
                                Real const& Eso4,       // SO4 effect
                                Real        pH)         // pH of solution
    {
        using std::pow;
        // compute [H+] concentration from pH
        Real hplus = pow(10.0,-pH);

        // HNO3
        Real Ehno3 = fact1_hno3/(1.0 + fact2_hno3*(1.0 + fact3_hno3/hplus));
        // SO2
        Real Eso2 = fact1_so2/(1.0 + fact2_so2*(1.0 + (fact3_so2/hplus) *(1.0 + fact4_so2/hplus)));
        // NH3 
        Real Enh3 = fact1_nh3/(1.0 + fact2_nh3*(1.0 + fact3_nh3*hplus)); 

        // Calculate charge balance
        Real nh4  = Enh3 * hplus;
        Real hso3 = Eso2 / hplus;
        Real so3  = hso3 * 2.0*fact4_so2/hplus;
        Real hco3 = Eco2 / hplus;
        Real oh   = Eh2o / hplus;
        Real no3  = Ehno3 / hplus;
        Real so4  = fact_so4*Eso4;
        Real pos  = hplus + nh4;
        Real neg  = oh + hco3 + no3 + hso3 + so3 + so4;

        return pos - neg;    
    } // compute_aqueous_charge





    /*---------------------------------------------------------------------------------------------
     *
     * The template parameter FUNC is used for interoperability with capturing lambda functions.
     * Functions (including lambdas) should have the signature:
     *     Real (*func)(Real)
     * 
     *-------------------------------------------------------------------------------------------*/
    template<typename FUNC>
    void find_root_bisection_method(FUNC func, Real x_low, Real x_high, 
                                    Real& x_root, bool& converged)
    {
        using std::abs;
        Real x, val, val_l, val_h;

        // Initialize bracketing values
        Real xl = x_low;
        Real xh = x_high;

        converged = false;
        int iter = 0;
        while (!converged)
        {
            iter++;

            // If iteration cap is exceeded, return midpoint
            if (iter > setsox_config::ITER_MAX)
            {
                x_root    = 0.5*(xl + xh);
                converged = false;
                return;
            }

            // First iteration: guess lower bound
            if (iter==1) {x = x_low;}
            // Second iteration: guess upper bound
            else if (iter==2) {x = x_high;}
            // Take mean value of guesses bracketing the root
            else {x = 0.5*(xl + xh);}

            // Get function value at current guess
            val = func(x);
 
            // Check convergence and update bracketing values
            if (abs(val) < setsox_config::TOLERANCE)
            {
                x_root    = x;
                converged = true;
                return;
            }
            // function value at 'x' is positive -> move upper bracket value
            else if (val > 0.0) 
            {
                if (iter==2) // x==x_low bounds do not bracket the root
                {
                    converged = false;
                    x_root    = x_low;
                    std::cerr << "WARNING: in find_root_bisection_method...\n"
                              << "Provided bounds do not bracket the root\n";
                    return;
                }
                xh = x;
            }
            // function value at 'x' is negative -> move lower bracket value
            else /*if (val <= 0.0)*/ 
            {
                if (iter==1) // x==x_high bounds do not bracket the root
                {
                    converged = false;
                    x_root    = x_high;
                    std::cerr << "WARNING: in find_root_bisection_method...\n"
                              << "Provided bounds do not bracket the root\n"; 
                    return;
                }
                xl = x;
            } // if (abs(val) < setsox_config::TOLERANCE)
        } // while (!converged)
    }





    #if 0 // wrong routine??

    /*---------------------------------------------------------------------------------------------
     * void sox_cldaero_update(...)
     *-------------------------------------------------------------------------------------------*/
    void sox_cldaero_update(//int ncol, 
                            //int lchnk, 
                            //int loffset, 
                            Real dtime, 
                            Real mbar, 
                            Real pdel, 
                            Real press, 
                            Real tfld, 
                            Real cldnum, 
                            Real cldfrc, 
                            Real cfact, 
                            Real xlwc, 
                            Real delso4_hprxn, 
                            Real xh2so4, 
                            Real xso4, 
                            Real xso4_init, 
                            Real nh3g, 
                            Real hno3g, 
                            Real xnh3, 
                            Real xhno3, 
                            Real xnh4c,  
                            Real xno3c, 
                            Real xmsa, 
                            Real xso2, 
                            Real xh2o2, 
                            //qcw, 
                            //qin
                            )
    {
        const int NSPECIES   = setsox_config::gas_pcnst;
        const int NTOT_AMODE = setsox_config::ntot_amode;
        
        // make sure dqdt is zero initially, for budgets
        Real dqdt_aqhprxn = 0.0
        Real dqdt_aqo3rxn = 0.0
        Real dqdt_aqso4[NSPECIES];
        Real dqdt_aqh2so4[NSPECIES];
        for (int i=0; i<NSPECIES; ++i)
        {
            dqdt_aqso4[i]   = 0;
            dqdt_aqh2so4[i] = 0;
        }

        Real qnum_c[NTOT_AMODE];
        Real faqgain_msa[NTOT_AMODE];
        Real faqgain_so4[NTOT_AMODE];

        Real xl = xlwc; // This is a relic from the Fortran version in which xlwc was an array
        if (xl <= setsox_config::SMALL_VALUE_LWC) // If cloud is present
        {
            Real delso4_o3rxn = xso4 - xso4_init;

            /*
            ORIGINAL CODE:
            
            if (id_nh3>0) then
               delnh3 = nh3g(i,k) - xnh3(i,k)
               delnh4 = - delnh3
            endif

            Authors seem to use id_nh3 (i.e. the index in the tracer array)
            to determine the presence of NH3. This sucks and doesn't really
            translate to EAMxx. It is assumed here that NH3 are always
            present.
            */
            Real delnh3 = nh3g - xnh3;
            Real delnh4 = -delnh3;

            /*-------------------------------------------------------------------------
             * compute factors for partitioning aerosol mass gains among modes
             * the factors are proportional to the activated particle MR for each
             * mode, which is the MR of cloud drops "associated with" the mode
             * thus we are assuming the cloud drop size is independent of the
             * associated aerosol mode properties (i.e., drops associated with
             * Aitken and coarse sea-salt particles are same size)
             *
             * qnum_c(n) = activated particle number MR for mode n (these are just
             * used for partitioning among modes, so don't need to divide by cldfrc)
             *-----------------------------------------------------------------------*/

            do n = 1, NTOT_AMODE
               qnum_c(n) = 0.0_r8
               l = numptrcw_amode(n) - loffset
               if (l > 0) qnum_c(n) = max( 0.0_r8, qcw(i,k,l) )
            end do

        }
        



    }
    #endif









} // namespace scream::gocart

#endif //SCREAM_GOCART_SETSOX_HPP