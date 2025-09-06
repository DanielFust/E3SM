module mo_gocart

    !---------------------------------------------------------------------
    !       ... Dry deposition velocity input data and code for netcdf input
    !---------------------------------------------------------------------
  
  !LKE (10/11/2010): added HCN, CH3CN, HCOOH
  
    use shr_kind_mod, only : r8 => shr_kind_r8, shr_kind_cl
    use chem_mods,    only : gas_pcnst
    use pmgrid,       only : plev, plevp
    use spmd_utils,   only : masterproc, iam
    use ppgrid,       only : pcols, pver, begchunk, endchunk
    use mo_tracname,  only : solsym
    use cam_abortutils,   only : endrun
    !use ioFileMod,    only : getfil
#ifdef SPMD
    use mpishorthand, only : mpicom, mpir8, mpiint, mpilog
#endif
    !use pio
    !use cam_pio_utils,only : cam_pio_openfile
    use cam_logfile,  only : iulog
    !use dyn_grid,     only : get_dyn_grid_parm, get_horiz_grid_d
    !use iop_data_mod, only:  single_column
  
    !use seq_drydep_mod, only : mapping
    use physconst,    only : karman
  
    implicit none
  
    save
  
  
    type gc_wetdep_inputs_t
       real(r8), pointer :: cldt(:,:) => null()  ! cloud fraction
       real(r8), pointer :: qme(:,:) => null()
       real(r8), pointer :: prain(:,:) => null()
       real(r8), pointer :: evapr(:,:) => null()
       real(r8), allocatable :: cldcu(:,:)  ! convective cloud fraction, currently empty
       real(r8), allocatable :: evapc(:,:)  ! Evaporation rate of convective precipitation
       real(r8), allocatable :: cmfdqr(:,:) ! convective production of rain
       real(r8), allocatable :: conicw(:,:) ! convective in-cloud water
       real(r8), allocatable :: totcond(:,:)! total condensate
       real(r8), allocatable :: cldv(:,:)   ! cloudy volume undergoing wet chem and scavenging
       real(r8), allocatable :: cldvcu(:,:) ! Convective precipitation area at the top interface of current layer
       real(r8), allocatable :: cldvst(:,:) ! Stratiform precipitation area at the top interface of current layer
    end type gc_wetdep_inputs_t
  
    real(r8), parameter :: cmftau = 3600._r8
    real(r8), parameter :: rhoh2o = 1000._r8            ! density of water
    real(r8), parameter :: molwta = 28.97_r8            ! molecular weight dry air gm/mole
    real(r8), parameter :: gravit = 9.80665_r8
    real(r8), parameter :: tmelt = 268.0_r8
  
    integer :: cld_idx             = 0
    integer :: qme_idx             = 0
    integer :: prain_idx           = 0
    integer :: nevapr_idx          = 0
  
    integer :: icwmrdp_idx         = 0
    integer :: icwmrsh_idx         = 0
    integer :: rprddp_idx          = 0
    integer :: rprdsh_idx          = 0
    integer :: sh_frac_idx         = 0
    integer :: dp_frac_idx         = 0
    integer :: nevapr_shcu_idx     = 0
    integer :: nevapr_dpcu_idx     = 0
    integer :: ixcldice, ixcldliq
  
    logical :: pergro_mods         = .false.
  
    integer,parameter :: nddvels=19
    character(len=20) :: drydep_list(nddvels)
    data drydep_list /'gcso4_a','gcsea_a1','gcsea_a2','gcsea_a3','gcsea_a4','gcsea_a5',&
                      'gcdst_a1','gcdst_a2','gcdst_a3','gcdst_a4','gcdst_a5','gcbco_a',&
                      'gcbci_a','gcoco_a','gcoci_a','gcnh3_a',&
                      'gcnh4_a','gcnit_a1','gcsoa_a'/
  
    real(r8)              :: dels
    real(r8), allocatable :: days(:)          ! day of year for soilw
    real(r8), allocatable :: dvel(:,:,:,:)    ! depvel array interpolated to model grid
    real(r8), allocatable :: dvel_interp(:,:,:) ! depvel array interpolated to grid and time
    integer :: last, next                     ! day indicies
    integer :: ndays                          !# of days in soilw file
    integer :: map(gas_pcnst)                 ! indices for drydep species
    integer :: nspecies                       ! number of depvel species in input file
  
    integer :: gcso4_ndx
    integer :: pan_ndx, mpan_ndx, no2_ndx, hno3_ndx, o3_ndx, &
               h2o2_ndx, onit_ndx, onitr_ndx, ch4_ndx, ch2o_ndx, &
               ch3ooh_ndx, pooh_ndx, ch3coooh_ndx, c2h5ooh_ndx, eooh_ndx, &
               c3h7ooh_ndx, rooh_ndx, ch3cocho_ndx, co_ndx, ch3coch3_ndx, &
               no_ndx, ho2no2_ndx, glyald_ndx, hyac_ndx, ch3oh_ndx, c2h5oh_ndx, &
               hydrald_ndx, h2_ndx, Pb_ndx, o3s_ndx, o3inert_ndx, macrooh_ndx, &
               xooh_ndx, ch3cho_ndx, isopooh_ndx
    integer :: alkooh_ndx, mekooh_ndx, tolooh_ndx, terpooh_ndx, ch3cooh_ndx
    integer :: soa_ndx, so4_ndx, cb1_ndx, cb2_ndx, oc1_ndx, oc2_ndx, nh3_ndx, nh4no3_ndx, &
               sa1_ndx, sa2_ndx, sa3_ndx, sa4_ndx, nh4_ndx
    integer :: soam_ndx, soai_ndx, soat_ndx, soab_ndx, soax_ndx, &
               sogm_ndx, sogi_ndx, sogt_ndx, sogb_ndx, sogx_ndx
  
    logical :: alkooh_dd, mekooh_dd, tolooh_dd, terpooh_dd, ch3cooh_dd
    logical :: soa_dd, so4_dd, cb1_dd, cb2_dd, oc1_dd, oc2_dd, nh3_dd, nh4no3_dd, &
               sa1_dd, sa2_dd, sa3_dd, sa4_dd, nh4_dd
    logical :: soam_dd, soai_dd, soat_dd, soab_dd, soax_dd, &
               sogm_dd, sogi_dd, sogt_dd, sogb_dd, sogx_dd
  
    logical :: pan_dd, mpan_dd, no2_dd, hno3_dd, o3_dd, isopooh_dd, ch4_dd,&
               h2o2_dd, onit_dd, onitr_dd, ch2o_dd, macrooh_dd, xooh_dd, &
               ch3ooh_dd, pooh_dd, ch3coooh_dd, c2h5ooh_dd, eooh_dd, ch3cho_dd, c2h5oh_dd, &
               c3h7ooh_dd, rooh_dd, ch3cocho_dd, co_dd, ch3coch3_dd, &
               glyald_dd, hyac_dd, ch3oh_dd, hydrald_dd, h2_dd, Pb_dd, o3s_dd, o3inert_dd
  
    integer :: so2_ndx
    integer :: ch3cn_ndx, hcn_ndx, hcooh_ndx
    logical :: ch3cn_dd,  hcn_dd, hcooh_dd
  
    integer :: o3a_ndx,xpan_ndx,xmpan_ndx,xno2_ndx,xhno3_ndx,xonit_ndx,xonitr_ndx,xno_ndx,xho2no2_ndx,xnh4no3_ndx
    logical :: o3a_dd, xpan_dd, xmpan_dd, xno2_dd, xhno3_dd, xonit_dd, xonitr_dd, xno_dd, xho2no2_dd, xnh4no3_dd
  
  ! chemUCI
    integer :: no3_ndx, n2o5_ndx
    logical :: no3_dd,  n2o5_dd
  
    integer :: cohc_ndx=-1, come_ndx=-1, co01_ndx=-1, co02_ndx=-1, co03_ndx=-1, co04_ndx=-1, co05_ndx=-1
    integer :: co06_ndx=-1, co07_ndx=-1, co08_ndx=-1, co09_ndx=-1, co10_ndx=-1
    integer :: co11_ndx=-1, co12_ndx=-1, co13_ndx=-1, co14_ndx=-1, co15_ndx=-1
    integer :: co16_ndx=-1, co17_ndx=-1, co18_ndx=-1, co19_ndx=-1, co20_ndx=-1
    integer :: co21_ndx=-1, co22_ndx=-1, co23_ndx=-1, co24_ndx=-1, co25_ndx=-1
    integer :: co26_ndx=-1, co27_ndx=-1, co28_ndx=-1, co29_ndx=-1, co30_ndx=-1
    integer :: co31_ndx=-1, co32_ndx=-1, co33_ndx=-1, co34_ndx=-1, co35_ndx=-1
    integer :: co36_ndx=-1, co37_ndx=-1, co38_ndx=-1, co39_ndx=-1, co40_ndx=-1
    integer :: co41_ndx=-1, co42_ndx=-1
  
  
    integer :: &
         o3_tab_ndx = -1, &
         h2o2_tab_ndx = -1, &
         ch3ooh_tab_ndx = -1, &
         co_tab_ndx = -1, &
         ch3cho_tab_ndx = -1
    logical :: &
         o3_in_tab = .false., &
         h2o2_in_tab = .false., &
         ch3ooh_in_tab = .false., &
         co_in_tab = .false., &
         ch3cho_in_tab = .false.
  
    real(r8), parameter    :: small_value = 1.e-36_r8
    real(r8), parameter    :: large_value = 1.e36_r8
    real(r8), parameter    :: diffm       = 1.789e-5_r8
    real(r8), parameter    :: diffk       = 1.461e-5_r8
    real(r8), parameter    :: difft       = 2.060e-5_r8
    real(r8), parameter    :: vonkar      = karman
    real(r8), parameter    :: ric         = 0.2_r8
    real(r8), parameter    :: r           = 287.04_r8
    real(r8), parameter    :: cp          = 1004._r8
    real(r8), parameter    :: grav        = 9.81_r8
    real(r8), parameter    :: p00         = 100000._r8
    real(r8), parameter    :: wh2o        = 18.0153_r8
    real(r8), parameter    :: ph          = 1.e-5_r8
    real(r8), parameter    :: ph_inv      = 1._r8/ph
    real(r8), parameter    :: rovcp = r/cp
  
    integer, pointer :: index_season_lai(:,:)
  
    logical, public :: has_dvel(gas_pcnst) = .false.
    integer         :: map_dvel(gas_pcnst) = 0
    real(r8) , allocatable            :: soilw_3d(:,:,:)
  
    logical, parameter :: dyn_soilw = .false.
  
    real(r8), allocatable  :: gc_fraction_landuse(:,:,:)
    real(r8), allocatable, dimension(:,:,:) :: dep_ra ! [s/m] aerodynamic resistance
    real(r8), allocatable, dimension(:,:,:) :: dep_rb ! [s/m] resistance across sublayer
    integer, parameter :: gc_n_land_type = 11
  
    integer, allocatable :: spc_ndx(:) ! nddvels
    real(r8), public :: crb 
  
    type lnd_dvel_type
       real(r8), pointer :: dvel(:,:)   ! deposition velocity over land (cm/s)
    end type lnd_dvel_type
  
    type(lnd_dvel_type), allocatable :: lnd(:)
    character(len=SHR_KIND_CL) :: gc_drydep_srf_file
  
  contains

    !---------------------------------------------------------------------------
    ! subroutine gc_dvel_inti_fromlnd 
    !
    ! Refactor Notes: 
    !   - Added 'use' statements for dependencies
    !   - ppgrid must be initialized beforehand to set 'begchunk' and 'endchunk'
    !
    ! FIX ME: Ideally we should get away from using the index list approach
    !         if possible
    !---------------------------------------------------------------------------
    subroutine gc_dvel_inti_fromlnd() 
      use mo_chem_utls,         only: get_spc_ndx
      use cam_abortutils,       only: endrun
      use chem_mods,            only: adv_mass
      use seq_drydep_mod,       only: dfoxd
      use ppgrid,               only: begchunk, endchunk
      implicit none
      !------- Local Variables ---------
      integer :: ispc, l
      !---------------------------------
  
      ! FIX ME: Unclear what 'begchunk' and 'endchunk' should be set to
      allocate(spc_ndx(nddvels))
      allocate( lnd(begchunk:endchunk) )
  
      do ispc = 1,nddvels
         spc_ndx(ispc) = get_spc_ndx(drydep_list(ispc))
         if (spc_ndx(ispc) < 1) then
            write(*,*) 'gc_drydep_inti: '//trim(drydep_list(ispc))//' is not included in species set'
            call endrun('gc_drydep_init: invalid dry deposition species')
         endif
      enddo
  
      crb = (difft/diffm)**(2._r8/3._r8) !.666666_r8
    end subroutine gc_dvel_inti_fromlnd


    !----------------------------------------------------------------------------------------------
    ! subroutine gc_drydep_fromlnd
    !
    ! Computes dry deposition from land??
    !
    !----------------------------------------------------------------------------------------------
    subroutine gc_drydep_fromlnd( ocnfrac, icefrac, ncdate, sfc_temp, pressure_sfc,  &
                                  wind_speed, spec_hum, air_temp, pressure_10m, rain, &
                                  snow, solar_flux, dvelocity, dflx, mmr, &
                                  tv, soilw, rh, ncol, lonndx, latndx, lchnk )
        
      !-------------------------------------------------------------------------------------
      ! combines the deposition velocities provided by the land model with deposition 
      ! velocities over ocean and sea ice 
      !-------------------------------------------------------------------------------------
      use ppgrid,         only : pcols
      use pmgrid,         only : plev, plevp
      use chem_mods,      only : gas_pcnst

#if (defined OFFLINE_DYN)
      use metdata, only: get_met_fields
#endif

      implicit none

      !-------------------------------------------------------------------------------------
      ! 	... dummy arguments
      !-------------------------------------------------------------------------------------

      real(r8), intent(in)    :: icefrac(pcols)            
      real(r8), intent(in)    :: ocnfrac(pcols)            

      integer, intent(in)   :: ncol
      integer, intent(in)   :: ncdate                   ! present date (yyyymmdd)
      real(r8), intent(in)      :: sfc_temp(pcols)          ! surface temperature (K)
      real(r8), intent(in)      :: pressure_sfc(pcols)      ! surface pressure (Pa)
      real(r8), intent(in)      :: wind_speed(pcols)        ! 10 meter wind speed (m/s)
      real(r8), intent(in)      :: spec_hum(pcols)          ! specific humidity (kg/kg)
      real(r8), intent(in)      :: rh(ncol,1)               ! relative humidity
      real(r8), intent(in)      :: air_temp(pcols)          ! surface air temperature (K)
      real(r8), intent(in)      :: pressure_10m(pcols)      ! 10 meter pressure (Pa)
      real(r8), intent(in)      :: rain(pcols)              
      real(r8), intent(in)      :: snow(pcols)              ! snow height (m)
      real(r8), intent(in)      :: soilw(pcols)             ! soil moisture fraction
      real(r8), intent(in)      :: solar_flux(pcols)        ! direct shortwave radiation at surface (W/m^2)
      real(r8), intent(in)      :: tv(pcols)                ! potential temperature
      real(r8), intent(in)      :: mmr(pcols,plev,gas_pcnst)    ! constituent concentration (kg/kg)
      real(r8), intent(out)     :: dvelocity(ncol,gas_pcnst)    ! deposition velocity (cm/s)
      real(r8), intent(inout)   :: dflx(pcols,gas_pcnst)        ! deposition flux (/cm^2/s)

      integer, intent(in)     ::   latndx(pcols)           ! chunk latitude indicies
      integer, intent(in)     ::   lonndx(pcols)           ! chunk longitude indicies
      integer, intent(in)     ::   lchnk                   ! chunk number

      !-------------------------------------------------------------------------------------
      ! 	... local variables
      !-------------------------------------------------------------------------------------
      real(r8) :: ocnice_dvel(ncol,gas_pcnst)
      real(r8) :: ocnice_dflx(pcols,gas_pcnst)

      real(r8), dimension(ncol) :: term    ! work array
      integer  :: ispec
      real(r8)  :: lndfrac(pcols)            
#if (defined OFFLINE_DYN)
      real(r8)  :: met_ocnfrac(pcols)
      real(r8)  :: met_icefrac(pcols)            
#endif

      lndfrac(:ncol) = 1._r8 - ocnfrac(:ncol) - icefrac(:ncol)

      where( lndfrac(:ncol) < 0._r8 ) 
      lndfrac(:ncol) = 0._r8 
      endwhere

#if (defined OFFLINE_DYN)
      call get_met_fields(lndfrac, met_ocnfrac, met_icefrac, lchnk, ncol)
#endif

      !-------------------------------------------------------------------------------------
      !   ... initialize
      !-------------------------------------------------------------------------------------
      dvelocity(:,:) = 0._r8

      !-------------------------------------------------------------------------------------
      !   ... compute the dep velocities over ocean and sea ice
      !       land type 7 is used for ocean
      !       land type 8 is used for sea ice
      !-------------------------------------------------------------------------------------
      call gc_drydep_xactive( ncdate, sfc_temp, pressure_sfc,  &
        wind_speed, spec_hum, air_temp, pressure_10m, rain, &
        snow, solar_flux, ocnice_dvel, ocnice_dflx, mmr, &
        tv, soilw, rh, ncol, lonndx, latndx, lchnk, &
#if (defined OFFLINE_DYN)
        ocnfrc=met_ocnfrac,icefrc=met_icefrac, beglandtype=7, endlandtype=8 )
#else
        ocnfrc=ocnfrac,icefrc=icefrac, beglandtype=7, endlandtype=8 )
#endif
      term(:ncol) = 1.e-2_r8 * pressure_10m(:ncol) / (r*tv(:ncol))

      species_loop3 : do ispec = 1,nddvels

        !-------------------------------------------------------------------------------------
        !        ... merge the land component with the non-land component
        !            ocn and ice already have fractions factored in
        !-------------------------------------------------------------------------------------
        dvelocity(:ncol,spc_ndx(ispec)) = lnd(lchnk)%dvel(:ncol,ispec)*lndfrac(:ncol) &
                  + ocnice_dvel(:ncol,spc_ndx(ispec))


        !-------------------------------------------------------------------------------------
        !        ... special adjustments
        !-------------------------------------------------------------------------------------
        if( spc_ndx(ispec) == mpan_ndx .or. spc_ndx(ispec) == xmpan_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,spc_ndx(ispec))/3._r8
        endif
        if( spc_ndx(ispec) == hcn_ndx .or. spc_ndx(ispec) == ch3cn_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = ocnice_dvel(:ncol,spc_ndx(ispec)) ! should be zero over land
        endif

        ! HCOOH, use CH3COOH dep.vel
        if( hcooh_ndx > 0 .and. ch3cooh_ndx > 0 ) then
          if( has_dvel(hcooh_ndx) ) then
            dvelocity(:ncol,hcooh_ndx) = dvelocity(:ncol,ch3cooh_ndx)
          end if
        end if

        !lke++
        !-------------------------------------------------------------------------------------
        !        ... assign CO tags to CO
        ! put this kludge in for now ...  
        !  -- should be able to set all these via the table mapping in seq_drydep_mod
        !-------------------------------------------------------------------------------------
        if( spc_ndx(ispec) == cohc_ndx  .or. spc_ndx(ispec) == come_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co01_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co02_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co03_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co04_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co05_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co06_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co07_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co08_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co09_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co10_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co11_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co12_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co13_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co14_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co15_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co16_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co17_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co18_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co19_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co20_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co21_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co22_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co23_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co24_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co25_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co26_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co27_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co28_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co29_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co30_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co31_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co32_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co33_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co34_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co35_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co36_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co37_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co38_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co39_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co40_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co41_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        if( spc_ndx(ispec) == co42_ndx ) then
          dvelocity(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,co_ndx)
        endif
        !lke--

        !-------------------------------------------------------------------------------------
        !        ... compute the deposition flux
        !-------------------------------------------------------------------------------------
        dflx(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,spc_ndx(ispec)) * term(:ncol) * mmr(:ncol,plev,spc_ndx(ispec))

      end do species_loop3

    end subroutine gc_drydep_fromlnd



    !-------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------
    subroutine gc_drydep_xactive( ncdate, sfc_temp, pressure_sfc,  &
                                  wind_speed, spec_hum, air_temp, pressure_10m, rain, &
                                  snow, solar_flux, dvel, dflx, mmr, &
                                  tv, soilw, rh, ncol, lonndx, latndx, lchnk, &
                                  ocnfrc, icefrc, beglandtype, endlandtype )
      !-------------------------------------------------------------------------------------
      !   code based on wesely (atmospheric environment, 1989, vol 23, p. 1293-1304) for
      !   calculation of r_c, and on walcek et. al. (atmospheric enviroment, 1986,
      !   vol. 20, p. 949-964) for calculation of r_a and r_b
      !
      !   as suggested in walcek (u_i)(u*_i) = (u_a)(u*_a)
      !   is kept constant where i represents a subgrid environment and a the
      !   grid average environment. thus the calculation proceeds as follows:
      !   va the grid averaged wind is calculated on dots
      !   z0(i) the grid averaged roughness coefficient is calculated
      !   ri(i) the grid averaged richardson number is calculated
      !   --> the grid averaged (u_a)(u*_a) is calculated
      !   --> subgrid scale u*_i is calculated assuming (u_i) given as above
      !   --> final deposotion velocity is weighted average of subgrid scale velocities
      !
      ! code written by P. Hess, rewritten in fortran 90 by JFL (August 2000)
      ! modified by JFL to be used in MOZART-2 (October 2002)
      !-------------------------------------------------------------------------------------

      use seq_drydep_mod, only: z0, rgso, rgss, h2_a, h2_b, h2_c, ri, rclo, rcls, rlu, rac
      use seq_drydep_mod, only: seq_drydep_setHCoeff, foxd, drat
      use physconst,      only: tmelt
      use seq_drydep_mod, only: drydep_method,  DD_XLND

      implicit none

      !-------------------------------------------------------------------------------------
      ! 	... dummy arguments
      !-------------------------------------------------------------------------------------
      integer, intent(in)   :: ncol
      integer, intent(in)   :: ncdate                   ! present date (yyyymmdd)
      real(r8), intent(in)      :: sfc_temp(pcols)          ! surface temperature (K)
      real(r8), intent(in)      :: pressure_sfc(pcols)      ! surface pressure (Pa)
      real(r8), intent(in)      :: wind_speed(pcols)        ! 10 meter wind speed (m/s)
      real(r8), intent(in)      :: spec_hum(pcols)          ! specific humidity (kg/kg)
      real(r8), intent(in)      :: rh(ncol,1)               ! relative humidity
      real(r8), intent(in)      :: air_temp(pcols)          ! surface air temperature (K)
      real(r8), intent(in)      :: pressure_10m(pcols)      ! 10 meter pressure (Pa)
      real(r8), intent(in)      :: rain(pcols)              
      real(r8), intent(in)      :: snow(pcols)              ! snow height (m)
      real(r8), intent(in)      :: soilw(pcols)             ! soil moisture fraction
      real(r8), intent(in)      :: solar_flux(pcols)        ! direct shortwave radiation at surface (W/m^2)
      real(r8), intent(in)      :: tv(pcols)                ! potential temperature
      real(r8), intent(in)      :: mmr(pcols,plev,gas_pcnst)    ! constituent concentration (kg/kg)
      real(r8), intent(out)     :: dvel(ncol,gas_pcnst)        ! deposition velocity (cm/s)
      real(r8), intent(inout)   :: dflx(pcols,gas_pcnst)        ! deposition flux (/cm^2/s)

      integer, intent(in)     ::   latndx(pcols)           ! chunk latitude indicies
      integer, intent(in)     ::   lonndx(pcols)           ! chunk longitude indicies
      integer, intent(in)     ::   lchnk                   ! chunk number

      integer, intent(in), optional     ::  beglandtype
      integer, intent(in), optional     ::  endlandtype

      real(r8), intent(in), optional      :: ocnfrc(pcols) 
      real(r8), intent(in), optional      :: icefrc(pcols) 

      !-------------------------------------------------------------------------------------
      ! 	... local variables
      !-------------------------------------------------------------------------------------
      real(r8), parameter :: scaling_to_cm_per_s = 100._r8
      real(r8), parameter :: rain_threshold      = 1.e-7_r8  ! of the order of 1cm/day expressed in m/s

      integer :: i, ispec, lt, m
      integer :: sndx
      integer :: month

      real(r8) :: slope = 0._r8
      real(r8) :: z0water ! revised z0 over water
      real(r8) :: p       ! pressure at midpoint first layer
      real(r8) :: pg      ! surface pressure
      real(r8) :: es      ! saturation vapor pressure
      real(r8) :: ws      ! saturation mixing ratio
      real(r8) :: hvar    ! constant to compute xmol
      real(r8) :: h       ! constant to compute xmol
      real(r8) :: psih    ! stability correction factor
      real(r8) :: rs      ! constant for calculating rsmx
      real(r8) :: rmx     ! resistance by vegetation
      real(r8) :: zovl    ! ratio of z to  m-o length
      real(r8) :: cvarb   ! cvar averaged over landtypes
      real(r8) :: bb      ! b averaged over landtypes
      real(r8) :: ustarb  ! ustar averaged over landtypes
      real(r8) :: tc(ncol)  ! temperature in celsius
      real(r8) :: cts(ncol) ! correction to rlu rcl and rgs for frost

      !-------------------------------------------------------------------------------------
      ! local arrays: dependent on location and species
      !-------------------------------------------------------------------------------------
      real(r8), dimension(ncol,nddvels) :: heff

      !-------------------------------------------------------------------------------------
      ! local arrays: dependent on location only
      !-------------------------------------------------------------------------------------
      integer                   :: index_season(ncol,gc_n_land_type)
      real(r8), dimension(ncol) :: tha     ! atmospheric virtual potential temperature
      real(r8), dimension(ncol) :: thg     ! ground virtual potential temperature
      real(r8), dimension(ncol) :: z       ! height of lowest level
      real(r8), dimension(ncol) :: va      ! magnitude of v on cross points
      real(r8), dimension(ncol) :: ribn    ! richardson number
      real(r8), dimension(ncol) :: qs      ! saturation specific humidity
      real(r8), dimension(ncol) :: crs     ! multiplier to calculate crs
      real(r8), dimension(ncol) :: rdc     ! part of lower canopy resistance
      real(r8), dimension(ncol) :: uustar  ! u*ustar (assumed constant over grid)
      real(r8), dimension(ncol) :: z0b     ! average roughness length over grid
      real(r8), dimension(ncol) :: wrk     ! work array
      real(r8), dimension(ncol) :: term    ! work array
      real(r8), dimension(ncol) :: resc    ! work array
      real(r8), dimension(ncol) :: lnd_frc ! work array
      logical,  dimension(ncol) :: unstable
      logical,  dimension(ncol) :: has_rain
      logical,  dimension(ncol) :: has_dew

      !-------------------------------------------------------------------------------------
      ! local arrays: dependent on location and landtype
      !-------------------------------------------------------------------------------------
      real(r8), dimension(ncol,gc_n_land_type) :: rds   ! resistance for deposition of sulfate
      real(r8), dimension(ncol,gc_n_land_type) :: b     ! buoyancy parameter for unstable conditions
      real(r8), dimension(ncol,gc_n_land_type) :: cvar  ! height parameter
      real(r8), dimension(ncol,gc_n_land_type) :: ustar ! friction velocity
      real(r8), dimension(ncol,gc_n_land_type) :: xmol  ! monin-obukhov length

      !-------------------------------------------------------------------------------------
      ! local arrays: dependent on location, landtype and species
      !-------------------------------------------------------------------------------------
      real(r8), dimension(ncol,gc_n_land_type,gas_pcnst) :: rsmx  ! vegetative resistance (plant mesophyll)
      real(r8), dimension(ncol,gc_n_land_type,gas_pcnst) :: rclx  ! lower canopy resistance
      real(r8), dimension(ncol,gc_n_land_type,gas_pcnst) :: rlux  ! vegetative resistance (upper canopy)
      real(r8), dimension(ncol,gc_n_land_type) :: rlux_o3  ! vegetative resistance (upper canopy)
      real(r8), dimension(ncol,gc_n_land_type,gas_pcnst) :: rgsx  ! ground resistance
      real(r8) :: pmid(ncol,1)                             ! for seasalt aerosols
      real(r8) :: tfld(ncol,1)                             ! for seasalt aerosols
      real(r8) :: fact, vds
      real(r8) :: rc                                    ! combined surface resistance
      real(r8) :: var_soilw, dv_soil_h2, fact_h2        ! h2 dvel wrking variables
      logical  :: fr_lnduse(ncol,gc_n_land_type)           ! wrking array
      real(r8) :: dewm                                  ! multiplier for rs when dew occurs

      real(r8) :: lcl_frc_landuse(ncol,gc_n_land_type) 

      integer :: beglt, endlt

      !-------------------------------------------------------------------------------------
      ! jfl : mods for PAN
      !-------------------------------------------------------------------------------------
      real(r8) :: dv_pan
      real(r8) :: c0_pan(11) = (/ 0.000_r8, 0.006_r8, 0.002_r8, 0.009_r8, 0.015_r8, &
                                  0.006_r8, 0.000_r8, 0.000_r8, 0.000_r8, 0.002_r8, 0.002_r8 /)
      real(r8) :: k_pan (11) = (/ 0.000_r8, 0.010_r8, 0.005_r8, 0.004_r8, 0.003_r8, &
                                  0.005_r8, 0.000_r8, 0.000_r8, 0.000_r8, 0.075_r8, 0.002_r8 /)

      if (present( beglandtype)) then
        beglt = beglandtype 
      else
        beglt = 1
      endif
      if (present( endlandtype)) then
        endlt = endlandtype 
      else
        endlt = gc_n_land_type
      endif

      !-------------------------------------------------------------------------------------
      ! initialize
      !-------------------------------------------------------------------------------------
      do m = gcso4_ndx,gas_pcnst
        dvel(:,m) = 0._r8
      end do

      if( all( .not. has_dvel(:) ) ) then
        return
      end if

      !-------------------------------------------------------------------------------------
      ! define species-dependent parameters (temperature dependent)
      !-------------------------------------------------------------------------------------
      call seq_drydep_setHCoeff( ncol, sfc_temp, heff )

      do lt = 1,gc_n_land_type
        dep_ra (:,lt,lchnk)   = 0._r8
        dep_rb (:,lt,lchnk)   = 0._r8
        rds(:,lt)   = 0._r8
      end do

      !-------------------------------------------------------------------------------------
      ! 	... set month
      !-------------------------------------------------------------------------------------
      month = mod( ncdate,10000 )/100

      !-------------------------------------------------------------------------------------
      ! define which season (relative to Northern hemisphere climate)
      !-------------------------------------------------------------------------------------

      !-------------------------------------------------------------------------------------
      ! define season index based on fixed LAI
      !-------------------------------------------------------------------------------------
      if ( drydep_method == DD_XLND ) then
        index_season = 4
      else
        do i = 1,ncol
          index_season(i,:) = index_season_lai(latndx(i),month)
        end do
      endif
      !-------------------------------------------------------------------------------------
      ! special case for snow covered terrain
      !-------------------------------------------------------------------------------------
      do i = 1,ncol
        if( snow(i) > .01_r8 ) then
          index_season(i,:) = 4
        end if
      end do
      !-------------------------------------------------------------------------------------
      ! scale rain and define logical arrays
      !-------------------------------------------------------------------------------------
      has_rain(:ncol) = rain(:ncol) > rain_threshold

      !-------------------------------------------------------------------------------------
      ! loop over longitude points
      !-------------------------------------------------------------------------------------
      col_loop :  do i = 1,ncol
        p   = pressure_10m(i)
        pg  = pressure_sfc(i)
        !-------------------------------------------------------------------------------------
        ! potential temperature
        !-------------------------------------------------------------------------------------
        tha(i) = air_temp(i) * (p00/p )**rovcp * (1._r8 + .61_r8*spec_hum(i))
        thg(i) = sfc_temp(i) * (p00/pg)**rovcp * (1._r8 + .61_r8*spec_hum(i))
        !-------------------------------------------------------------------------------------
        ! height of 1st level
        !-------------------------------------------------------------------------------------
        z(i) = - r/grav * air_temp(i) * (1._r8 + .61_r8*spec_hum(i)) * log(p/pg)
        !-------------------------------------------------------------------------------------
        ! wind speed
        !-------------------------------------------------------------------------------------
        va(i) = max( .01_r8,wind_speed(i) )
        !-------------------------------------------------------------------------------------
        ! Richardson number
        !-------------------------------------------------------------------------------------
        ribn(i) = z(i) * grav * (tha(i) - thg(i))/thg(i) / (va(i)*va(i))
        ribn(i) = min( ribn(i),ric )
        unstable(i) = ribn(i) < 0._r8
        !-------------------------------------------------------------------------------------
        ! saturation vapor pressure (Pascals)
        ! saturation mixing ratio
        ! saturation specific humidity
        !-------------------------------------------------------------------------------------
        es    = 611._r8*exp( 5414.77_r8*(sfc_temp(i) - tmelt)/(tmelt*sfc_temp(i)) )
        ws    = .622_r8*es/(pg - es)
        qs(i) = ws/(1._r8 + ws)
        has_dew(i) = .false.
        if( qs(i) <= spec_hum(i) ) then
          has_dew(i) = .true.
        end if
        if( sfc_temp(i) < tmelt ) then
          has_dew(i) = .false.
        end if
        !-------------------------------------------------------------------------------------
        ! constant in determining rs
        !-------------------------------------------------------------------------------------
        tc(i) = sfc_temp(i) - tmelt
        if( sfc_temp(i) > tmelt .and. sfc_temp(i) < 313.15_r8 ) then
          crs(i) = (1._r8 + (200._r8/(solar_flux(i) + .1_r8))**2) * (400._r8/(tc(i)*(40._r8 - tc(i))))
        else
          crs(i) = large_value
        end if
        !-------------------------------------------------------------------------------------
        ! rdc (lower canopy res)
        !-------------------------------------------------------------------------------------
        rdc(i) = 100._r8*(1._r8 + 1000._r8/(solar_flux(i) + 10._r8))/(1._r8 + 1000._r8*slope)
      end do col_loop

      !-------------------------------------------------------------------------------------
      ! 	... form working arrays
      !-------------------------------------------------------------------------------------
      do lt = 1,gc_n_land_type
        do i=1,ncol
          if ( drydep_method == DD_XLND ) then
            lcl_frc_landuse(i,lt) = 0._r8
          else
            lcl_frc_landuse(i,lt) = gc_fraction_landuse(i,lt,lchnk)
          endif
        enddo
      end do
      if ( present(ocnfrc) .and. present(icefrc) ) then
        do i=1,ncol
          ! land type 7 is used for ocean
          ! land type 8 is used for sea ice
          lcl_frc_landuse(i,7) = ocnfrc(i)
          lcl_frc_landuse(i,8) = icefrc(i)
        enddo
      endif
      do lt = 1,gc_n_land_type
        do i=1,ncol
          fr_lnduse(i,lt) = lcl_frc_landuse(i,lt) > 0._r8
        enddo
      end do

      !-------------------------------------------------------------------------------------
      ! find grid averaged z0: z0bar (the roughness length) z_o=exp[S(f_i*ln(z_oi))]
      ! this is calculated so as to find u_i, assuming u*u=u_i*u_i
      !-------------------------------------------------------------------------------------
      z0b(:) = 0._r8
      do lt = 1,gc_n_land_type
        do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
            z0b(i) = z0b(i) + lcl_frc_landuse(i,lt) * log( z0(index_season(i,lt),lt) )
          end if
        end do
      end do

      !-------------------------------------------------------------------------------------
      ! find the constant velocity uu*=(u_i)(u*_i)
      !-------------------------------------------------------------------------------------
      do i = 1,ncol
        z0b(i) = exp( z0b(i) )
        cvarb  = vonkar/log( z(i)/z0b(i) )
        !-------------------------------------------------------------------------------------
        ! unstable and stable cases
        !-------------------------------------------------------------------------------------
        if( unstable(i) ) then
          bb = 9.4_r8*(cvarb**2)*sqrt( abs(ribn(i))*z(i)/z0b(i) )
          ustarb = cvarb * va(i) * sqrt( 1._r8 - (9.4_r8*ribn(i)/(1._r8 + 7.4_r8*bb)) )
        else
          ustarb = cvarb * va(i)/(1._r8 + 4.7_r8*ribn(i))
        end if
        uustar(i) = va(i)*ustarb
      end do

      !-------------------------------------------------------------------------------------
      ! calculate the friction velocity for each land type u_i=uustar/u*_i
      !-------------------------------------------------------------------------------------
      do lt = beglt,endlt
        do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
            if( unstable(i) ) then
              cvar(i,lt)  = vonkar/log( z(i)/z0(index_season(i,lt),lt) )
              b(i,lt)     = 9.4_r8*(cvar(i,lt)**2)* sqrt( abs(ribn(i))*z(i)/z0(index_season(i,lt),lt) )
              ustar(i,lt) = sqrt( cvar(i,lt)*uustar(i)*sqrt( 1._r8 - (9.4_r8*ribn(i)/(1._r8 + 7.4_r8*b(i,lt))) ) )
            else
              cvar(i,lt)  = vonkar/log( z(i)/z0(index_season(i,lt),lt) )
              ustar(i,lt) = sqrt( cvar(i,lt)*uustar(i)/(1._r8 + 4.7_r8*ribn(i)) )
            end if
          end if
        end do
      end do

      !-------------------------------------------------------------------------------------
      ! revise calculation of friction velocity and z0 over water
      !-------------------------------------------------------------------------------------
      lt = 7    
      do i = 1,ncol
        if( fr_lnduse(i,lt) ) then
          if( unstable(i) ) then
            z0water     = (.016_r8*(ustar(i,lt)**2)/grav) + diffk/(9.1_r8*ustar(i,lt))
            cvar(i,lt)  = vonkar/(log( z(i)/z0water ))
            b(i,lt)     = 9.4_r8*(cvar(i,lt)**2)*sqrt( abs(ribn(i))*z(i)/z0water )
            ustar(i,lt) = sqrt( cvar(i,lt)*uustar(i)* sqrt( 1._r8 - (9.4_r8*ribn(i)/(1._r8+ 7.4_r8*b(i,lt))) ) )
          else
            z0water     = (.016_r8*(ustar(i,lt)**2)/grav) + diffk/(9.1_r8*ustar(i,lt))
            cvar(i,lt)  = vonkar/(log(z(i)/z0water))
            ustar(i,lt) = sqrt( cvar(i,lt)*uustar(i)/(1._r8 + 4.7_r8*ribn(i)) )
          end if
        end if
      end do

      !-------------------------------------------------------------------------------------
      ! compute monin-obukhov length for unstable and stable conditions/ sublayer resistance
      !-------------------------------------------------------------------------------------
      do lt = beglt,endlt
        do i = 1,ncol
        if( fr_lnduse(i,lt) ) then
          hvar = (va(i)/0.74_r8) * (tha(i) - thg(i)) * (cvar(i,lt)**2)
          if( unstable(i) ) then                      ! unstable
            h = hvar*(1._r8 - (9.4_r8*ribn(i)/(1._r8 + 5.3_r8*b(i,lt))))
          else
            h = hvar/((1._r8+4.7_r8*ribn(i))**2)
          end if
          xmol(i,lt) = thg(i) * ustar(i,lt) * ustar(i,lt) / (vonkar * grav * h)
        end if
        end do
      end do

      !-------------------------------------------------------------------------------------
      ! psih
      !-------------------------------------------------------------------------------------
      do lt = beglt,endlt
        do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
            if( xmol(i,lt) < 0._r8 ) then
              zovl = z(i)/xmol(i,lt)
              zovl = max( -1._r8,zovl )
              psih = exp( .598_r8 + .39_r8*log( -zovl ) - .09_r8*(log( -zovl ))**2 )
              vds  = 2.e-3_r8*ustar(i,lt) * (1._r8 + (300/(-xmol(i,lt)))**0.666_r8)
            else
              zovl = z(i)/xmol(i,lt)
              zovl = min( 1._r8,zovl )
              psih = -5._r8 * zovl
              vds  = 2.e-3_r8*ustar(i,lt)
            end if
            dep_ra (i,lt,lchnk) = (vonkar - psih*cvar(i,lt))/(ustar(i,lt)*vonkar*cvar(i,lt))
            dep_rb (i,lt,lchnk) = (2._r8/(vonkar*ustar(i,lt))) * crb
            rds(i,lt) = 1._r8/vds
          end if
        end do
      end do

      !-------------------------------------------------------------------------------------
      ! surface resistance : depends on both land type and species
      ! land types are computed seperately, then resistance is computed as average of values
      ! following wesely rc=(1/(rs+rm) + 1/rlu +1/(rdc+rcl) + 1/(rac+rgs))**-1
      !
      ! compute rsmx = 1/(rs+rm) : multiply by 3 if surface is wet
      !-------------------------------------------------------------------------------------
      rlux = 0._r8
      species_loop1 :  do ispec = gcso4_ndx,gas_pcnst
        if( has_dvel(ispec) ) then
          m = map_dvel(ispec)
          do lt = beglt,endlt
            do i = 1,ncol
              if( fr_lnduse(i,lt) ) then
                sndx = index_season(i,lt)
                if( ispec == o3_ndx .or. ispec == o3a_ndx .or. ispec == so2_ndx ) then
                  rmx = 0._r8
                else
                  rmx = 1._r8/(heff(i,m)/3000._r8 + 100._r8*foxd(m))
                end if
                cts(i) = 1000._r8*exp( - tc(i) - 4._r8 )                 ! correction for frost
                rgsx(i,lt,ispec) = cts(i) + 1._r8/((heff(i,m)/(1.e5_r8*rgss(sndx,lt))) + (foxd(m)/rgso(sndx,lt)))
                !-------------------------------------------------------------------------------------
                ! special case for H2 and CO;; CH4 is set ot a fraction of dv(H2)
                !-------------------------------------------------------------------------------------
                if( ispec == h2_ndx .or. ispec == co_ndx .or. ispec == ch4_ndx ) then
                  if( ispec == co_ndx ) then
                    fact_h2 = 1.0_r8
                  elseif ( ispec == h2_ndx ) then
                    fact_h2 = 0.5_r8
                  elseif ( ispec == ch4_ndx ) then
                    fact_h2 = 50.0_r8
                  end if
                  !-------------------------------------------------------------------------------------
                  ! no deposition on snow, ice, desert, and water
                  !-------------------------------------------------------------------------------------
                  if( lt == 1 .or. lt == 7 .or. lt == 8 .or. sndx == 4 ) then
                    rgsx(i,lt,ispec) = large_value
                  else
                    var_soilw = max( .1_r8,min( soilw(i),.3_r8 ) )
                    if( lt == 3 ) then
                      var_soilw = log( var_soilw )
                    end if
                    dv_soil_h2 = h2_c(lt) + var_soilw*(h2_b(lt) + var_soilw*h2_a(lt))
                    if( dv_soil_h2 > 0._r8 ) then
                      rgsx(i,lt,ispec) = fact_h2/(dv_soil_h2*1.e-4_r8)
                    end if
                  end if
                end if
                if( lt == 7 ) then
                  rclx(i,lt,ispec) = large_value
                  rsmx(i,lt,ispec) = large_value
                  rlux(i,lt,ispec) = large_value
                else
                  rs = ri(sndx,lt)*crs(i)
                  if ( has_dew(i) .or. has_rain(i) ) then
                    dewm = 3._r8
                  else
                    dewm = 1._r8
                  end if
                  rsmx(i,lt,ispec) = (dewm*rs*drat(m) + rmx)
                  !-------------------------------------------------------------------------------------
                  ! jfl : special case for PAN
                  !-------------------------------------------------------------------------------------
                  if( ispec == pan_ndx .or. ispec == xpan_ndx ) then
                    dv_pan =  c0_pan(lt) * (1._r8 - exp( -k_pan(lt)*(dewm*rs*drat(m))*1.e-2_r8 ))
                    if( dv_pan > 0._r8 .and. sndx /= 4 ) then
                      rsmx(i,lt,ispec) = ( 1._r8/dv_pan )
                    end if
                  end if
                  rclx(i,lt,ispec) = cts(i) + 1._r8/((heff(i,m)/(1.e5_r8*rcls(sndx,lt))) + (foxd(m)/rclo(sndx,lt)))
                  rlux(i,lt,ispec) = cts(i) + rlu(sndx,lt)/(1.e-5_r8*heff(i,m) + foxd(m))
                end if
              end if
            end do
          end do
        end if
      end do species_loop1

      do lt = beglt,endlt
        if( lt /= 7 ) then
          do i = 1,ncol
            if( fr_lnduse(i,lt) ) then
              sndx = index_season(i,lt)
              !-------------------------------------------------------------------------------------
              ! 	... no effect if sfc_temp < O C
              !-------------------------------------------------------------------------------------
              if( sfc_temp(i) > tmelt ) then
                if( has_dew(i) ) then
                  rlux_o3(i,lt)     = 3000._r8*rlu(sndx,lt)/(1000._r8 + rlu(sndx,lt))
                  if( o3_ndx > 0 ) then
                    rlux(i,lt,o3_ndx) = rlux_o3(i,lt)
                  endif
                  if( o3a_ndx > 0 ) then
                    rlux(i,lt,o3a_ndx) = rlux_o3(i,lt)
                  endif
                end if
                if( has_rain(i) ) then
                  ! rlux(i,lt,o3_ndx) = 1./(1.e-3 + (1./(3.*rlu(sndx,lt))))
                  rlux_o3(i,lt)     = 3000._r8*rlu(sndx,lt)/(1000._r8 + 3._r8*rlu(sndx,lt))
                  if( o3_ndx > 0 ) then
                    rlux(i,lt,o3_ndx) = rlux_o3(i,lt)
                  endif
                  if( o3a_ndx > 0 ) then
                    rlux(i,lt,o3a_ndx) = rlux_o3(i,lt)
                  endif
                end if
              end if

              if ( o3_ndx > 0 ) then
                rclx(i,lt,o3_ndx) = cts(i) + rclo(index_season(i,lt),lt)
                rlux(i,lt,o3_ndx) = cts(i) + rlux(i,lt,o3_ndx)
              end if
              if ( o3a_ndx > 0 ) then
                rclx(i,lt,o3a_ndx) = cts(i) + rclo(index_season(i,lt),lt)
                rlux(i,lt,o3a_ndx) = cts(i) + rlux(i,lt,o3a_ndx)
              end if

            end if
          end do
        end if
      end do

      species_loop2 : do ispec = gcso4_ndx,gas_pcnst
        m = map_dvel(ispec)
        if( has_dvel(ispec) ) then
          if( ispec /= o3_ndx .and. ispec /= o3a_ndx .and. ispec /= so2_ndx ) then
            do lt = beglt,endlt
              if( lt /= 7 ) then
                do i = 1,ncol
                  if( fr_lnduse(i,lt) ) then
                    !-------------------------------------------------------------------------------------
                    ! no effect if sfc_temp < O C
                    !-------------------------------------------------------------------------------------
                    if( sfc_temp(i) > tmelt ) then
                      if( has_dew(i) ) then
                          rlux(i,lt,ispec) = 1._r8/((1._r8/(3._r8*rlux(i,lt,ispec))) &
                              + 1.e-7_r8*heff(i,m) + foxd(m)/rlux_o3(i,lt))
                      end if
                    end if

                  end if
                end do
              end if
            end do
            else if( ispec == so2_ndx ) then
              do lt = beglt,endlt
                if( lt /= 7 ) then
                  do i = 1,ncol
                    if( fr_lnduse(i,lt) ) then
                      !-------------------------------------------------------------------------------------
                      ! no effect if sfc_temp < O C
                      !-------------------------------------------------------------------------------------
                      if( sfc_temp(i) > tmelt ) then
                        if( qs(i) <= spec_hum(i) ) then
                          rlux(i,lt,ispec) = 100._r8
                        end if
                        if( has_rain(i) ) then
                          !                               rlux(i,lt,ispec) = 1./(2.e-4 + (1./(3.*rlu(index_season(i,lt),lt))))
                          rlux(i,lt,ispec) = 15._r8*rlu(index_season(i,lt),lt)/(5._r8 + 3.e-3_r8*rlu(index_season(i,lt),lt))
                        end if
                      end if
                      rclx(i,lt,ispec) = cts(i) + rcls(index_season(i,lt),lt)
                      rlux(i,lt,ispec) = cts(i) + rlux(i,lt,ispec)

                  end if
                end do
              end if
            end do
            do i = 1,ncol
              if( fr_lnduse(i,1) .and. (has_dew(i) .or. has_rain(i)) ) then
                rlux(i,1,ispec) = 50._r8
              end if
            end do
          end if
        end if
      end do species_loop2

      !-------------------------------------------------------------------------------------
      ! compute rc
      !-------------------------------------------------------------------------------------
      term(:ncol) = 1.e-2_r8 * pressure_10m(:ncol) / (r*tv(:ncol))
      species_loop3 : do ispec = gcso4_ndx,gas_pcnst
        if( has_dvel(ispec) ) then
          wrk(:) = 0._r8
          lt_loop: do lt = beglt,endlt
          do i = 1,ncol
            if (fr_lnduse(i,lt)) then
              resc(i) = 1._r8/( 1._r8/rsmx(i,lt,ispec) + 1._r8/rlux(i,lt,ispec) &
                          + 1._r8/(rdc(i) + rclx(i,lt,ispec)) &
                          + 1._r8/(rac(index_season(i,lt),lt) + rgsx(i,lt,ispec)))

              resc(i) = max( 10._r8,resc(i) )

              lnd_frc(i) = lcl_frc_landuse(i,lt)
            endif
          enddo
          !-------------------------------------------------------------------------------------
          ! 	... compute average deposition velocity
          !-------------------------------------------------------------------------------------
          select case( solsym(ispec) )
            case( 'SO2' )
              if( lt == 7 ) then
                where( fr_lnduse(:ncol,lt) )
                  ! assume no surface resistance for SO2 over water`
                  wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt,lchnk) + dep_rb(:ncol,lt,lchnk)) 
                endwhere
              else
                where( fr_lnduse(:ncol,lt) )
                  wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt,lchnk) + dep_rb(:ncol,lt,lchnk) + resc(:))
                endwhere
              end if
            case( 'SO4' )
              where( fr_lnduse(:ncol,lt) )
                wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt,lchnk) + rds(:,lt))
              endwhere
            case( 'gcso4_a' )
              where( fr_lnduse(:ncol,lt) )
                wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt,lchnk) + rds(:,lt))
              endwhere
            case( 'NH4', 'NH4NO3', 'XNH4NO3' )
              where( fr_lnduse(:ncol,lt) )
                wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt,lchnk) + 0.5_r8*rds(:,lt))
              endwhere

            !-------------------------------------------------------------------------------------
            !  ... special case for Pb (for consistency with offline code)
            !-------------------------------------------------------------------------------------
            case( 'Pb' )
              if( lt == 7 ) then
                where( fr_lnduse(:ncol,lt) )
                  wrk(:) = wrk(:) + lnd_frc(:) * 0.05e-2_r8
                endwhere
              else
                where( fr_lnduse(:ncol,lt) )
                  wrk(:ncol) = wrk(:ncol) + lnd_frc(:ncol) * 0.2e-2_r8
                endwhere
              end if

              !-------------------------------------------------------------------------------------
              !  ... special case for carbon aerosols
              !-------------------------------------------------------------------------------------
              case( 'CB1', 'CB2', 'OC1', 'OC2', 'SOAM', 'SOAI', 'SOAT', 'SOAB','SOAX' )
                if ( drydep_method == DD_XLND ) then
                  where( fr_lnduse(:ncol,lt) )
                    wrk(:ncol) = wrk(:ncol) + lnd_frc(:ncol) * 0.10e-2_r8
                  endwhere
                else
                  wrk(:ncol) = 0.10e-2_r8
                endif

              !-------------------------------------------------------------------------------------
              ! deposition over ocean for HCN, CH3CN
              !    velocity estimated from aircraft measurements (E.Apel, INTEX-B)
              !-------------------------------------------------------------------------------------
              case( 'HCN','CH3CN' )
                if( lt == 7 ) then ! over ocean only
                  where( fr_lnduse(:ncol,lt) .and. snow(:ncol) < 0.01_r8  )
                    wrk(:ncol) = wrk(:ncol) + lnd_frc(:ncol) * 0.2e-2_r8
                  endwhere
                end if
              case default
                where( fr_lnduse(:ncol,lt) )
                  wrk(:ncol) = wrk(:ncol) + lnd_frc(:ncol)/(dep_ra(:ncol,lt,lchnk) + dep_rb(:ncol,lt,lchnk) + resc(:ncol))
                endwhere
            end select
          end do lt_loop
          dvel(:ncol,ispec) = wrk(:ncol) * scaling_to_cm_per_s
          dflx(:ncol,ispec) = term(:ncol) * dvel(:ncol,ispec) * mmr(:ncol,plev,ispec)
        end if

      end do species_loop3

      if ( beglt > 1 ) return

      !-------------------------------------------------------------------------------------
      ! 	... special adjustments
      !-------------------------------------------------------------------------------------
      if( mpan_ndx > 0 ) then
        if( has_dvel(mpan_ndx) ) then
          dvel(:ncol,mpan_ndx) = dvel(:ncol,mpan_ndx)/3._r8
          dflx(:ncol,mpan_ndx) = term(:ncol) * dvel(:ncol,mpan_ndx) * mmr(:ncol,plev,mpan_ndx)
        end if
      end if
      if( xmpan_ndx > 0 ) then
        if( has_dvel(xmpan_ndx) ) then
          dvel(:ncol,xmpan_ndx) = dvel(:ncol,xmpan_ndx)/3._r8
          dflx(:ncol,xmpan_ndx) = term(:ncol) * dvel(:ncol,xmpan_ndx) * mmr(:ncol,plev,xmpan_ndx)
        end if
      end if

      ! HCOOH, use CH3COOH dep.vel
      if( hcooh_ndx > 0) then
        if( has_dvel(hcooh_ndx) ) then
          dvel(:ncol,hcooh_ndx) = dvel(:ncol,ch3cooh_ndx)
          dflx(:ncol,hcooh_ndx) = term(:ncol) * dvel(:ncol,hcooh_ndx) * mmr(:ncol,plev,hcooh_ndx)
        end if
      end if
      !
      ! SOG species
      !
      if( sogm_ndx > 0) then
        if( has_dvel(sogm_ndx) ) then
          dvel(:ncol,sogm_ndx) = dvel(:ncol,ch3cooh_ndx)
          dflx(:ncol,sogm_ndx) = term(:ncol) * dvel(:ncol,sogm_ndx) * mmr(:ncol,plev,sogm_ndx)
        end if
      end if
      if( sogi_ndx > 0) then
        if( has_dvel(sogi_ndx) ) then
          dvel(:ncol,sogi_ndx) = dvel(:ncol,ch3cooh_ndx)
          dflx(:ncol,sogi_ndx) = term(:ncol) * dvel(:ncol,sogi_ndx) * mmr(:ncol,plev,sogi_ndx)
        end if
      end if
      if( sogt_ndx > 0) then
        if( has_dvel(sogt_ndx) ) then
          dvel(:ncol,sogt_ndx) = dvel(:ncol,ch3cooh_ndx)
          dflx(:ncol,sogt_ndx) = term(:ncol) * dvel(:ncol,sogt_ndx) * mmr(:ncol,plev,sogt_ndx)
        end if
      end if
      if( sogb_ndx > 0) then
        if( has_dvel(sogb_ndx) ) then
          dvel(:ncol,sogb_ndx) = dvel(:ncol,ch3cooh_ndx)
          dflx(:ncol,sogb_ndx) = term(:ncol) * dvel(:ncol,sogb_ndx) * mmr(:ncol,plev,sogb_ndx)
        end if
      end if
      if( sogx_ndx > 0) then
        if( has_dvel(sogx_ndx) ) then
          dvel(:ncol,sogx_ndx) = dvel(:ncol,ch3cooh_ndx)
          dflx(:ncol,sogx_ndx) = term(:ncol) * dvel(:ncol,sogx_ndx) * mmr(:ncol,plev,sogx_ndx)
        end if
      end if
      !
    end subroutine gc_drydep_xactive



  
    
  !------------------------------------------------------------------------------
  end module mo_gocart
  