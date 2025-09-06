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
    use ioFileMod,    only : getfil
  #ifdef SPMD
    use mpishorthand, only : mpicom, mpir8, mpiint, mpilog
  #endif
    use pio
    use cam_pio_utils,only : cam_pio_openfile
    use cam_logfile,  only : iulog
    use dyn_grid,     only : get_dyn_grid_parm, get_horiz_grid_d
    use iop_data_mod, only:  single_column
  
    use seq_drydep_mod, only : mapping
    use physconst,    only : karman
  
    implicit none
  
    save
  
    interface gc_drydep_inti
       module procedure gc_dvel_inti_table
       module procedure gc_dvel_inti_xactive
       module procedure gc_dvel_inti_fromlnd
    end interface
  
    interface gc_drydep
       module procedure gc_drydep_table
       module procedure gc_drydep_xactive
       module procedure gc_drydep_fromlnd
    end interface
  
    private
    public :: gc_drydep_inti, gc_drydep, set_soilw, chk_soilw, gc_has_drydep
    public :: gc_n_land_type, gc_fraction_landuse, gc_drydep_srf_file
    public :: gc_wetdep,gc_clddiag
  
    public :: gc_wetdep_inputs_t
    public :: gc_wetdep_init
    public :: gc_wetdep_inputs_set
    public :: gc_wetdep_inputs_unset
  
    public :: aerosol_depvel_compute
  
    public :: NIthermo
  
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
    integer :: ndays                          ! # of days in soilw file
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
    !---------------------------------------------------------------------------
    subroutine gc_dvel_inti_fromlnd 
      use mo_chem_utls,         only : get_spc_ndx
      use cam_abortutils,           only : endrun
      use chem_mods,            only : adv_mass
      use seq_drydep_mod,       only : dfoxd
  
      implicit none
  
      integer :: ispc, l
  
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
  
    endsubroutine gc_dvel_inti_fromlnd
  
    !-------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------
    subroutine gc_drydep_fromlnd( ocnfrac, icefrac, ncdate, sfc_temp, pressure_sfc,  &
                               wind_speed, spec_hum, air_temp, pressure_10m, rain, &
                               snow, solar_flux, dvelocity, dflx, mmr, &
                               tv, soilw, rh, ncol, lonndx, latndx, lchnk )
                            
      !-------------------------------------------------------------------------------------
      ! combines the deposition velocities provided by the land model with deposition 
      ! velocities over ocean and sea ice 
      !-------------------------------------------------------------------------------------
  
      use ppgrid,         only : pcols
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
  
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    subroutine gc_dvel_inti_table( depvel_file )
      !---------------------------------------------------------------------------
      !       ... Initialize module, depvel arrays, and a few other variables.
      !           The depvel fields will be linearly interpolated to the correct time
      !---------------------------------------------------------------------------
  
      use mo_constants,  only : d2r, r2d
      use ioFileMod,     only : getfil
      use string_utils,  only : to_lower, GLC
      use mo_chem_utls,  only : get_spc_ndx
      use constituents,  only : pcnst
      use interpolate_data, only : lininterp_init, lininterp, lininterp_finish,interp_type
      use mo_constants,     only : pi
      use phys_grid, only : get_ncols_p, get_rlat_all_p, get_rlon_all_p
  
      implicit none
  
      character(len=*), intent(in) :: depvel_file
  
      !---------------------------------------------------------------------------
      !       ... Local variables
      !---------------------------------------------------------------------------
      integer :: nlat, nlon, nmonth, ndims
      integer :: dimid_lat, dimid_lon, dimid_species, dimid_time
      integer :: dimid(4), count(4), start(4)
      integer :: m, ispecies, nchar, ierr
      real(r8)    :: scale_factor
  
      real(r8), allocatable :: dvel_lats(:), dvel_lons(:)
      real(r8), allocatable :: dvel_in(:,:,:,:)                          ! input depvel array
      character(len=50) :: units
      character(len=20), allocatable :: species_names(:)             ! names of depvel species
      logical :: found
      type(file_desc_t) :: piofile
      type(var_desc_t) :: vid, vid_dvel
  
      character(len=shr_kind_cl) :: locfn
      integer :: mm,n
      integer :: plat, plon
  
      integer :: i, c, ncols
      real(r8) :: to_lats(pcols), to_lons(pcols)
      type(interp_type) :: lon_wgts, lat_wgts
      real(r8), parameter :: zero=0._r8, twopi=2._r8*pi
  
      mm = 1
      do m = 1,pcnst
         if ( len_trim(drydep_list(m))==0 ) exit
         n = get_spc_ndx(drydep_list(m))
         if ( n < 1 ) then
            write(iulog,*) 'gc_drydep_inti: '//drydep_list(m)//' is not included in species set'
            call endrun('gc_drydep_init: invalid dry deposition species')
         endif
      enddo
  
      if( masterproc ) then
         write(iulog,*) 'gc_drydep_inti: following species have dry deposition'
         do i=1,nddvels
            if( len_trim(drydep_list(i)) > 0 ) then
               write(iulog,*) 'gc_drydep_inti: '//trim(drydep_list(i))//' is requested to have dry dep'
            endif
         enddo
         write(iulog,*) 'gc_drydep_inti:'
      endif
  
      if ( nddvels < 1 ) return
  
      plat = get_dyn_grid_parm('plat')
      plon = get_dyn_grid_parm('plon')
  
      !---------------------------------------------------------------------------
      !       ... Setup species maps
      !---------------------------------------------------------------------------
      o3a_ndx   = get_spc_ndx( 'O3A')
      xpan_ndx  = get_spc_ndx( 'XPAN')
      xmpan_ndx = get_spc_ndx( 'XMPAN')
      xno2_ndx  = get_spc_ndx( 'XNO2')
      xhno3_ndx = get_spc_ndx( 'XHNO3')
      xonit_ndx     = get_spc_ndx( 'XONIT')
      xonitr_ndx    = get_spc_ndx( 'XONITR')
      xno_ndx       = get_spc_ndx( 'XNO')
      xho2no2_ndx   = get_spc_ndx( 'XHO2NO2')
      o3a_dd   = gc_has_drydep( 'O3A')
      xpan_dd  = gc_has_drydep( 'XPAN')
      xmpan_dd = gc_has_drydep( 'XMPAN')
      xno2_dd  = gc_has_drydep( 'XNO2')
      xhno3_dd = gc_has_drydep( 'XHNO3')
      xonit_dd     = gc_has_drydep( 'XONIT')
      xonitr_dd    = gc_has_drydep( 'XONITR')
      xno_dd       = gc_has_drydep( 'XNO')
      xho2no2_dd   = gc_has_drydep( 'XHO2NO2')
  
      pan_ndx  = get_spc_ndx( 'PAN')
      mpan_ndx = get_spc_ndx( 'MPAN')
      no2_ndx  = get_spc_ndx( 'NO2')
      hno3_ndx = get_spc_ndx( 'HNO3')
      co_ndx   = get_spc_ndx( 'CO')
      o3_ndx   = get_spc_ndx( 'O3')
      if( o3_ndx < 1 ) then
         o3_ndx = get_spc_ndx( 'OX')
      end if
      h2o2_ndx     = get_spc_ndx( 'H2O2')
      onit_ndx     = get_spc_ndx( 'ONIT')
      onitr_ndx    = get_spc_ndx( 'ONITR')
      ch4_ndx      = get_spc_ndx( 'CH4')
      ch2o_ndx     = get_spc_ndx( 'CH2O')
      ch3ooh_ndx   = get_spc_ndx( 'CH3OOH')
      ch3cho_ndx   = get_spc_ndx( 'CH3CHO')
      ch3cocho_ndx = get_spc_ndx( 'CH3COCHO')
      pooh_ndx     = get_spc_ndx( 'POOH')
      ch3coooh_ndx = get_spc_ndx( 'CH3COOOH')
      c2h5ooh_ndx  = get_spc_ndx( 'C2H5OOH')
      eooh_ndx     = get_spc_ndx( 'EOOH')
      c3h7ooh_ndx  = get_spc_ndx( 'C3H7OOH')
      rooh_ndx     = get_spc_ndx( 'ROOH')
      ch3coch3_ndx = get_spc_ndx( 'CH3COCH3')
      no_ndx       = get_spc_ndx( 'NO')
      ho2no2_ndx   = get_spc_ndx( 'HO2NO2')
      glyald_ndx   = get_spc_ndx( 'GLYALD')
      hyac_ndx     = get_spc_ndx( 'HYAC')
      ch3oh_ndx    = get_spc_ndx( 'CH3OH')
      c2h5oh_ndx   = get_spc_ndx( 'C2H5OH')
      macrooh_ndx  = get_spc_ndx( 'MACROOH')
      isopooh_ndx  = get_spc_ndx( 'ISOPOOH')
      xooh_ndx     = get_spc_ndx( 'XOOH')
      hydrald_ndx  = get_spc_ndx( 'HYDRALD')
      h2_ndx       = get_spc_ndx( 'H2')
      Pb_ndx       = get_spc_ndx( 'Pb')
      o3s_ndx      = get_spc_ndx( 'O3S')
      o3inert_ndx  = get_spc_ndx( 'O3INERT')
      alkooh_ndx  = get_spc_ndx( 'ALKOOH')
      mekooh_ndx  = get_spc_ndx( 'MEKOOH')
      tolooh_ndx  = get_spc_ndx( 'TOLOOH')
      terpooh_ndx = get_spc_ndx( 'TERPOOH')
      ch3cooh_ndx = get_spc_ndx( 'CH3COOH')
      soam_ndx    = get_spc_ndx( 'SOAM' )
      soai_ndx    = get_spc_ndx( 'SOAI' )
      soat_ndx    = get_spc_ndx( 'SOAT' )
      soab_ndx    = get_spc_ndx( 'SOAB' )
      soax_ndx    = get_spc_ndx( 'SOAX' )
      sogm_ndx    = get_spc_ndx( 'SOGM' )
      sogi_ndx    = get_spc_ndx( 'SOGI' )
      sogt_ndx    = get_spc_ndx( 'SOGT' )
      sogb_ndx    = get_spc_ndx( 'SOGB' )
      sogx_ndx    = get_spc_ndx( 'SOGX' )
      soa_ndx     = get_spc_ndx( 'SOA' )
      so4_ndx     = get_spc_ndx( 'SO4' )
      cb1_ndx     = get_spc_ndx( 'CB1' )
      cb2_ndx     = get_spc_ndx( 'CB2' )
      oc1_ndx     = get_spc_ndx( 'OC1' )
      oc2_ndx     = get_spc_ndx( 'OC2' )
      nh3_ndx     = get_spc_ndx( 'NH3' )
      nh4no3_ndx  = get_spc_ndx( 'NH4NO3' )
      xnh4no3_ndx  = get_spc_ndx( 'XNH4NO3' )
      sa1_ndx     = get_spc_ndx( 'SA1' )
      sa2_ndx     = get_spc_ndx( 'SA2' )
      sa3_ndx     = get_spc_ndx( 'SA3' )
      sa4_ndx     = get_spc_ndx( 'SA4' )
      nh4_ndx     = get_spc_ndx( 'NH4' )
      alkooh_dd  = gc_has_drydep( 'ALKOOH')
      mekooh_dd  = gc_has_drydep( 'MEKOOH')
      tolooh_dd  = gc_has_drydep( 'TOLOOH')
      terpooh_dd = gc_has_drydep( 'TERPOOH')
      ch3cooh_dd = gc_has_drydep( 'CH3COOH')
      soam_dd    = gc_has_drydep( 'SOAM' )
      soai_dd    = gc_has_drydep( 'SOAI' )
      soat_dd    = gc_has_drydep( 'SOAT' )
      soab_dd    = gc_has_drydep( 'SOAB' )
      soax_dd    = gc_has_drydep( 'SOAX' )
      sogm_dd    = gc_has_drydep( 'SOGM' )
      sogi_dd    = gc_has_drydep( 'SOGI' )
      sogt_dd    = gc_has_drydep( 'SOGT' )
      sogb_dd    = gc_has_drydep( 'SOGB' )
      sogx_dd    = gc_has_drydep( 'SOGX' )
      soa_dd     = gc_has_drydep( 'SOA' )
      so4_dd     = gc_has_drydep( 'SO4' )
      cb1_dd     = gc_has_drydep( 'CB1' )
      cb2_dd     = gc_has_drydep( 'CB2' )
      oc1_dd     = gc_has_drydep( 'OC1' )
      oc2_dd     = gc_has_drydep( 'OC2' )
      nh3_dd     = gc_has_drydep( 'NH3' )
      nh4no3_dd  = gc_has_drydep( 'NH4NO3' )
      xnh4no3_dd = gc_has_drydep( 'XNH4NO3' )
      sa1_dd     = gc_has_drydep( 'SA1' ) 
      sa2_dd     = gc_has_drydep( 'SA2' )
      sa3_dd     = gc_has_drydep( 'SA3' ) 
      sa4_dd     = gc_has_drydep( 'SA4' )
      nh4_dd     = gc_has_drydep( 'NH4' ) 
      pan_dd  = gc_has_drydep( 'PAN')
      mpan_dd = gc_has_drydep( 'MPAN')
      no2_dd  = gc_has_drydep( 'NO2')
      hno3_dd = gc_has_drydep( 'HNO3')
      co_dd   = gc_has_drydep( 'CO')
      o3_dd   = gc_has_drydep( 'O3')
      if( .not. o3_dd ) then
         o3_dd = gc_has_drydep( 'OX')
      end if
      h2o2_dd     = gc_has_drydep( 'H2O2')
      onit_dd     = gc_has_drydep( 'ONIT')
      onitr_dd    = gc_has_drydep( 'ONITR')
      ch4_dd      = gc_has_drydep( 'CH4')
      ch2o_dd     = gc_has_drydep( 'CH2O')
      ch3ooh_dd   = gc_has_drydep( 'CH3OOH')
      ch3cho_dd   = gc_has_drydep( 'CH3CHO')
      c2h5oh_dd   = gc_has_drydep( 'C2H5OH')
      eooh_dd     = gc_has_drydep( 'EOOH')
      ch3cocho_dd = gc_has_drydep( 'CH3COCHO')
      pooh_dd     = gc_has_drydep( 'POOH')
      ch3coooh_dd = gc_has_drydep( 'CH3COOOH')
      c2h5ooh_dd  = gc_has_drydep( 'C2H5OOH')
      c3h7ooh_dd  = gc_has_drydep( 'C3H7OOH')
      rooh_dd     = gc_has_drydep( 'ROOH')
      ch3coch3_dd = gc_has_drydep( 'CH3COCH3')
      glyald_dd   = gc_has_drydep( 'GLYALD')
      hyac_dd     = gc_has_drydep( 'HYAC')
      ch3oh_dd    = gc_has_drydep( 'CH3OH')
      macrooh_dd  = gc_has_drydep( 'MACROOH')
      isopooh_dd  = gc_has_drydep( 'ISOPOOH')
      xooh_dd     = gc_has_drydep( 'XOOH')
      hydrald_dd  = gc_has_drydep( 'HYDRALD')
      h2_dd       = gc_has_drydep( 'H2')
      Pb_dd       = gc_has_drydep( 'Pb')
      o3s_dd      = gc_has_drydep( 'O3S')
      o3inert_dd  = gc_has_drydep( 'O3INERT')
      ch3cn_dd    = gc_has_drydep( 'CH3CN')
      hcn_dd      = gc_has_drydep( 'HCN')
      hcooh_dd    = gc_has_drydep( 'HCOOH')
      ch3cn_ndx   = get_spc_ndx( 'CH3CN')
      hcn_ndx     = get_spc_ndx( 'HCN')
      hcooh_ndx   = get_spc_ndx( 'HCOOH' )
  
      if( masterproc ) then
         write(iulog,*) 'dvel_inti: diagnostics'
         write(iulog,'(10i5)') pan_ndx, mpan_ndx, no2_ndx, hno3_ndx, o3_ndx, &
              h2o2_ndx, onit_ndx, onitr_ndx, ch4_ndx, ch2o_ndx, &
              ch3ooh_ndx, pooh_ndx, ch3coooh_ndx, c2h5ooh_ndx, eooh_ndx, &
              c3h7ooh_ndx, rooh_ndx, ch3cocho_ndx, co_ndx, ch3coch3_ndx, &
              no_ndx, ho2no2_ndx, glyald_ndx, hyac_ndx, ch3oh_ndx, c2h5oh_ndx, &
              hydrald_ndx, h2_ndx, Pb_ndx, o3s_ndx, o3inert_ndx, macrooh_ndx, &
              xooh_ndx, ch3cho_ndx, isopooh_ndx
         write(iulog,*) pan_dd, mpan_dd, no2_dd, hno3_dd, o3_dd, isopooh_dd, ch4_dd,&
              h2o2_dd, onit_dd, onitr_dd, ch2o_dd, macrooh_dd, xooh_dd, &
              ch3ooh_dd, pooh_dd, ch3coooh_dd, c2h5ooh_dd, eooh_dd, ch3cho_dd, c2h5oh_dd, &
              c3h7ooh_dd, rooh_dd, ch3cocho_dd, co_dd, ch3coch3_dd, &
              glyald_dd, hyac_dd, ch3oh_dd, hydrald_dd, h2_dd, Pb_dd, o3s_dd, o3inert_dd
      endif
      !---------------------------------------------------------------------------
      !       ... Open NetCDF file
      !---------------------------------------------------------------------------
      call getfil (depvel_file, locfn, 0)
      call cam_pio_openfile (piofile, trim(locfn), PIO_NOWRITE)
  
      !---------------------------------------------------------------------------
      !       ... Get variable ID for dep vel array
      !---------------------------------------------------------------------------
      ierr = pio_inq_varid( piofile, 'dvel', vid_dvel )
  
      !---------------------------------------------------------------------------
      !       ... Inquire about dimensions
      !---------------------------------------------------------------------------
      ierr = pio_inq_dimid( piofile, 'lon', dimid_lon )
      ierr = pio_inq_dimlen( piofile, dimid_lon, nlon )
      ierr = pio_inq_dimid( piofile, 'lat', dimid_lat )
      ierr = pio_inq_dimlen( piofile, dimid_lat, nlat )
      ierr = pio_inq_dimid( piofile, 'species', dimid_species )
      ierr = pio_inq_dimlen( piofile, dimid_species, nspecies )
      ierr = pio_inq_dimid( piofile, 'time', dimid_time )
      ierr = pio_inq_dimlen( piofile, dimid_time, nmonth )
      if(masterproc) write(iulog,*) 'dvel_inti: dimensions (nlon,nlat,nspecies,nmonth) = ',nlon,nlat,nspecies,nmonth
  
      !---------------------------------------------------------------------------
      !       ... Check dimensions of dvel variable. Must be (lon, lat, species, month).
      !---------------------------------------------------------------------------
      ierr = pio_inq_varndims( piofile, vid_dvel, ndims )
  
      if( masterproc .and. ndims /= 4 ) then
         write(iulog,*) 'dvel_inti: dvel has ',ndims,' dimensions. Expecting 4.'
         call endrun
      end if
      ierr = pio_inq_vardimid( piofile, vid_dvel, dimid )
  
      if( dimid(1) /= dimid_lon .or. dimid(2) /= dimid_lat .or. &
           dimid(3) /= dimid_species .or. dimid(4) /= dimid_time ) then
         write(iulog,*) 'dvel_inti: Dimensions in wrong order for dvel'
         write(iulog,*) '...      Expecting (lon, lat, species, month)'
         call endrun
      end if
  
      !---------------------------------------------------------------------------
      !       ... Allocate depvel lats, lons and read
      !---------------------------------------------------------------------------
      allocate( dvel_lats(nlat), stat=ierr )
      if( ierr /= 0 ) then
         write(iulog,*) 'dvel_inti: Failed to allocate dvel_lats vector'
         call endrun
      end if
      allocate( dvel_lons(nlon), stat=ierr )
      if( ierr /= 0 ) then
         write(iulog,*) 'dvel_inti: Failed to allocate dvel_lons vector'
         call endrun
      end if
  
      ierr = pio_inq_varid( piofile, 'lat', vid )
      ierr = pio_get_var( piofile, vid, dvel_lats )
      ierr = pio_inq_varid( piofile, 'lon', vid )
      ierr = pio_get_var( piofile, vid, dvel_lons )
  
      !---------------------------------------------------------------------------
      !       ... Set the transform from inputs lats to simulation lats
      !---------------------------------------------------------------------------
      dvel_lats(:nlat) = d2r * dvel_lats(:nlat)
      dvel_lons(:nlon) = d2r * dvel_lons(:nlon)
  
      !---------------------------------------------------------------------------
      !     	... Allocate dvel and read data from file
      !---------------------------------------------------------------------------
      allocate( dvel_in(nlon, nlat ,nspecies, nmonth), stat=ierr )
      if( ierr /= 0 ) then
         write(iulog,*) 'dvel_inti: Failed to allocate dvel_in'
         call endrun
      end if
      start = (/ 1, 1, 1, 1 /)
      count = (/ nlon, nlat, nspecies, nmonth /)
  
      ierr = pio_get_var( piofile, vid_dvel, start, count, dvel_in )
  
  
      !---------------------------------------------------------------------------
      !     	... Check units of deposition velocity. If necessary, convert to cm/s.
      !---------------------------------------------------------------------------
      units(:) = ' '
      ierr = pio_get_att( piofile, vid_dvel, 'units', units )
      if( to_lower(trim(units(:GLC(units)))) == 'm/s' ) then
  #ifdef DEBUG
         if(masterproc)  write(iulog,*) 'dvel_inti: depvel units = m/s. Converting to cm/s'
  #endif
         scale_factor = 100._r8
      elseif( to_lower(trim(units(:GLC(units)))) == 'cm/s' ) then
  #ifdef DEBUG
         if(masterproc)  write(iulog,*) 'dvel_inti: depvel units = cm/s'
  #endif
         scale_factor = 1._r8
      else
  #ifdef DEBUG
         if(masterproc) then
            write(iulog,*) 'dvel_inti: Warning! depvel units unknown = ', to_lower(trim(units)) 
            write(iulog,*) '           ...      proceeding with scale_factor=1'
         end if
  #endif
         scale_factor = 1._r8
      end if
  
      dvel_in(:,:,:,:) = scale_factor*dvel_in(:,:,:,:)
  
      !---------------------------------------------------------------------------
      !     	... Regrid deposition velocities
      !---------------------------------------------------------------------------
      allocate( dvel(pcols,begchunk:endchunk,nspecies,nmonth),stat=ierr )
      if( ierr /= 0 ) then
         write(iulog,*) 'dvel_inti: Failed to allocate dvel'
         call endrun
      end if
  
      do c=begchunk,endchunk
         ncols = get_ncols_p(c)
         call get_rlat_all_p(c, pcols, to_lats)
         call get_rlon_all_p(c, pcols, to_lons)
         call lininterp_init(dvel_lons, nlon, to_lons, ncols, 2, lon_wgts, zero, twopi)
         call lininterp_init(dvel_lats, nlat, to_lats, ncols, 1, lat_wgts)
  
         do ispecies = 1,nspecies
            do m = 1,12
               call lininterp( dvel_in( :,:,ispecies,m ), nlon, nlat, dvel(:,c,ispecies,m), ncols,lon_wgts,lat_wgts)
            end do
         end do
  
         call lininterp_finish(lat_wgts)
         call lininterp_finish(lon_wgts)
      end do
  
      deallocate( dvel_in )
      deallocate( dvel_lats, dvel_lons )
  
      !---------------------------------------------------------------------------
      !     	... Read in species names and determine mapping to tracer numbers
      !---------------------------------------------------------------------------
      allocate( species_names(nspecies), stat=ierr )
      if( ierr /= 0 ) then
         write(iulog,*) 'dvel_inti: species_names allocation error = ',ierr
         call endrun
      end if
      ierr = pio_inq_varid( piofile, 'species_name', vid )
      ierr = pio_inq_varndims( piofile, vid, ndims )
  
      ierr = pio_inq_vardimid( piofile, vid, dimid )
  
      ierr = pio_inq_dimlen( piofile, dimid(1), nchar )
      map(:) = 0
      do ispecies = 1,nspecies
         start(:2) = (/ 1, ispecies /)
         count(:2) = (/ nchar, 1 /)
         species_names(ispecies)(:) = ' '
         ierr = pio_get_var( piofile, vid, start(1:2), count(1:2), species_names(ispecies:ispecies) )
         if( species_names(ispecies) == 'O3' ) then
            o3_in_tab  = .true.
            o3_tab_ndx = ispecies
         else if( species_names(ispecies) == 'H2O2' ) then
            h2o2_in_tab  = .true.
            h2o2_tab_ndx = ispecies
         else if( species_names(ispecies) == 'CH3OOH' ) then
            ch3ooh_in_tab  = .true.
            ch3ooh_tab_ndx = ispecies
         else if( species_names(ispecies) == 'CO' ) then
            co_in_tab  = .true.
            co_tab_ndx = ispecies
         else if( species_names(ispecies) == 'CH3CHO' ) then
            ch3cho_in_tab  = .true.
            ch3cho_tab_ndx = ispecies
         end if
         found = .false.
         do m = gcso4_ndx,gas_pcnst
            if( species_names(ispecies) == solsym(m) .or. &
                 (species_names(ispecies) == 'O3' .and. solsym(m) == 'OX') .or. &
                 (species_names(ispecies) == 'HNO4' .and. solsym(m) == 'HO2NO2') ) then
               if ( gc_has_drydep( solsym(m) ) ) then
                  map(m) = ispecies
                  found = .true.
  #ifdef DEBUG
                  if( masterproc ) then
                     write(iulog,*) 'dvel_inti: ispecies, m, tracnam = ',ispecies,m,trim(solsym(m))
                  end if
  #endif
                  exit
               end if
            end if
         end do
         if( .not. found ) then
            write(iulog,*) 'dvel_inti: Warning! DVEL species ',trim(species_names(ispecies)),' not found'
         endif
      end do
      deallocate( species_names )
  
      call pio_closefile( piofile )
  
      !---------------------------------------------------------------------------
      !     	... Allocate dvel_interp array
      !---------------------------------------------------------------------------
      allocate( dvel_interp(pcols,begchunk:endchunk,nspecies),stat=ierr )
      if( ierr /= 0 ) then
         write(iulog,*) 'dvel_inti: Failed to allocate dvel_interp; error = ',ierr
         call endrun
      end if
  
    end subroutine gc_dvel_inti_table
  
    !-------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------
    subroutine interpdvel( calday, ncol, lchnk )
      !---------------------------------------------------------------------------
      ! 	... Interpolate the fields whose values are required at the
      !           begining of a timestep.
      !---------------------------------------------------------------------------
  
      use time_manager,  only : get_calday
  
      implicit none
  
      !---------------------------------------------------------------------------
      ! 	... Dummy arguments
      !---------------------------------------------------------------------------
      real(r8), intent(in) :: calday   ! Interpolate the input data to calday
      integer, intent(in) :: ncol, lchnk
  
      !---------------------------------------------------------------------------
      ! 	... Local variables
      !---------------------------------------------------------------------------
      integer :: m, last, next
      integer  ::  dates(12) = (/ 116, 214, 316, 415,  516,  615, &
                                  716, 816, 915, 1016, 1115, 1216 /)
      real(r8) :: calday_loc, last_days, next_days
      real(r8), save ::  dys(12)
      logical, save  ::  entered = .false.
  
      if( .not. entered ) then
         do m = 1,12
            dys(m) = get_calday( dates(m), 0 )
         end do
         entered = .true.
      end if
  
      if( calday < dys(1) ) then
         next = 1
         last = 12
      else if( calday >= dys(12) ) then
         next = 1
         last = 12
      else
         do m = 11,1,-1
            if( calday >= dys(m) ) then
               exit
            end if
         end do
         last = m
         next = m + 1
      end if
  
      last_days  = dys( last )
      next_days  = dys( next )
      calday_loc = calday
  
      if( next_days < last_days ) then
         next_days = next_days + 365._r8
      end if
      if( calday_loc < last_days ) then
         calday_loc = calday_loc + 365._r8
      end if
  
      do m = 1,nspecies
         call intp2d( last_days, next_days, calday_loc, ncol, lchnk, &
                      dvel(:,lchnk,m,last), &
                      dvel(:,lchnk,m,next), &
                      dvel_interp(:,lchnk,m) )
      end do
  
    end subroutine interpdvel
  
    !-------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------
    subroutine intp2d( t1, t2, tint, ncol, lchnk, f1, f2, fint )
      !-----------------------------------------------------------------------
      ! 	... Linearly interpolate between f1(t1) and f2(t2) to fint(tint).
      !-----------------------------------------------------------------------
  
      implicit none
  
      !-----------------------------------------------------------------------
      ! 	... Dummy arguments
      !-----------------------------------------------------------------------
      real(r8), intent(in) :: &
           t1, &            ! time level of f1
           t2, &            ! time level of f2
           tint             ! interpolant time
      real(r8), dimension(pcols), intent(in) :: &
           f1, &            ! field at time t1
           f2               ! field at time t2
  
      integer, intent(in) :: ncol, lchnk
  
      real(r8), intent(out) :: &
           fint(pcols) ! field at time tint
  
  
      !-----------------------------------------------------------------------
      ! 	... Local variables
      !-----------------------------------------------------------------------
      integer  :: j, plat
      real(r8) :: factor
  
      factor = (tint - t1)/(t2 - t1)
  
      fint(:ncol) = f1(:ncol) + (f2(:ncol) - f1(:ncol))*factor
  
    end subroutine intp2d
  
    !-------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------
    subroutine gc_drydep_table( calday, tsurf, zen_angle, &
                             depvel, dflx, q, p, &
                             tv, ncol, icefrac, ocnfrac, lchnk )
      !--------------------------------------------------------
      !       ... Form the deposition velocities for this
      !           latitude slice
      !--------------------------------------------------------
  
      use physconst,     only : rair,pi
      use dycore,        only : dycore_is
  
      implicit none
  
      !--------------------------------------------------------
      !       ... Dummy arguments
      !--------------------------------------------------------
      integer, intent(in)     ::   ncol                    ! columns in chunk
      real(r8), intent(in)    ::  q(pcols,plev,gas_pcnst) ! tracer mmr (kg/kg)
      real(r8), intent(in)    ::  p(pcols)            ! midpoint pressure in surface layer (Pa)
      real(r8), intent(in)    ::  tv(pcols)           ! virtual temperature in surface layer (K)
      real(r8), intent(in)    ::  calday              ! time of year in days
      real(r8), intent(in)    ::  tsurf(pcols)        ! surface temperature (K)
      real(r8), intent(in)    ::  zen_angle(ncol)    ! zenith angle (radians)
      real(r8), intent(inout) ::  dflx(pcols,gas_pcnst)   ! flux due to dry deposition (kg/m^2/sec)
      real(r8), intent(out)   ::  depvel(ncol,gas_pcnst) ! deposition vel (cm/s)
  
      real(r8), intent(in) :: icefrac(pcols)          ! sea-ice areal fraction
      real(r8), intent(in) :: ocnfrac(pcols)          ! ocean areal fraction
      
      integer, intent(in)     :: lchnk
      !-----------------------------------------------------------------------
      ! 	... Local variables
      !-----------------------------------------------------------------------
      integer :: m, spc_ndx, tmp_ndx, i
      real(r8), dimension(ncol) :: vel, glace, temp_fac, wrk, tmp
      real(r8), dimension(ncol) :: o3_tab_dvel
      real(r8), dimension(ncol) :: ocean 
  
      real(r8), parameter :: pid2 = .5_r8 * pi
  
      if(dycore_is('UNSTRUCTURED')) then
         call endrun( 'Option not supported for unstructured atmosphere grids ')
      end if
  
      !-----------------------------------------------------------------------
      !       ... Note the factor 1.e-2 in the wrk array calculation is
      !           to transform the incoming dep vel from cm/s to m/s
      !-----------------------------------------------------------------------
      wrk(:ncol) =  1.e-2_r8 * p(:ncol) / (rair * tv(:ncol))
  
      !--------------------------------------------------------
      !       ... Initialize all deposition velocities to zero
      !--------------------------------------------------------
      depvel(:,:) = 0._r8
  
      !--------------------------------------------------------
      !       ... Time interpolate primary depvel array
      !           (also seaice and npp)
      !--------------------------------------------------------
      call interpdvel( calday, ncol, lchnk )
  
      if( o3_in_tab ) then
         do i=1,ncol
            o3_tab_dvel(i) = dvel_interp(i,lchnk,o3_tab_ndx)
         enddo
      end if
  
      !--------------------------------------------------------
      !       ... Set deposition velocities
      !--------------------------------------------------------
      do m = gcso4_ndx,gas_pcnst
         if( map(m) /= 0 ) then
            do i = 1,ncol
               depvel(i,m) = dvel_interp(i,lchnk,map(m))
               dflx(i,m)   = wrk(i) * depvel(i,m) * q(i,plev,m)
            enddo
         end if
      end do
  
      !--------------------------------------------------------
      !       ... Set some variables needed for some dvel calculations
      !--------------------------------------------------------
      temp_fac(:ncol)   = min( 1._r8, max( 0._r8, (tsurf(:ncol) - 268._r8) / 5._r8 ) )
      ocean(:ncol)  = icefrac(:ncol)+ocnfrac(:ncol)
      glace(:ncol)  = icefrac(:ncol) + (1._r8 - ocean(:ncol)) * (1._r8 - temp_fac(:ncol))
      glace(:ncol)  = min( 1._r8,glace(:ncol) )
  
      !--------------------------------------------------------
      !       ... Set pan & mpan
      !--------------------------------------------------------
      if( o3_in_tab ) then
         tmp(:ncol) = o3_tab_dvel(:ncol) / 3._r8
      else
         tmp(:) = 0._r8
      end if
      if( pan_dd ) then
         if( map(pan_ndx) == 0 ) then
            depvel(:ncol,pan_ndx) = tmp(:ncol)
            dflx(:ncol,pan_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,pan_ndx)
         end if
      end if
      if( mpan_dd ) then
         if( map(mpan_ndx) == 0 ) then
            depvel(:ncol,mpan_ndx) = tmp(:ncol)
            dflx(:ncol,mpan_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,mpan_ndx)
         end if
      end if
  
      !--------------------------------------------------------
      !       ... Set no2 dvel
      !--------------------------------------------------------
      if( no2_dd ) then
         if( map(no2_ndx) == 0 .and. o3_in_tab ) then
            depvel(:ncol,no2_ndx) = (.6_r8*o3_tab_dvel(:ncol) + .055_r8*ocean(:ncol)) * .9_r8
            dflx(:ncol,no2_ndx)   = wrk(:) * depvel(:ncol,no2_ndx) * q(:ncol,plev,no2_ndx)
         end if
      end if
  
      !--------------------------------------------------------
      !       ... Set hno3 dvel
      !--------------------------------------------------------
      tmp(:ncol) = (2._r8 - ocnfrac(:ncol)) * (1._r8 - glace(:ncol)) + .05_r8 * glace(:ncol)
      if( hno3_dd ) then
         if( map(hno3_ndx) == 0 ) then
            depvel(:ncol,hno3_ndx) = tmp(:ncol)
            dflx(:ncol,hno3_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,hno3_ndx)
         else
            tmp(:ncol) = depvel(:ncol,hno3_ndx)
         end if
      end if
      if( onitr_dd ) then
         if( map(onitr_ndx) == 0 ) then
            depvel(:ncol,onitr_ndx) = tmp(:ncol)
            dflx(:ncol,onitr_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,onitr_ndx)
         end if
      end if
      if( isopooh_dd ) then
         if( map(isopooh_ndx) == 0 ) then
            depvel(:ncol,isopooh_ndx) = tmp(:ncol)
            dflx(:ncol,isopooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,isopooh_ndx)
         end if
      end if
  
      !--------------------------------------------------------
      !       ... Set h2o2 dvel
      !--------------------------------------------------------
      if( .not. h2o2_in_tab ) then
         if( o3_in_tab ) then
            tmp(:ncol) = .05_r8*glace(:ncol) + ocean(:ncol) - icefrac(:ncol) &
                 + (1._r8 - (glace(:) + ocean(:ncol)) + icefrac(:ncol)) &
                 *max( 1._r8,1._r8/(.5_r8 + 1._r8/(6._r8*o3_tab_dvel(:ncol))) )
         else
            tmp(:ncol) = 0._r8
         end if
      else
         do i=1,ncol
            tmp(i) = dvel_interp(i,lchnk,h2o2_tab_ndx)
         enddo
      end if
      if( h2o2_dd ) then
         if( map(h2o2_ndx) == 0 ) then
            depvel(:ncol,h2o2_ndx) = tmp(:ncol)
            dflx(:ncol,h2o2_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,h2o2_ndx)
         end if
      end if
      !--------------------------------------------------------
      !       ... Set hcn dvel
      !--------------------------------------------------------
      if( hcn_dd ) then
         if( map(hcn_ndx) == 0 ) then
            depvel(:ncol,hcn_ndx) = ocnfrac(:ncol)*0.2_r8
         endif
      endif
      !--------------------------------------------------------
      !       ... Set ch3cn dvel
      !--------------------------------------------------------
      if( ch3cn_dd ) then
         if( map(ch3cn_ndx) == 0 ) then
            depvel(:,ch3cn_ndx) = ocnfrac(:ncol)*0.2_r8
         endif
      endif
      !--------------------------------------------------------
      !       ... Set onit
      !--------------------------------------------------------
      if( onit_dd ) then
         if( map(onit_ndx) == 0 ) then
            depvel(:ncol,onit_ndx) = tmp(:ncol)
            dflx(:ncol,onit_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,onit_ndx)
         end if
      end if
      if( ch3cocho_dd ) then
         if( map(ch3cocho_ndx) == 0 ) then
            depvel(:ncol,ch3cocho_ndx) = tmp(:ncol)
            dflx(:ncol,ch3cocho_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,ch3cocho_ndx)
         end if
      end if
      if( ch3ooh_in_tab ) then
         do i=1,ncol
            tmp(i) = dvel_interp(i,lchnk,ch3ooh_tab_ndx)
         enddo
      else
         tmp(:ncol) = .5_r8 * tmp(:ncol)
      end if
      if( ch3ooh_dd ) then
         if( map(ch3ooh_ndx) == 0 ) then
            depvel(:ncol,ch3ooh_ndx) = tmp(:ncol)
            dflx(:ncol,ch3ooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,ch3ooh_ndx)
         end if
      end if
      if( pooh_dd ) then
         if( map(pooh_ndx) == 0 ) then
            depvel(:ncol,pooh_ndx) = tmp(:ncol)
            dflx(:ncol,pooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,pooh_ndx)
         end if
      end if
      if( ch3coooh_dd ) then
         if( map(ch3coooh_ndx) == 0 ) then
            depvel(:ncol,ch3coooh_ndx) = tmp(:ncol)
            dflx(:ncol,ch3coooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,ch3coooh_ndx)
         end if
      end if
      if( c2h5ooh_dd ) then
         if( map(c2h5ooh_ndx) == 0 ) then
            depvel(:ncol,c2h5ooh_ndx) = tmp(:ncol)
            dflx(:ncol,c2h5ooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,c2h5ooh_ndx)
         end if
      end if
      if( c3h7ooh_dd ) then
         if( map(c3h7ooh_ndx) == 0 ) then
            depvel(:ncol,c3h7ooh_ndx) = tmp(:ncol)
            dflx(:ncol,c3h7ooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,c3h7ooh_ndx)
         end if
      end if
      if( rooh_dd ) then
         if( map(rooh_ndx) == 0 ) then
            depvel(:ncol,rooh_ndx) = tmp(:ncol)
            dflx(:ncol,rooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,rooh_ndx)
         end if
      end if
      if( macrooh_dd ) then
         if( map(macrooh_ndx) == 0 ) then
            depvel(:ncol,macrooh_ndx) = tmp(:ncol)
            dflx(:ncol,macrooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,macrooh_ndx)
         end if
      end if
      if( xooh_dd ) then
         if( map(xooh_ndx) == 0 ) then
            depvel(:ncol,xooh_ndx) = tmp(:ncol)
            dflx(:ncol,xooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,xooh_ndx)
         end if
      end if
      if( ch3oh_dd ) then
         if( map(ch3oh_ndx) == 0 ) then
            depvel(:ncol,ch3oh_ndx) = tmp(:ncol)
            dflx(:ncol,ch3oh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,ch3oh_ndx)
         end if
      end if
      if( c2h5oh_dd ) then
         if( map(c2h5oh_ndx) == 0 ) then
            depvel(:ncol,c2h5oh_ndx) = tmp(:ncol)
            dflx(:ncol,c2h5oh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,c2h5oh_ndx)
         end if
      end if
      if( alkooh_dd ) then
         if( map(alkooh_ndx) == 0 ) then
            depvel(:ncol,alkooh_ndx) = tmp(:ncol)
            dflx(:ncol,alkooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,alkooh_ndx)
         end if
      end if
      if( mekooh_dd ) then
         if( map(mekooh_ndx) == 0 ) then
            depvel(:ncol,mekooh_ndx) = tmp(:ncol)
            dflx(:ncol,mekooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,mekooh_ndx)
         end if
      end if
      if( tolooh_dd ) then
         if( map(tolooh_ndx) == 0 ) then
            depvel(:ncol,tolooh_ndx) = tmp(:ncol)
            dflx(:ncol,tolooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,tolooh_ndx)
         end if
      end if
      if( terpooh_dd ) then
         if( map(terpooh_ndx) == 0 ) then
            depvel(:ncol,terpooh_ndx) = tmp(:ncol)
            dflx(:ncol,terpooh_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,terpooh_ndx)
         end if
      end if
  
      if( o3_in_tab ) then
         tmp(:ncol) = o3_tab_dvel(:ncol)
      else
         tmp(:ncol) = 0._r8
      end if
      if( ch2o_dd ) then
         if( map(ch2o_ndx) == 0 ) then
            depvel(:ncol,ch2o_ndx) = tmp(:ncol)
            dflx(:ncol,ch2o_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,ch2o_ndx)
         end if
      end if
  
      if( hydrald_dd ) then
         if( map(hydrald_ndx) == 0 ) then
            depvel(:ncol,hydrald_ndx) = tmp(:ncol)
            dflx(:ncol,hydrald_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,hydrald_ndx)
         end if
      end if
      if( ch3cooh_dd  ) then
         if( map(ch3cooh_ndx) == 0 ) then
            depvel(:ncol,ch3cooh_ndx) = depvel(:ncol,ch2o_ndx)
            dflx(:ncol,ch3cooh_ndx) = wrk(:ncol) * depvel(:ncol,ch3cooh_ndx) * q(:ncol,plev,ch3cooh_ndx)
         end if
      end if
      if( eooh_dd ) then
         if( map(eooh_ndx) == 0 ) then
            depvel(:ncol,eooh_ndx) = depvel(:ncol,ch2o_ndx)
            dflx(:ncol,eooh_ndx) = wrk(:ncol) * depvel(:ncol,eooh_ndx) * q(:ncol,plev,eooh_ndx)
         end if
      end if
      ! HCOOH - set to CH3COOH
      if( hcooh_dd  ) then
         if( map(hcooh_ndx) == 0 ) then
            depvel(:ncol,hcooh_ndx) = depvel(:ncol,ch2o_ndx)
            dflx(:ncol,hcooh_ndx) = wrk(:ncol) * depvel(:ncol,hcooh_ndx) * q(:ncol,plev,hcooh_ndx)
         end if
      end if
  
      !--------------------------------------------------------
      !       ... Set co and related species dep vel
      !--------------------------------------------------------
      if( co_in_tab ) then
         do i=1,ncol
            tmp(i) = dvel_interp(i,lchnk,co_tab_ndx)
         enddo
      else
         tmp(:) = 0._r8
      end if
      if( co_dd ) then
         if( map(co_ndx) == 0 ) then
            depvel(:ncol,co_ndx) = tmp(:ncol)
            dflx(:ncol,co_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,co_ndx)
         end if
      end if
      if( ch3coch3_dd ) then
         if( map(ch3coch3_ndx) == 0 ) then
            depvel(:ncol,ch3coch3_ndx) = tmp(:ncol)
            dflx(:ncol,ch3coch3_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,ch3coch3_ndx)
         end if
      end if
      if( hyac_dd ) then
         if( map(hyac_ndx) == 0 ) then
            depvel(:ncol,hyac_ndx) = tmp(:ncol)
            dflx(:ncol,hyac_ndx)   = wrk(:ncol) * tmp(:ncol) * q(:ncol,plev,hyac_ndx)
         end if
      end if
      if( h2_dd ) then
         if( map(h2_ndx) == 0 ) then
            depvel(:ncol,h2_ndx) = tmp(:ncol) * 1.5_r8                ! Hough(1991)
            dflx(:ncol,h2_ndx)   = wrk(:ncol) * depvel(:ncol,h2_ndx) * q(:ncol,plev,h2_ndx)
         end if
      end if
  
      !--------------------------------------------------------
      !       ... Set glyald
      !--------------------------------------------------------
      if( glyald_dd ) then
         if( map(glyald_ndx) == 0 ) then
            if( ch3cho_dd ) then
               depvel(:ncol,glyald_ndx) = depvel(:ncol,ch3cho_ndx)
            else if( ch3cho_in_tab ) then
               do i=1,ncol
                  depvel(i,glyald_ndx) = dvel_interp(i,lchnk,ch3cho_tab_ndx)
               enddo
            else
               depvel(:ncol,glyald_ndx) = 0._r8
            end if
            dflx(:ncol,glyald_ndx)   = wrk(:ncol) * depvel(:ncol,glyald_ndx) * q(:ncol,plev,glyald_ndx)
         end if
      end if
  
      !--------------------------------------------------------
      !       ... Lead deposition
      !--------------------------------------------------------
      if( Pb_dd ) then
         if( map(Pb_ndx) == 0 ) then
            depvel(:ncol,Pb_ndx) = ocean(:ncol)  * .05_r8 + (1._r8 - ocean(:ncol)) * .2_r8
            dflx(:ncol,Pb_ndx)   = wrk(:ncol) * depvel(:ncol,Pb_ndx) * q(:ncol,plev,Pb_ndx)
         end if
      end if
  
      !--------------------------------------------------------
      !       ... diurnal dependence for OX dvel
      !--------------------------------------------------------
      if( o3_dd .or. o3s_dd .or. o3inert_dd ) then
         if( o3_dd .or. o3_in_tab ) then
            if( o3_dd ) then
               tmp(:ncol) = max( 1._r8,sqrt( (depvel(:ncol,o3_ndx) - .2_r8)**3/.27_r8 + 4._r8*depvel(:ncol,o3_ndx) + .67_r8 ) )
               vel(:ncol) = depvel(:ncol,o3_ndx)
            else if( o3_in_tab ) then
               tmp(:ncol) = max( 1._r8,sqrt( (o3_tab_dvel(:ncol) - .2_r8)**3/.27_r8 + 4._r8*o3_tab_dvel(:ncol) + .67_r8 ) )
               vel(:ncol) = o3_tab_dvel(:ncol)
            end if
            where( abs( zen_angle(:) ) > pid2 )
               vel(:) = vel(:) / tmp(:)
            elsewhere
               vel(:) = vel(:) * tmp(:)
            endwhere
  
         else
            vel(:ncol) = 0._r8
         end if
         if( o3_dd ) then
            depvel(:ncol,o3_ndx) = vel(:ncol)
            dflx(:ncol,o3_ndx)   = wrk(:ncol) * vel(:ncol) * q(:ncol,plev,o3_ndx)
         end if
         !--------------------------------------------------------
         !       ... Set stratospheric O3 deposition
         !--------------------------------------------------------
         if( o3s_dd ) then
            depvel(:ncol,o3s_ndx) = vel(:ncol)
            dflx(:ncol,o3s_ndx)   = wrk(:ncol) * vel(:ncol) * q(:ncol,plev,o3s_ndx)
         end if
         if( o3inert_dd ) then
            depvel(:ncol,o3inert_ndx) = vel(:ncol)
            dflx(:ncol,o3inert_ndx)   = wrk(:ncol) * vel(:ncol) * q(:ncol,plev,o3inert_ndx)
         end if
      end if
  
      if( xno2_dd ) then 
         if( map(xno2_ndx) == 0 ) then
            depvel(:ncol,xno2_ndx) = depvel(:ncol,no2_ndx)
            dflx(:ncol,xno2_ndx)   = wrk(:ncol) * depvel(:ncol,xno2_ndx) * q(:ncol,plev,xno2_ndx)
         end if
      endif
      if( o3a_dd ) then 
         if( map(o3a_ndx) == 0 ) then
            depvel(:ncol,o3a_ndx) = depvel(:ncol,o3_ndx)
            dflx(:ncol,o3a_ndx)   = wrk(:ncol) * depvel(:ncol,o3a_ndx) * q(:ncol,plev,o3a_ndx)
         end if
      endif
      if( xhno3_dd ) then 
         if( map(xhno3_ndx) == 0 ) then
            depvel(:ncol,xhno3_ndx) = depvel(:ncol,hno3_ndx)
            dflx(:ncol,xhno3_ndx)   = wrk(:ncol) * depvel(:ncol,xhno3_ndx) * q(:ncol,plev,xhno3_ndx)
         end if
      endif
      if( xnh4no3_dd ) then 
         if( map(xnh4no3_ndx) == 0 ) then
            depvel(:ncol,xnh4no3_ndx) = depvel(:ncol,nh4no3_ndx)
            dflx(:ncol,xnh4no3_ndx)   = wrk(:ncol) * depvel(:ncol,xnh4no3_ndx) * q(:ncol,plev,xnh4no3_ndx)
         end if
      endif
      if( xpan_dd ) then 
         if( map(xpan_ndx) == 0 ) then
            depvel(:ncol,xpan_ndx) = depvel(:ncol,pan_ndx)
            dflx(:ncol,xpan_ndx)   = wrk(:ncol) * depvel(:ncol,xpan_ndx) * q(:ncol,plev,xpan_ndx)
         end if
      endif
      if( xmpan_dd ) then 
         if( map(xmpan_ndx) == 0 ) then
            depvel(:ncol,xmpan_ndx) = depvel(:ncol,mpan_ndx)
            dflx(:ncol,xmpan_ndx)   = wrk(:ncol) * depvel(:ncol,xmpan_ndx) * q(:ncol,plev,xmpan_ndx)
         end if
      endif
      if( xonit_dd ) then 
         if( map(xonit_ndx) == 0 ) then
            depvel(:ncol,xonit_ndx) = depvel(:ncol,onit_ndx)
            dflx(:ncol,xonit_ndx)   = wrk(:ncol) * depvel(:ncol,xonit_ndx) * q(:ncol,plev,xonit_ndx)
         end if
      endif
      if( xonitr_dd ) then 
         if( map(xonitr_ndx) == 0 ) then
            depvel(:ncol,xonitr_ndx) = depvel(:ncol,onitr_ndx)
            dflx(:ncol,xonitr_ndx)   = wrk(:ncol) * depvel(:ncol,xonitr_ndx) * q(:ncol,plev,xonitr_ndx)
         end if
      endif
      if( xno_dd ) then 
         if( map(xno_ndx) == 0 ) then
            depvel(:ncol,xno_ndx) = depvel(:ncol,no_ndx)
            dflx(:ncol,xno_ndx)   = wrk(:ncol) * depvel(:ncol,xno_ndx) * q(:ncol,plev,xno_ndx)
         end if
      endif
      if( xho2no2_dd ) then 
         if( map(xho2no2_ndx) == 0 ) then
            depvel(:ncol,xho2no2_ndx) = depvel(:ncol,ho2no2_ndx)
            dflx(:ncol,xho2no2_ndx)   = wrk(:ncol) * depvel(:ncol,xho2no2_ndx) * q(:ncol,plev,xho2no2_ndx)
         end if
      endif
  
    end subroutine gc_drydep_table
  
    !-------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------
    subroutine gc_dvel_inti_xactive( depvel_lnd_file, clim_soilw_file, season_wes_file )
      !-------------------------------------------------------------------------------------
      ! 	... intialize interactive drydep
      !-------------------------------------------------------------------------------------
      use dycore,        only : dycore_is
      use mo_constants,  only : r2d
      use chem_mods,     only : adv_mass
      use mo_chem_utls,  only : get_spc_ndx
      use seq_drydep_mod,only : drydep_method, DD_XATM, DD_XLND
      use phys_control,  only : phys_getopts
  
      implicit none
  
      !-------------------------------------------------------------------------------------
      ! 	... dummy arguments
      !-------------------------------------------------------------------------------------
      character(len=*), intent(in) :: depvel_lnd_file, clim_soilw_file, season_wes_file 
  
      !-------------------------------------------------------------------------------------
      ! 	... local variables
      !-------------------------------------------------------------------------------------
      integer :: i, j, ii, jj, jl, ju
      integer :: nlon_veg, nlat_veg, npft_veg
      integer :: nlat_lai, npft_lai, pos_min, imin
      integer :: dimid
      integer :: m, n, l, id
      integer :: length1, astat
      integer, allocatable :: wk_lai(:,:,:)
      integer, allocatable :: index_season_lai_j(:,:)
      integer :: k, num_max, k_max
      integer :: num_seas(5)
      integer :: plon, plat
      integer :: ierr
  
      real(r8)              :: spc_mass
      real(r8)              :: diff_min, target_lat
      real(r8), allocatable :: vegetation_map(:,:,:)
      real(r8), pointer     :: soilw_map(:,:,:)
      real(r8), allocatable :: work(:,:)
      real(r8), allocatable :: landmask(:,:)
      real(r8), allocatable :: urban(:,:)
      real(r8), allocatable :: lake(:,:)
      real(r8), allocatable :: wetland(:,:)
      real(r8), allocatable :: lon_veg(:)
      real(r8), allocatable :: lon_veg_edge(:)
      real(r8), allocatable :: lat_veg(:)
      real(r8), allocatable :: lat_veg_edge(:)
      real(r8), allocatable :: lat_lai(:)
      real(r8), allocatable :: clat(:)
      character(len=32) :: test_name
      type(file_desc_t) :: piofile
      type(var_desc_t) :: vid
      logical :: do_soilw
  
      character(len=shr_kind_cl) :: locfn
      logical :: prog_modal_aero
  
      ! determine if modal aerosols are active so that gc_fraction_landuse array is initialized for modal aerosal dry dep
      call phys_getopts(prog_modal_aero_out=prog_modal_aero)
  
      call gc_dvel_inti_fromlnd()
  
      if( masterproc ) then
         write(iulog,*) 'gc_drydep_inti: following species have dry deposition'
         do i=1,nddvels
            if( len_trim(drydep_list(i)) > 0 ) then
               write(iulog,*) 'gc_drydep_inti: '//trim(drydep_list(i))//' is requested to have dry dep'
            endif
         enddo
         write(iulog,*) 'gc_drydep_inti:'
      endif
  
      !-------------------------------------------------------------------------------------
      ! 	... get species indices
      !-------------------------------------------------------------------------------------
      gcso4_ndx     = get_spc_ndx( 'gcso4_a' )
      xpan_ndx      = get_spc_ndx( 'XPAN' )
      xmpan_ndx     = get_spc_ndx( 'XMPAN' )
      o3a_ndx       = get_spc_ndx( 'O3A' )
  
      ch4_ndx      = get_spc_ndx( 'CH4' )
      h2_ndx       = get_spc_ndx( 'H2' )
      co_ndx       = get_spc_ndx( 'CO' )
      Pb_ndx       = get_spc_ndx( 'Pb' )
      pan_ndx      = get_spc_ndx( 'PAN' )
      mpan_ndx     = get_spc_ndx( 'MPAN' )
      o3_ndx       = get_spc_ndx( 'OX' )
      if( o3_ndx < 0 ) then
         o3_ndx  = get_spc_ndx( 'O3' )
      end if
      so2_ndx     = get_spc_ndx( 'SO2' )
      alkooh_ndx  = get_spc_ndx( 'ALKOOH')
      mekooh_ndx  = get_spc_ndx( 'MEKOOH')
      tolooh_ndx  = get_spc_ndx( 'TOLOOH')
      terpooh_ndx = get_spc_ndx( 'TERPOOH')
      ch3cooh_ndx = get_spc_ndx( 'CH3COOH')
      soa_ndx     = get_spc_ndx( 'SOA' )
      so4_ndx     = get_spc_ndx( 'SO4' )
      cb1_ndx     = get_spc_ndx( 'CB1' )
      cb2_ndx     = get_spc_ndx( 'CB2' )
      oc1_ndx     = get_spc_ndx( 'OC1' )
      oc2_ndx     = get_spc_ndx( 'OC2' )
      nh3_ndx     = get_spc_ndx( 'NH3' )
      nh4no3_ndx  = get_spc_ndx( 'NH4NO3' )
      sa1_ndx     = get_spc_ndx( 'SA1' )
      sa2_ndx     = get_spc_ndx( 'SA2' )
      sa3_ndx     = get_spc_ndx( 'SA3' )
      sa4_ndx     = get_spc_ndx( 'SA4' )
      nh4_ndx     = get_spc_ndx( 'NH4' )
      alkooh_dd  = gc_has_drydep( 'ALKOOH')
      mekooh_dd  = gc_has_drydep( 'MEKOOH')
      tolooh_dd  = gc_has_drydep( 'TOLOOH')
      terpooh_dd = gc_has_drydep( 'TERPOOH')
      ch3cooh_dd = gc_has_drydep( 'CH3COOH')
      soa_dd     = gc_has_drydep( 'SOA' )
      so4_dd     = gc_has_drydep( 'SO4' )
      cb1_dd     = gc_has_drydep( 'CB1' )
      cb2_dd     = gc_has_drydep( 'CB2' )
      oc1_dd     = gc_has_drydep( 'OC1' )
      oc2_dd     = gc_has_drydep( 'OC2' )
      nh3_dd     = gc_has_drydep( 'NH3' )
      nh4no3_dd  = gc_has_drydep( 'NH4NO3' )
      sa1_dd     = gc_has_drydep( 'SA1' ) 
      sa2_dd     = gc_has_drydep( 'SA2' )
      sa3_dd     = gc_has_drydep( 'SA3' ) 
      sa4_dd     = gc_has_drydep( 'SA4' )
      nh4_dd     = gc_has_drydep( 'NH4' ) 
  !
      soam_ndx   = get_spc_ndx( 'SOAM' )
      soai_ndx   = get_spc_ndx( 'SOAI' )
      soat_ndx   = get_spc_ndx( 'SOAT' )
      soab_ndx   = get_spc_ndx( 'SOAB' )
      soax_ndx   = get_spc_ndx( 'SOAX' )
      sogm_ndx   = get_spc_ndx( 'SOGM' )
      sogi_ndx   = get_spc_ndx( 'SOGI' )
      sogt_ndx   = get_spc_ndx( 'SOGT' )
      sogb_ndx   = get_spc_ndx( 'SOGB' )
      sogx_ndx   = get_spc_ndx( 'SOGX' )
      soam_dd    = gc_has_drydep ( 'SOAM' )
      soai_dd    = gc_has_drydep ( 'SOAI' )
      soat_dd    = gc_has_drydep ( 'SOAT' )
      soab_dd    = gc_has_drydep ( 'SOAB' )
      soax_dd    = gc_has_drydep ( 'SOAX' )
      sogm_dd    = gc_has_drydep ( 'SOGM' )
      sogi_dd    = gc_has_drydep ( 'SOGI' )
      sogt_dd    = gc_has_drydep ( 'SOGT' )
      sogb_dd    = gc_has_drydep ( 'SOGB' )
      sogx_dd    = gc_has_drydep ( 'SOGX' )
  !
      hcn_ndx     = get_spc_ndx( 'HCN')
      ch3cn_ndx   = get_spc_ndx( 'CH3CN')
  
  ! chemUCI
      no3_ndx     = get_spc_ndx( 'NO3')
      n2o5_ndx    = get_spc_ndx( 'N2O5')
      no3_dd      = gc_has_drydep( 'NO3' )
      n2o5_dd     = gc_has_drydep( 'N2O5' )
  
    ! for the cotags kludge..
        cohc_ndx     = get_spc_ndx( 'COhc' )
        come_ndx     = get_spc_ndx( 'COme' )
        co01_ndx     = get_spc_ndx( 'CO01' )
        co02_ndx     = get_spc_ndx( 'CO02' )
        co03_ndx     = get_spc_ndx( 'CO03' )
        co04_ndx     = get_spc_ndx( 'CO04' )
        co05_ndx     = get_spc_ndx( 'CO05' )
        co06_ndx     = get_spc_ndx( 'CO06' )
        co07_ndx     = get_spc_ndx( 'CO07' )
        co08_ndx     = get_spc_ndx( 'CO08' )
        co09_ndx     = get_spc_ndx( 'CO09' )
        co10_ndx     = get_spc_ndx( 'CO10' )
        co11_ndx     = get_spc_ndx( 'CO11' )
        co12_ndx     = get_spc_ndx( 'CO12' )
        co13_ndx     = get_spc_ndx( 'CO13' )
        co14_ndx     = get_spc_ndx( 'CO14' )
        co15_ndx     = get_spc_ndx( 'CO15' )
        co16_ndx     = get_spc_ndx( 'CO16' )
        co17_ndx     = get_spc_ndx( 'CO17' )
        co18_ndx     = get_spc_ndx( 'CO18' )
        co19_ndx     = get_spc_ndx( 'CO19' )
        co20_ndx     = get_spc_ndx( 'CO20' )
        co21_ndx     = get_spc_ndx( 'CO21' )
        co22_ndx     = get_spc_ndx( 'CO22' )
        co23_ndx     = get_spc_ndx( 'CO23' )
        co24_ndx     = get_spc_ndx( 'CO24' )
        co25_ndx     = get_spc_ndx( 'CO25' )
        co26_ndx     = get_spc_ndx( 'CO26' )
        co27_ndx     = get_spc_ndx( 'CO27' )
        co28_ndx     = get_spc_ndx( 'CO28' )
        co29_ndx     = get_spc_ndx( 'CO29' )
        co30_ndx     = get_spc_ndx( 'CO30' )
        co31_ndx     = get_spc_ndx( 'CO31' )
        co32_ndx     = get_spc_ndx( 'CO32' )
        co33_ndx     = get_spc_ndx( 'CO33' )
        co34_ndx     = get_spc_ndx( 'CO34' )
        co35_ndx     = get_spc_ndx( 'CO35' )
        co36_ndx     = get_spc_ndx( 'CO36' )
        co37_ndx     = get_spc_ndx( 'CO37' )
        co38_ndx     = get_spc_ndx( 'CO38' )
        co39_ndx     = get_spc_ndx( 'CO39' )
        co40_ndx     = get_spc_ndx( 'CO40' )
        co41_ndx     = get_spc_ndx( 'CO41' )
        co42_ndx     = get_spc_ndx( 'CO42' )
  
      do i=1,nddvels
         test_name = drydep_list(i)
         m = get_spc_ndx( test_name )
         has_dvel(m) = .true.
         map_dvel(m) = i
      enddo
  
      if( all( .not. has_dvel(:) ) ) then
         return
      end if
  
      !---------------------------------------------------------------------------
      ! 	... allocate module variables
      !---------------------------------------------------------------------------
      allocate( dep_ra(pcols,gc_n_land_type,begchunk:endchunk),stat=astat )
      if( astat /= 0 ) then
         write(iulog,*) 'dvel_inti: failed to allocate dep_ra; error = ',astat
         call endrun
      end if
      allocate( dep_rb(pcols,gc_n_land_type,begchunk:endchunk),stat=astat )
      if( astat /= 0 ) then
         write(iulog,*) 'dvel_inti: failed to allocate dep_rb; error = ',astat
         call endrun
      end if
  
      if (drydep_method == DD_XLND .and. (.not.prog_modal_aero)) then
         return
      endif
  
      do_soilw = .not. dyn_soilw .and. (gc_has_drydep( 'H2' ) .or. gc_has_drydep( 'CO' ))
      allocate( gc_fraction_landuse(pcols,gc_n_land_type, begchunk:endchunk),stat=astat )
      if( astat /= 0 ) then
         write(iulog,*) 'dvel_inti: failed to allocate gc_fraction_landuse; error = ',astat
         call endrun
      end if
      if(do_soilw) then
         allocate(soilw_3d(pcols,12,begchunk:endchunk),stat=astat)
         if( astat /= 0 ) then
            write(iulog,*) 'dvel_inti: failed to allocate soilw_3d error = ',astat
            call endrun
         end if
      end if
  
      plon = get_dyn_grid_parm('plon')
      plat = get_dyn_grid_parm('plat')
      allocate( index_season_lai_j(gc_n_land_type,12),stat=astat )
      if( astat /= 0 ) then
         write(iulog,*) 'dvel_inti: failed to allocate index_season_lai_j; error = ',astat
         call endrun
      end if
      if(dycore_is('UNSTRUCTURED') ) then
         call get_landuse_and_soilw_from_file(do_soilw)
         allocate( index_season_lai(plon,12),stat=astat )
         if( astat /= 0 ) then
            write(iulog,*) 'dvel_inti: failed to allocate index_season_lai; error = ',astat
            call endrun
         end if
      else
         allocate( index_season_lai(plat,12),stat=astat )
         if( astat /= 0 ) then
            write(iulog,*) 'dvel_inti: failed to allocate index_season_lai; error = ',astat
            call endrun
         end if
         !---------------------------------------------------------------------------
         ! 	... read landuse map
         !---------------------------------------------------------------------------
         call getfil (depvel_lnd_file, locfn, 0)
         call cam_pio_openfile (piofile, trim(locfn), PIO_NOWRITE)
         !---------------------------------------------------------------------------
         ! 	... get the dimensions
         !---------------------------------------------------------------------------
         ierr = pio_inq_dimid( piofile, 'lon', dimid )
         ierr = pio_inq_dimlen( piofile, dimid, nlon_veg )
         ierr = pio_inq_dimid( piofile, 'lat', dimid )
         ierr = pio_inq_dimlen( piofile, dimid, nlat_veg )
         ierr = pio_inq_dimid( piofile, 'pft', dimid )
         ierr = pio_inq_dimlen( piofile, dimid, npft_veg )
         !---------------------------------------------------------------------------
         ! 	... allocate arrays
         !---------------------------------------------------------------------------
         allocate( vegetation_map(nlon_veg,nlat_veg,npft_veg), work(nlon_veg,nlat_veg), stat=astat )
         if( astat /= 0 ) then
            write(iulog,*) 'dvel_inti: failed to allocate vegation_map; error = ',astat
            call endrun
         end if
         allocate( urban(nlon_veg,nlat_veg), lake(nlon_veg,nlat_veg), &
              landmask(nlon_veg,nlat_veg), wetland(nlon_veg,nlat_veg), stat=astat )
         if( astat /= 0 ) then
            write(iulog,*) 'dvel_inti: failed to allocate vegation_map; error = ',astat
            call endrun
         end if
         allocate( lon_veg(nlon_veg), lat_veg(nlat_veg), &
              lon_veg_edge(nlon_veg+1), lat_veg_edge(nlat_veg+1), stat=astat )
         if( astat /= 0 ) then
            write(iulog,*) 'dvel_inti: failed to allocate vegation lon, lat arrays; error = ',astat
            call endrun
         end if
         !---------------------------------------------------------------------------
         ! 	... read the vegetation map and landmask
         !---------------------------------------------------------------------------
         ierr = pio_inq_varid( piofile, 'PCT_PFT', vid )
         ierr = pio_get_var( piofile, vid, vegetation_map )
  
         ierr = pio_inq_varid( piofile, 'LANDMASK', vid )
         ierr = pio_get_var( piofile, vid, landmask )
  
         ierr = pio_inq_varid( piofile, 'PCT_URBAN', vid )
         ierr = pio_get_var( piofile, vid, urban )
  
         ierr = pio_inq_varid( piofile, 'PCT_LAKE', vid )
         ierr = pio_get_var( piofile, vid, lake )
  
         ierr = pio_inq_varid( piofile, 'PCT_WETLAND', vid )
         ierr = pio_get_var( piofile, vid, wetland )
  
         call pio_closefile( piofile )
  
         !---------------------------------------------------------------------------
         ! scale vegetation, urban, lake, and wetland to fraction
         !---------------------------------------------------------------------------
         vegetation_map(:,:,:) = .01_r8 * vegetation_map(:,:,:)
         wetland(:,:)          = .01_r8 * wetland(:,:)
         lake(:,:)             = .01_r8 * lake(:,:)
         urban(:,:)            = .01_r8 * urban(:,:)
  #ifdef DEBUG
         if(masterproc) then
            write(iulog,*) 'minmax vegetation_map ',minval(vegetation_map),maxval(vegetation_map)
            write(iulog,*) 'minmax wetland        ',minval(wetland),maxval(wetland)
            write(iulog,*) 'minmax landmask       ',minval(landmask),maxval(landmask)
         end if
  #endif
         !---------------------------------------------------------------------------
         ! 	... define lat-lon of vegetation map (1x1)
         !---------------------------------------------------------------------------
         lat_veg(:)      = (/ (-89.5_r8 + (i-1),i=1,nlat_veg  ) /)
         lon_veg(:)      = (/ (  0.5_r8 + (i-1),i=1,nlon_veg  ) /)
         lat_veg_edge(:) = (/ (-90.0_r8 + (i-1),i=1,nlat_veg+1) /)
         lon_veg_edge(:) = (/ (  0.0_r8 + (i-1),i=1,nlon_veg+1) /)
         !---------------------------------------------------------------------------
         ! 	... read soilw table if necessary
         !---------------------------------------------------------------------------
  
         if( do_soilw ) then
            call soilw_inti( clim_soilw_file, nlon_veg, nlat_veg, soilw_map )
         end if
  
         !---------------------------------------------------------------------------
         ! 	... regrid to model grid
         !---------------------------------------------------------------------------
  
         call interp_map( plon, plat, nlon_veg, nlat_veg, npft_veg, lat_veg, lat_veg_edge, &
              lon_veg, lon_veg_edge, landmask, urban, lake, &
              wetland, vegetation_map, soilw_map, do_soilw )
  
         deallocate( vegetation_map, work, stat=astat )
         deallocate( lon_veg, lat_veg, lon_veg_edge, lat_veg_edge, stat=astat )
         deallocate( landmask, urban, lake, wetland, stat=astat )
         if( do_soilw ) then
            deallocate( soilw_map, stat=astat )
         end if
      endif  ! Unstructured grid
  
      if (drydep_method == DD_XLND) then
         return
      endif
  
      !---------------------------------------------------------------------------
      ! 	... read LAI based season indeces
      !---------------------------------------------------------------------------
      call getfil (season_wes_file, locfn, 0)
      call cam_pio_openfile (piofile, trim(locfn), PIO_NOWRITE)
      !---------------------------------------------------------------------------
      ! 	... get the dimensions
      !---------------------------------------------------------------------------
      ierr = pio_inq_dimid( piofile, 'lat', dimid )
      ierr = pio_inq_dimlen( piofile, dimid, nlat_lai )
      ierr = pio_inq_dimid( piofile, 'pft', dimid )
      ierr = pio_inq_dimlen( piofile, dimid, npft_lai )
      !---------------------------------------------------------------------------
      ! 	... allocate arrays
      !---------------------------------------------------------------------------
      allocate( lat_lai(nlat_lai), wk_lai(nlat_lai,npft_lai,12), stat=astat )
      if( astat /= 0 ) then
         write(iulog,*) 'dvel_inti: failed to allocate vegation_map; error = ',astat
         call endrun
      end if
      !---------------------------------------------------------------------------
      ! 	... read the latitude and the season indicies
      !---------------------------------------------------------------------------
      ierr = pio_inq_varid( piofile, 'lat', vid )
      ierr = pio_get_var( piofile, vid, lat_lai )
  
      ierr = pio_inq_varid( piofile, 'season_wes', vid )
      ierr = pio_get_var( piofile, vid, wk_lai )
  
      call pio_closefile( piofile )
  
  
      if(dycore_is('UNSTRUCTURED') ) then
         ! For unstructured grids plon is the 1d horizontal grid size and plat=1
         ! So this code averages at the latitude of each grid point - not an ideal solution
         allocate(clat(plon))
         call get_horiz_grid_d(plon, clat_d_out=clat)
         jl = 1
         ju = plon
      else
         allocate(clat(plat))
         call get_horiz_grid_d(plat, clat_d_out=clat)
         jl = 1
         ju = plat
      end if
      imin = 1
      do j = 1,ju
         diff_min = 10._r8
         pos_min  = -99
         target_lat = clat(j)*r2d
         do i = imin,nlat_lai
            if( abs(lat_lai(i) - target_lat) < diff_min ) then
               diff_min = abs(lat_lai(i) - target_lat)
               pos_min  = i
            end if
         end do
         if( pos_min < 0 ) then
            write(iulog,*) 'dvel_inti: cannot find ',target_lat,' at j,pos_min,diff_min = ',j,pos_min,diff_min
            write(iulog,*) 'dvel_inti: imin,nlat_lai = ',imin,nlat_lai
            write(iulog,*) 'dvel_inti: lat_lai'
            write(iulog,'(1p,10g12.5)') lat_lai(:)
            call endrun
         end if
         if(dycore_is('UNSTRUCTURED') ) then
            imin=1
         else
            imin = pos_min
         end if
         index_season_lai_j(:,:) = wk_lai(pos_min,:,:)
  
         !---------------------------------------------------------------------------
         ! specify the season as the most frequent in the 11 vegetation classes
         ! this was done to remove a banding problem in dvel (JFL Oct 04)
         !---------------------------------------------------------------------------
         do m = 1,12
            num_seas = 0
            do l = 1,11
               do k = 1,5
                  if( index_season_lai_j(l,m) == k ) then
                     num_seas(k) = num_seas(k) + 1
                     exit
                  end if
               end do
            end do
  
            num_max = -1
            do k = 1,5
               if( num_seas(k) > num_max ) then
                  num_max = num_seas(k)
                  k_max = k
               endif
            end do
  
            index_season_lai(j,m) = k_max
         end do
      end do
  
      deallocate( lat_lai, wk_lai, clat, index_season_lai_j)
  
    end subroutine gc_dvel_inti_xactive
  
    !-------------------------------------------------------------------------------------
    subroutine get_landuse_and_soilw_from_file(do_soilw)
      use cam_pio_utils, only : cam_pio_openfile
      use ncdio_atm, only : infld
      use cam_control_mod, only: aqua_planet
      use mo_drydep,       only: drydep_srf_file
      logical, intent(in) :: do_soilw
      logical :: readvar
      
      type(file_desc_t) :: piofile
      character(len=shr_kind_cl) :: locfn
      logical :: lexist
      
      if (aqua_planet) then
        gc_fraction_landuse = 0.
      else
  
        gc_drydep_srf_file = drydep_srf_file
  
        call getfil (gc_drydep_srf_file, locfn, 1, lexist)
        if(lexist) then
           call cam_pio_openfile(piofile, locfn, PIO_NOWRITE)
  
           call infld('gc_fraction_landuse', piofile, 'ncol','class',' ',1,pcols,1,gc_n_land_type, begchunk,endchunk, &
                gc_fraction_landuse, readvar, gridname='physgrid')
  
           if(do_soilw) then
              call infld('soilw', piofile, 'ncol','month',' ',1,pcols,1,12, begchunk,endchunk, &
                   soilw_3d, readvar, gridname='physgrid')
           end if
  
           call pio_closefile(piofile)
        else
           call endrun('Unstructured grids require gc_drydep_srf_file ')
        end if
      
      end if ! aqua_planet
    end subroutine get_landuse_and_soilw_from_file
  
    !-------------------------------------------------------------------------------------
    subroutine interp_map( plon, plat, nlon_veg, nlat_veg, npft_veg, lat_veg, lat_veg_edge, &
                           lon_veg, lon_veg_edge, landmask, urban, lake, &
                           wetland, vegetation_map, soilw_map, do_soilw )
  
      use mo_constants, only : r2d
      use iop_data_mod, only : latiop,loniop,scmlat,scmlon,use_replay
      use shr_scam_mod  , only: shr_scam_getCloseLatLon  ! Standardized system subroutines
      use filenames, only: ncdata
      use dycore, only : dycore_is
      use phys_grid, only : scatter_field_to_chunk
      implicit none
  
      !-------------------------------------------------------------------------------------
      ! 	... dummy arguments
      !-------------------------------------------------------------------------------------
      integer, intent(in)      ::  plon, plat, nlon_veg, nlat_veg, npft_veg
      real(r8), pointer            :: soilw_map(:,:,:)
      real(r8), intent(in)         :: landmask(nlon_veg,nlat_veg)
      real(r8), intent(in)         :: urban(nlon_veg,nlat_veg)
      real(r8), intent(in)         :: lake(nlon_veg,nlat_veg)
      real(r8), intent(in)         :: wetland(nlon_veg,nlat_veg)
      real(r8), intent(in)         :: vegetation_map(nlon_veg,nlat_veg,npft_veg)
      real(r8), intent(in)         :: lon_veg(nlon_veg)
      real(r8), intent(in)         :: lon_veg_edge(nlon_veg+1)
      real(r8), intent(in)         :: lat_veg(nlat_veg)
      real(r8), intent(in)         :: lat_veg_edge(nlat_veg+1)
      logical,  intent(in)         :: do_soilw
  
      !-------------------------------------------------------------------------------------
      ! 	... local variables
      !-------------------------------------------------------------------------------------
      real(r8) :: closelat,closelon
      integer :: latidx,lonidx
  
      integer, parameter           :: veg_ext = 20
      type(file_desc_t)            :: piofile
      integer                      :: i, j, ii, jj, jl, ju, i_ndx, n
      integer, dimension(plon+1)   :: ind_lon
      integer, dimension(plat+1)  :: ind_lat
      real(r8)                         :: total_land
      real(r8), dimension(plon+1)      :: lon_edge
      real(r8), dimension(plat+1)     :: lat_edge
      real(r8)                         :: lat1, lat2, lon1, lon2
      real(r8)                         :: x1, x2, y1, y2, dx, dy
      real(r8)                         :: area, total_area
      real(r8), dimension(npft_veg+3)  :: fraction
      real(r8)                         :: total_soilw_area
      real(r8)                         :: fraction_soilw
      real(r8)                         :: total_soilw(12)
      
      real(r8),    dimension(-veg_ext:nlon_veg+veg_ext) :: lon_veg_edge_ext
      integer, dimension(-veg_ext:nlon_veg+veg_ext) :: mapping_ext
  
      real(r8), allocatable :: lam(:), phi(:), garea(:)
  
      logical, parameter :: has_npole = .true.
      integer :: ploniop,platiop
      character(len=shr_kind_cl) :: ncdata_loc
      real(r8) :: tmp_frac_lu(plon,gc_n_land_type,plat), tmp_soilw_3d(plon,12,plat)
  
      allocate(lam(plon), phi(plat))
      call get_horiz_grid_d(plon, clon_d_out=lam)
      call get_horiz_grid_d(plat, clat_d_out=phi)
  
  
  
      jl = 1
      ju = plon
  
      if (single_column) then
         if (use_replay) then
            call getfil (ncdata, ncdata_loc)
            call cam_pio_openfile (piofile, trim(ncdata_loc), PIO_NOWRITE)
            call shr_scam_getCloseLatLon(piofile,scmlat,scmlon,closelat,closelon,latidx,lonidx)
            call pio_closefile ( piofile)
            ploniop=size(loniop)
            platiop=size(latiop)
         else 
            latidx=1
            lonidx=1
            ploniop=1
            platiop=1
         end if
        
         lon_edge(1) = loniop(lonidx) * r2d - .5_r8*(loniop(2) - loniop(1)) * r2d
  
         if (lonidx.lt.ploniop) then
            lon_edge(2) = loniop(lonidx+1) * r2d - .5_r8*(loniop(2) - loniop(1)) * r2d
         else
            lon_edge(2) = lon_edge(1) + (loniop(2) - loniop(1)) * r2d
         end if
  
         lat_edge(1) = latiop(latidx) * r2d - .5_r8*(latiop(2) - latiop(1)) * r2d
  
         if (latidx.lt.platiop) then
            lat_edge(2) = latiop(latidx+1) * r2d - .5_r8*(latiop(2) - latiop(1)) * r2d
         else
            lat_edge(2) = lat_edge(1) + (latiop(2) - latiop(1)) * r2d
         end if       
      else
         do i = 1,plon
            lon_edge(i) = lam(i) * r2d - .5_r8*(lam(2) - lam(1)) * r2d
         end do
         lon_edge(plon+1) = lon_edge(plon) + (lam(2) - lam(1)) * r2d
         if( .not. has_npole ) then
            do j = 1,plat+1
               lat_edge(j) = phi(j) * r2d - .5_r8*(phi(2) - phi(1)) * r2d
            end do
         else
            do j = 1,plat
               lat_edge(j) = phi(j) * r2d - .5_r8*(phi(2) - phi(1)) * r2d
            end do
            lat_edge(plat+1) = lat_edge(plat) + (phi(2) - phi(1)) * r2d
         end if
      end if
      do j = 1,plat+1
         lat_edge(j) = min( lat_edge(j), 90._r8 )
         lat_edge(j) = max( lat_edge(j),-90._r8 )
      end do
  
      !-------------------------------------------------------------------------------------
      ! wrap around the longitudes
      !-------------------------------------------------------------------------------------
      do i = -veg_ext,0
         lon_veg_edge_ext(i) = lon_veg_edge(nlon_veg+i) - 360._r8
         mapping_ext     (i) =              nlon_veg+i
      end do
      do i = 1,nlon_veg
         lon_veg_edge_ext(i) = lon_veg_edge(i)
         mapping_ext     (i) =              i
      end do
      do i = nlon_veg+1,nlon_veg+veg_ext
         lon_veg_edge_ext(i) = lon_veg_edge(i-nlon_veg) + 360._r8
         mapping_ext     (i) =              i-nlon_veg
      end do
  #ifdef DEBUG
      write(iulog,*) 'interp_map : lon_edge ',lon_edge
      write(iulog,*) 'interp_map : lat_edge ',lat_edge
      write(iulog,*) 'interp_map : mapping_ext ',mapping_ext
  #endif
      do j = 1,plon+1
         lon1 = lon_edge(j) 
         do i = -veg_ext,nlon_veg+veg_ext
            dx = lon_veg_edge_ext(i  ) - lon1
            dy = lon_veg_edge_ext(i+1) - lon1
            if( dx*dy <= 0._r8 ) then
               ind_lon(j) = i
               exit
            end if
         end do
      end do
  
      do j = 1,plat+1
         lat1 = lat_edge(j)
         do i = 1,nlat_veg
            dx = lat_veg_edge(i  ) - lat1
            dy = lat_veg_edge(i+1) - lat1
            if( dx*dy <= 0._r8 ) then
               ind_lat(j) = i
               exit
            end if
         end do
      end do
  #ifdef DEBUG
      write(iulog,*) 'interp_map : ind_lon ',ind_lon
      write(iulog,*) 'interp_map : ind_lat ',ind_lat
  #endif
      lat_loop : do j = 1,plat
         lon_loop : do i = 1,plon
            total_area       = 0._r8
            fraction         = 0._r8
            total_soilw(:)   = 0._r8
            total_soilw_area = 0._r8
            do jj = ind_lat(j),ind_lat(j+1)
               y1 = max( lat_edge(j),lat_veg_edge(jj) )
               y2 = min( lat_edge(j+1),lat_veg_edge(jj+1) ) 
               dy = (y2 - y1)/(lat_veg_edge(jj+1) - lat_veg_edge(jj))
               do ii =ind_lon(i),ind_lon(i+1)
                  i_ndx = mapping_ext(ii)
                  x1 = max( lon_edge(i),lon_veg_edge_ext(ii) )
                  x2 = min( lon_edge(i+1),lon_veg_edge_ext(ii+1) ) 
                  dx = (x2 - x1)/(lon_veg_edge_ext(ii+1) - lon_veg_edge_ext(ii))
                  area = dx * dy
                  total_area = total_area + area
                  !-----------------------------------------------------------------
                  ! 	... special case for ocean grid point 
                  !-----------------------------------------------------------------
                  if( nint(landmask(i_ndx,jj)) == 0 ) then
                     fraction(npft_veg+1) = fraction(npft_veg+1) + area
                  else
                     do n = 1,npft_veg
                        fraction(n) = fraction(n) + vegetation_map(i_ndx,jj,n) * area
                     end do
                     fraction(npft_veg+1) = fraction(npft_veg+1) + area * lake   (i_ndx,jj)
                     fraction(npft_veg+2) = fraction(npft_veg+2) + area * wetland(i_ndx,jj)
                     fraction(npft_veg+3) = fraction(npft_veg+3) + area * urban  (i_ndx,jj)
                     !-----------------------------------------------------------------
                     ! 	... check if land accounts for the whole area.
                     !           If not, the remaining area is in the ocean
                     !-----------------------------------------------------------------
                     total_land = sum(vegetation_map(i_ndx,jj,:)) &
                                + urban  (i_ndx,jj) &
                                + lake   (i_ndx,jj) &
                                + wetland(i_ndx,jj)
                     if( total_land < 1._r8 ) then
                        fraction(npft_veg+1) = fraction(npft_veg+1) + (1._r8 - total_land) * area
                     end if
                     !-------------------------------------------------------------------------------------
                     ! 	... compute weighted average of soilw over grid (non-water only)
                     !-------------------------------------------------------------------------------------
                     if( do_soilw ) then
                        fraction_soilw = total_land  - (lake(i_ndx,jj) + wetland(i_ndx,jj))
                        total_soilw_area = total_soilw_area + fraction_soilw * area
                        total_soilw(:)   = total_soilw(:) + fraction_soilw * area * soilw_map(i_ndx,jj,:)
                     end if
                  end if
               end do
            end do
            !-------------------------------------------------------------------------------------
            ! 	... divide by total area of grid box
            !-------------------------------------------------------------------------------------
            fraction(:) = fraction(:)/total_area
            !-------------------------------------------------------------------------------------
            ! 	... make sure we don't have too much or too little
            !-------------------------------------------------------------------------------------
            if( abs( sum(fraction) - 1._r8) > .001_r8 ) then
               fraction(:) = fraction(:)/sum(fraction)
            end if
            !-------------------------------------------------------------------------------------
            ! 	... map to Wesely land classification
            !-------------------------------------------------------------------------------------
  
            
  
  
            tmp_frac_lu(i, 1, j) =     fraction(20)
            tmp_frac_lu(i, 2, j) = sum(fraction(16:17))
            tmp_frac_lu(i, 3, j) = sum(fraction(13:15))
            tmp_frac_lu(i, 4, j) = sum(fraction( 5: 9))
            tmp_frac_lu(i, 5, j) = sum(fraction( 2: 4))
            tmp_frac_lu(i, 6, j) =     fraction(19)
            tmp_frac_lu(i, 7, j) =     fraction(18)
            tmp_frac_lu(i, 8, j) =     fraction( 1)
            tmp_frac_lu(i, 9, j) = 0._r8
            tmp_frac_lu(i,10, j) = 0._r8
            tmp_frac_lu(i,11, j) = sum(fraction(10:12))
            if( do_soilw ) then
               if( total_soilw_area > 0._r8 ) then
                  tmp_soilw_3d(i,:,j) = total_soilw(:)/total_soilw_area
               else
                  tmp_soilw_3d(i,:,j) = -99._r8
               end if
            end if
         end do lon_loop
      end do lat_loop
      !-------------------------------------------------------------------------------------
      ! 	... reshape according to lat-lon blocks
      !-------------------------------------------------------------------------------------
      call scatter_field_to_chunk(1,gc_n_land_type,1,plon,tmp_frac_lu,gc_fraction_landuse)
      if(do_soilw) call scatter_field_to_chunk(1,12,1,plon,tmp_soilw_3d,soilw_3d)
      !-------------------------------------------------------------------------------------
      ! 	... make sure there are no out of range values
      !-------------------------------------------------------------------------------------
      where (gc_fraction_landuse < 0._r8) gc_fraction_landuse = 0._r8
      where (gc_fraction_landuse > 1._r8) gc_fraction_landuse = 1._r8
  
    end subroutine interp_map
    
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
      integer                :: index_season(ncol,gc_n_land_type)
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
  
    !-------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------
    subroutine soilw_inti( ncfile, nlon_veg, nlat_veg, soilw_map )
      !------------------------------------------------------------------
      !	... read primary soil moisture table
      !------------------------------------------------------------------
  
      use time_manager,  only : get_calday
  
      implicit none
  
      !------------------------------------------------------------------
      !	... dummy args
      !------------------------------------------------------------------
      integer, intent(in) :: &
           nlon_veg, &
           nlat_veg
      real(r8), pointer :: soilw_map(:,:,:)
      character(len=*), intent(in) :: ncfile ! file name of netcdf file containing data
  
      !------------------------------------------------------------------
      !	... local variables
      !------------------------------------------------------------------
      integer :: gndx = 0
      integer :: nlat, &             ! # of lats in soilw file
                 nlon                ! # of lons in soilw file
      integer :: i, ip, k, m
      integer :: j, jl, ju
      integer :: lev, day, ierr
      type(file_desc_t) :: piofile
      type(var_desc_t) :: vid
      
      integer :: dimid_lat, dimid_lon, dimid_time
      integer :: dates(12) = (/ 116, 214, 316, 415,  516,  615, &
                                716, 816, 915, 1016, 1115, 1216 /)
  
      character(len=shr_kind_cl) :: locfn
  
      !-----------------------------------------------------------------------
      !       ... open netcdf file
      !-----------------------------------------------------------------------
      call getfil (ncfile, locfn, 0)
      call cam_pio_openfile (piofile, trim(locfn), PIO_NOWRITE)
  
      !-----------------------------------------------------------------------
      !       ... get longitudes
      !-----------------------------------------------------------------------
      ierr = pio_inq_dimid( piofile, 'lon', dimid_lon )
      ierr = pio_inq_dimlen( piofile, dimid_lon, nlon )
      if( nlon /= nlon_veg ) then
         write(iulog,*) 'soilw_inti: soil and vegetation lons differ; ',nlon, nlon_veg
         call endrun
      end if
      !-----------------------------------------------------------------------
      !       ... get latitudes
      !-----------------------------------------------------------------------
      ierr = pio_inq_dimid( piofile, 'lat', dimid_lat )
      ierr = pio_inq_dimlen( piofile, dimid_lat, nlat )
      if( nlat /= nlat_veg ) then
         write(iulog,*) 'soilw_inti: soil and vegetation lats differ; ',nlat, nlat_veg
         call endrun
      end if
      !-----------------------------------------------------------------------
      !       ... set times (days of year)
      !-----------------------------------------------------------------------
      ierr = pio_inq_dimid( piofile, 'time', dimid_time )
      ierr = pio_inq_dimlen( piofile, dimid_time, ndays )
      if( ndays /= 12 ) then
         write(iulog,*) 'soilw_inti: dataset not a cyclical year'
         call endrun
      end if
      allocate( days(ndays),stat=ierr )
      if( ierr /= 0 ) then
         write(iulog,*) 'soilw_inti: days allocation error = ',ierr
         call endrun
      end if
      do m = 1,min(12,ndays)
         days(m) = get_calday( dates(m), 0 )
      end do
  
      !------------------------------------------------------------------
      !	... allocate arrays
      !------------------------------------------------------------------
      allocate( soilw_map(nlon,nlat,ndays), stat=ierr )
      if( ierr /= 0 ) then
         write(iulog,*) 'soilw_inti: soilw_map allocation error = ',ierr
         call endrun
      end if
  
      !------------------------------------------------------------------
      !	... read in the soil moisture
      !------------------------------------------------------------------
      ierr = pio_inq_varid( piofile, 'SOILW', vid )
      ierr = pio_get_var( piofile, vid, soilw_map )
      !------------------------------------------------------------------
      !	... close file
      !------------------------------------------------------------------
      call pio_closefile( piofile )
  
    end subroutine soilw_inti
    
    !-------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------
    subroutine chk_soilw( calday )
      !--------------------------------------------------------------------
      !	... check timing for ub values
      !--------------------------------------------------------------------
  
      use mo_constants, only : dayspy
  
      implicit none
  
      !--------------------------------------------------------------------
      !	... dummy args
      !--------------------------------------------------------------------
      real(r8), intent(in)    :: calday
  
      !--------------------------------------------------------------------
      !	... local variables
      !--------------------------------------------------------------------
      integer  ::  m, upper
      real(r8)     ::  numer, denom
  
      !--------------------------------------------------------
      !	... setup the time interpolation
      !--------------------------------------------------------
      if( calday < days(1) ) then
         next = 1
         last = ndays
      else
         if( days(ndays) < dayspy ) then
            upper = ndays
         else
            upper = ndays - 1
         end if
         do m = upper,1,-1
            if( calday >= days(m) ) then
               exit
            end if
         end do
         last = m
         next = mod( m,ndays ) + 1
      end if
      numer = calday - days(last)
      denom = days(next) - days(last)
      if( numer < 0._r8 ) then
         numer = dayspy + numer
      end if
      if( denom < 0._r8 ) then
         denom = dayspy + denom
      end if
      dels = max( min( 1._r8,numer/denom ),0._r8 )
  
    end subroutine chk_soilw
  
    !-------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------
    subroutine set_soilw( soilw, lchnk, calday )
      !--------------------------------------------------------------------
      !	... set the soil moisture
      !--------------------------------------------------------------------
  
      implicit none
  
      !--------------------------------------------------------------------
      !	... dummy args
      !--------------------------------------------------------------------
      real(r8), intent(inout) :: soilw(pcols)
      integer,  intent(in)    :: lchnk           ! chunk indice
      real(r8), intent(in)    :: calday
  
  
      integer :: i, ilon,ilat
  
      call chk_soilw( calday )
  
      soilw(:) = soilw_3d(:,last,lchnk) + dels *( soilw_3d(:,next,lchnk) - soilw_3d(:,last,lchnk))
  
    end subroutine set_soilw
  
    !-------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------
    function gc_has_drydep( name )
  
      implicit none
  
      character(len=*), intent(in) :: name
  
      logical :: gc_has_drydep
      integer :: i
  
      gc_has_drydep = .false.
  
      do i=1,nddvels
         if ( trim(name) == trim(drydep_list(i)) ) then
           gc_has_drydep = .true.
           exit
         endif
      enddo
  
    endfunction gc_has_drydep
  
  !----------------------------------------------------------------------------------------
  
  !==============================================================================
  !==============================================================================
  subroutine gc_wetdep_init()
    use physics_buffer, only: pbuf_get_index
    use constituents,   only: cnst_get_ind
    use phys_control,   only: phys_getopts
  
    cld_idx             = pbuf_get_index('CLD')    
    qme_idx             = pbuf_get_index('QME')    
    prain_idx           = pbuf_get_index('PRAIN')  
    nevapr_idx          = pbuf_get_index('NEVAPR') 
  
    icwmrdp_idx         = pbuf_get_index('ICWMRDP') 
    rprddp_idx          = pbuf_get_index('RPRDDP')  
    icwmrsh_idx         = pbuf_get_index('ICWMRSH') 
    rprdsh_idx          = pbuf_get_index('RPRDSH')  
    sh_frac_idx         = pbuf_get_index('SH_FRAC' )
    dp_frac_idx         = pbuf_get_index('DP_FRAC') 
    nevapr_shcu_idx     = pbuf_get_index('NEVAPR_SHCU') 
    nevapr_dpcu_idx     = pbuf_get_index('NEVAPR_DPCU') 
  
    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)
    call phys_getopts(pergro_mods_out = pergro_mods)
  
  endsubroutine gc_wetdep_init
  
  !==============================================================================
  ! gathers up the inputs needed for the gc_wetdepa routines
  !==============================================================================
  subroutine gc_wetdep_inputs_set( ncol, tfld, pmid, pdel, qliq, qice, pbuf, inputs )
    use phys_control,   only: cam_physpkg_is
    use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx
  
    ! args
    real(r8) :: tfld(pcols,pver),pmid(pcols,pver),pdel(pcols,pver)
    real(r8) :: qliq(ncol,pver),qice(ncol,pver)
    type(physics_buffer_desc), pointer :: pbuf(:)         !! physics buffer
    type(gc_wetdep_inputs_t), intent(out) :: inputs          !! collection of gc_wetdepa inputs
  
    ! local vars
  
    real(r8), pointer :: icwmrdp(:,:)    ! in cloud water mixing ratio, deep convection
    real(r8), pointer :: rprddp(:,:)     ! rain production, deep convection
    real(r8), pointer :: icwmrsh(:,:)    ! in cloud water mixing ratio, deep convection
    real(r8), pointer :: rprdsh(:,:)     ! rain production, deep convection
    real(r8), pointer :: sh_frac(:,:)    ! Shallow convective cloud fraction
    real(r8), pointer :: dp_frac(:,:)    ! Deep convective cloud fraction
    real(r8), pointer :: evapcsh(:,:)    ! Evaporation rate of shallow convective precipitation >=0.
    real(r8), pointer :: evapcdp(:,:)    ! Evaporation rate of deep    convective precipitation >=0.
  
    real(r8) :: rainmr(pcols,pver)       ! mixing ratio of rain within cloud volume
    real(r8) :: cldst(pcols,pver)        ! Stratiform cloud fraction
  
    integer :: itim, ncol
    integer :: ierror
  
    itim = pbuf_old_tim_idx()
  
    allocate (inputs%cldcu(pcols,pver), stat=ierror)
    if ( ierror /= 0 ) call endrun('WETDEP_INPUTS_SET error: allocation error cldcu')
  
    allocate (inputs%evapc(pcols,pver), stat=ierror)
    if ( ierror /= 0 ) call endrun('WETDEP_INPUTS_SET error: allocation error evapc')
  
    allocate (inputs%cmfdqr(pcols,pver), stat=ierror)
    if ( ierror /= 0 ) call endrun('WETDEP_INPUTS_SET error: allocation error cmfdqr')
  
    allocate (inputs%conicw(pcols,pver), stat=ierror)
    if ( ierror /= 0 ) call endrun('WETDEP_INPUTS_SET error: allocation error conicw')
  
    allocate (inputs%totcond(pcols,pver), stat=ierror)
    if ( ierror /= 0 ) call endrun('WETDEP_INPUTS_SET error: allocation error totcond')
  
    allocate (inputs%cldv(pcols,pver), stat=ierror)
    if ( ierror /= 0 ) call endrun('WETDEP_INPUTS_SET error: allocation error cldv')
  
    allocate (inputs%cldvcu(pcols,pver), stat=ierror)
    if ( ierror /= 0 ) call endrun('WETDEP_INPUTS_SET error: allocation error cldvcu')
  
    allocate (inputs%cldvst(pcols,pver), stat=ierror)
    if ( ierror /= 0 ) call endrun('WETDEP_INPUTS_SET error: allocation error cldvst')
  
    call pbuf_get_field(pbuf, cld_idx,         inputs%cldt, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
    call pbuf_get_field(pbuf, qme_idx,         inputs%qme     )
    call pbuf_get_field(pbuf, prain_idx,       inputs%prain   )
    call pbuf_get_field(pbuf, nevapr_idx,      inputs%evapr   )
    call pbuf_get_field(pbuf, icwmrdp_idx,     icwmrdp )
    call pbuf_get_field(pbuf, icwmrsh_idx,     icwmrsh )
    call pbuf_get_field(pbuf, rprddp_idx,      rprddp  )
    call pbuf_get_field(pbuf, rprdsh_idx,      rprdsh  )
    call pbuf_get_field(pbuf, sh_frac_idx,     sh_frac )
    call pbuf_get_field(pbuf, dp_frac_idx,     dp_frac )
    call pbuf_get_field(pbuf, nevapr_shcu_idx, evapcsh )
    call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp )
  
    inputs%cldcu(:ncol,:)  = dp_frac(:ncol,:) + sh_frac(:ncol,:)
    cldst(:ncol,:)          = inputs%cldt(:ncol,:) - inputs%cldcu(:ncol,:)       ! Stratiform cloud fraction
    inputs%evapc(:ncol,:)  = evapcsh(:ncol,:) + evapcdp(:ncol,:)
    inputs%cmfdqr(:ncol,:) = rprddp(:ncol,:)  + rprdsh(:ncol,:)
  
    ! sum deep and shallow convection contributions
    if (cam_physpkg_is('default')) then
       ! Dec.29.2009. Sungsu
       inputs%conicw(:ncol,:) = (icwmrdp(:ncol,:)*dp_frac(:ncol,:) + icwmrsh(:ncol,:)*sh_frac(:ncol,:))/ &
                                max(0.01_r8, sh_frac(:ncol,:) + dp_frac(:ncol,:))
    else
       inputs%conicw(:ncol,:) = icwmrdp(:ncol,:) + icwmrsh(:ncol,:)
    end if
  
    inputs%totcond(:ncol,:) = qliq(:ncol,:) + qice(:ncol,:)
  
    call gc_clddiag( tfld,     pmid,               pdel,   inputs%cmfdqr, inputs%evapc, &
                 inputs%cldt,  inputs%cldcu,       cldst,  inputs%qme,    inputs%evapr, &
                 inputs%prain, inputs%cldv, inputs%cldvcu, inputs%cldvst,       rainmr, &
                 ncol )
  
  end subroutine gc_wetdep_inputs_set
  
  !==============================================================================
  ! deallocate storage assoicated with gc_wetdep_inputs_t type variable
  !==============================================================================
  subroutine gc_wetdep_inputs_unset(inputs)
  
    ! args
    type(gc_wetdep_inputs_t), intent(inout) :: inputs          !! collection of gc_wetdepa inputs
  
    deallocate(inputs%cldcu)
    deallocate(inputs%evapc)
    deallocate(inputs%cmfdqr)
    deallocate(inputs%conicw)
    deallocate(inputs%totcond)
    deallocate(inputs%cldv)
    deallocate(inputs%cldvcu)
    deallocate(inputs%cldvst)
  
  end subroutine gc_wetdep_inputs_unset
  
  subroutine gc_clddiag(t, pmid, pdel, cmfdqr, evapc, &
                     cldt, cldcu, cldst, cme, evapr, &
                     prain, cldv, cldvcu, cldvst, rain, &
                     ncol)
  
      use physconst,     only : rair,pi
  
     ! ------------------------------------------------------------------------------------ 
     ! Estimate the cloudy volume which is occupied by rain or cloud water as
     ! the max between the local cloud amount or the
     ! sum above of (cloud*positive precip production)      sum total precip from above
     !              ----------------------------------   x ------------------------
     ! sum above of     (positive precip           )        sum positive precip from above
     ! Author: P. Rasch
     !         Sungsu Park. Mar.2010 
     ! ------------------------------------------------------------------------------------
  
     ! Input arguments:
     real(r8), intent(in) :: t(pcols,pver)        ! temperature (K)
     real(r8), intent(in) :: pmid(pcols,pver)     ! pressure at layer midpoints
     real(r8), intent(in) :: pdel(pcols,pver)     ! pressure difference across layers
     real(r8), intent(in) :: cmfdqr(pcols,pver)   ! dq/dt due to convective rainout 
     real(r8), intent(in) :: evapc(pcols,pver)    ! Evaporation rate of convective precipitation ( >= 0 ) 
     real(r8), intent(in) :: cldt(pcols,pver)    ! total cloud fraction
     real(r8), intent(in) :: cldcu(pcols,pver)    ! Cumulus cloud fraction
     real(r8), intent(in) :: cldst(pcols,pver)    ! Stratus cloud fraction
     real(r8), intent(in) :: cme(pcols,pver)      ! rate of cond-evap within the cloud
     real(r8), intent(in) :: evapr(pcols,pver)    ! rate of evaporation of falling precipitation (kg/kg/s)
     real(r8), intent(in) :: prain(pcols,pver)    ! rate of conversion of condensate to precipitation (kg/kg/s)
     integer, intent(in) :: ncol
  
     ! Output arguments:
     real(r8), intent(out) :: cldv(pcols,pver)     ! fraction occupied by rain or cloud water 
     real(r8), intent(out) :: cldvcu(pcols,pver)   ! Convective precipitation volume
     real(r8), intent(out) :: cldvst(pcols,pver)   ! Stratiform precipitation volume
     real(r8), intent(out) :: rain(pcols,pver)     ! mixing ratio of rain (kg/kg)
  
     ! Local variables:
     integer  i, k
     real(r8) convfw         ! used in fallspeed calculation; taken from findmcnew
     real(r8) sumppr(pcols)        ! precipitation rate (kg/m2-s)
     real(r8) sumpppr(pcols)       ! sum of positive precips from above
     real(r8) cldv1(pcols)         ! precip weighted cloud fraction from above
     real(r8) lprec                ! local production rate of precip (kg/m2/s)
     real(r8) lprecp               ! local production rate of precip (kg/m2/s) if positive
     real(r8) rho                  ! air density
     real(r8) vfall
     real(r8) sumppr_cu(pcols)     ! Convective precipitation rate (kg/m2-s)
     real(r8) sumpppr_cu(pcols)    ! Sum of positive convective precips from above
     real(r8) cldv1_cu(pcols)      ! Convective precip weighted convective cloud fraction from above
     real(r8) lprec_cu             ! Local production rate of convective precip (kg/m2/s)
     real(r8) lprecp_cu            ! Local production rate of convective precip (kg/m2/s) if positive
     real(r8) sumppr_st(pcols)     ! Stratiform precipitation rate (kg/m2-s)
     real(r8) sumpppr_st(pcols)    ! Sum of positive stratiform precips from above
     real(r8) cldv1_st(pcols)      ! Stratiform precip weighted stratiform cloud fraction from above
     real(r8) lprec_st             ! Local production rate of stratiform precip (kg/m2/s)
     real(r8) lprecp_st            ! Local production rate of stratiform precip (kg/m2/s) if positive
     ! -----------------------------------------------------------------------
  
     convfw = 1.94_r8*2.13_r8*sqrt(rhoh2o*gravit*2.7e-4_r8)
     do i=1,ncol
        sumppr(i) = 0._r8
        cldv1(i) = 0._r8
        sumpppr(i) = 1.e-36_r8
        sumppr_cu(i)  = 0._r8
        cldv1_cu(i)   = 0._r8
        sumpppr_cu(i) = 1.e-36_r8
        sumppr_st(i)  = 0._r8
        cldv1_st(i)   = 0._r8
        sumpppr_st(i) = 1.e-36_r8
     end do
  
     do k = 1,pver
        do i = 1,ncol
           cldv(i,k) = &
              max(min(1._r8, &
              cldv1(i)/sumpppr(i) &
              )*sumppr(i)/sumpppr(i), &
              cldt(i,k) &
              )
           lprec = pdel(i,k)/gravit &
              *(prain(i,k)+cmfdqr(i,k)-evapr(i,k))
           lprecp = max(lprec,1.e-30_r8)
           cldv1(i) = cldv1(i)  + cldt(i,k)*lprecp
           sumppr(i) = sumppr(i) + lprec
           sumpppr(i) = sumpppr(i) + lprecp
  
           ! For convective precipitation volume at the top interface of each layer. Neglect the current layer.
           cldvcu(i,k)   = max(min(1._r8,cldv1_cu(i)/sumpppr_cu(i))*(sumppr_cu(i)/sumpppr_cu(i)),0._r8)
           lprec_cu      = (pdel(i,k)/gravit)*(cmfdqr(i,k)-evapc(i,k))
           lprecp_cu     = max(lprec_cu,1.e-30_r8)
           cldv1_cu(i)   = cldv1_cu(i) + cldcu(i,k)*lprecp_cu
           sumppr_cu(i)  = sumppr_cu(i) + lprec_cu
           sumpppr_cu(i) = sumpppr_cu(i) + lprecp_cu
  
           ! For stratiform precipitation volume at the top interface of each layer. Neglect the current layer.
           cldvst(i,k)   = max(min(1._r8,cldv1_st(i)/sumpppr_st(i))*(sumppr_st(i)/sumpppr_st(i)),0._r8)
           lprec_st      = (pdel(i,k)/gravit)*(prain(i,k)-evapr(i,k))
           lprecp_st     = max(lprec_st,1.e-30_r8)
           cldv1_st(i)   = cldv1_st(i) + cldst(i,k)*lprecp_st
           sumppr_st(i)  = sumppr_st(i) + lprec_st
           sumpppr_st(i) = sumpppr_st(i) + lprecp_st
  
           rain(i,k) = 0._r8
           if(t(i,k) .gt. tmelt) then
              rho = pmid(i,k)/(rair*t(i,k))
              vfall = convfw/sqrt(rho)
              rain(i,k) = sumppr(i)/(rho*vfall)
              if (rain(i,k).lt.1.e-14_r8) rain(i,k) = 0._r8
           endif
        end do
     end do
  
  end subroutine gc_clddiag
  
  subroutine gc_wetdep( ncol, deltat, &
                         t, p, q, zi, pdel, &
                         cmfdqr, evapc, dlf, conicw, &
                         precs, conds, evaps, cwat, &
                         cldt, cldc, cldv, cldvcu, cldvst, &
                         tracer, scavt )
  
        !----------------------------------------------------------------------- 
        ! Purpose: 
        ! scavenging code for very soluble aerosols
        ! 
        ! Author: P. Rasch
        ! Modified by T. Bond 3/2003 to track different removals
        ! Sungsu Park. Mar.2010 : Impose consistencies with a few changes in physics.
        !-----------------------------------------------------------------------
  
        use phys_control, only: phys_getopts
        use physconst,    only: gravit, rair, tmelt
  
        implicit none
  
        integer, intent(in) :: ncol
  
        real(r8), intent(in) ::&
           deltat,               &! time step
           t(pcols,pver),        &! temperature
           p(pcols,pver),        &! pressure
           q(pcols,pver),        &! moisture
           zi(pcols,pver+1),     &! HGT
           pdel(pcols,pver),     &! pressure thikness
           cmfdqr(pcols,pver),   &! rate of production of convective precip
  ! Sungsu
           evapc(pcols,pver),    &! Evaporation rate of convective precipitation
           dlf(pcols,pver),      &! Detrainment of convective condensate [kg/kg/s]
  ! Sungsu
           conicw(pcols,pver),   &! convective cloud water
           precs(pcols,pver),    &! rate of production of stratiform precip
           conds(pcols,pver),    &! rate of production of condensate
           evaps(pcols,pver),    &! rate of evaporation of precip
           cwat(pcols,pver),     &! cloud water amount 
           cldt(pcols,pver),     &! total cloud fraction
           cldc(pcols,pver),     &! convective cloud fraction
           cldv(pcols,pver),     &! total cloud fraction
  ! Sungsu
           cldvcu(pcols,pver),   &! Convective precipitation area at the top interface of each layer
           cldvst(pcols,pver),   &! Stratiform precipitation area at the top interface of each layer
  ! Sungsu
           tracer(pcols,pver,19)     ! trace species
        ! If subroutine is called with just sol_fact:
              ! sol_fact is used for both in- and below-cloud scavenging
        ! If subroutine is called with optional argument sol_facti_in:
              ! sol_fact  is used for below cloud scavenging
              ! sol_facti is used for in cloud scavenging
  
        real(r8), intent(out) :: scavt(pcols,pver,19)   ! scavenging tend 
  
        ! local variables
  
        integer i                 ! x index
        integer k                 ! z index
        integer m                 ! spc index
  
        real(r8) adjfac               ! factor stolen from cmfmca
        real(r8) aqfrac               ! fraction of tracer in aqueous phase
        real(r8) cwatc                ! local convective total water amount 
        real(r8) cwats                ! local stratiform total water amount 
        real(r8) cwatp                ! local water amount falling from above precip
        real(r8) fracp                ! fraction of cloud water converted to precip
        real(r8) gafrac               ! fraction of tracer in gas phasea
        real(r8) hconst               ! henry's law solubility constant when equation is expressed
                                  ! in terms of mixing ratios
        real(r8) mpla                 ! moles / liter H2O entering the layer from above
        real(r8) mplb                 ! moles / liter H2O leaving the layer below
        real(r8) omsm                 ! 1 - (a small number)
        real(r8) part                 !  partial pressure of tracer in atmospheres
        real(r8) patm                 ! total pressure in atmospheres
        real(r8) pdog                 ! work variable (pdel/gravit)
        real(r8) precabc(pcols)       ! conv precip from above (work array)
        real(r8) precabs(pcols)       ! strat precip from above (work array)
        real(r8) precbl               ! precip falling out of level (work array)
        real(r8) precmin              ! minimum convective precip causing scavenging
        real(r8) rat(pcols)           ! ratio of amount available to amount removed
        real(r8) scavab(19)           ! scavenged tracer flux from above (work array)
        real(r8) scavabc(19)          ! scavenged tracer flux from above (work array)
        real(r8) srcc                 ! tend for convective rain
        real(r8) srcs                 ! tend for stratiform rain
        real(r8) srct(pcols)          ! work variable
        real(r8) tracab(pcols)        ! column integrated tracer amount
  !      real(r8) vfall                ! fall speed of precip
        real(r8) fins                 ! fraction of rem. rate by strat rain
        real(r8) finc                 ! fraction of rem. rate by conv. rain
        real(r8) srcs1                ! work variable
        real(r8) srcs2                ! work variable
        real(r8) tc                   ! temp in celcius
        real(r8) weight               ! fraction of condensate which is ice
        real(r8) cldmabs(pcols)       ! maximum cloud at or above this level
        real(r8) cldmabc(pcols)       ! maximum cloud at or above this level
        real(r8) odds                 ! limit on removal rate (proportional to prec)
        real(r8) dblchek(pcols)
        logical :: found
  
      ! Jan.16.2009. Sungsu for wet scavenging below clouds.
      ! real(r8) cldovr_cu(pcols)     ! Convective precipitation area at the base of each layer
      ! real(r8) cldovr_st(pcols)     ! Stratiform precipitation area at the base of each layer
  
        real(r8) tracer_incu
        real(r8) tracer_mean
  
      ! End by Sungsu
  
  !     real(r8) sol_facti,  sol_factb  ! in cloud and below cloud fraction of aerosol scavenged
  !     real(r8) sol_factii, sol_factbi ! in cloud and below cloud fraction of aerosol scavenged by ice
  !     real(r8) sol_factic(pcols,pver)             ! sol_facti for convective clouds
  !     real(r8) sol_factiic            ! sol_factii for convective clouds
        ! sol_factic & solfact_iic added for MODAL_AERO.  
        ! For stratiform cloud, cloudborne aerosol is treated explicitly,
        !    and sol_facti is 1.0 for cloudborne, 0.0 for interstitial.
        ! For convective cloud, cloudborne aerosol is not treated explicitly,
        !    and sol_factic is 1.0 for both cloudborne and interstitial.
  
        integer  jstrcnv
  
        real(r8), parameter :: prec_smallaa = 1.0e-30_r8  ! 1e-30 kg/m2/s (or mm/s) = 3.2e-23 mm/yr
        real(r8), parameter :: x_smallaa = 1.0e-30_r8
  
        real(r8) arainx
        real(r8) evapx
        real(r8) pprdx
        real(r8) precabc_base(pcols)  ! conv precip at an effective cloud base for calculations in a particular layer
        real(r8) precabs_base(pcols)  ! strat precip at an effective cloud base for calculations in a particular layer
        real(r8) precabx_old, precabx_tmp, precabx_new
        real(r8) precabx_base_old, precabx_base_tmp, precabx_base_new
        real(r8) precnums_base(pcols)  ! stratiform precip number flux at the bottom of a particular layer
        real(r8) precnumc_base(pcols)  ! convective precip number flux at the bottom of a particular layer
        real(r8) precnumx_base_old, precnumx_base_tmp, precnumx_base_new
        real(r8) resusp_c          ! aerosol mass re-suspension in a particular layer from convective rain
        real(r8) resusp_s          ! aerosol mass re-suspension in a particular layer from stratiform rain
  
        real(r8) resusp_x
        real(r8) resusp_c_sv(pcols)
        real(r8) resusp_s_sv(pcols)
        real(r8) scavabx_old, scavabx_tmp, scavabx_new
        real(r8) srcx
        real(r8) tmpa, tmpb
        real(r8) u_old, u_tmp
        real(r8) x_old, x_tmp, x_ratio
  
        real(r8) fracev,evrate,rrate,blr,wrate(4)
        real(r8) cair
        real(r8) evouts(pcols,pver,19)
        real(r8) routs(pcols,pver,19)
        real(r8) wouts(pcols,pver,19)
        real(r8) evoutc(pcols,pver,19)
        real(r8) routc(pcols,pver,19)
        real(r8) woutc(pcols,pver,19)
        real(r8) reff(19,3),reff2(19),reff3(19)
        integer  it,wtype(19)
  
        data reff(:,1) /1.0_r8, &
                        1.0_r8,1.0_r8,1.0_r8,1.0_r8,1.0_r8, &
                        1.0_r8,1.0_r8,1.0_r8,1.0_r8,1.0_r8, &
                        0.0_r8,1.0_r8,0.0_r8,1.0_r8, &
                        0.0_r8,1.0_r8,1.0_r8,0.8_r8/
        data reff(:,2) /0.0_r8, &
                        0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8, &
                        1.0_r8,1.0_r8,1.0_r8,1.0_r8,1.0_r8, &
                        1.0_r8,0.0_r8,1.0_r8,0.0_r8, &
                        0.0_r8,0.0_r8,0.0_r8,0.0_r8/
        data reff(:,3) /0.4_r8, &
                        0.4_r8,0.4_r8,0.4_r8,0.4_r8,0.4_r8, &
                        1.0_r8,1.0_r8,1.0_r8,1.0_r8,1.0_r8, &
                        0.5_r8,0.4_r8,0.0_r8,0.4_r8, &
                        0.0_r8,0.4_r8,0.4_r8,0.36_r8/
  
        data wtype /2, 2,2,3,3,3, 1,1,3,3,3, &
                    1,2,1,2, 4,2,2,2/
  
  ! ------------------------------------------------------------------------
  !      omsm = 1.-1.e-10          ! used to prevent roundoff errors below zero
        omsm = 1._r8-2*epsilon(1._r8) ! used to prevent roundoff errors below zero
        precmin =  0.1_r8/8.64e4_r8      ! set critical value to 0.1 mm/day in kg/m2/s
  
        adjfac = deltat/(max(deltat,cmftau)) ! adjustment factor from hack scheme
  
        scavt=0.0_r8
        do i = 1,ncol
           precabs(i) = 0.0_r8
           precabc(i) = 0.0_r8
           cldmabs(i) = 0.0_r8
           cldmabc(i) = 0.0_r8
        enddo
  
        evouts=0.0_r8
        routs=0.0_r8
        wouts=0.0_r8
        evoutc=0.0_r8
        routc=0.0_r8
        woutc=0.0_r8
  
        do i=1,ncol
  
          scavab(:) = 0.0_r8
          scavabc(:) = 0.0_r8
  
        do k=1,pver
          cair=p(i,k)/(t(i,k)*287.05_r8)*(zi(i,k)-zi(i,k+1))
  
          if(evaps(i,k)>1.0e-20_r8)then
            fracev = max(0._r8,min(1._r8,evaps(i,k)*pdel(i,k)/gravit &
                    /max(1.e-12_r8,precabs(i))))
            evrate = 0.5_r8*fracev
            do m=1,19
              evouts(i,k,m)=scavab(m)*evrate/cair
              scavab(m)=scavab(m)*(1.0_r8-evrate)
            enddo
          endif
  
          if(precs(i,k)>1.0e-20_r8)then
            fracp = precs(i,k)*deltat/max(cwat(i,k)+precs(i,k)*deltat,1.e-12_r8)
            rrate = max(0.0_r8,(1._r8-exp(-fracp)))*max(0._r8,(cldt(i,k)-cldc(i,k)))
            if(t(i,k)>258.0_r8)then
              it=1
            else if(t(i,k)>237.0_r8)then
              it=2
            else
              it=3
            endif
            do m=1,19
              routs(i,k,m)=tracer(i,k,m)*rrate*reff(m,it)
              scavab(m)=scavab(m)+routs(i,k,m)*cair
            enddo
          endif
  
          cldmabs(i) = cldvst(i,k)
  
          precabs(i) = precabs(i) + (precs(i,k)- evaps(i,k))*pdel(i,k)/gravit
          if(precabs(i)*cldmabs(i)>1.0e-20_r8)then
            blr=precabs(i)/max(cldmabs(i),1.e-5_r8)
            if(t(i,k)>268.0_r8)then
              wrate(1)=cldmabs(i)*(1._r8-exp(-5.0e-7_r8*(blr*3.6e4_r8 )**0.7_r8*deltat))
              wrate(2)=cldmabs(i)*(1._r8-exp(-1.0e-5_r8*(blr*3.6e4_r8 )**0.7_r8*deltat))
              wrate(3)=cldmabs(i)*(1._r8-exp(-2.0e-4_r8*(blr*3.6e4_r8 )**0.85_r8*deltat))
            else if(t(i,k)>248.0_r8)then
              wrate(1)=cldmabs(i)*(1._r8-exp(-1.0e-5_r8*(blr*3.6e4_r8 )**0.66_r8*deltat))
              wrate(2)=cldmabs(i)*(1._r8-exp(-2.0e-4_r8*(blr*3.6e4_r8 )**0.66_r8*deltat))
              wrate(3)=cldmabs(i)*(1._r8-exp(-2.0e-3_r8*(blr*3.6e4_r8 )**0.7_r8*deltat))
            else
              wrate(1)=cldmabs(i)*(1._r8-exp(-2.0e-6_r8*(blr*3.6e4_r8 )**0.66_r8*deltat))
              wrate(2)=cldmabs(i)*(1._r8-exp(-4.0e-5_r8*(blr*3.6e4_r8 )**0.66_r8*deltat))
              wrate(3)=cldmabs(i)*(1._r8-exp(-4.0e-4_r8*(blr*3.6e4_r8 )**0.7_r8*deltat))
            endif
            wrate(4)=0.0_r8
            do m=1,19
              wouts(i,k,m)=tracer(i,k,m)*wrate(wtype(m))
              scavab(m)=scavab(m)+wouts(i,k,m)*cair
            enddo
          endif
  
          if(evapc(i,k)>1.0e-20_r8)then
            fracev = max(0._r8,min(1._r8,evapc(i,k)*pdel(i,k)/gravit &
                    /max(1.e-12_r8,precabc(i))))
            evrate = 0.5_r8*fracev
            do m=1,19
              evoutc(i,k,m)=scavabc(m)*evrate/cair
              scavabc(m)=scavabc(m)*(1.0_r8-evrate)
            enddo
          endif
  
          if(cmfdqr(i,k)>1.0e-20_r8)then
            fracp = cmfdqr(i,k)*deltat/max(cldc(i,k)*conicw(i,k)+(cmfdqr(i,k)+dlf(i,k))*deltat,1.e-12_r8)
            rrate = max(0.0_r8,(1._r8-exp(-fracp)))*max(0._r8,cldc(i,k))
            if(t(i,k)>258.0_r8)then
              it=1
            else if(t(i,k)>237.0_r8)then
              it=2
            else
              it=3
            endif
            do m=1,19
              routc(i,k,m)=tracer(i,k,m)*rrate*reff(m,it)
              scavabc(m)=scavabc(m)+routc(i,k,m)*cair
            enddo
          endif
  
          cldmabc(i) = cldvcu(i,k)
  
          precabc(i) = precabc(i) + (cmfdqr(i,k) - evapc(i,k))*pdel(i,k)/gravit
          if(precabc(i)>1.0e-20_r8)then
            blr=precabc(i)/max(cldmabc(i),1.e-5_r8)
            if(t(i,k)>268.0_r8)then
              wrate(1)=cldmabc(i)*(1._r8-exp(-5.0e-7_r8*(blr*3.6e4_r8 )**0.7_r8*deltat))
              wrate(2)=cldmabc(i)*(1._r8-exp(-1.0e-5_r8*(blr*3.6e4_r8 )**0.7_r8*deltat))
              wrate(3)=cldmabc(i)*(1._r8-exp(-2.0e-4_r8*(blr*3.6e4_r8 )**0.85_r8*deltat))
            else if(t(i,k)>248.0_r8)then
              wrate(1)=cldmabc(i)*(1._r8-exp(-1.0e-5_r8*(blr*3.6e4_r8 )**0.66_r8*deltat))
              wrate(2)=cldmabc(i)*(1._r8-exp(-2.0e-4_r8*(blr*3.6e4_r8 )**0.66_r8*deltat))
              wrate(3)=cldmabc(i)*(1._r8-exp(-2.0e-3_r8*(blr*3.6e4_r8 )**0.7_r8*deltat))
            else
              wrate(1)=cldmabc(i)*(1._r8-exp(-2.0e-6_r8*(blr*3.6e4_r8 )**0.66_r8*deltat))
              wrate(2)=cldmabc(i)*(1._r8-exp(-4.0e-5_r8*(blr*3.6e4_r8 )**0.66_r8*deltat))
              wrate(3)=cldmabc(i)*(1._r8-exp(-4.0e-4_r8*(blr*3.6e4_r8 )**0.7_r8*deltat))
            endif
            do m=1,19
              woutc(i,k,m)=tracer(i,k,m)*wrate(wtype(m))
              scavabc(m)=scavabc(m)+woutc(i,k,m)*cair
            enddo
          endif
  
          do m=1,19
            scavt(i,k,m)=max((-tracer(i,k,m)+1.e-20_r8), &
                        (evouts(i,k,m)-routs(i,k,m)-wouts(i,k,m)+ &
                         evoutc(i,k,m)-routc(i,k,m)-woutc(i,k,m)))/deltat
          enddo
        enddo
  
        enddo
  
     end subroutine gc_wetdep
  
  !==============================================================================
        function flux_precnum_vs_flux_prec_mpln( flux_prec, jstrcnv )
        real(r8) :: flux_precnum_vs_flux_prec_mpln
        real(r8), intent(in) :: flux_prec
        integer,  intent(in) :: jstrcnv
  
        if (jstrcnv <= 1) then
           flux_precnum_vs_flux_prec_mpln = flux_precnum_vs_flux_prec_mp( flux_prec )
        else
           flux_precnum_vs_flux_prec_mpln = flux_precnum_vs_flux_prec_ln( flux_prec )
        end if
  
        return
        end function flux_precnum_vs_flux_prec_mpln
  
  
  !==============================================================================
        function faer_resusp_vs_fprec_evap_mpln( fprec_evap, jstrcnv )
        real(r8) :: faer_resusp_vs_fprec_evap_mpln
        real(r8), intent(in) :: fprec_evap
        integer,  intent(in) :: jstrcnv
  
        if (jstrcnv <= 1) then
           faer_resusp_vs_fprec_evap_mpln = faer_resusp_vs_fprec_evap_mp( fprec_evap )
        else
           faer_resusp_vs_fprec_evap_mpln = faer_resusp_vs_fprec_evap_ln( fprec_evap )
        end if
  
        return
        end function faer_resusp_vs_fprec_evap_mpln
  
  
  !==============================================================================
        function fprecn_resusp_vs_fprec_evap_mpln( fprec_evap, jstrcnv )
        real(r8) :: fprecn_resusp_vs_fprec_evap_mpln
        real(r8), intent(in) :: fprec_evap
        integer,  intent(in) :: jstrcnv
  
        if (jstrcnv <= 1) then
           fprecn_resusp_vs_fprec_evap_mpln = fprecn_resusp_vs_fprec_evap_mp( fprec_evap )
        else
           fprecn_resusp_vs_fprec_evap_mpln = fprecn_resusp_vs_fprec_evap_ln( fprec_evap )
        end if
  
        return
        end function fprecn_resusp_vs_fprec_evap_mpln
  
  
  !==============================================================================
        function flux_precnum_vs_flux_prec_mp( flux_prec )
  !
  !  flux_prec = precipitation mass flux at the cloud base (kg/m^2/s)
  !  flux_precnum_vs_flux_prec_mp = precipitation number flux
  !     at the cloud base (drops/m^2/s), assuming marshall-palmer raindrop size distribution
  !
  !
        real(r8) :: flux_precnum_vs_flux_prec_mp
        real(r8), intent(in) :: flux_prec
  
        real(r8), parameter :: a0 =  1.0885896550304022E+01_r8
        real(r8), parameter :: a1 =  4.3660645528167907E-01_r8
  
        real(r8) :: x, y
     
        if (flux_prec >= 1.0e-36_r8) then
           x = log( flux_prec )
           y = exp( a0 + a1*x )    
        else
           y = 0.0_r8
        end if
        flux_precnum_vs_flux_prec_mp = y
  
        return
        end function flux_precnum_vs_flux_prec_mp
  
  
  !==============================================================================
        function flux_precnum_vs_flux_prec_ln( flux_prec )
  !
  !  flux_prec = precipitation mass flux at the cloud base (kg/m^2/s)
  !  flux_precnum_vs_flux_prec_ln = precipitation number flux
  !     at the cloud base (drops/m^2/s), assuming log-normal raindrop size distribution
  !
  !
        real(r8) :: flux_precnum_vs_flux_prec_ln
        real(r8), intent(in) :: flux_prec
  
        real(r8), parameter :: a0 =  9.9067806476181524E+00_r8
        real(r8), parameter :: a1 =  4.2690709912134056E-01_r8
  
        real(r8) :: x, y
     
        if (flux_prec >= 1.0e-36_r8) then
           x = log( flux_prec )
           y = exp( a0 + a1*x )    
        else
           y = 0.0_r8
        end if
        flux_precnum_vs_flux_prec_ln = y
  
        return
        end function flux_precnum_vs_flux_prec_ln
  
  
  !==============================================================================
        function faer_resusp_vs_fprec_evap_mp( fprec_evap )
  !
  !  fprec_evap = fraction of precipitation flux that has evaporated (below cloud base)
  !  faer_resusp_vs_fprec_evap_mp = corresponding fraction of precipitation-borne aerosol
  !     flux that is resuspended, assuming marshall-palmer raindrop size distribution
  !
  !  note that these fractions are relative to the cloud-base fluxes,
  !      and not to the layer immediately above fluxes
  !
        real(r8) :: faer_resusp_vs_fprec_evap_mp
        real(r8), intent(in) :: fprec_evap
  
        real(r8), parameter :: a01 =  8.6591133737322856E-02_r8
        real(r8), parameter :: a02 = -1.7389168499601941E+00_r8
        real(r8), parameter :: a03 =  2.7401882373663732E+01_r8
        real(r8), parameter :: a04 = -1.5861714653209464E+02_r8
        real(r8), parameter :: a05 =  5.1338179363011193E+02_r8
        real(r8), parameter :: a06 = -9.6835933124501412E+02_r8
        real(r8), parameter :: a07 =  1.0588489932213311E+03_r8
        real(r8), parameter :: a08 = -6.2184513459217271E+02_r8
        real(r8), parameter :: a09 =  1.5184126886039758E+02_r8
        real(r8), parameter :: x_lox_lin =  5.0000000000000003E-02_r8
        real(r8), parameter :: y_lox_lin =  2.5622471203221014E-03_r8
  
        real(r8) :: x, y
  
        x = max( 0.0_r8, min( 1.0_r8, fprec_evap ) )
        if (x < x_lox_lin) then
           y = y_lox_lin * (x/x_lox_lin)
        else
           y = x*( a01 + x*( a02 + x*( a03 + x*( a04 + x*( a05 &
             + x*( a06 + x*( a07 + x*( a08 + x*a09 ))))))))
        end if
        faer_resusp_vs_fprec_evap_mp = y
  
        return
        end function faer_resusp_vs_fprec_evap_mp
  
  
  !==============================================================================
        function faer_resusp_vs_fprec_evap_ln( fprec_evap )
  !
  !  fprec_evap = fraction of precipitation flux that has evaporated (below cloud base)
  !  faer_resusp_vs_fprec_evap_ln = corresponding fraction of precipitation-borne aerosol
  !     flux that is resuspended, assuming log-normal raindrop size distribution
  !
  !  note that these fractions are relative to the cloud-base fluxes,
  !      and not to the layer immediately above fluxes
  !
        real(r8) :: faer_resusp_vs_fprec_evap_ln
        real(r8), intent(in) :: fprec_evap
  
        real(r8), parameter :: a01 =  6.1944215103685640E-02_r8
        real(r8), parameter :: a02 = -2.0095166685965378E+00_r8
        real(r8), parameter :: a03 =  2.3882460251821236E+01_r8
        real(r8), parameter :: a04 = -1.2695611774753374E+02_r8
        real(r8), parameter :: a05 =  4.0086943562320101E+02_r8
        real(r8), parameter :: a06 = -7.4954272875943707E+02_r8
        real(r8), parameter :: a07 =  8.1701055892023624E+02_r8
        real(r8), parameter :: a08 = -4.7941894659538502E+02_r8
        real(r8), parameter :: a09 =  1.1710291076059025E+02_r8
        real(r8), parameter :: x_lox_lin =  1.0000000000000001E-01_r8
        real(r8), parameter :: y_lox_lin =  6.2227889828044350E-04_r8
  
        real(r8) :: x, y
  
        x = max( 0.0_r8, min( 1.0_r8, fprec_evap ) )
        if (x < x_lox_lin) then
           y = y_lox_lin * (x/x_lox_lin)
        else
           y = x*( a01 + x*( a02 + x*( a03 + x*( a04 + x*( a05 &
             + x*( a06 + x*( a07 + x*( a08 + x*a09 ))))))))
        end if
        faer_resusp_vs_fprec_evap_ln = y
  
        return
        end function faer_resusp_vs_fprec_evap_ln
  
  
  !==============================================================================
        function fprecn_resusp_vs_fprec_evap_mp( fprec_evap )
  !
  !  fprec_evap = fraction of precipitation flux that has evaporated (below cloud base)
  !  fprecn_resusp_vs_fprec_evap_mp = Rain number evaporation fraction, 
  !                                  assuming marshall-palmer raindrop size distribution
  !
  !  note that these fractions are relative to the cloud-base fluxes,
  !      and not to the layer immediately above fluxes
  !
        real(r8) :: fprecn_resusp_vs_fprec_evap_mp
        real(r8), intent(in) :: fprec_evap
  
        real(r8), parameter :: a01 =  4.5461070198414655E+00_r8
        real(r8), parameter :: a02 = -3.0381753620077529E+01_r8
        real(r8), parameter :: a03 =  1.7959619926085665E+02_r8
        real(r8), parameter :: a04 = -6.7152282193785618E+02_r8
        real(r8), parameter :: a05 =  1.5651931323557126E+03_r8
        real(r8), parameter :: a06 = -2.2743927701175126E+03_r8
        real(r8), parameter :: a07 =  2.0004645897056735E+03_r8
        real(r8), parameter :: a08 = -9.7351466279626209E+02_r8
        real(r8), parameter :: a09 =  2.0101198012962413E+02_r8
        real(r8), parameter :: x_lox_lin =  5.0000000000000003E-02_r8
        real(r8), parameter :: y_lox_lin =  1.7005858490684875E-01_r8
  
        real(r8) :: x, y
  
        x = max( 0.0_r8, min( 1.0_r8, fprec_evap ) )
        if (x < x_lox_lin) then
           y = y_lox_lin * (x/x_lox_lin)
        else
           y = x*( a01 + x*( a02 + x*( a03 + x*( a04 + x*( a05 &
             + x*( a06 + x*( a07 + x*( a08 + x*a09 ))))))))
        end if
        fprecn_resusp_vs_fprec_evap_mp = y
  
        return
        end function fprecn_resusp_vs_fprec_evap_mp
  
  
  !==============================================================================
        function fprecn_resusp_vs_fprec_evap_ln( fprec_evap )
  !
  !  fprec_evap = fraction of precipitation flux that has evaporated (below cloud base)
  !  fprecn_resusp_vs_fprec_evap_ln = Rain number evaporation fraction, 
  !                                  assuming log-normal raindrop size distribution
  !
  !  note that these fractions are relative to the cloud-base fluxes,
  !      and not to the layer immediately above fluxes
  !
        real(r8) :: fprecn_resusp_vs_fprec_evap_ln
        real(r8), intent(in) :: fprec_evap
  
        real(r8), parameter :: a01 = -5.2335291116884175E-02_r8
        real(r8), parameter :: a02 =  2.7203158069178226E+00_r8
        real(r8), parameter :: a03 =  9.4730878152409375E+00_r8
        real(r8), parameter :: a04 = -5.0573187592544798E+01_r8
        real(r8), parameter :: a05 =  9.4732631441282862E+01_r8
        real(r8), parameter :: a06 = -8.8265926556465814E+01_r8
        real(r8), parameter :: a07 =  3.5247835268269142E+01_r8
        real(r8), parameter :: a08 =  1.5404586576716444E+00_r8
        real(r8), parameter :: a09 = -3.8228795492549068E+00_r8
        real(r8), parameter :: x_lox_lin =  1.0000000000000001E-01_r8
        real(r8), parameter :: y_lox_lin =  2.7247994766566485E-02_r8
  
        real(r8) :: x, y
  
        x = max( 0.0_r8, min( 1.0_r8, fprec_evap ) )
        if (x < x_lox_lin) then
           y = y_lox_lin * (x/x_lox_lin)
        else
           y = x*( a01 + x*( a02 + x*( a03 + x*( a04 + x*( a05 &
             + x*( a06 + x*( a07 + x*( a08 + x*a09 ))))))))
        end if
        fprecn_resusp_vs_fprec_evap_ln = y
  
        return
        end function fprecn_resusp_vs_fprec_evap_ln
  
  !--------------------------------------------------------------------------------
  ! settling velocity
  !--------------------------------------------------------------------------------
  subroutine aerosol_depvel_compute( ncol, nlev, naero, t, pmid, ram1, fv, diam, stk_crc, dns_aer, &
                                     vlc_dry, vlc_trb, vlc_grv )
  
    use shr_kind_mod, only: r8 => shr_kind_r8
    use physconst,    only: pi, gravit, rair, boltz
  
    ! !ARGUMENTS:
    !
    implicit none
    !
    integer,  intent(in) :: ncol,nlev
    integer,  intent(in) :: naero
    real(r8), intent(in) :: t(:,:)          !atm temperature (K)
    real(r8), intent(in) :: pmid(:,:)       !atm pressure (Pa)
    real(r8), intent(in) :: fv(:)           !friction velocity (m/s)
    real(r8), intent(in) :: ram1(:)         !aerodynamical resistance (s/m)
    real(r8), intent(in) :: diam(:,:,:)
    real(r8), intent(in) :: stk_crc(:)
    real(r8), intent(in) :: dns_aer(:,:,:)
  
    real(r8), intent(out) :: vlc_trb(:,:)    !Turbulent deposn velocity (m/s)
    real(r8), intent(out) :: vlc_grv(:,:,:)  !grav deposn velocity (m/s)
    real(r8), intent(out) :: vlc_dry(:,:,:)  !dry deposn velocity (m/s)
  
    !------------------------------------------------------------------------
    ! Local Variables
    integer  :: m,i,k          !indices
    real(r8) :: vsc_dyn_atm(ncol,nlev)   ![kg m-1 s-1] Dynamic viscosity of air
    real(r8) :: vsc_knm_atm(ncol,nlev)   ![m2 s-1] Kinematic viscosity of atmosphere
    real(r8) :: shm_nbr_xpn   ![frc] Sfc-dep exponent for aerosol-diffusion dependence on Schmidt number
    real(r8) :: shm_nbr       ![frc] Schmidt number
    real(r8) :: stk_nbr       ![frc] Stokes number
    real(r8) :: mfp_atm(ncol,nlev)       ![m] Mean free path of air
    real(r8) :: dff_aer       ![m2 s-1] Brownian diffusivity of particle
    real(r8) :: rss_trb       ![s m-1] Resistance to turbulent deposition
    real(r8) :: slp_crc(ncol,nlev,naero) ![frc] Slip correction factor
    real(r8) :: rss_lmn(naero) ![s m-1] Quasi-laminar layer resistance
    real(r8) :: tmp          !temporary 
  
    ! constants
    real(r8),parameter::shm_nbr_xpn_lnd=-2._r8/3._r8 ![frc] shm_nbr_xpn over land
    real(r8),parameter::shm_nbr_xpn_ocn=-1._r8/2._r8 ![frc] shm_nbr_xpn over ccean
  
    real(r8) :: rho                 !atm density (kg/m**3)
  
    ! needs fv and ram1 passed in from lnd model
  
    !------------------------------------------------------------------------
  
    do k=1,nlev
       do i=1,ncol
          rho = pmid(i,k)/rair/t(i,k)
          ! from subroutine dst_dps_dry (consider adding sanity checks from line 212)
          ! when code asks to use midlayer density, pressure, temperature,
          ! I use the data coming in from the atmosphere, ie t(i,k), pmid(i,k)
  
          ! Quasi-laminar layer resistance: call rss_lmn_get
          ! Size-independent thermokinetic properties
          vsc_dyn_atm(i,k) = 1.72e-5_r8 * ((t(i,k)/273.0_r8)**1.5_r8) * 393.0_r8 / &
               (t(i,k)+120.0_r8)      ![kg m-1 s-1] RoY94 p. 102
          mfp_atm(i,k) = 2.0_r8 * vsc_dyn_atm(i,k) / &   ![m] SeP97 p. 455
               (pmid(i,k)*sqrt(8.0_r8/(pi*rair*t(i,k))))
          vsc_knm_atm(i,k) = vsc_dyn_atm(i,k) / rho ![m2 s-1] Kinematic viscosity of air
  
          do m = 1, naero
             slp_crc(i,k,m) = 1.0_r8 + 2.0_r8 * mfp_atm(i,k) * &
                  (1.257_r8+0.4_r8*exp(-1.1_r8*diam(i,k,m)/(2.0_r8*mfp_atm(i,k)))) / &
                  diam(i,k,m)   ![frc] Slip correction factor SeP97 p. 464
             vlc_grv(i,k,m) = (1.0_r8/18.0_r8) * diam(i,k,m) * diam(i,k,m) * dns_aer(i,k,m) * &
                  gravit * slp_crc(i,k,m) / vsc_dyn_atm(i,k) ![m s-1] Stokes' settling velocity SeP97 p. 466
             vlc_grv(i,k,m) = vlc_grv(i,k,m) * stk_crc(m)         ![m s-1] Correction to Stokes settling velocity
             vlc_dry(i,k,m)=vlc_grv(i,k,m)
          end do
  
       enddo
    enddo
    k=nlev  ! only look at bottom level for next part
    do m = 1, naero
       do i=1,ncol
          stk_nbr = vlc_grv(i,k,m) * fv(i) * fv(i) / (gravit*vsc_knm_atm(i,k))    ![frc] SeP97 p.965
          dff_aer = boltz * t(i,k) * slp_crc(i,k,m) / &    ![m2 s-1]
               (3.0_r8*pi*vsc_dyn_atm(i,k)*diam(i,k,m)) !SeP97 p.474
          shm_nbr = vsc_knm_atm(i,k) / dff_aer                        ![frc] SeP97 p.972
          shm_nbr_xpn = shm_nbr_xpn_lnd                          ![frc]
          !           if(ocnfrac.gt.0.5) shm_nbr_xpn=shm_nbr_xpn_ocn
          ! fxm: Turning this on dramatically reduces
          ! deposition velocity in low wind regimes
          ! Schmidt number exponent is -2/3 over solid surfaces and
          ! -1/2 over liquid surfaces SlS80 p. 1014
          ! if (oro(i)==0.0) shm_nbr_xpn=shm_nbr_xpn_ocn else shm_nbr_xpn=shm_nbr_xpn_lnd
          ! [frc] Surface-dependent exponent for aerosol-diffusion dependence on Schmidt # 
          tmp = shm_nbr**shm_nbr_xpn + 10.0_r8**(-3.0_r8/stk_nbr)
          rss_lmn(m) = 1.0_r8 / (tmp*fv(i)) ![s m-1] SeP97 p.972,965
  
          rss_trb = ram1(i) + rss_lmn(m) + ram1(i)*rss_lmn(m)*vlc_grv(i,k,m) ![s m-1]
          vlc_trb(i,m) = 1.0_r8 / rss_trb                            ![m s-1]
          vlc_dry(i,k,m) = vlc_trb(i,m)  +vlc_grv(i,k,m)
       end do !ncol
    end do
  
  end subroutine aerosol_depvel_compute
  
  !==================================================================================
  !BOP
  ! !IROUTINE: NIthermo
  
     subroutine NIthermo (ncol, km, klid, grav, pmid, tmpu, rh, &
                          SO4, NH3, NH4a, NO3an1, HNO3)
  
     use physconst,     only : rair,pi
  
  ! !USES:
     implicit NONE
  
  ! !INPUT PARAMETERS:
     integer, intent(in) :: ncol  ! total model ncols
     integer, intent(in) :: km    ! total model levels
     integer, intent(in) :: klid   ! index for pressure lid
     real*8, intent(in)    :: grav  ! gravity [m s-2]
     real*8, dimension(:,:), intent(in)  :: pmid   ! pressure [Pa]
     real*8, dimension(:,:), intent(in)  :: tmpu   ! Layer temperature [K]
     real*8, dimension(:,:), intent(in)  :: rh     ! relative humidity [0-1]
  
  ! !INOUTPUT PARAMETERS:
     real*8, dimension(:,:), intent(inout)  :: SO4    ! Sulphate aerosol [kg kg-1]
     real*8, dimension(:,:), intent(inout)  :: NH3    ! Ammonia (NH3, gas phase) [kg kg-1]
     real*8, dimension(:,:), intent(inout)  :: NO3an1 ! Nitrate size bin 001 [kg kg-1]
     real*8, dimension(:,:), intent(inout)  :: NH4a   ! Ammonium ion (NH4+, aerosol phase) [kg kg-1]
     real*8, dimension(:,:), intent(inout)  :: HNO3  ! buffer for NITRATE_HNO3 [kg kg-1]
  
  ! !DESCRIPTION: Prepares variables and calls the RPMARES (thermodynamics module)
  !
  ! !REVISION HISTORY:
  !
  ! Aug2020 E.Sherman - Refactored for process library
  !
  
  ! !Local Variables
     real   :: fmmr_to_conc, rhoa
     real(kind=r8) :: SO4_, GNO3, GNH3, RH_, TEMP, ASO4, AHSO4, AH2O, ANO3, ANH4
     integer :: k, j, i
  
     integer :: status
  
  !EOP
  !-------------------------------------------------------------------------
  !  Begin...
  
     do k = klid, km
      do i = 1, ncol
  
  !     if(tmpu(i,k)>258.0)then
        rhoa = pmid(i,k)/(rair*tmpu(i,k))
  !     Conversion of mass mixing ratio to concentration (ug m-3)
        fmmr_to_conc = 1.e9 * rhoa
  
  !     Unit conversion for input to thermodynamic module
  !     Per grid box call to RPMARES thermodynamic module
  !     We do not presently treat chemistry of sulfate completely,
  !     hence we ignore terms for ASO4, AHSO4, AH2O, and we do
  !     not update SO4 on output from RPMARES.
  !     At present we are importing HNO3 from offline file, so we
  !     do not update on return.
        SO4_  = max(1.d-32,SO4(i,k) * fmmr_to_conc)*96.0/29.0
        GNO3  = max(1.d-32,HNO3(i,k) * fmmr_to_conc)*63.0/29.0
        GNH3  = max(1.d-32,NH3(i,k)  * fmmr_to_conc)*17.0/29.0
        RH_   = rh(i,k)
        TEMP  = tmpu(i,k)
        ASO4  = 1.d-32
        AHSO4 = 1.d-32
        ANO3  = max(1.d-32,NO3an1(i,k) * fmmr_to_conc)*62.0/29.0
        AH2O  = 1.d-32
        ANH4  = max(1.d-32,NH4a(i,k) * fmmr_to_conc)*18.0/29.0
  
        call RPMARES (  SO4_, GNO3,  GNH3, RH_,  TEMP, &
                        ASO4, AHSO4, ANO3, AH2O, ANH4 )
  
  !     Unit conversion back on return from thermodynamic module
        NH3(i,k)    = max(1.d-32, GNH3 / fmmr_to_conc)*29.0/17.0
        NO3an1(i,k) = max(1.d-32, ANO3 / fmmr_to_conc)*29.0/62.0
        NH4a(i,k)   = max(1.d-32, ANH4 / fmmr_to_conc)*29.0/18.0
        HNO3(i,k) = max(1.d-32, GNO3 / fmmr_to_conc)*29.0/63.0
  
      enddo
     enddo
  
     end subroutine NIthermo
  
  !==================================================================================
  !BOP
  ! !IROUTINE: RPMARES
  
     subroutine RPMARES( SO4,  GNO3,  GNH3, RH,   TEMP, &
                         ASO4, AHSO4, ANO3, AH2O, ANH4 )
  
  ! !USES:
     implicit NONE
  
  ! !INPUT PARAMETERS:
     real(kind=r8) :: SO4              ! Total sulfate in micrograms / m**3
     real(kind=r8) :: GNO3             ! Gas-phase nitric acid in micrograms / m**3
     real(kind=r8) :: GNH3             ! Gas-phase ammonia in micrograms / m**3
     real(kind=r8) :: RH               ! Fractional relative humidity
     real(kind=r8) :: TEMP             ! Temperature in Kelvin
     real(kind=r8) :: ASO4             ! Aerosol sulfate in micrograms / m**3
     real(kind=r8) :: AHSO4            ! Aerosol bisulfate in micrograms / m**3
     real(kind=r8) :: ANO3             ! Aerosol nitrate in micrograms / m**3
     real(kind=r8) :: AH2O             ! Aerosol liquid water content water in
                                       !   micrograms / m**3
     real(kind=r8) :: ANH4             ! Aerosol ammonium in micrograms / m**3
  
  ! !DESCRIPTION:
  !   ARES calculates the chemical composition of a sulfate/nitrate/
  !   ammonium/water aerosol based on equilibrium thermodynamics.
  !
  !   This code considers two regimes depending upon the molar ratio
  !   of ammonium to sulfate.
  !
  !   For values of this ratio less than 2,the code solves a cubic for
  !   hydrogen ion molality, H+,  and if enough ammonium and liquid
  !   water are present calculates the dissolved nitric acid. For molal
  !   ionic strengths greater than 50, nitrate is assumed not to be present.
  !
  !   For values of the molar ratio of 2 or greater, all sulfate is assumed
  !   to be ammonium sulfate and a calculation is made for the presence of
  !   ammonium nitrate.
  !
  !   The Pitzer multicomponent approach is used in subroutine ACTCOF to
  !   obtain the activity coefficients. Abandoned -7/30/97 FSB
  !
  !   The Bromley method of calculating the multicomponent activity coefficients
  !    is used in this version 7/30/97 SJR/FSB
  !
  !   The calculation of liquid water
  !   is done in subroutine water. Details for both calculations are given
  !   in the respective subroutines.
  !
  !   Based upon MARS due to
  !   P. Saxena, A.B. Hudischewskyj, C. Seigneur, and J.H. Seinfeld,
  !   Atmos. Environ., vol. 20, Number 7, Pages 1471-1483, 1986.
  !
  !   and SCAPE due to
  !   Kim, Seinfeld, and Saxeena, Aerosol Sience and Technology,
  !   Vol 19, number 2, pages 157-181 and pages 182-198, 1993.
  !
  ! NOTE: All concentrations supplied to this subroutine are TOTAL
  !       over gas and aerosol phases
  
  !
  ! !REVISION HISTORY:
  !
  !      Who       When        Detailed description of changes
  !   ---------   --------  -------------------------------------------
  !   S.Roselle   11/10/87  Received the first version of the MARS code
  !   S.Roselle   12/30/87  Restructured code
  !   S.Roselle   2/12/88   Made correction to compute liquid-phase
  !                         concentration of H2O2.
  !   S.Roselle   5/26/88   Made correction as advised by SAI, for
  !                         computing H+ concentration.
  !   S.Roselle   3/1/89    Modified to operate with EM2
  !   S.Roselle   5/19/89   Changed the maximum ionic strength from
  !                         100 to 20, for numerical stability.
  !   F.Binkowski 3/3/91    Incorporate new method for ammonia rich case
  !                         using equations for nitrate budget.
  !   F.Binkowski 6/18/91   New ammonia poor case which
  !                         omits letovicite.
  !   F.Binkowski 7/25/91   Rearranged entire code, restructured
  !                         ammonia poor case.
  !   F.Binkowski 9/9/91    Reconciled all cases of ASO4 to be output
  !                         as SO4--
  !   F.Binkowski 12/6/91   Changed the ammonia defficient case so that
  !                         there is only neutralized sulfate (ammonium
  !                         sulfate) and sulfuric acid.
  !   F.Binkowski 3/5/92    Set RH bound on AWAS to 37 % to be in agreement
  !                          with the Cohen et al. (1987)  maximum molality
  !                          of 36.2 in Table III.( J. Phys Chem (91) page
  !                          4569, and Table IV p 4587.)
  !   F.Binkowski 3/9/92    Redid logic for ammonia defficient case to remove
  !                         possibility for denomenator becoming zero;
  !                         this involved solving for H+ first.
  !                         Note that for a relative humidity
  !                          less than 50%, the model assumes that there is no
  !                          aerosol nitrate.
  !   F.Binkowski 4/17/95   Code renamed  ARES (AeRosol Equilibrium System)
  !                          Redid logic as follows
  !                         1. Water algorithm now follows Spann & Richardson
  !                         2. Pitzer Multicomponent method used
  !                         3. Multicomponent practical osmotic coefficient
  !                            use to close iterations.
  !                         4. The model now assumes that for a water
  !                            mass fraction WFRAC less than 50% there is
  !                            no aerosol nitrate.
  !   F.Binkowski 7/20/95   Changed how nitrate is calculated in ammonia poor
  !                         case, and changed the WFRAC criterion to 40%.
  !                         For ammonium to sulfate ratio less than 1.0
  !                         all ammonium is aerosol and no nitrate aerosol
  !                         exists.
  !   F.Binkowski 7/21/95   Changed ammonia-ammonium in ammonia poor case to
  !                         allow gas-phase ammonia to exist.
  !   F.Binkowski 7/26/95   Changed equilibrium constants to values from
  !                         Kim et al. (1993)
  !   F.Binkowski 6/27/96   Changed to new water format
  !   F.Binkowski 7/30/97   Changed to Bromley method for multicomponent
  !                         activity coefficients. The binary activity
  !                         coefficients
  !                         are the same as the previous version
  !   F.Binkowski 8/1/97    Changed minimum sulfate from 0.0 to 1.0e-6 i.e.
  !                         1 picogram per cubic meter
  !   F.Binkowski 2/23/98   Changes to code made by Ingmar Ackermann to
  !                         deal with precision problems on workstations
  !                         incorporated in to this version.  Also included
  !                         are his improved descriptions of variables.
  !  F. Binkowski 8/28/98   changed logic as follows:
  !                         If iterations fail, initial values of nitrate
  !                          are retained.
  !                         Total mass budgets are changed to account for gas
  !                         phase returned.
  !  F.Binkowski 10/01/98   Removed setting RATIO to 5 for low to
  !                         to zero sulfate sulfate case.
  !  F.Binkowski 01/10/2000 reconcile versions
  !
  !  F.Binkowski 05/17/2000 change to logic for calculating RATIO
  !  F.Binkowski 04/09/2001 change for very low values of RATIO,
  !                         RATIO < 0.5, no iterative calculations are done
  !                         in low ammonia case a MAX(1.0e-10, MSO4) IS
  !                         applied, and the iteration count is
  !                         reduced to fifty for each iteration loop.
  !  R. Yantosca 09/25/2002 Bundled into "rpmares_mod.f".  Declared all REALs
  !                          as REAL*8's.  Cleaned up comments.  Also now force
  !                          double precision explicitly with "D" exponents.
  !  P. Le Sager and        Bug fix for low ammonia case -- prevent floating
  !  R. Yantosca 04/10/2008  point underflow and NaN's.
  !  S. Steenrod 04/15/2010 Modified to include into GMI model
  !  E. Sherman  08/06/2020 Moved to GOCART2G process library
  
  ! !Local Variables
    !=================================================================
    ! PARAMETERS and their descriptions:
    !=================================================================
    ! Molecular weights
     real(kind=r8), PARAMETER :: MWNO3  = 62.0049d0                ! NO3
     real(kind=r8), PARAMETER :: MWHNO3 = 63.01287d0               ! HNO3
     real(kind=r8), PARAMETER :: MWSO4  = 96.0576d0                ! SO4
     real(kind=r8), PARAMETER :: MWNH3  = 17.03061d0               ! NH3
     real(kind=r8), PARAMETER :: MWNH4  = 18.03858d0               ! NH4
  
     ! Minimum value of sulfate aerosol concentration
     real(kind=r8), PARAMETER :: MINSO4 = 1.0d-6 / MWSO4
  
     ! Minimum total nitrate cncentration
     real(kind=r8), PARAMETER :: MINNO3 = 1.0d-6 / MWNO3
  
     ! Force a minimum concentration
     real(kind=r8), PARAMETER :: FLOOR  = 1.0d-30
  
     ! Tolerances for convergence test.  NOTE: We now have made these
     ! parameters so they don't lose their values (phs, bmy, 4/10/08)
     real(kind=r8), PARAMETER :: TOLER1 = 0.00001d0
     real(kind=r8), PARAMETER :: TOLER2 = 0.001d0
  
     ! Limit to test for zero ionic activity (phs, bmy, 4/10/08)
     real(kind=r8), PARAMETER :: EPS    = 1.0d-30
  
     !=================================================================
     ! SCRATCH LOCAL VARIABLES and their descriptions:
     !=================================================================
  
     INTEGER :: IRH              ! Index set to percent relative humidity
     INTEGER :: NITR             ! Number of iterations for activity
                                 !   coefficients
     INTEGER :: NNN              ! Loop index for iterations
     INTEGER :: NR               ! Number of roots to cubic equation for
                                 ! H+ ciaprecision
     real(kind=r8)  :: A0        ! Coefficients and roots of
     real(kind=r8)  :: A1        ! Coefficients and roots of
     real(kind=r8)  :: A2        ! Coefficients and roots of
     REAL    :: AA               ! Coefficients and discriminant for
                                 ! quadratic equation for ammonium nitrate
     real(kind=r8)  :: BAL       ! internal variables ( high ammonia case)
     real(kind=r8)  :: BB        ! Coefficients and discriminant for
                                 !   quadratic equation for ammonium nitrate
     real(kind=r8)  :: BHAT      ! Variables used for ammonia solubility
     real(kind=r8)  :: CC        ! Coefficients and discriminant for
                                 !   quadratic equation for ammonium nitrate
     real(kind=r8)  :: CONVT     ! Factor for conversion of units
     real(kind=r8)  :: DD        ! Coefficients and discriminant for
                                 !   quadratic equation for ammonium nitrate
     real(kind=r8)  :: DISC      ! Coefficients and discriminant for
                                 !   quadratic equation for ammonium nitrate
     real(kind=r8)  :: EROR      ! Relative error used for convergence test
     real(kind=r8)  :: FNH3      ! "Free ammonia concentration", that
                                 !   which exceeds TWOSO4
     real(kind=r8)  :: GAMAAB    ! Activity Coefficient for (NH4+,
                                 !   HSO4-)GAMS( 2,3 )
     real(kind=r8)  :: GAMAAN    ! Activity coefficient for (NH4+, NO3-)
                                 !   GAMS( 2,2 )
     real(kind=r8)  :: GAMAHAT   ! Variables used for ammonia solubility
     real(kind=r8)  :: GAMANA    ! Activity coefficient for (H+ ,NO3-)
                                 !   GAMS( 1,2 )
     real(kind=r8)  :: GAMAS1    ! Activity coefficient for (2H+, SO4--)
                                 !   GAMS( 1,1 )
     real(kind=r8)  :: GAMAS2    ! Activity coefficient for (H+, HSO4-)
                                 !   GAMS( 1,3 )
     real(kind=r8)  :: GAMOLD    ! used for convergence of iteration
     real(kind=r8)  :: GASQD     ! internal variables ( high ammonia case)
     real(kind=r8)  :: HPLUS     ! Hydrogen ion (low ammonia case) (moles
                                 !   / kg water)
     real(kind=r8)  :: K1A       ! Equilibrium constant for ammonia to
                                 !   ammonium
     real(kind=r8)  :: K2SA      ! Equilibrium constant for
                                 !   sulfate-bisulfate (aqueous)
     real(kind=r8)  :: K3        ! Dissociation constant for ammonium
                                 !   nitrate
     real(kind=r8)  :: KAN       ! Equilibrium constant for ammonium
                                 !   nitrate (aqueous)
     real(kind=r8)  :: KHAT      ! Variables used for ammonia solubility
     real(kind=r8)  :: KNA       ! Equilibrium constant for nitric acid
                                 !   (aqueous)
     real(kind=r8)  :: KPH       ! Henry's Law Constant for ammonia
     real(kind=r8)  :: KW        ! Equilibrium constant for water
                                 !  dissociation
     real(kind=r8)  :: KW2       ! Internal variable using KAN
     real(kind=r8)  :: MAN       ! Nitrate (high ammonia case) (moles /
                                 !   kg water)
     real(kind=r8)  :: MAS       ! Sulfate (high ammonia case) (moles /
                                 !   kg water)
     real(kind=r8)  :: MHSO4     ! Bisulfate (low ammonia case) (moles /
                                 !   kg water)
     real(kind=r8)  :: MNA       ! Nitrate (low ammonia case) (moles / kg
                                 !   water)
     real(kind=r8)  :: MNH4      ! Ammonium (moles / kg water)
     real(kind=r8)  :: MOLNU     ! Total number of moles of all ions
     real(kind=r8)  :: MSO4      ! Sulfate (low ammonia case) (moles / kg
                                 !   water)
     real(kind=r8)  :: PHIBAR    ! Practical osmotic coefficient
     real(kind=r8)  :: PHIOLD    ! Previous value of practical osmotic
                                 !   coefficient used for convergence of
                                 !   iteration
     real(kind=r8)  :: RATIO     ! Molar ratio of ammonium to sulfate
     real(kind=r8)  :: RK2SA     ! Internal variable using K2SA
     real(kind=r8)  :: RKNA      ! Internal variables using KNA
     real(kind=r8)  :: RKNWET    ! Internal variables using KNA
     real(kind=r8)  :: RR1
     real(kind=r8)  :: RR2
     real(kind=r8)  :: STION     ! Ionic strength
     real(kind=r8)  :: T1        ! Internal variables for temperature
                                 !   corrections
     real(kind=r8)  :: T2        ! Internal variables for temperature
                                 !   corrections
     real(kind=r8)  :: T21       ! Internal variables of convenience (low
                                 !   ammonia case)
     real(kind=r8)  :: T221      ! Internal variables of convenience (low
                                 !   ammonia case)
     real(kind=r8)  :: T3        ! Internal variables for temperature
                                 !   corrections
     real(kind=r8)  :: T4        ! Internal variables for temperature
                                 !   corrections
     real(kind=r8)  :: T6        ! Internal variables for temperature
                                 !   corrections
     real(kind=r8)  :: TNH4      ! Total ammonia and ammonium in
                                 !   micromoles / meter ** 3
     real(kind=r8)  :: TNO3      ! Total nitrate in micromoles / meter ** 3
     !-----------------------------------------------------------------------
     ! Prior to 4/10/08:
     ! Now make these PARAMETERS instead of variables (bmy, 4/10/08)
     !real(kind=r8)  :: TOLER1           ! Tolerances for convergence test
     !real(kind=r8)  :: TOLER2           ! Tolerances for convergence test
     !-----------------------------------------------------------------------
     real(kind=r8)  :: TSO4      ! Total sulfate in micromoles / meter ** 3
     real(kind=r8)  :: TWOSO4    ! 2.0 * TSO4  (high ammonia case) (moles
                                 !   / kg water)
     real(kind=r8)  :: WFRAC     ! Water mass fraction
     real(kind=r8)  :: WH2O      ! Aerosol liquid water content (internally)
                                 !   micrograms / meter **3 on output
                                 !   internally it is 10 ** (-6) kg (water)
                                 !   / meter ** 3
                                 !   the conversion factor (1000 g = 1 kg)
                                 !   is applied for AH2O output
     real(kind=r8)  :: WSQD      ! internal variables ( high ammonia case)
     real(kind=r8)  :: XNO3      ! Nitrate aerosol concentration in
                                 ! micromoles / meter ** 3
     real(kind=r8)  :: XXQ       ! Variable used in quadratic solution
     real(kind=r8)  :: YNH4      ! Ammonium aerosol concentration in
                                 !  micromoles / meter** 3
     real(kind=r8)  :: ZH2O      ! Water variable saved in case ionic
                                 !  strength too high.
     real(kind=r8)  :: ZSO4      ! Total sulfate molality - mso4 + mhso4
                                 !  (low ammonia case) (moles / kg water)
     real(kind=r8)  :: CAT( 2 )  ! Array for cations (1, H+); (2, NH4+)
                                 !  (moles / kg water)
     real(kind=r8)  :: AN ( 3 )  ! Array for anions (1, SO4--); (2,
                                 !   NO3-); (3, HSO4-)  (moles / kg water)
     real(kind=r8)  :: CRUTES( 3 )      ! Coefficients and roots of
     real(kind=r8)  :: GAMS( 2, 3 )     ! Array of activity coefficients
     real(kind=r8)  :: TMASSHNO3        ! Total nitrate (vapor and particle)
     real(kind=r8)  :: GNO3_IN, ANO3_IN
     character (len=75) :: err_msg
  
     integer :: status
  
  !EOP
  !-------------------------------------------------------------------------
  !  Begin...
  
        ! For extremely low relative humidity ( less than 1% ) set the
        ! water content to a minimum and skip the calculation.
        IF ( RH .LT. 0.01 ) THEN
           AH2O = FLOOR
           RETURN
        ENDIF
  
        ! total sulfate concentration
        TSO4 = MAX( FLOOR, SO4 / MWSO4  )
        ASO4 = SO4
  
        !Cia models3 merge NH3/NH4 , HNO3,NO3 here
        !c *** recommended by Dr. Ingmar Ackermann
  
        ! total nitrate
        TNO3      = MAX( 0.0d0, ( ANO3 / MWNO3 + GNO3 / MWHNO3 ) )
  
        ! total ammonia
        TNH4      = MAX( 0.0d0, ( GNH3 / MWNH3 + ANH4 / MWNH4 )  )
  
        GNO3_IN   = GNO3
        ANO3_IN   = ANO3
        TMASSHNO3 = MAX( 0.0d0, GNO3 + ANO3 )
  
        ! set the  molar ratio of ammonium to sulfate
        RATIO = TNH4 / TSO4
  
        ! validity check for negative concentration
        IF ( TSO4 < 0.0d0 .OR. TNO3 < 0.0d0 .OR. TNH4 < 0.0d0 ) THEN
            PRINT*, 'TSO4 : ', TSO4
            PRINT*, 'TNO3 : ', TNO3
            PRINT*, 'TNH4 : ', TNH4
        ENDIF
  
        ! now set humidity index IRH as a percent
        IRH = NINT( 100.0 * RH )
  
        ! now set humidity index IRH as a percent
        IRH = MAX(  1, IRH )
        IRH = MIN( 99, IRH )
  
        !=================================================================
        ! Specify the equilibrium constants at  correct temperature.
        ! Also change units from ATM to MICROMOLE/M**3 (for KAN, KPH, and K3 )
        ! Values from Kim et al. (1993) except as noted.
        ! Equilibrium constant in Kim et al. (1993)
        !   K = K0 exp[ a(T0/T -1) + b(1+log(T0/T)-T0/T) ], T0 = 298.15 K
        !   K = K0 EXP[ a T3 + b T4 ] in the code here.
        !=================================================================
        CONVT = 1.0d0 / ( 0.082d0 * TEMP )
        T6    = 0.082d-9 *  TEMP
        T1    = 298.0d0 / TEMP
        T2    = LOG( T1 )
        T3    = T1 - 1.0d0
        T4    = 1.0d0 + T2 - T1
  
        !=================================================================
        ! Equilibrium Relation
        !
        ! HSO4-(aq)         = H+(aq)   + SO4--(aq)  ; K2SA
        ! NH3(g)            = NH3(aq)               ; KPH
        ! NH3(aq) + H2O(aq) = NH4+(aq) + OH-(aq)    ; K1A
        ! HNO3(g)           = H+(aq)   + NO3-(aq)   ; KNA
        ! NH3(g) + HNO3(g)  = NH4NO3(s)             ; K3
        ! H2O(aq)           = H+(aq)   + OH-(aq)    ; KW
        !=================================================================
        KNA  = 2.511d+06 *  EXP(  29.17d0 * T3 + 16.83d0 * T4 ) * T6
        K1A  = 1.805d-05 *  EXP(  -1.50d0 * T3 + 26.92d0 * T4 )
        K2SA = 1.015d-02 *  EXP(   8.85d0 * T3 + 25.14d0 * T4 )
        KW   = 1.010d-14 *  EXP( -22.52d0 * T3 + 26.92d0 * T4 )
        KPH  = 57.639d0  *  EXP(  13.79d0 * T3 -  5.39d0 * T4 ) * T6
        !K3   =  5.746E-17 * EXP( -74.38 * T3 + 6.12  * T4 ) * T6 * T6
        KHAT =  KPH * K1A / KW
        KAN  =  KNA * KHAT
  
        ! Compute temperature dependent equilibrium constant for NH4NO3
        ! (from Mozurkewich, 1993)
        K3 = EXP( 118.87d0  - 24084.0d0 / TEMP -  6.025d0  * LOG( TEMP ) )
  
        ! Convert to (micromoles/m**3) **2
        K3     = K3 * CONVT * CONVT
  
        WH2O   = 0.0d0
        STION  = 0.0d0
  !.sds      AH2O   = 0.0d0
        AH2O   = FLOOR
  
        MAS    = 0.0d0
        MAN    = 0.0d0
        HPLUS  = 0.0d0
        !--------------------------------------------------------------
        ! Prior to 4/10/08:
        ! Now make these parameters so that they won't lose their
        ! values. (phs, bmy, 4/10/08)
        !TOLER1 = 0.00001d0
        !TOLER2 = 0.001d0
        !--------------------------------------------------------------
        NITR   = 0
        NR     = 0
        GAMAAN = 1.0d0
        GAMOLD = 1.0d0
  
        ! If there is very little sulfate and  nitrate
        ! set concentrations to a very small value and return.
        IF ( ( TSO4 .LT. MINSO4 ) .AND. ( TNO3 .LT. MINNO3 ) ) THEN
           ASO4  = MAX( FLOOR, ASO4  )
           AHSO4 = MAX( FLOOR, AHSO4 ) ! [rjp, 12/12/01]
           ANO3  = MAX( FLOOR, ANO3  )
           ANH4  = MAX( FLOOR, ANH4  )
           WH2O  = FLOOR
           AH2O  = FLOOR
           GNH3  = MAX( FLOOR, GNH3  )
           GNO3  = MAX( FLOOR, GNO3  )
  
           RETURN
        ENDIF
        !=================================================================
        ! High Ammonia Case
        !=================================================================
        IF ( RATIO .GT. 2.0d0 ) THEN
  
           GAMAAN = 0.1d0
  
           ! Set up twice the sulfate for future use.
           TWOSO4 = 2.0d0 * TSO4
           XNO3   = 0.0d0
           YNH4   = TWOSO4
  
           ! Treat different regimes of relative humidity
           !
           ! ZSR relationship is used to set water levels. Units are
           !  10**(-6) kg water/ (cubic meter of air)
           !  start with ammomium sulfate solution without nitrate
  
           CALL AWATER( IRH, TSO4, YNH4, TNO3, AH2O ) !**** note TNO3
           WH2O = 1.0d-3 * AH2O
  
           ASO4 = TSO4   * MWSO4
           ! In sulfate poor case, Sulfate ion is preferred
           ! Set bisulfate equal to zero [rjp, 12/12/01]
           AHSO4 = 0.0d0
           ANO3  = 0.0d0
           ANH4  = YNH4 * MWNH4
           WFRAC = AH2O / ( ASO4 + ANH4 +  AH2O )
  
          !IF ( WFRAC .EQ. 0.0 )  RETURN   ! No water
          IF ( WFRAC .LT. 0.2d0 ) THEN
  
             ! "dry" ammonium sulfate and ammonium nitrate
             ! compute free ammonia
             FNH3 = TNH4 - TWOSO4
             CC   = TNO3 * FNH3 - K3
  
             ! check for not enough to support aerosol
             IF ( CC .LE. 0.0d0 ) THEN
                XNO3 = 0.0d0
             ELSE
                AA   = 1.0d0
                BB   = -( TNO3 + FNH3 )
                DISC = BB * BB - 4.0d0 * CC
  
                ! Check for complex roots of the quadratic
                ! set retain initial values of nitrate and RETURN
                ! if complex roots are found
                IF ( DISC .LT. 0.0d0 ) THEN
                   XNO3  = 0.0d0
                   AH2O  = 1000.0d0 * WH2O
                   YNH4  = TWOSO4
                   ASO4  = TSO4 * MWSO4
                   AHSO4 = 0.0d0
                   ANH4  = YNH4 * MWNH4
                   GNH3  = MWNH3 * MAX( FLOOR, ( TNH4 - YNH4 ) )
                   GNO3  = GNO3_IN
                   ANO3  = ANO3_IN
                   RETURN
                ENDIF
  
                ! to get here, BB .lt. 0.0, CC .gt. 0.0 always
                DD  = SQRT( DISC )
                XXQ = -0.5d0 * ( BB + SIGN ( 1.0d0, BB ) * DD )
  
  
                ! Since both roots are positive, select smaller root.
                XNO3 = MIN( XXQ / AA, CC / XXQ )
  
             ENDIF                ! CC .LE. 0.0
  
             AH2O  = 1000.0d0 * WH2O
             YNH4  = TWOSO4 + XNO3
             ASO4  = TSO4 * MWSO4
             AHSO4 = FLOOR
             ANO3  = XNO3 * MWNO3
             ANH4  = YNH4 * MWNH4
             GNH3  = MWNH3 * MAX( FLOOR, ( TNH4 - YNH4 )  )
             GNO3  = MAX( FLOOR, ( TMASSHNO3 - ANO3 ) )
             RETURN
          ENDIF                  ! WFRAC .LT. 0.2
  
          ! liquid phase containing completely neutralized sulfate and
          ! some nitrate.  Solve for composition and quantity.
          MAS    = TSO4 / WH2O
          MAN    = 0.0d0
          XNO3   = 0.0d0
          YNH4   = TWOSO4
          PHIOLD = 1.0d0
  
          !===============================================================
          ! Start loop for iteration
          !
          ! The assumption here is that all sulfate is ammonium sulfate,
          ! and is supersaturated at lower relative humidities.
          !===============================================================
          DO NNN = 1, 50 ! loop count reduced 0409/2001 by FSB
  
             NITR  = NNN
             GASQD = GAMAAN * GAMAAN
             WSQD  = WH2O * WH2O
             KW2   = KAN * WSQD / GASQD
             AA    = 1.0 - KW2
             BB    = TWOSO4 + KW2 * ( TNO3 + TNH4 - TWOSO4 )
             CC    = -KW2 * TNO3 * ( TNH4 - TWOSO4 )
  
             ! This is a quadratic for XNO3 [MICROMOLES / M**3]
             ! of nitrate in solution
             DISC = BB * BB - 4.0d0 * AA * CC
  
             ! Check for complex roots, retain inital values and RETURN
             IF ( DISC .LT. 0.0 ) THEN
                XNO3  = 0.0d0
                AH2O  = 1000.0d0 * WH2O
                YNH4  = TWOSO4
                ASO4  = TSO4 * MWSO4
                AHSO4 = FLOOR     ! [rjp, 12/12/01]
                ANH4  = YNH4 * MWNH4
                GNH3  = MWNH3 * MAX( FLOOR, (TNH4 - YNH4 ) )
                GNO3  = GNO3_IN
                ANO3  = ANO3_IN
                RETURN
             ENDIF
  
             ! Deal with degenerate case (yoj)
             IF ( AA .NE. 0.0d0 ) THEN
                DD  = SQRT( DISC )
                XXQ = -0.5d0 * ( BB + SIGN( 1.0d0, BB ) * DD )
                RR1 = XXQ / AA
                RR2 = CC / XXQ
  
                ! choose minimum positve root
                IF ( ( RR1 * RR2 ) .LT. 0.0d0 ) THEN
                   XNO3 = MAX( RR1, RR2 )
                ELSE
                   XNO3 = MIN( RR1, RR2 )
                ENDIF
             ELSE
                XNO3 = - CC / BB  ! AA equals zero here.
             ENDIF
  
             XNO3 = MIN( XNO3, TNO3 )
  
             ! This version assumes no solid sulfate forms (supersaturated )
             ! Now update water
             CALL AWATER ( IRH, TSO4, YNH4, XNO3, AH2O )
  
             ! ZSR relationship is used to set water levels. Units are
             ! 10**(-6) kg water/ (cubic meter of air).  The conversion
             ! from micromoles to moles is done by the units of WH2O.
             WH2O = 1.0d-3 * AH2O
  
             ! Ionic balance determines the ammonium in solution.
             MAN  = XNO3 / WH2O
             MAS  = TSO4 / WH2O
             MNH4 = 2.0d0 * MAS + MAN
             YNH4 = MNH4 * WH2O
  
             ! MAS, MAN and MNH4 are the aqueous concentrations of sulfate,
             ! nitrate, and ammonium in molal units (moles/(kg water) ).
             STION    = 3.0d0 * MAS + MAN
             CAT( 1 ) = 0.0d0
             CAT( 2 ) = MNH4
             AN ( 1 ) = MAS
             AN ( 2 ) = MAN
             AN ( 3 ) = 0.0d0
  !           CALL ACTCOF ( CAT, AN, GAMS, MOLNU, PHIBAR )
             CALL ACTCOF ( CAT, AN, GAMS )
             GAMAAN = GAMS( 2, 2 )
  
             ! Use GAMAAN for convergence control
             EROR   = ABS( GAMOLD - GAMAAN ) / GAMOLD
             GAMOLD = GAMAAN
  
             ! Check to see if we have a solution
             IF ( EROR .LE. TOLER1 ) THEN
                ASO4  = TSO4 * MWSO4
                AHSO4 = 0.0d0       ! [rjp, 12/12/01]
                ANO3  = XNO3 * MWNO3
                ANH4  = YNH4 * MWNH4
                GNO3  = MAX( FLOOR, ( TMASSHNO3  - ANO3 ) )
                GNH3  = MWNH3 * MAX( FLOOR, ( TNH4 - YNH4 ) )
                AH2O  = 1000.0d0 * WH2O
                RETURN
             ENDIF
  
          ENDDO
  
          ! If after NITR iterations no solution is found, then:
          ! FSB retain the initial values of nitrate particle and vapor
          ! note whether or not convert all bisulfate to sulfate
          ASO4  = TSO4 * MWSO4
          AHSO4 = FLOOR
          XNO3  = TNO3 / MWNO3
          YNH4  = TWOSO4
          ANH4  = YNH4 * MWNH4
  
          CALL AWATER ( IRH, TSO4, YNH4, XNO3, AH2O )
  
          GNO3  = GNO3_IN
          ANO3  = ANO3_IN
          GNH3  = MAX( FLOOR, MWNH3 * (TNH4 - YNH4 ) )
          RETURN
  
        !================================================================
        ! Low Ammonia Case
        !
        ! Coded by Dr. Francis S. Binkowski 12/8/91.(4/26/95)
        ! modified 8/28/98
        ! modified 04/09/2001
        !
        ! All cases covered by this logic
        !=================================================================
        ELSE
  
           WH2O = 0.0d0
           CALL AWATER ( IRH, TSO4, TNH4, TNO3, AH2O )
           WH2O = 1.0d-3 * AH2O
           ZH2O = AH2O
  
           ! convert 10**(-6) kg water/(cubic meter of air) to micrograms
           ! of water per cubic meter of air (1000 g = 1 kg)
           ! in sulfate rich case, preferred form is HSO4-
           !ASO4 = TSO4 * MWSO4
           ASO4  = FLOOR          ![rjp, 12/12/01]
           AHSO4 = TSO4 * MWSO4   ![rjp, 12/12/01]
           ANH4  = TNH4 * MWNH4
           ANO3  = ANO3_IN
           GNO3  = TMASSHNO3 - ANO3
           GNH3  = FLOOR
  
           !==============================================================
           ! *** Examine special cases and return if necessary.
           !
           ! FSB For values of RATIO less than 0.5 do no further
           ! calculations.  The code will cycle and still predict the
           ! same amount of ASO4, ANH4, ANO3, AH2O so terminate early
           ! to swame computation
           !==============================================================
           IF ( RATIO .LT. 0.5d0 ) RETURN ! FSB 04/09/2001
  
           ! Check for zero water.
           IF ( WH2O .EQ. 0.0d0 ) RETURN
           ZSO4 = TSO4 / WH2O
  
           ! ZSO4 is the molality of total sulfate i.e. MSO4 + MHSO4
           ! do not solve for aerosol nitrate for total sulfate molality
           ! greater than 11.0 because the model parameters break down
           !### IF ( ZSO4 .GT. 11.0 ) THEN
           !IF ( ZSO4 .GT. 9.0 ) THEN ! 18 June 97
           !IF ( ZSO4 .GT. 9.d0 ) THEN ! H. Bian 24 June 2015
           IF ( ZSO4 .GT. 9.00 ) THEN ! H. Bian 24 June 2015
              RETURN
           ENDIF
           IF ( ZSO4 .GT. 0.1d0 .and. TEMP .le. 220.d0) THEN ! H. Bian 24 June 2015
              RETURN
           ENDIF
  
           ! *** Calculation may now proceed.
           !
           ! First solve with activity coeffs of 1.0, then iterate.
           PHIOLD = 1.0d0
           GAMANA = 1.0d0
           GAMAS1 = 1.0d0
           GAMAS2 = 1.0d0
           GAMAAB = 1.0d0
           GAMOLD = 1.0d0
  
           ! All ammonia is considered to be aerosol ammonium.
           MNH4 = TNH4 / WH2O
  
           ! MNH4 is the molality of ammonium ion.
           YNH4 = TNH4
  
           ! loop for iteration
           DO NNN = 1, 50    ! loop count reduced 04/09/2001 by FSB
              NITR = NNN
  
              ! set up equilibrium constants including activities
              ! solve the system for hplus first then sulfate & nitrate
              RK2SA  = K2SA * GAMAS2 * GAMAS2 / (GAMAS1 * GAMAS1 * GAMAS1)
              RKNA   = KNA / ( GAMANA * GAMANA )
              RKNWET = RKNA * WH2O
              T21    = ZSO4 - MNH4
              T221   = ZSO4 + T21
  
              ! set up coefficients for cubic
              A2 = RK2SA + RKNWET - T21
              A1 = RK2SA * RKNWET - T21 * ( RK2SA + RKNWET ) &
       &           - RK2SA * ZSO4 - RKNA * TNO3
              A0 = - (T21 * RK2SA * RKNWET &
       &           + RK2SA * RKNWET * ZSO4 + RK2SA * RKNA * TNO3 )
  
              CALL CUBIC ( A2, A1, A0, NR, CRUTES )
  
              ! Code assumes the smallest positive root is in CRUTES(1)
              HPLUS = CRUTES( 1 )
              BAL   = HPLUS **3 + A2 * HPLUS**2 + A1 * HPLUS + A0
  
              ! molality of sulfate ion
              MSO4  = RK2SA * ZSO4 / ( HPLUS + RK2SA )
  
              ! molality of bisulfate ion
              ! MAX added 04/09/2001 by FSB
              MHSO4 = MAX( 1.0d-10, ZSO4 - MSO4 )
  
              ! molality of nitrate ion
              MNA   = RKNA * TNO3 / ( HPLUS + RKNWET )
              MNA   = MAX( 0.0d0, MNA )
              MNA   = MIN( MNA, TNO3 / WH2O )
              XNO3  = MNA * WH2O
              ANO3  = MNA * WH2O * MWNO3
              GNO3  = MAX( FLOOR, TMASSHNO3 - ANO3 )
              ASO4  = MSO4 * WH2O * MWSO4 ![rjp, 12/12/01]
              AHSO4 = MHSO4 * WH2O * MWSO4 ![rjp, 12/12/01]
  
              ! Calculate ionic strength
              STION = 0.5d0 * ( HPLUS + MNA + MNH4 + MHSO4 + 4.0d0 * MSO4)
              ! Update water
              CALL AWATER ( IRH, TSO4, YNH4, XNO3, AH2O )
  
              ! Convert 10**(-6) kg water/(cubic meter of air) to micrograms
              ! of water per cubic meter of air (1000 g = 1 kg)
              WH2O     = 1.0d-3 * AH2O
              CAT( 1 ) = HPLUS
              CAT( 2 ) = MNH4
              AN ( 1 ) = MSO4
              AN ( 2 ) = MNA
              AN ( 3 ) = MHSO4
  
              CALL ACTCOF ( CAT, AN, GAMS )
  
              GAMANA = GAMS( 1, 2 )
              GAMAS1 = GAMS( 1, 1 )
              GAMAS2 = GAMS( 1, 3 )
              GAMAAN = GAMS( 2, 2 )
  
              !------------------------------------------------------------
              ! Add robustness: now check if GAMANA or GAMAS1 is too small
              ! for the division in RKNA and RK2SA. If they are, return w/
              ! original values: basically replicate the procedure used
              ! after the current DO-loop in case of no-convergence
              ! (phs, bmy, rjp, 4/10/08)
              !--------------------------------------------------------------
              IF ( ( ABS( GAMANA ) < EPS ) .OR. ( ABS( GAMAS1 ) < EPS ) ) THEN
  
                 ! Reset to original values
                 ANH4  = TNH4 * MWNH4
                 GNH3  = FLOOR
                 GNO3  = GNO3_IN
                 ANO3  = ANO3_IN
                 ASO4  = TSO4 * MWSO4
                 AHSO4 = FLOOR
  
                 ! Update water
                 CALL AWATER ( IRH, TSO4, TNH4, TNO3, AH2O )
  
                 ! Exit this subroutine
                 RETURN
              ENDIF
  
              GAMAHAT = ( GAMAS2 * GAMAS2 / ( GAMAAB * GAMAAB ) )
              BHAT = KHAT * GAMAHAT
              !### EROR = ABS ( ( PHIOLD - PHIBAR ) / PHIOLD )
              !### PHIOLD = PHIBAR
              EROR = ABS ( GAMOLD - GAMAHAT ) / GAMOLD
              GAMOLD = GAMAHAT
              ! return with good solution
              IF ( EROR .LE. TOLER2 ) THEN
                 RETURN
              ENDIF
  
           ENDDO
  
           ! after NITR iterations, failure to solve the system
           ! convert all ammonia to aerosol ammonium and return input
           ! values of NO3 and HNO3
           ANH4 = TNH4 * MWNH4
           GNH3 = FLOOR
           GNO3 = GNO3_IN
           ANO3 = ANO3_IN
           ASO4 = TSO4 * MWSO4    ! [rjp, 12/17/01]
           AHSO4= FLOOR           ! [rjp, 12/17/01]
  
           CALL AWATER ( IRH, TSO4, TNH4, TNO3, AH2O )
  
           RETURN
  
        ENDIF                     ! ratio .gt. 2.0
  
        ! Return to calling program
  
     end subroutine RPMARES
  
  !------------------------------------------------------------------------------
  
        SUBROUTINE AWATER( IRHX, MSO4, MNH4, MNO3, WH2O )
  !
  !******************************************************************************
  ! NOTE!!! wh2o is returned in micrograms / cubic meter
  !         mso4,mnh4,mno3 are in microMOLES / cubic meter
  !
  !  This  version uses polynomials rather than tables, and uses empirical
  ! polynomials for the mass fraction of solute (mfs) as a function of water
  ! activity
  !   where:
  !
  !            mfs = ms / ( ms + mw)
  !             ms is the mass of solute
  !             mw is the mass of water.
  !
  !  Define y = mw/ ms
  !
  !  then  mfs = 1 / (1 + y)
  !
  !    y can then be obtained from the values of mfs as
  !
  !             y = (1 - mfs) / mfs
  !
  !
  !     the aerosol is assumed to be in a metastable state if the rh is
  !     is below the rh of deliquescence, but above the rh of crystallization.
  !
  !     ZSR interpolation is used for sulfates with x ( the molar ratio of
  !     ammonium to sulfate in eh range 0 <= x <= 2, by sections.
  !     section 1: 0 <= x < 1
  !     section 2: 1 <= x < 1.5
  !     section 3: 1.5 <= x < 2.0
  !     section 4: 2 <= x
  !     In sections 1 through 3, only the sulfates can affect the amount of water
  !     on the particles.
  !     In section 4, we have fully neutralized sulfate, and extra ammonium which
  !     allows more nitrate to be present. Thus, the ammount of water is
  !     calculated
  !     using ZSR for ammonium sulfate and ammonium nitrate. Crystallization is
  !     assumed to occur in sections 2,3,and 4. See detailed discussion below.
  !
  ! definitions:
  !     mso4, mnh4, and mno3 are the number of micromoles/(cubic meter of air)
  !      for sulfate, ammonium, and nitrate respectively
  !     irhx is the relative humidity (%)
  !     wh2o is the returned water amount in micrograms / cubic meter of air
  !     x is the molar ratio of ammonium to sulfate
  !     y0,y1,y1.5, y2 are the water contents in mass of water/mass of solute
  !     for pure aqueous solutions with x equal 1, 1.5, and 2 respectively.
  !     y3 is the value of the mass ratio of water to solute for
  !     a pure ammonium nitrate  solution.
  !
  !
  !     coded by Dr. Francis S. Binkowski, 4/8/96.
  !
  ! *** modified 05/30/2000
  !     The use of two values of mfs at an ammonium to sulfate ratio
  !     representative of ammonium sulfate led to an minor inconsistancy
  !     in nitrate behavior as the ratio went from a value less than two
  !     to a value greater than two and vice versa with either ammonium
  !     held constant and sulfate changing, or sulfate held constant and
  !     ammonium changing. the value of Chan et al. (1992) is the only value
  !     now used.
  !
  ! *** Modified 09/25/2002
  !     Ported into "rpmares_mod.f".  Now declare all variables with REAL*8.
  !     Also cleaned up comments and made cosmetic changes.  Force double
  !     precision explicitly with "D" exponents.
  !******************************************************************************
  !
        ! Arguments
        INTEGER           :: IRHX
        REAL*8            :: MSO4, MNH4, MNO3, WH2O
  
        ! Local variables
        INTEGER           :: IRH
        REAL*8            :: TSO4,  TNH4,  TNO3,  X,      AW,     AWC
        REAL*8            :: MFS0,  MFS1,  MFS15, Y
        REAL*8            :: Y0,    Y1,    Y15,   Y2,     Y3,     Y40
        REAL*8            :: Y140,  Y1540, YC,    MFSSO4, MFSNO3
  
        ! Molecular weight parameters
        REAL*8, PARAMETER :: MWSO4  = 96.0636d0
        REAL*8, PARAMETER :: MWNH4  = 18.0985d0
        REAL*8, PARAMETER :: MWNO3  = 62.0649d0
        REAL*8, PARAMETER :: MW2    = MWSO4 + 2.0d0 * MWNH4
        REAL*8, PARAMETER :: MWANO3 = MWNO3 + MWNH4
  
        !=================================================================
        ! The polynomials use data for aw as a function of mfs from Tang
        ! and Munkelwitz, JGR 99: 18801-18808, 1994.  The polynomials were
        ! fit to Tang's values of water activity as a function of mfs.
        !
        ! *** coefficients of polynomials fit to Tang and Munkelwitz data
        !     now give mfs as a function of water activity.
        !=================================================================
        REAL*8 :: C1(4)  = (/ 0.9995178d0,  -0.7952896d0, &
       &                      0.99683673d0, -1.143874d0 /)
  
        REAL*8 :: C15(4) = (/ 1.697092d0, -4.045936d0, &
       &                      5.833688d0, -3.463783d0 /)
  
        !=================================================================
        ! The following coefficients are a fit to the data in Table 1 of
        !    Nair & Vohra, J. Aerosol Sci., 6: 265-271, 1975
        !      data c0/0.8258941, -1.899205, 3.296905, -2.214749 /
        !
        ! New data fit to data from
        !       Nair and Vohra J. Aerosol Sci., 6: 265-271, 1975
        !       Giaque et al. J.Am. Chem. Soc., 82: 62-70, 1960
        !       Zeleznik J. Phys. Chem. Ref. Data, 20: 157-1200
        !=================================================================
        REAL*8 :: C0(4)  =  (/ 0.798079d0, -1.574367d0, &
       &                       2.536686d0, -1.735297d0 /)
  
        !=================================================================
        ! Polynomials for ammonium nitrate and ammonium sulfate are from:
        ! Chan et al.1992, Atmospheric Environment (26A): 1661-1673.
        !=================================================================
        REAL*8 :: KNO3(6) = (/  0.2906d0,   6.83665d0, -26.9093d0, &
       &                       46.6983d0, -38.803d0,    11.8837d0 /)
  
        REAL*8 :: KSO4(6) = (/   2.27515d0, -11.147d0,   36.3369d0, &
       &                       -64.2134d0,   56.8341d0, -20.0953d0 /)
  
        !=================================================================
        ! AWATER begins here!
        !=================================================================
  
        ! Check range of per cent relative humidity
        IRH  = IRHX
        IRH  = MAX( 1, IRH )
        IRH  = MIN( IRH, 100 )
  
        ! Water activity = fractional relative humidity
        AW   = DBLE( IRH ) / 100.0d0
        TSO4 = MAX( MSO4 , 0.0d0 )
        TNH4 = MAX( MNH4 , 0.0d0 )
        TNO3 = MAX( MNO3 , 0.0d0 )
        X    = 0.0d0
  
        ! If there is non-zero sulfate calculate the molar ratio
        ! otherwise check for non-zero nitrate and ammonium
        IF ( TSO4 .GT. 0.0d0 ) THEN
           X = TNH4 / TSO4
        ELSE
           IF ( TNO3 .GT. 0.0d0 .AND. TNH4 .GT. 0.0d0 ) X = 10.0d0
        ENDIF
  
        ! *** begin screen on x for calculating wh2o
        IF ( X .LT. 1.0d0 ) THEN
           MFS0 = nh3_POLY4( C0, AW )
           MFS1 = nh3_POLY4( C1, AW )
           Y0   = ( 1.0d0 - MFS0 ) / MFS0
           Y1   = ( 1.0d0 - MFS1 ) / MFS1
           Y    = ( 1.0d0 - X    ) * Y0 + X * Y1
  
        ELSE IF ( X .LT. 1.5d0 ) THEN
  
           IF ( IRH .GE. 40 ) THEN
              MFS1  = nh3_POLY4( C1,  AW )
              MFS15 = nh3_POLY4( C15, AW )
              Y1    = ( 1.0d0 - MFS1  ) / MFS1
              Y15   = ( 1.0d0 - MFS15 ) / MFS15
              Y     = 2.0d0 * ( Y1 * ( 1.5d0 - X ) + Y15 *( X - 1.0d0 ) )
  
           !==============================================================
           ! Set up for crystalization
           !
           ! Crystallization is done as follows:
           !
           ! For 1.5 <= x, crystallization is assumed to occur
           ! at rh = 0.4
           !
           ! For x <= 1.0, crystallization is assumed to occur at an
           ! rh < 0.01, and since the code does not allow ar rh < 0.01,
           ! crystallization is assumed not to occur in this range.
           !
           ! For 1.0 <= x <= 1.5 the crystallization curve is a straignt
           ! line from a value of y15 at rh = 0.4 to a value of zero at
           ! y1. From point B to point A in the diagram.  The algorithm
           ! does a double interpolation to calculate the amount of
           ! water.
           !
           !        y1(0.40)               y15(0.40)
           !         +                     + Point B
           !
           !
           !
           !
           !         +--------------------+
           !       x=1                   x=1.5
           !      Point A
           !==============================================================
           ELSE
  
              ! rh along the crystallization curve.
              AWC = 0.80d0 * ( X - 1.0d0 )
              Y   = 0.0d0
  
              ! interpolate using crystalization curve
              IF ( AW .GE. AWC ) THEN
                 MFS1  = nh3_POLY4( C1,  0.40d0 )
                 MFS15 = nh3_POLY4( C15, 0.40d0 )
                 Y140  = ( 1.0d0 - MFS1  ) / MFS1
                 Y1540 = ( 1.0d0 - MFS15 ) / MFS15
                 Y40   = 2.0d0 * ( Y140  * ( 1.5d0 - X ) + &
       &                           Y1540 * ( X - 1.0d0 ) )
  
                 ! Y along crystallization curve
                 YC   = 2.0d0 * Y1540 * ( X - 1.0d0 )
                 Y    = Y40 - (Y40 - YC) * (0.40d0 - AW) / (0.40d0 - AWC)
              ENDIF
           ENDIF
  
        ELSE IF ( X .LT. 2.0d0 ) then               ! changed 12/11/2000 by FSB
           Y = 0.0D0
  
           IF ( IRH .GE. 40 ) THEN
              MFS15  = nh3_POLY4( C15, AW )
              !MFS2  = nh3_POLY4( C2,  AW )
              Y15    = ( 1.0d0 - MFS15 ) / MFS15
              !y2    = ( 1.0d0 - MFS2  ) / MFS2
              MFSSO4 = nh3_POLY6( KSO4, AW )             ! Changed 05/30/2000 by FSB
              Y2     = ( 1.0d0 - MFSSO4 ) / MFSSO4
              Y      = 2.0d0 * (Y15 * (2.0d0 - X) + Y2 * (X - 1.5d0) )
           ENDIF
  
        ELSE                                 ! 2.0 <= x changed 12/11/2000 by FSB
  
           !==============================================================
           ! Regime where ammonium sulfate and ammonium nitrate are
           ! in solution.
           !
           ! following cf&s for both ammonium sulfate and ammonium nitrate
           ! check for crystallization here. their data indicate a 40%
           ! value is appropriate.
           !==============================================================
           Y2 = 0.0d0
           Y3 = 0.0d0
  
           IF ( IRH .GE. 40 ) THEN
              MFSSO4 = nh3_POLY6( KSO4, AW )
              MFSNO3 = nh3_POLY6( KNO3, AW )
              Y2     = ( 1.0d0 - MFSSO4 ) / MFSSO4
              Y3     = ( 1.0d0 - MFSNO3 ) / MFSNO3
  
           ENDIF
  
        ENDIF                     ! end of checking on x
  
        !=================================================================
        ! Now set up output of WH2O
        ! WH2O units are micrograms (liquid water) / cubic meter of air
        !=================================================================
        IF ( X .LT. 2.0D0 ) THEN  ! changed 12/11/2000 by FSB
  
           WH2O =  Y * ( TSO4 * MWSO4 + MWNH4 * TNH4 )
  
        ELSE
  
           ! this is the case that all the sulfate is ammonium sulfate
           ! and the excess ammonium forms ammonum nitrate
           WH2O =   Y2 * TSO4 * MW2 + Y3 * TNO3 * MWANO3
  
        ENDIF
  
        ! Return to calling program
        END SUBROUTINE AWATER
  
  !------------------------------------------------------------------------------
  
        FUNCTION nh3_POLY4( A, X ) RESULT( Y )
  
        ! Arguments
        REAL*8, INTENT(IN) :: A(4), X
  
        ! Return value
        REAL*8             :: Y
  
        !=================================================================
        ! nh3_POLY4 begins here!
        !=================================================================
        Y = A(1) + X * ( A(2) + X * ( A(3) + X * ( A(4) )))
  
        ! Return to calling program
        END FUNCTION nh3_POLY4
  
  !------------------------------------------------------------------------------
  
        FUNCTION nh3_POLY6( A, X ) RESULT( Y )
  
        ! Arguments
        REAL*8, INTENT(IN) :: A(6), X
  
        ! Return value
        REAL*8             :: Y
  
        !=================================================================
        ! nh3_POLY6 begins here!
        !=================================================================
        Y = A(1) + X * ( A(2) + X * ( A(3) + X * ( A(4) +  &
       &           X * ( A(5) + X * ( A(6)  )))))
  
        ! Return to calling program
        END FUNCTION nh3_POLY6
  
  !------------------------------------------------------------------------------
  
        SUBROUTINE CUBIC( A2, A1, A0, NR, CRUTES )
  
  !
  !******************************************************************************
  ! Subroutine to find the roots of a cubic equation / 3rd order polynomial
  ! Formulae can be found in numer. recip.  on page 145
  !   kiran  developed  this version on 25/4/1990
  !   Dr. Francis S. Binkowski modified the routine on 6/24/91, 8/7/97
  ! ***
  ! *** modified 2/23/98 by fsb to incorporate Dr. Ingmar Ackermann's
  !     recommendations for setting a0, a1,a2 as real*8 variables.
  !
  ! Modified by Bob Yantosca (10/15/02)
  ! - Now use upper case / white space
  ! - force double precision with "D" exponents
  ! - updated comments / cosmetic changes
  ! - now call ERROR_STOP from "error_mod.f" to stop the run safely
  !******************************************************************************
  !
        ! Arguments
        INTEGER           :: NR
        REAL*8            :: A2, A1, A0
        REAL*8            :: CRUTES(3)
  
        ! Local variables
        REAL*8            :: QQ,    RR,    A2SQ,  THETA, DUM1, DUM2
        REAL*8            :: PART1, PART2, PART3, RRSQ,  PHI,  YY1
        REAL*8            :: YY2,   YY3,   COSTH, SINTH
        REAL*8, PARAMETER :: ONE    = 1.0d0
        REAL*8, PARAMETER :: SQRT3  = 1.732050808d0
        REAL*8, PARAMETER :: ONE3RD = 0.333333333d0
        ! !LOCAL VARIABLES:
        character (len=75) :: err_msg
  
        integer :: status
  
        !=================================================================
        ! CUBIC begins here!
        !=================================================================
        A2SQ = A2 * A2
        QQ   = ( A2SQ - 3.d0*A1 ) / 9.d0
        RR   = ( A2*( 2.d0*A2SQ - 9.d0*A1 ) + 27.d0*A0 ) / 54.d0
  
        ! CASE 1 THREE REAL ROOTS or  CASE 2 ONLY ONE REAL ROOT
        DUM1 = QQ * QQ * QQ
        RRSQ = RR * RR
        DUM2 = DUM1 - RRSQ
  
        IF ( DUM2 .GE. 0.d0 ) THEN
  
           ! Now we have three real roots
           PHI = SQRT( DUM1 )
  
           IF ( ABS( PHI ) .LT. 1.d-20 ) THEN
              CRUTES(1) = 0.0d0
              CRUTES(2) = 0.0d0
              CRUTES(3) = 0.0d0
              NR        = 0
              print *,'PHI < 1d-20 in  CUBIC (rpmares_mod.f)'
           ENDIF
  
           THETA = ACOS( RR / PHI ) / 3.0d0
           COSTH = COS( THETA )
           SINTH = SIN( THETA )
  
           ! Use trig identities to simplify the expressions
           ! Binkowski's modification
           PART1     = SQRT( QQ )
           YY1       = PART1 * COSTH
           YY2       = YY1 - A2/3.0d0
           YY3       = SQRT3 * PART1 * SINTH
           CRUTES(3) = -2.0d0*YY1 - A2/3.0d0
           CRUTES(2) = YY2 + YY3
           CRUTES(1) = YY2 - YY3
  
           ! Set negative roots to a large positive value
           IF ( CRUTES(1) .LT. 0.0d0 ) CRUTES(1) = 1.0d9
           IF ( CRUTES(2) .LT. 0.0d0 ) CRUTES(2) = 1.0d9
           IF ( CRUTES(3) .LT. 0.0d0 ) CRUTES(3) = 1.0d9
  
           ! Put smallest positive root in crutes(1)
           CRUTES(1) = MIN( CRUTES(1), CRUTES(2), CRUTES(3) )
           NR        = 3
  
        ELSE
  
           ! Now here we have only one real root
           PART1     = SQRT( RRSQ - DUM1 )
           PART2     = ABS( RR )
           PART3     = ( PART1 + PART2 )**ONE3RD
           CRUTES(1) = -SIGN(ONE,RR) * ( PART3 + (QQ/PART3) ) - A2/3.D0
           CRUTES(2) = 0.D0
           CRUTES(3) = 0.D0
           NR        = 1
  
        ENDIF
  
        ! Return to calling program
        END SUBROUTINE CUBIC
  
  !------------------------------------------------------------------------------
  
         SUBROUTINE ACTCOF( CAT, AN, GAMA, MOLNU, PHIMULT )
  !
  !******************************************************************************
  !
  ! DESCRIPTION:
  !
  !  This subroutine computes the activity coefficients of (2NH4+,SO4--),
  !  (NH4+,NO3-),(2H+,SO4--),(H+,NO3-),AND (H+,HSO4-) in aqueous
  !  multicomponent solution, using Bromley's model and Pitzer's method.
  !
  ! REFERENCES:
  !
  !   Bromley, L.A. (1973) Thermodynamic properties of strong electrolytes
  !     in aqueous solutions.  AIChE J. 19, 313-320.
  !
  !   Chan, C.K. R.C. Flagen, & J.H.  Seinfeld (1992) Water Activities of
  !     NH4NO3 / (NH4)2SO4 solutions, Atmos. Environ. (26A): 1661-1673.
  !
  !   Clegg, S.L. & P. Brimblecombe (1988) Equilibrium partial pressures
  !     of strong acids over saline solutions - I HNO3,
  !     Atmos. Environ. (22): 91-100
  !
  !   Clegg, S.L. & P. Brimblecombe (1990) Equilibrium partial pressures
  !     and mean activity and osmotic coefficients of 0-100% nitric acid
  !     as a function of temperature,   J. Phys. Chem (94): 5369 - 5380
  !
  !   Pilinis, C. and J.H. Seinfeld (1987) Continued development of a
  !     general equilibrium model for inorganic multicomponent atmospheric
  !     aerosols.  Atmos. Environ. 21(11), 2453-2466.
  !
  !
  !
  !
  ! ARGUMENT DESCRIPTION:
  !
  !     CAT(1) : conc. of H+    (moles/kg)
  !     CAT(2) : conc. of NH4+  (moles/kg)
  !     AN(1)  : conc. of SO4-- (moles/kg)
  !     AN(2)  : conc. of NO3-  (moles/kg)
  !     AN(3)  : conc. of HSO4- (moles/kg)
  !     GAMA(2,1)    : mean molal ionic activity coeff for (2NH4+,SO4--)
  !     GAMA(2,2)    :  "    "     "       "       "    "  (NH4+,NO3-)
  !     GAMA(2,3)    :  "    "     "       "       "    "  (NH4+. HSO4-)
  !     GAMA(1,1)    :  "    "     "       "       "    "  (2H+,SO4--)
  !     GAMA(1,2)    :  "    "     "       "       "    "  (H+,NO3-)
  !     GAMA(1,3)    :  "    "     "       "       "    "  (H+,HSO4-)
  !     MOLNU   : the total number of moles of all ions.
  !     PHIMULT : the multicomponent paractical osmotic coefficient.
  !
  ! REVISION HISTORY:
  !      Who       When        Detailed description of changes
  !   ---------   --------  -------------------------------------------
  !   S.Roselle   7/26/89   Copied parts of routine BROMLY, and began this
  !                         new routine using a method described by Pilinis
  !                         and Seinfeld 1987, Atmos. Envirn. 21 pp2453-2466.
  !   S.Roselle   7/30/97   Modified for use in Models-3
  !   F.Binkowski 8/7/97    Modified coefficients BETA0, BETA1, CGAMA
  !   R.Yantosca  9/25/02   Ported into "rpmares_mod.f" for GEOS-CHEM.  Cleaned
  !                         up comments, etc.  Also force double precision by
  !                         declaring REALs as REAL*8 and by using "D" exponents.
  !******************************************************************************
        ! Error codes
  
  
  
        !=================================================================
        ! PARAMETERS and their descriptions:
        !=================================================================
        INTEGER, PARAMETER :: NCAT = 2         ! number of cation
        INTEGER, PARAMETER :: NAN  = 3         ! number of anions
        REAL*8,  PARAMETER :: XSTAT0 = 0       ! Normal, successful completion
        REAL*8,  PARAMETER :: XSTAT1 = 1       ! File I/O error
        REAL*8,  PARAMETER :: XSTAT2 = 2       ! Execution error
        REAL*8,  PARAMETER :: XSTAT3 = 3       ! Special  error
  
        !=================================================================
        ! ARGUMENTS and their descriptions
        !=================================================================
        REAL*8, optional   :: MOLNU            ! tot # moles of all ions
        REAL*8, optional   :: PHIMULT          ! multicomponent paractical
                                               !   osmotic coef
        REAL*8             :: CAT(NCAT)        ! cation conc in moles/kg (input)
        REAL*8             :: AN(NAN)          ! anion conc in moles/kg (input)
        REAL*8             :: GAMA(NCAT,NAN)   ! mean molal ionic activity coefs
  
        !=================================================================
        ! SCRATCH LOCAL VARIABLES and their descriptions:
        !=================================================================
        INTEGER            :: IAN              ! anion indX
        INTEGER            :: ICAT             ! cation indX
        REAL*8             :: FGAMA            !
        REAL*8             :: I                ! ionic strength
        REAL*8             :: R                !
        REAL*8             :: S                !
        REAL*8             :: TA               !
        REAL*8             :: TB               !
        REAL*8             :: TC               !
        REAL*8             :: TEXPV            !
        REAL*8             :: TRM              !
        REAL*8             :: TWOI             ! 2*ionic strength
        REAL*8             :: TWOSRI           ! 2*sqrt of ionic strength
        REAL*8             :: ZBAR             !
        REAL*8             :: ZBAR2            !
        REAL*8             :: ZOT1             !
        REAL*8             :: SRI              ! square root of ionic strength
        REAL*8             :: F2(NCAT)         !
        REAL*8             :: F1(NAN)          !
        REAL*8             :: BGAMA (NCAT,NAN) !
        REAL*8             :: X     (NCAT,NAN) !
        REAL*8             :: M     (NCAT,NAN) ! molality of each electrolyte
        REAL*8             :: LGAMA0(NCAT,NAN) ! binary activity coefficients
        REAL*8             :: Y     (NAN,NCAT) !
        REAL*8             :: BETA0 (NCAT,NAN) ! binary activity coef parameter
        REAL*8             :: BETA1 (NCAT,NAN) ! binary activity coef parameter
        REAL*8             :: CGAMA (NCAT,NAN) ! binary activity coef parameter
        REAL*8             :: V1    (NCAT,NAN) ! # of cations in electrolyte
                                               !   formula
        REAL*8             :: V2    (NCAT,NAN) ! # of anions in electrolyte
                                               !   formula
        ! absolute value of charges of cation
        REAL*8             :: ZP(NCAT) = (/ 1.0d0, 1.0d0 /)
  
        ! absolute value of charges of anion
        REAL*8             :: ZM(NAN)  = (/ 2.0d0, 1.0d0, 1.0d0 /)
  
        ! Character values.
        CHARACTER(LEN=120)      :: XMSG  = ' '
  !      CHARACTER(LEN=16), SAVE :: PNAME = ' driver program name'
  
        !================================================================
        ! *** Sources for the coefficients BETA0, BETA1, CGAMA
        ! (1,1);(1,3)  - Clegg & Brimblecombe (1988)
        ! (2,3)        - Pilinis & Seinfeld (1987), cgama different
        ! (1,2)        - Clegg & Brimblecombe (1990)
        ! (2,1);(2,2)  - Chan, Flagen & Seinfeld (1992)
        !================================================================
  
        ! now set the basic constants, BETA0, BETA1, CGAMA
        DATA BETA0(1,1) /2.98d-2/,      BETA1(1,1) / 0.0d0/,  &
       &     CGAMA(1,1) /4.38d-2/                                 ! 2H+SO4-
  
        DATA BETA0(1,2) /  1.2556d-1/,  BETA1(1,2) / 2.8778d-1/,  &
       &     CGAMA(1,2) / -5.59d-3/                               ! HNO3
  
        DATA BETA0(1,3) / 2.0651d-1/,   BETA1(1,3) / 5.556d-1/,  &
       &     CGAMA(1,3) /0.0d0/                                   ! H+HSO4-
  
        DATA BETA0(2,1) / 4.6465d-2/,   BETA1(2,1) /-0.54196d0/,  &
       &     CGAMA(2,1) /-1.2683d-3/                              ! (NH4)2SO4
  
        DATA BETA0(2,2) /-7.26224d-3/,  BETA1(2,2) /-1.168858d0/,  &
       &     CGAMA(2,2) / 3.51217d-5/                             ! NH4NO3
  
        DATA BETA0(2,3) / 4.494d-2/,    BETA1(2,3) / 2.3594d-1/,  &
       &     CGAMA(2,3) /-2.962d-3/                               ! NH4HSO4
  
        DATA V1(1,1), V2(1,1) / 2.0d0, 1.0d0 /     ! 2H+SO4-
        DATA V1(2,1), V2(2,1) / 2.0d0, 1.0d0 /     ! (NH4)2SO4
        DATA V1(1,2), V2(1,2) / 1.0d0, 1.0d0 /     ! HNO3
        DATA V1(2,2), V2(2,2) / 1.0d0, 1.0d0 /     ! NH4NO3
        DATA V1(1,3), V2(1,3) / 1.0d0, 1.0d0 /     ! H+HSO4-
        DATA V1(2,3), V2(2,3) / 1.0d0, 1.0d0 /     ! NH4HSO4
  
        !=================================================================
        ! ACTCOF begins here!
        !=================================================================
  
        ! Compute ionic strength
        I = 0.0d0
        DO ICAT = 1, NCAT
           I = I + CAT( ICAT ) * ZP( ICAT ) * ZP( ICAT )
        ENDDO
  
        DO IAN = 1, NAN
           I = I + AN( IAN ) * ZM( IAN ) * ZM( IAN )
        ENDDO
  
        I = 0.5d0 * I
  
        ! check for problems in the ionic strength
        IF ( I .EQ. 0.0d0 ) THEN
  
           DO IAN  = 1, NAN
           DO ICAT = 1, NCAT
              GAMA( ICAT, IAN ) = 0.0d0
           ENDDO
           ENDDO
  
           XMSG = 'Ionic strength is zero...returning zero activities'
           !CALL M3WARN ( PNAME, 0, 0, XMSG )
           RETURN
  
        ELSE IF ( I .LT. 0.0d0 ) THEN
           XMSG = 'Ionic strength below zero...negative concentrations'
           !CALL M3EXIT ( PNAME, 0, 0, XMSG, XSTAT1 )
        ENDIF
  
        ! Compute some essential expressions
        SRI    = SQRT( I )
        TWOSRI = 2.0d0 * SRI
        TWOI   = 2.0d0 * I
        TEXPV  = 1.0d0 - EXP( -TWOSRI ) * ( 1.0d0 + TWOSRI - TWOI )
        R      = 1.0d0 + 0.75d0 * I
        S      = 1.0d0 + 1.5d0  * I
        ZOT1   = 0.511d0 * SRI / ( 1.0d0 + SRI )
  
        ! Compute binary activity coeffs
        FGAMA = -0.392d0 * ( ( SRI / ( 1.0d0 + 1.2d0 * SRI )  &
       &      + ( 2.0d0 / 1.2d0 ) * LOG( 1.0d0 + 1.2d0 * SRI ) ) )
  
        DO ICAT = 1, NCAT
        DO IAN  = 1, NAN
  
           BGAMA( ICAT, IAN ) = 2.0d0 * BETA0( ICAT, IAN )  &
       &        + ( 2.0d0 * BETA1( ICAT, IAN ) / ( 4.0d0 * I ) )  &
       &        * TEXPV
  
           ! Compute the molality of each electrolyte for given ionic strength
           M( ICAT, IAN ) = ( CAT( ICAT )**V1( ICAT, IAN )  &
       &                   *   AN( IAN )**V2( ICAT, IAN ) )**( 1.0d0  &
       &                   / ( V1( ICAT, IAN ) + V2( ICAT, IAN ) ) )
  
           ! Calculate the binary activity coefficients
           LGAMA0( ICAT, IAN ) = ( ZP( ICAT ) * ZM( IAN ) * FGAMA  &
       &        + M( ICAT, IAN )  &
       &        * ( 2.0d0 * V1( ICAT, IAN ) * V2( ICAT, IAN )  &
       &        / ( V1( ICAT, IAN ) + V2( ICAT, IAN ) )  &
       &        * BGAMA( ICAT, IAN ) )  &
       &        + M( ICAT, IAN ) * M( ICAT, IAN )  &
       &        * ( 2.0d0 * ( V1( ICAT, IAN )  &
       &        * V2( ICAT, IAN ) )**1.5d0  &
       &        / ( V1( ICAT, IAN ) + V2( ICAT, IAN ) )  &
       &        * CGAMA( ICAT, IAN ) ) ) / 2.302585093d0
  
        ENDDO
        ENDDO
  
        ! prepare variables for computing the multicomponent activity coeffs
        DO IAN = 1, NAN
        DO ICAT = 1, NCAT
           ZBAR           = ( ZP( ICAT ) + ZM( IAN ) ) * 0.5d0
           ZBAR2          = ZBAR * ZBAR
           Y( IAN, ICAT ) = ZBAR2 * AN( IAN ) / I
           X( ICAT, IAN ) = ZBAR2 * CAT( ICAT ) / I
        ENDDO
        ENDDO
  
        DO IAN = 1, NAN
           F1( IAN ) = 0.0d0
           DO ICAT = 1, NCAT
              F1( IAN ) = F1( IAN ) + X( ICAT, IAN ) * LGAMA0( ICAT, IAN )  &
       &                + ZOT1 * ZP( ICAT ) * ZM( IAN ) * X( ICAT, IAN )
           ENDDO
        ENDDO
  
        DO ICAT = 1, NCAT
           F2( ICAT ) = 0.0d0
           DO IAN = 1, NAN
              F2( ICAT ) = F2( ICAT ) + Y( IAN, ICAT ) * LGAMA0(ICAT, IAN)  &
       &                 + ZOT1 * ZP( ICAT ) * ZM( IAN ) * Y( IAN, ICAT )
           ENDDO
        ENDDO
  
        ! now calculate the multicomponent activity coefficients
        DO IAN  = 1, NAN
        DO ICAT = 1, NCAT
  
           TA  = -ZOT1 * ZP( ICAT ) * ZM( IAN )
           TB  = ZP( ICAT ) * ZM( IAN ) / ( ZP( ICAT ) + ZM( IAN ) )
           TC  = ( F2( ICAT ) / ZP( ICAT ) + F1( IAN ) / ZM( IAN ) )
           TRM = TA + TB * TC
  
           IF ( TRM .GT. 30.0d0 ) THEN
              GAMA( ICAT, IAN ) = 1.0d+30
           ELSE
              GAMA( ICAT, IAN ) = 10.0d0**TRM
           ENDIF
  
        ENDDO
        ENDDO
  
        ! Return to calling program
        END SUBROUTINE ACTCOF
  
  !============================================================================
  !BOP
  !
  ! !IROUTINE: HNO3_reaction_rate
  !
  ! !INTERFACE:
     subroutine HNO3_reaction_rate(i1, i2, j1, j2, km, klid, rmed, fnum, rhoa, temp, rh, q, kan, &
                                   AVOGAD, AIRMW, PI, RUNIV, fMassHNO3)
  
  ! !DESCRIPTION:
  
  ! !USES:
     implicit NONE
  
  ! !INPUT PARAMETERS:
     integer, intent(in)  :: i1, i2, j1, j2        ! grid dimension
     integer, intent(in)                 :: km     ! model levels
     integer, intent(in)                 :: klid   ! index for pressure lid
     real, intent(in)                    :: rmed   ! aerosol radius [um]
     real, intent(in)                    :: fnum   ! number of aerosol particles per kg mass
     real, dimension(:,:,:), intent(in)  :: rhoa   ! Layer air density [kg/m^3]
     real, dimension(:,:,:), intent(in)  :: temp   ! Layer temperature [K]
     real, dimension(:,:,:), intent(in)  :: rh     ! relative humidity [1]
     real, dimension(:,:,:), intent(in)  :: q      ! aerosol
     real, intent(in)                    :: AVOGAD ! Avogadro's number [1/kmol]
     real, intent(in)                    :: AIRMW  ! molecular weight of air [kg/kmol]
     real, intent(in)                    :: PI     ! pi constant
     real, intent(in)                    :: RUNIV  ! ideal gas constant [J/(Kmole*K)]
     real, intent(in)                    :: fMassHNO3      ! gram molecular weight
  
  ! !OUTPUT PARAMETERS:
     real, dimension(:,:,:), intent(out) :: kan
  
  ! !Local Variables
     integer :: i, j, k
  
     real :: f_sad
     real :: f_ad
     real :: radius
     real :: ad
     real :: sad
  
  !EOP
  !------------------------------------------------------------------------------------
  !  Begin..
  
     f_ad = 1.e-6 * AVOGAD / AIRMW  ! air number density # cm-3 per unit air density
  
     ! surface area per unit air density and unit aerosol mass mixing ratio
     f_sad  = 0.01 * 4 * PI * rmed**2 * fnum
  
     ! radius in 'cm'
     radius = 100 * rmed
  
     do k = klid, km
        do j = j1, j2
           do i = i1, i2
              ad   = f_ad  * rhoa(i,j,k)             ! air number density # cm-3
              sad  = f_sad * rhoa(i,j,k) * q(i,j,k)  ! surface area density cm2 cm-3
  
              kan(i,j,k) = sktrs_hno3(temp(i,j,k), rh(i,j,k), sad, ad, radius, PI, &
                                      RUNIV, fMassHNO3)
           end do
        end do
     end do
  
     end subroutine HNO3_reaction_rate
  
  !============================================================================
  !BOP
  !
  ! !IROUTINE: SSLT_reaction_rate
  !
  ! !INTERFACE:
     subroutine SSLT_reaction_rate(i1, i2, j1, j2, km, klid, rmed, fnum, rhoa, temp, rh, q, kan, &
                                   AVOGAD, AIRMW, PI, RUNIV, fMassHNO3)
  
  ! !DESCRIPTION:
  
  ! !USES:
     implicit NONE
  
  ! !INPUT PARAMETERS:
     integer, intent(in)  :: i1, i2, j1, j2        ! grid dimension
     integer, intent(in)                 :: km     ! model levels
     integer, intent(in)                 :: klid   ! index for pressure lid
     real, intent(in)                    :: rmed   ! aerosol radius [um]
     real, intent(in)                    :: fnum   ! number of aerosol particles per kg mass
     real, dimension(:,:,:), intent(in)  :: rhoa   ! Layer air density [kg/m^3]
     real, dimension(:,:,:), intent(in)  :: temp   ! Layer temperature [K]
     real, dimension(:,:,:), intent(in)  :: rh     ! relative humidity [1]
     real, dimension(:,:,:), intent(in)  :: q      ! aerosol
     real, intent(in)                    :: AVOGAD ! Avogadro's number [1/kmol]
     real, intent(in)                    :: AIRMW  ! molecular weight of air [kg/kmol]
     real, intent(in)                    :: PI     ! pi constant
     real, intent(in)                    :: RUNIV  ! ideal gas constant [J/(Kmole*K)]
     real, intent(in)                    :: fMassHNO3      ! gram molecular weight
  
  ! !OUTPUT PARAMETERS:
     real, dimension(:,:,:), intent(out) :: kan
  
  ! !Local Variables
     integer :: i, j, k
  
     real :: f_sad
     real :: f_ad
     real :: radius
     real :: ad
     real :: sad
  
  !EOP
  !------------------------------------------------------------------------------------
  !  Begin..
  
        f_ad = 1.e-6 * AVOGAD / AIRMW  ! air number density # cm-3 per unit air density
  
        ! surface area per unit air density and unit aerosol mass mixing ratio
        f_sad  = 0.01 * 4 * PI * rmed**2 * fnum
  
        ! radius in 'cm'
        radius = 100 * rmed
  
        do k = klid, km
         do j = j1, j2
           do i = i1, i2
            ad   = f_ad  * rhoa(i,j,k)             ! air number density # cm-3
            sad  = f_sad * rhoa(i,j,k) * q(i,j,k)  ! surface area density cm2 cm-3
  
            kan(i,j,k) = sktrs_sslt(temp(i,j,k), sad, ad, radius, PI, RUNIV, fMassHNO3)
           end do
         end do
        end do
  
     end subroutine SSLT_reaction_rate
  
  !============================================================================
  !BOP
  !
  ! !IROUTINE: apportion_reaction_rate
  !
  ! !INTERFACE:
     subroutine apportion_reaction_rate (i1, i2, j1, j2, km, kan, kan_total)
  
  ! !DESCRIPTION:
  
  ! !USES:
     implicit NONE
  
  
     integer, intent(in) :: i1, i2, j1, j2, km
  
     real, dimension(i1:i2,j1:j2,km), intent(inout) :: kan
     real, dimension(i1:i2,j1:j2,km), intent(in)    :: kan_total
  !EOP
  !------------------------------------------------------------------------------------
  !  Begin..
  
     where (kan_total > tiny(kan_total))
         kan = kan / kan_total
     else where
         kan = 0.0
     end where
  
     end subroutine apportion_reaction_rate
  
  !============================================================================
  !BOP
  !
  ! !IROUTINE: sktrs_hno3
  !
  ! !INTERFACE:
     function sktrs_hno3 ( tk, rh, sad, ad, radA, pi, rgas, fMassHNO3 )
  
  ! !DESCRIPTION:
  ! Below are the series of heterogeneous reactions
  ! The reactions sktrs_hno3n1, sktrs_hno3n2, and sktrs_hno3n3 are provided
  ! as given by Huisheng Bian.  As written they depend on knowing the GOCART
  ! structure and operate per column but the functions themselves are
  ! repetitive.  I cook up a single sktrs_hno3 function which is called per
  ! grid box per species with an optional parameter gamma being passed.
  ! Following is objective:
  ! loss rate (k = 1/s) of species on aerosol surfaces
  !
  ! k = sad * [ radA/Dg +4/(vL) ]^(-1)
  !
  ! where
  ! Dg = gas phase diffusion coefficient (cm2/s)
  ! L = sticking coefficient (unitless)  = gamma
  ! v = mean molecular speed (cm/s) = [ 8RT / (pi*M) ]^1/2
  !
  ! radA/Dg = uptake by gas-phase diffusion to the particle surface
  ! 4/(vL) = uptake by free molecular collisions of gas molecules with the surface
  
     implicit none
  
  ! !INPUT PARAMETERS:
     real, intent(in) ::  tk   ! temperature [K]
     real, intent(in) ::  rh   ! fractional relative humidity [0 - 1]
     real, intent(in) ::  sad  ! aerosol surface area density [cm2 cm-3]
     real, intent(in) ::  ad   ! air number concentration [# cm-3]
     real, intent(in) ::  radA ! aerosol radius [cm]
  
     real  :: pi   ! pi constant
     real  :: rgas ! ideal gas constant [J/(K*mol)]
     real  :: fMassHNO3 ! gram molecular weight of HNO3
     real :: sktrs_hno3
  
  ! !Local Variables
     real, parameter :: fmassHNO3_hno3 = 63.013
  
  !   REAL,  PARAMETER :: GAMMA_HNO3 = 0.1
     REAL,  PARAMETER :: GAMMA_HNO3 = 1.0e-3
  !   REAL,  PARAMETER :: GAMMA_HNO3 = 5.0e-4
  
     real :: dfkg
     real :: avgvel
     real :: gamma
     real :: f_rh
     real :: sqrt_tk
  !   real(kind=r8) :: pi_dp = pi
  !   real(kind=r8) :: rgas_dp = rgas
  
  !   real, parameter :: p_dfkg   = sqrt(3.472e-2 + 1.0/fmassHNO3)
  !   real, parameter :: p_avgvel = sqrt(8.0 * rgas_dp * 1000.0 / (pi_dp * fmassHNO3))
  
     real(kind=r8) :: pi_dp
     real(kind=r8) :: rgas_dp
  
     real :: p_dfkg
     real :: p_avgvel
  
  !EOP
  !------------------------------------------------------------------------------------
  !  Begin..
  
     pi_dp = pi
     rgas_dp = rgas
     p_dfkg   = sqrt(3.472e-2 + 1.0/fmassHNO3)
     p_avgvel = sqrt(8.0 * rgas_dp * 1000.0 / (pi_dp * fmassHNO3))
  
        ! RH factor - Figure 1 in Duncan et al. (2010)
        f_rh = 0.03
  
        if (rh >= 0.1 .and. rh < 0.3)       then
           f_rh = 0.03 + 0.8  * (rh - 0.1)
        else if (rh >= 0.3 .and. rh < 0.5 ) then
           f_rh = 0.19 + 2.55 * (rh - 0.3)
        else if (rh >= 0.5 .and. rh < 0.6)  then
           f_rh = 0.7  + 3.0  * (rh - 0.5)
        else if (rh >= 0.6 .and. rh < 0.7)  then
           f_rh = 1.0  + 3.0  * (rh - 0.6)
        else if (rh >= 0.7 .and. rh < 0.8)  then
           f_rh = 1.3  + 7.0  * (rh - 0.7)
        else if (rh >= 0.8 )                then
           f_rh = 2.0
        end if
  
  !     Following uptake coefficients of Liu et al.(2007)
        gamma = gamma_hno3 * f_rh
  
        sqrt_tk = sqrt(tk)
  
  !     calculate gas phase diffusion coefficient (cm2/s)
        dfkg = 9.45e17 / ad * sqrt_tk * p_dfkg
  
  !     calculate mean molecular speed (cm/s)
        avgvel = 100.0 * sqrt_tk * p_avgvel
  
  !     calculate rate coefficient
        sktrs_hno3 = sad / ( 4.0 / (gamma * avgvel) + radA / dfkg )
  
        END FUNCTION sktrs_hno3
  
  !============================================================================
  !BOP
  !
  ! !IROUTINE: sktrs_sslt
  !
  ! !INTERFACE:
     function sktrs_sslt ( tk, sad, ad, radA, pi, rgas, fMassHNO3 )
  
  ! !DESCRIPTION:
  ! Below are the series of heterogeneous reactions
  ! The reactions sktrs_hno3n1, sktrs_hno3n2, and sktrs_hno3n3 are provided
  ! as given by Huisheng Bian.  As written they depend on knowing the GOCART
  ! structure and operate per column but the functions themselves are
  ! repetitive.  I cook up a single sktrs_hno3 function which is called per
  ! grid box per species with an optional parameter gamma being passed.
  ! Following is objective:
  ! loss rate (k = 1/s) of species on aerosol surfaces
  !
  ! k = sad * [ radA/Dg +4/(vL) ]^(-1)
  !
  ! where
  ! Dg = gas phase diffusion coefficient (cm2/s)
  ! L = sticking coefficient (unitless)  = gamma
  ! v = mean molecular speed (cm/s) = [ 8RT / (pi*M) ]^1/2
  !
  ! radA/Dg = uptake by gas-phase diffusion to the particle surface
  ! 4/(vL) = uptake by free molecular collisions of gas molecules with the surface
  
     implicit none
  
  ! !INPUT PARAMETERS:
     real  :: tk   ! temperature [K]
     real  :: sad  ! aerosol surface area density [cm2 cm-3]
     real  :: ad   ! air number concentration [# cm-3]
     real  :: radA ! aerosol radius [cm]
     real  :: sktrs_sslt
     real  :: pi   ! pi constant
     real  :: rgas ! ideal gas constant [J/(K*mol)]
     real  :: fMassHNO3 ! gram molecular weight of HNO3
  !   real(kind=r8), optional  :: gammaInp ! optional uptake coefficient (e.g., 0.2 for SS, else calculated)
  
  !  Locals
     REAL,  PARAMETER :: GAMMA_SSLT = 0.1e0
  
     real :: dfkg
     real :: avgvel
     real :: sqrt_tk
  
     real(kind=r8) :: pi_dp
     real(kind=r8) :: rgas_dp
  
     real :: p_dfkg
     real :: p_avgvel
  
  !EOP
  !------------------------------------------------------------------------------------
  !  Begin..
     pi_dp = pi
     rgas_dp = rgas
  
     p_dfkg   = sqrt(3.472e-2 + 1.0/fmassHNO3)
     p_avgvel = sqrt(8.0 * rgas_dp * 1000.0 / (pi_dp * fmassHNO3))
  
  !  Initialize
     sqrt_tk = sqrt(tk)
  
  !     calculate gas phase diffusion coefficient (cm2/s)
        dfkg = 9.45e17 / ad * sqrt_tk * p_dfkg
  
  !     calculate mean molecular speed (cm/s)
        avgvel = 100.0 * sqrt_tk * p_avgvel
  
  !     calculate rate coefficient
        sktrs_sslt = sad / ( 4.0 / (gamma_sslt * avgvel) + radA / dfkg )
  
     end function sktrs_sslt
  
  !------------------------------------------------------------------------------
  end module mo_gocart
  