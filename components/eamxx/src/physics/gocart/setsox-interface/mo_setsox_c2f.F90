!==================================================================================================
! module mo_setsox_c2f
!
! This module serves as an interface between the mo_setsox module, adapted from EAM, and EAMxx
!==================================================================================================
module mo_setsox_c2f
  implicit none
  
  contains

    !----------------------------------------------------------------------------------------------
    ! subroutine sox_inti_c2f(...)
    !
    ! Calls mo_setsox initialization routine. Also initializes sox_cldaero_mod module that is 
    ! required for aqueous chemistry calculation.
    ! Note: Based on the use of the pp_chemuci_linozv3_mam5, it is likely that 
    !       `use_modal_aerosols` should be set to `.true.` however it is less clear what `use_mmf`
    !       and `use_ecpp` should be set to.
    !----------------------------------------------------------------------------------------------
    subroutine sox_inti_c2f(use_modal_aerosols, use_mmf, use_ecpp) bind(c,name="sox_inti_c2f")
      use iso_c_binding,   only: c_bool
      use mo_setsox,       only: sox_inti 
      use sox_cldaero_mod, only: sox_cldaero_init 
      implicit none
      !--------------- Dummy Arguments -----------------
      logical(c_bool), intent(in), value :: use_modal_aerosols
      logical(c_bool), intent(in), value :: use_mmf  ! (???)
      logical(c_bool), intent(in), value :: use_ecpp ! explicit cloud parameterized pollutants
      !-------------------------------------------------
      write(*,*) "In sox_inti_c2f..."
      call sox_cldaero_init()
      call sox_inti( logical(use_modal_aerosols), &
                     logical(use_mmf),            &
                     logical(use_ecpp) )
    end subroutine sox_inti_c2f


    !----------------------------------------------------------------------------------------------
    ! subroutine setsox_c2f
    !
    ! Bridge to primary `SETSOX` subprogram of the `mo_setsox` module for cloud aqueous 
    ! computations.
    !
    ! Notes:
    !   - `invariants` has been removed as it is appears unused in `SETSOX` when the 
    !     `pp_chemuci_linozv3_mam5` aerosol model and due to the difficuly in correctly 
    !     initializing it in an EAMxx bridge.
    !
    !   - In the bridge, each processor is prescribed as a single chunk so `lchnk` is simply
    !     set to unity and excluded as an argument for clarity and convenience. In `mo_setsox`
    !     it appears that the chunk ID is only used for writing history data, which has been
    !     excluded from the bridge, so the `lchnk` may be effectively unused regardless, and
    !     unimportant to the function of `SETSOX`
    !
    !   - The 'q' arrays of species/tracers contain only the tracers specified in the aerosol 
    !     model, therefore the offset `loffset` is set to zero.
    !----------------------------------------------------------------------------------------------
    subroutine setsox_c2f(ncol, nlev, &! lchnk, loffset, 
                          dtime,  press,  pdel, &
                          tfld,   mbar,   lwc,   &
                          cldfrc, cldnum, xhnm,   &! invariants, &
                          qcw, qin) bind(c,name="setsox_c2f")
      use iso_c_binding, only: c_int, c_double
      use mo_setsox,     only: SETSOX
      use chem_mods,     only: gas_pcnst
      !--------------- Dummy Arguments ------------------
      integer(c_int),         intent(in), value :: ncol              ! num of columns in chunk
      integer(c_int),         intent(in), value :: nlev              ! num of elevation levels
      !integer(c_int),         intent(in), value :: lchnk             ! chunk id
      !integer(c_int),         intent(in), value    :: loffset           ! offset of chem tracers in the advected tracers array
      real(c_double),         intent(in), value :: dtime             ! time step (sec)
      real(c_double),         intent(in)        :: press(ncol,nlev)!(:,:)        ! midpoint pressure ( Pa )
      real(c_double),         intent(in)        :: pdel(ncol,nlev)!(:,:)         ! pressure thickness of levels (Pa)
      real(c_double),         intent(in)        :: tfld(ncol,nlev)!(:,:)         ! temperature
      real(c_double),         intent(in)        :: mbar(ncol,nlev)!(:,:)         ! mean wet atmospheric mass ( amu )
      !real(c_double), target, intent(in)        :: lwc(:,:)          ! cloud liquid water content (kg/kg)
      !real(c_double), target, intent(in)        :: cldfrc(:,:)       ! cloud fraction
      real(c_double),         intent(in)        :: lwc(ncol,nlev)!(:,:)          ! cloud liquid water content (kg/kg)
      real(c_double),         intent(in)        :: cldfrc(ncol,nlev)!(:,:)       ! cloud fraction
      real(c_double),         intent(in)        :: cldnum(ncol,nlev)!(:,:)       ! droplet number concentration (#/kg)
      real(c_double),         intent(in)        :: xhnm(ncol,nlev)!(:,:)         ! total atms density ( /cm**3)
      !real(c_double),         intent(in)    :: invariants(:,:,:)
      !real(c_double), target, intent(inout)     :: qcw(:,:,:)        ! cloud-borne aerosol (vmr)
      real(c_double), intent(inout)     :: qcw(ncol,nlev,gas_pcnst)!(:,:,:)        ! cloud-borne aerosol (vmr)
      real(c_double),         intent(inout)     :: qin(ncol,nlev,gas_pcnst)!(:,:,:)        ! transported species ( vmr )
      !------------------ Local -------------------------
      integer, parameter :: lchnk=1
      integer, parameter :: loffset=0
      !--------------------------------------------------
      write(*,*) 'In setsox_c2f...'
      call SETSOX(ncol,   lchnk, loffset, dtime, press,  &
                  pdel,   tfld,  mbar,    lwc,   cldfrc, &
                  cldnum, xhnm,   &!       invariants, &
                  qcw,    qin)
    end subroutine setsox_c2f
end module mo_setsox_c2f