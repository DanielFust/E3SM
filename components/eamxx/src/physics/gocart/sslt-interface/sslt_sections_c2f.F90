!--------------------------------------------------------------------------------------------------
! module sslt_sections_c2f
!
! Interface module for the EAM 'sslt_sections' Fortran module for seasalt aerosol emissions
!--------------------------------------------------------------------------------------------------
module sslt_sections_c2f
  use iso_c_binding  
  use sslt_sections, only: sslt_sections_init, fluxes, nsections, Dg, rdry
  implicit none  
  private

  !---- Exposed module procedures ----
  public :: sslt_sections_c2f_init
  public :: sslt_fluxes_c2f
  !public :: nsections_c2f
  !public :: Dg_c2f
  !public :: rdry_c2f

  contains

    !----------------------------------------------------------------------------------------------
    ! subroutine sslt_section_c2f_init()
    !
    ! Calls initilzation procedure of 'sslt_section' to set module variables
    !----------------------------------------------------------------------------------------------
    subroutine sslt_sections_c2f_init() bind(c, name="sslt_sections_init")
      implicit none  
      write(*,*) 'In sslt_sections_c2f_init...'
      call sslt_sections_init() 
    end subroutine sslt_sections_c2f_init

    !----------------------------------------------------------------------------------------------
    ! subroutine sslt_fluxes(...)
    !
    ! Estimates the flux of (???) due to seasalt.
    ! 
    ! Outputs:
    !   fi(ncol,nsections)......flux of (???) into the atmosphere
    ! Inputs:
    !   sst(ncol)...............Sea surface temperature (K)
    !   u10cubed(ncol)..........10m wind speed "cubed" (3.41) power according to according to Gong 
    !                           et al., 1997
    !   ncol....................Number of columns in the model grid
    !----------------------------------------------------------------------------------------------
    subroutine sslt_fluxes_c2f(ncol, nsec, fi, sst, u10cubed) bind(c, name="sslt_fluxes_c2f")
      implicit none
      !--------------- Arguments ----------------
      real(c_double), intent(out)       :: fi(ncol,nsec) !(:,:)
      real(c_double), intent(in)        :: sst(ncol)     !(:)
      real(c_double), intent(in)        :: u10cubed(ncol)!(:)
      integer(c_int), intent(in), value :: ncol
      integer(c_int), intent(in), value :: nsec
      !------------------------------------------
      !write(*,*) "In sslt_fluxes_c2f..."
      fi = fluxes(sst,u10cubed,ncol)
    end subroutine sslt_fluxes_c2f


    !==================== Getters for module variables ========================

    !---- I don't know what 'Dg' is ----
    !subroutine sslt_sections_Dg(Dg_vals) bind(c, name="sslt_sections_Dg")
    !  implicit none
    !  real(c_double), intent(out) :: Dg_vals(nsections)
    !  Dg_vals(1:nsections) = Dg(1:nsections)
    !subroutine sslt_sections_Dg

    !---- I don't know what 'rdry' is ----
    !subroutine sslt_sections_rdry(rdry_vals) bind(c, name="sslt_sections_rdry")
    !  implicit none
    !  real(c_double), intent(out) :: rdry_vals(nsections)
    !  rdry_vals(1:nsections) = rdry(1:nsections)
    !subroutine sslt_sections_rdry

    !---- Number of Sections ----
    function sslt_sections_nsections() result(nsec) bind(c, name="sslt_sections_nsections")
      implicit none
      integer(c_int) :: nsec
      nsec = nsections
    end function sslt_sections_nsections  

    !---- Pointer to rdry array (this appears to be set in 'sslt_sections_init') 
    !     but 'rdry' is not set as a parameter, so it may be intended to be modified
    !     externally
    ! This isn't possible without adding the 'target' attribute to 'rdry' in the
    ! 'sslt_sections' module.
    !function sslt_sections_rdry() result(rdry_c_ptr) bind("c", name="sslt_sections_rdry")
    !  real(c_ptr)         :: rdry_c_ptr
    !  real(r8),   pointer :: rdry_f_ptr(:)
    !
    !  rdry_f_ptr => rdry 
    !
    !  call c_f_pointer(rdry_c_ptr, rdry_f_ptr, [nsections])
    !function sslt_sections_rdry






end module sslt_sections_c2f