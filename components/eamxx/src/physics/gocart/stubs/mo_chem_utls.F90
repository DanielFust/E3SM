!==================================================================================================
! module mo_chem_utls
!
! This module is a stub for the EAM 'mo_chem_utls' module.
!
! Direct dependees:
!   - mo_gocart:
!        gc_dvel_inti_fromlnd:   get_spc_ndx       
!
!       Daniel Fust (Aug 2025)
!==================================================================================================
module mo_chem_utls
  implicit none


  public :: get_spc_ndx

  contains

    !------------------------------------------------------------------------------------------------
    ! function get_spc_ndx
    !
    ! Returns the index of a species in the GOCART chemical scheme.
    !------------------------------------------------------------------------------------------------
    function get_spc_ndx(spc_name) result(spc_ndx)
      character(len=*), intent(in) :: spc_name
      integer :: spc_ndx

      ! For now, we return a dummy value. This should be replaced with actual logic to find the index.
      spc_ndx = 0
    end function get_spc_ndx
end module mo_chem_utls