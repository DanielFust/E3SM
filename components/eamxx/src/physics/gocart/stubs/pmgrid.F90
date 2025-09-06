!==================================================================================================
! module pmgrid (stub)
!
! Refactor Notes:
!   - It is unclear whether any of the compile-time data can be set in EAMxx. The 'parameter'
!     attibute has been removed.
!   - An initialization subroutine has been added to set module variables
!   - There are no lattitude pairs in EAMxx. Corresponding variables have been removed.
!==================================================================================================
module pmgrid
  !----------------------------------------------------------------------- 
  ! 
  ! Purpose: Parameters and variables related to the dynamics grid
  ! 
  ! Author: 
  ! 
  !-----------------------------------------------------------------------
    
  implicit none
    
  public 
    
  !integer, parameter :: plon   = PLON                     ! number of longitudes
  !integer, parameter :: plev   = PLEV                     ! number of vertical levels
  !integer, parameter :: plat   = PLAT                     ! number of latitudes
  !integer, parameter :: plevp  = plev + 1                 ! plev + 1
  !integer, parameter :: plnlv  = plon*plev                ! Length of multilevel field slice
    
  integer :: plon   = -1    ! number of longitudes
  integer :: plev   = -1    ! number of vertical levels
  integer :: plat   = -1    ! number of latitudes
  integer :: plevp  = -1    ! plev + 1
  integer :: plnlv  = -1    ! Length of multilevel field slice
    
       !-- These no longer exist in EAMxx --
       !integer :: beglat     ! beg. index for latitudes owned by a given proc
       !integer :: endlat     ! end. index for latitudes owned by a given proc
       !integer :: begirow    ! beg. index for latitude pairs owned by a given proc
       !integer :: endirow    ! end. index for latitude pairs owned by a given proc
  integer :: numlats                 ! number of latitudes owned by a given proc
  logical :: dyndecomp_set = .false. ! flag indicates dynamics grid has been set for history
    
!#if ( ! defined SPMD )
!       parameter (beglat   = 1)
!       parameter (endlat   = plat)
!       parameter (begirow  = 1)
!       parameter (endirow  = plat/2)
!       parameter (numlats  = plat)
!#endif

  contains

    !----------------------------------------------------------------------------------------------
    ! subroutine init_pmgrid
    !
    ! Initializes pmgrid data
    !----------------------------------------------------------------------------------------------
    subroutine init_pmgrid(nlon,nlat,nlev)
      use mpi
#ifdef SPMD
      use mpishorthand,    only: mpicom
#endif
      implicit none
      !------ Dummy Arguments -------
      integer, intent(in) :: nlon    ! number of longitudes on proc
      integer, intent(in) :: nlat    ! number of latitudes on proc
      integer, intent(in) :: nlev    ! number of elevations
      !------ Local Variables -------
      integer :: ierr
      !------------------------------

#ifdef SPMD
      call MPI_ALLREDUCE(nlon,plon,1,MPI_INTEGER, MPI_MAX, mpicom, ierr)
      call MPI_ALLREDUCE(nlat,plat,1,MPI_INTEGER, MPI_MAX, mpicom, ierr)
#else
      call MPI_ALLREDUCE(nlon,plon,1,MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(nlat,plat,1,MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
#endif
      plev  = nlev
      plnlv = plon*plev  

      numlats = nlat
    end subroutine init_pmgrid
end module pmgrid
    
    