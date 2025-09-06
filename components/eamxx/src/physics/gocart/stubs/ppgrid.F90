!==================================================================================================
! module ppgrid
!
! This module is a stub for the EAM 'ppgrid' module. For GOCART and many of its dependencies,
! it is primarily used to retrieve grid size and related information for sizing arrays. In
! EAM, this data is known at compile time, but this is presently not garaunteed for EAMxx. 
! An initialization routine must now be called during the 'AtmosphereProcess' 
! 'set_grids' or 'initialize_impl' implementation before downstream processes can use the grid 
! size.
!
! Direct dependees:
!   - physconst:      pcols, pver, pverp, begchunk, endchunk 
!
!       Daniel Fust (Aug 2025)
!==================================================================================================
module ppgrid
  implicit none

  !-- module variables are publically accessible and saved --
  public
  save

  !-- module variables --
  integer :: pcols    = -1
  integer :: pver     = -1
  integer :: pverp    = -1
  integer :: begchunk = -1
  integer :: endchunk = -1

  !-- procedure access --
  public :: init_ppgrid

  contains

    !------------------------------------------------------------------------------------------------
    ! subroutine init_ppgrid
    !
    ! Sets the module varibles
    !------------------------------------------------------------------------------------------------
    subroutine init_ppgrid(ncols, nver)
      use mpi
#ifdef SPMD
      use mpishorthand,    only: mpicom
#endif
      implicit none
      !---------- Dummy Arguments -------------
      integer, intent(in) :: ncols
      integer, intent(in) :: nver
      !integer, intent(in) :: verp
      !integer, intent(in) :: beg_chunk
      !integer, intent(in) :: end_chunk
      !----------- Local variables ------------
      integer :: ierr
      !----------------------------------------

      pver  = nver
      pverp = pver+1
      !-- largest number of columns taken as pcols --
#ifdef SPMD
      call MPI_ALLREDUCE(ncols,pcols,1,MPI_INTEGER, MPI_MAX, mpicom, ierr)
#else
      call MPI_ALLREDUCE(ncols,pcols,1,MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
#endif 

      !-- Assuming there is just one chunk... --
      begchunk = 1
      endchunk = ncols !<-- not 'pcols' in case there are fewer columns on this processor
    end subroutine init_ppgrid
end module ppgrid  