!==================================================================================================
! module mo_gocart_c2f
!
! This module serves as the Fortran-side interface from the 'mo_gocart' module used to drive
! GOCART-2G and the GOCART AtmosphereProcess in EAMxx
!
! Daniel Fust (Aug 2025)
!==================================================================================================
module mo_gocart_c2f
  use iso_c_binding  
  implicit none  

  contains

    !----------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------
    subroutine init_mo_gocart_c2f(cols,levs,nlon,nlat,eam_nlfile,nl_len) bind(c, name="init_mo_gocart_c2f")
      use iso_c_binding,   only: c_char, c_int
      use mpi
      use spmd_utils,      only: spmdinit
      use gas_wetdep_opts, only: gas_wetdep_readnl, gas_wetdep_cnt, gas_wetdep_list
      use mo_tracname, only: solsym
      use chem_mods,   only: gas_pcnst
      use mo_sim_dat,  only: set_sim_dat
      use ppgrid,      only: init_ppgrid
      use pmgrid,      only: init_pmgrid
      !use mo_gocart,   only: 
#ifdef SPMD
      use mpishorthand,    only: mpicom
#endif
      implicit none
      !---------------- Arguments ------------------
      integer,           intent(in), value        :: cols       ! number of columns on processor
      integer,           intent(in), value        :: levs       ! number of elevation strata on processor
      integer,           intent(in), value        :: nlon       ! number of longitudes on processor
      integer,           intent(in), value        :: nlat       ! number of lattitudes on processor
      character(c_char), intent(in), dimension(*) :: eam_nlfile ! filepath for file containing EAM namelist input
      integer(c_int),    intent(in), value        :: nl_len     ! length of namelist string
      !-------------- Local Variables --------------
      character(len=nl_len) :: f90_nlfile  ! Fortran-style string for EAM namelist input
      integer :: i

      integer :: comm, rank, ierr
      !---------------------------------------------
      
      !-- We must first initialize all EAM infrastructure required by mo_gocart --



      !-- mo_gocart deposition velocity initalization procedure --
      !call gc_dvel_inti_fromlnd() 
      

    end subroutine init_mo_gocart_c2f
end module mo_gocart_c2f