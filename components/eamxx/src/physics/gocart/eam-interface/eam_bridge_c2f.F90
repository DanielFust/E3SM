!--------------------------------------------------------------------------------------------------
! module eam_bridge_c2f
!
! This module provides a C/Fortran bridge to infrastructure in EAM potentially required by multiple
! bridged modules.
!--------------------------------------------------------------------------------------------------
module eam_bridge_c2f
  implicit none

  contains

    !----------------------------------------------------------------------------------------------
    ! subroutine eam_bridge_init(...)
    !
    ! Calls initiazation routines for EAM data structures and processes
    !----------------------------------------------------------------------------------------------
    subroutine eam_bridge_init_c2f(cols,levs,nlon,nlat,eam_nlfile,nl_len) bind(c,name="eam_bridge_init_c2f")
      use mpi
      use iso_c_binding,   only: c_int, c_char
      use spmd_utils,      only: spmdinit
      use gas_wetdep_opts, only: gas_wetdep_readnl, gas_wetdep_cnt, gas_wetdep_list
      use mo_tracname,     only: solsym
      use chem_mods,       only: gas_pcnst
      use mo_sim_dat,      only: set_sim_dat
      use ppgrid,          only: init_ppgrid
      use pmgrid,          only: init_pmgrid
      use modal_aero_initialize_data, only: modal_aero_initialize
      use modal_aero_data, only: lptr_so4_cw_amode, lptr_msa_cw_amode, lptr_nh4_cw_amode, ntot_amode
#ifdef SPMD
      use mpishorthand,    only: mpicom
#endif
      implicit none
      !---------------- Arguments ------------------
      integer(c_int),    intent(in), value        :: cols       ! number of columns on processor
      integer(c_int),    intent(in), value        :: levs       ! number of elevation strata on processor
      integer(c_int),    intent(in), value        :: nlon       ! number of longitudes on processor
      integer(c_int),    intent(in), value        :: nlat       ! number of lattitudes on processor
      character(c_char), intent(in), dimension(*) :: eam_nlfile ! filepath for file containing EAM namelist input
      integer(c_int),    intent(in), value        :: nl_len     ! length of namelist string
      !-------------- Local Variables --------------
      character(len=nl_len) :: f90_nlfile  ! Fortran-style string for EAM namelist input
      integer :: i
      !---------------------------------------------

            ! must first call this routine from spmd_utils to initialize MPI data stuff
#ifdef SPDM
      call spmdinit(mpicom)
#else
      call spmdinit(MPI_COMM_WORLD)
#endif

      !-- set grid size data --
      call init_ppgrid(cols,levs)
      call init_pmgrid(nlon,nlat,levs)

      !-- sets tracer names and other data --
      call set_sim_dat()

      !-- translate c-string to Fortan string 
      !   c-strings are terminated by '\0' whereas Fortran strings a fixed length 
      !   and padded with spaces
      do i=1,nl_len
        f90_nlfile(i:i) = eam_nlfile(i)
      enddo 

      !-- call to initalize wet dep list and method from input file
      call gas_wetdep_readnl(f90_nlfile)

      !-- checking to see if tracer list is initialized --
      !write(*,*) 'tracers:'
      !do i=1,gas_pcnst
      !  write(*,*) solsym(i)
      !enddo

      !-- initialize modal aero data --
      !   NOTE: This is currently only factored to support running `mo_setsox`
      !         cloud aqueous chemistry. If other infrastructure is needed
      !         to support wet/dry deposition, intermodal interactions of
      !         aerosols, and seasalt interactions, then this stub will need
      !         to be revised and others created.
      call modal_aero_initialize() !<-- FIX ME: not working to set IDs
      !--- sketchy workaround for setting indices for sox_cldaero_mod ----
      !    NOTE: the indices work differently for 7 modes, so this
      !          would need revision if other chemistry modal models 
      !          are used.
      lptr_msa_cw_amode(1:ntot_amode) = -1 !<-- not used in pp_chemuci_linozv3_mam5?
      lptr_so4_cw_amode(1:ntot_amode) = -1
      lptr_nh4_cw_amode(1:ntot_amode) = -1 !<-- not used in pp_chemuci_linozv3_mam5?
      do i=1,gas_pcnst
        if (trim(solsym(i)) == "so4_a1") then
          lptr_so4_cw_amode(1) = i
          write(*,*) "so4_a1 index:", i
        endif
        if (trim(solsym(i)) == "so4_a2") then
          lptr_so4_cw_amode(2) = i
          write(*,*) "so4_a2 index:", i
        endif
        if (trim(solsym(i)) == "so4_a3") then 
          lptr_so4_cw_amode(3) = i
          write(*,*) "so4_a3 index:", i
        endif
        if (trim(solsym(i)) == "so4_a4") then 
          lptr_so4_cw_amode(4) = i
          write(*,*) "so4_a4 index:", i
        endif
#if (defined MODAL_AERO_5MODE)
        if (trim(solsym(i)) == "so4_a5") then
          lptr_so4_cw_amode(5) = i !<-- 5 not 4. Not sure what 4 is supposed to be
          write(*,*) "so4_a5 index:", i
        endif
#endif        
      enddo
      
    end subroutine eam_bridge_init_c2f



    !----------------------------------------------------------------------------------------------
    ! function gas_pcnst()
    !
    ! Provides read access to `gas_pcnst` module parameter for number of advected gas species 
    ! expected by bridged EAM code.
    !----------------------------------------------------------------------------------------------
    function gas_pcnst_c2f() result(num_gas_species) bind(c,name="gas_pcnst")
      use iso_c_binding, only: c_int
      use chem_mods,     only: gas_pcnst
      implicit none
      !---------------------------------
      integer(c_int) :: num_gas_species
      !---------------------------------
      num_gas_species = gas_pcnst
    end function gas_pcnst_c2f






    !----------------------------------------------------------------------------------------------
    ! function get_spc_ndx_c2f
    !
    ! Bridge to `get_spc_ndx` function from `mo_chem_utls` module to retrieve the index of a 
    ! chemical species within the contiguous array expected by EAM.
    !----------------------------------------------------------------------------------------------  
    function get_spc_ndx_c2f( spc_name, len ) result(index) bind(c,name="get_spc_ndx_c2f")
      use iso_c_binding, only: c_int, c_char
      use mo_chem_utls,  only: get_spc_ndx
      implicit none
      !----------- Dummy Arguments --------------
      character(c_char), intent(in), dimension(*) :: spc_name ! species name
      integer(c_int),    intent(in), value        :: len      ! length of name string
      integer(c_int)                              :: index
      !----------- Local Variables --------------
      character(len=len) :: str_f90  ! Fortran-style string for species name
      integer            :: i
      !------------------------------------------
      do i=1,len
        str_f90(i:i) = spc_name(i)
      enddo
      index = get_spc_ndx(trim(str_f90))
    end function get_spc_ndx_c2f





    !----------------------------------------------------------------------------------------------
    ! function solsym_c2f(...) 
    !
    ! Provides a read-only bridge to the `solsym` list of tracer names in the `mo_tracname` 
    ! module. Returns the name of an EAM advected gas tracer species of the corresponding index 
    ! (1 <= index <= gas_pcnst) as a C-style string.
    !----------------------------------------------------------------------------------------------
    subroutine solsym_c2f(index,spc_name) bind(c,name="solsym_c2f")
      use iso_c_binding, only: c_int, c_char
      use mo_tracname,   only: solsym
      implicit none
      !------------------- Dummy Arguments --------------------------
      integer(c_int),          intent(in), value        :: index
      character(kind=c_char), intent(out), dimension(*) :: spc_name
      !------------------- Local Variables --------------------------
      integer           :: i, n
      character(len=16) :: spc_name_f90
      !--------------------------------------------------------------

      !-- The Fortran string is fixed length (16 characters) but we must add
      !   a "\0" termination suffix for our C-style string. The C-side interface
      !   is garaunteed to have a char[18] buffer to recieve the string
      spc_name_f90 = trim(solsym(index)) !<-- in case there are leading spaces
      n = len_trim(spc_name_f90)
      do i=1,n
        spc_name(i:i) = spc_name_f90(i:i)
      enddo
      spc_name(n+1:n+1) = char(0) !<-- null terminator
    end subroutine solsym_c2f



    !----------------------------------------------------------------------------------------------
    ! subroutine mmr2vmr_c2f(...)
    !
    ! Interface for EAM procedure for converting mass mixing ratio (mmr), "units" expected for 
    ! tracer advection, to volume mixing ratio (vmr), "units" expected for many chemical processes.
    !----------------------------------------------------------------------------------------------
    subroutine mmr2vmr_c2f(mmr, vmr, mbar, ncol) bind(c,name="mmr2vmr_c2f")
      use iso_c_binding,  only: c_int, c_double
      use ppgrid,         only: pcols, pver
      use chem_mods,      only: gas_pcnst
      use mo_mass_xforms, only: mmr2vmr
      implicit none
      !----------------- Dummy Arguments ---------------------
      real(c_double), intent(in)          :: mmr(pcols,pver,gas_pcnst) !mass mixing ratio of tracers
      real(c_double), intent(inout)       :: vmr(ncol,pver,gas_pcnst)  !volume mixing ratio of tracers
      real(c_double), intent(in)          :: mbar(ncol,pver)           !mean wet atmospheric mass (amu)
      integer(c_int), intent(in),   value :: ncol                      !number of columns
      !-------------------------------------------------------
      call mmr2vmr(mmr,vmr,mbar,ncol)
    end subroutine mmr2vmr_c2f


    !----------------------------------------------------------------------------------------------
    ! subroutine vmr2mmr_c2f(...)
    !
    ! Interface for EAM procedure for converting volume mixing ratio (vmr), "units" expected for 
    ! many chemical processes, to mass mixing ratio (mmr), "units" expected for tracer advection.
    !----------------------------------------------------------------------------------------------
    subroutine vmr2mmr_c2f(vmr, mmr, mbar, ncol) bind(c,name="vmr2mmr_c2f")
      use iso_c_binding,  only: c_int, c_double
      use ppgrid,         only: pcols, pver
      use chem_mods,      only: gas_pcnst
      use mo_mass_xforms, only: vmr2mmr
      implicit none
      !----------------- Dummy Arguments ---------------------
      real(c_double), intent(in)          :: vmr(ncol,pver,gas_pcnst)  !volume mixing ratio of tracers
      real(c_double), intent(inout)       :: mmr(pcols,pver,gas_pcnst) !mass mixing ratio of tracers
      real(c_double), intent(in)          :: mbar(ncol,pver)           !mean wet atmospheric mass (amu)
      integer(c_int), intent(in),   value :: ncol                      !number of columns
      !-------------------------------------------------------
      call vmr2mmr(vmr,mmr,mbar,ncol)
    end subroutine vmr2mmr_c2f
end module eam_bridge_c2f