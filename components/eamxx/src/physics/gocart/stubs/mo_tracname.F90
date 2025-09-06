!==================================================================================================
! module mo_tracname
!
! This is a complete copy of the EAM 'mo_tracname' module. This module contains an array of strings
! of the names of the tracer species which enables a mapping between a species name and its index.
!
! The name array 'solsym' is set at compile time by the subroutine: 'set_sim_dat' from the
! 'mo_sim_dat' module from: ${EAM}/src/chemistry/pp_chemuci_linozv3_mam5_vbs
!
! Daniel Fust (Aug 2025)
!==================================================================================================
module mo_tracname
  !-----------------------------------------------------------
  ! ... List of advected and non-advected trace species, and
  !     surface fluxes for the advected species.
  !-----------------------------------------------------------
  use chem_mods, only : grpcnt, gas_pcnst
  implicit none
    
  character(len=16) :: solsym(gas_pcnst)   ! species names
    
end module mo_tracname
    