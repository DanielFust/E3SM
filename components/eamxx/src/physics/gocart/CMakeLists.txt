# Source files for GOCART process
set(GOCART_SRCS
  eamxx_gocart_process_interface.cpp
  surface-emission/eamxx_gocart_surface_emission.cpp
  eam-interface/eam_bridge_c2f.F90
  sslt-interface/sslt_sections_c2f.F90
  setsox-interface/mo_setsox_c2f.F90
  gocart-interface/mo_gocart_c2f.F90
)

# E3SM root directory
set(E3SM_ROOT ${SCREAM_BASE_DIR}/../../)

# EAM root directory
set(EAM_ROOT ${SCREAM_BASE_DIR}/../eam/)


# Share library for low level utils
# Note: This library is (hopefully) deprecated.
#       The utilities in this standalone library shoud in principle 
#       be available from the 'csm_share' library and by proxy
#       via 'scream_share' 
set(SHR_UTILS_SRC
    ${E3SM_ROOT}/share/util/shr_file_mod.F90
    ${E3SM_ROOT}/share/util/shr_kind_mod.F90
    ${E3SM_ROOT}/share/util/shr_const_mod.F90
    ${E3SM_ROOT}/share/util/shr_sys_mod.F90
    ${E3SM_ROOT}/share/util/shr_log_mod.F90
    ${E3SM_ROOT}/share/util/shr_strconvert_mod.F90
    ${E3SM_ROOT}/share/util/shr_timer_mod.F90
    ${E3SM_ROOT}/share/util/shr_abort_mod.F90
    ${E3SM_ROOT}/share/util/shr_mpi_mod.F90
    ${E3SM_ROOT}/share/util/shr_nl_mod.F90
    ${E3SM_ROOT}/share/util/shr_string_mod.F90
    ${E3SM_ROOT}/share/util/shr_taskmap_mod.F90
    stubs/share/shr_infnan_mod.F90
    stubs/share/shr_assert_mod.F90
)
add_library(shr_utils ${SHR_UTILS_SRC})
target_include_directories(shr_utils
    PRIVATE ${E3SM_ROOT}/share/include stubs/share
)


# Build 'sslt_sections' external library
set(SSLT_EXTERNAL_SRC
    ${EAM_ROOT}/src/chemistry/aerosol/sslt_sections.F90
)
add_library(sslt_sections ${SSLT_EXTERNAL_SRC})
# ${E3SM_ROOT}/share/util/shr_kind_mod.F90  # unnecessary to include


# ESMF timing library: Not clear if this is built by default by some other name. If so, this may be replacable.
set(ESMF_SRC
    ${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF_AlarmClockMod.F90
    ${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF_AlarmMod.F90
    ${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF_BaseMod.F90
    ${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF_BaseTimeMod.F90
    ${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF_CalendarMod.F90
    ${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF_ClockMod.F90
    ${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF_FractionMod.F90
    ${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF_ShrTimeMod.F90
    ${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF_Stubs.F90
    ${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF_TimeIntervalMod.F90
    ${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF_TimeMod.F90
    ${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF.F90
    #${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF_TimeMgr.inc
    #${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF_Macros.inc
    ${E3SM_ROOT}/share/esmf_wrf_timemgr/MeatMod.F90
    ${E3SM_ROOT}/share/esmf_wrf_timemgr/wrf_error_fatal.F90
    ${E3SM_ROOT}/share/esmf_wrf_timemgr/wrf_message.F90
)
add_library(esmf ${ESMF_SRC})
target_include_directories(esmf
    PRIVATE ${E3SM_ROOT}/share/esmf_wrf_timemgr
)
target_link_libraries(esmf csm_share)




# Build 'eam_share' library for EAM infrastruture required by multiple bridged items
set(EAM_SHARE_SRC
    ${EAM_ROOT}/src/control/mpishorthand.F
    ${EAM_ROOT}/src/physics/crm/samxx/test/perf_mod.F90
    ${EAM_ROOT}/src/utils/namelist_utils.F90
    ${EAM_ROOT}/src/control/cam_logfile.F90
    ${EAM_ROOT}/src/utils/cam_abortutils.F90
    ${EAM_ROOT}/src/control/units.F90
    ${EAM_ROOT}/src/control/physconst.F90
    ${EAM_ROOT}/src/utils/spmd_utils.F90
    stubs/ppgrid.F90
    stubs/pmgrid.F90
    stubs/mo_tracname.F90
    stubs/constituents.F90
    ${EAM_ROOT}/src/chemistry/pp_chemuci_linozv3_mam5_vbs/chem_mods.F90
    ${EAM_ROOT}/src/chemistry/pp_chemuci_linozv3_mam5_vbs/mo_sim_dat.F90
    ${EAM_ROOT}/src/chemistry/pp_chemuci_linozv3_mam5_vbs/m_spc_id.F90
    stubs/modal_aero_data.F90              #<-- requires preprocessor defines MODAL_AERO_*MODE and VBS_SOA
    ${EAM_ROOT}/src/physics/rrtmg/radconstants.F90
    ${EAM_ROOT}/src/chemistry/utils/mo_constants.F90
    ${EAM_ROOT}/src/chemistry/mozart/mo_chem_utls.F90
    ${EAM_ROOT}/src/chemistry/mozart/mo_mass_xforms.F90
    stubs/modal_aero_data.F90
    stubs/modal_aero_initialize_data.F90
)
add_library(eam_share ${EAM_SHARE_SRC})
target_include_directories(eam_share 
    INTERFACE ${E3SM_ROOT}/share/include
)
target_compile_definitions(eam_share
    PUBLIC MODAL_AERO_5MODE
    INTERFACE VBS_SOA
) 
target_link_libraries(eam_share csm_share)





# Build 'mo_setsox' library
set(MO_SETSOX_SRC
    ${EAM_ROOT}/src/chemistry/aerosol/cldaero_mod.F90
    stubs/sox_cldaero_mod.F90
    setsox-interface/mo_setsox.F90
)
add_library(mo_setsox ${MO_SETSOX_SRC})
target_link_libraries(mo_setsox eam_share)



# Build 'mo_gocart' library
set(MO_GOCART_EXTERNAL_SRC
    #stubs/mo_sim_dat.F90
    ${EAM_ROOT}/src/chemistry/mozart/gas_wetdep_opts.F90
    stubs/seq_drydep_mod.F90
    gocart-interface/mo_gocart_refactor.F90
)
add_library(mo_gocart ${MO_GOCART_EXTERNAL_SRC} ${SHR_UTILS_SRC})
#target_include_directories(mo_gocart
#    PRIVATE ${E3SM_ROOT}/share/include stubs/share
#)
target_link_libraries(mo_gocart eam_share)


#target_include_directories(mo_gocart
#    PRIVATE ${E3SM_ROOT}/share/esmf_wrf_timemgr
#)


# Create library for GOCART
set(GOCART_LIBS "gocart")
add_library(gocart ${GOCART_SRCS})
target_compile_definitions(gocart PUBLIC EAMXX_HAS_GOCART)
target_link_libraries(gocart 
                      physics_share 
                      scream_share 
                      sslt_sections 
                      mo_setsox
                      mo_gocart)
target_compile_options(gocart PUBLIC)

include(ScreamUtils)

# Add this to eamxx_physics
target_link_libraries(eamxx_physics INTERFACE gocart)














    #stubs/share/shr_infnan_mod.F90
    #stubs/share/shr_assert_mod.F90
     #${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF_AlarmClockMod.F90
    #${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF_AlarmMod.F90
    #${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF_BaseMod.F90
    #${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF_BaseTimeMod.F90
    #${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF_CalendarMod.F90
    #${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF_ClockMod.F90
    #${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF_FractionMod.F90
    #${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF_ShrTimeMod.F90
    #${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF_Stubs.F90
    #${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF_TimeIntervalMod.F90
    #${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF_TimeMod.F90
    #${E3SM_ROOT}/share/esmf_wrf_timemgr/ESMF.F90
    #${E3SM_ROOT}/share/esmf_wrf_timemgr/MeatMod.F90
     #${E3SM_ROOT}/share/util/shr_file_mod.F90
    #${E3SM_ROOT}/share/util/shr_kind_mod.F90
    #${E3SM_ROOT}/share/util/shr_const_mod.F90
    #${E3SM_ROOT}/share/util/shr_sys_mod.F90
    #${E3SM_ROOT}/share/util/shr_log_mod.F90
    #${E3SM_ROOT}/share/util/shr_strconvert_mod.F90
    #${E3SM_ROOT}/share/util/shr_timer_mod.F90
    #${E3SM_ROOT}/share/util/shr_abort_mod.F90
    #${E3SM_ROOT}/share/util/shr_mpi_mod.F90
    #${E3SM_ROOT}/share/util/shr_nl_mod.F90
    #stubs/share/shr_string_mod.F90
    #${E3SM_ROOT}/share/util/shr_taskmap_mod.F90