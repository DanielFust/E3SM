# Overview
This directory contains source files for a bridge to the **EAM** Fortran implementation of the **GOCART-2G** aerosol model and constituent processes, recast as an **EAMxx** `AtmosphereProcess` class. At the request of the collaborators and authors of the **EAM** implementation, the `GOCART` `AtmosphereProcess` aims to either bridge or reimplement the following processes and modules:
  - Surface and elevated emissions from prescribed file I/O
  - Seasalt online emissions from `sslt_sections.F90`
  - Cloud aqueous chemistry from `mo_setsox.F90`
  - Wet and dry deposition from `mo_gocart.F90`
This process is still, as of 9/5/2025, in the development stage and is intended to be run through the newly created unit-test, which is located at `${SCREAM_BASE_DIR}/ctest-build/full_debug/tests/single-process/gocart`.

## Code structure
The base gocart directory contains the process interface header and implementation files in addition to the `CMakeLists.txt` file and reference files such as this README and `atm_in` copied from the provided example run of **GOCART** in **EAM** upon which this bridge was based. The `default-data` subdirectory contains user input from the example case that is required to correctly initialize some of the bridged **EAM** infrastructure. Generally speaking, each major class resides in its own subdirectory containing at a minimum a header and implementation file, and optionally helper functions and classes that are not required by external classes or functions. Currently, the only example of this organizational pattern is the `surface-emission` directory containing the `SurfaceEmission` class for prescribed I/O of emission data. Shared small helper `struct`s are intended to be stored in the `eamxx_gocart_data_structs.hpp` header, though this may be refactored as a directory of headers if the number of `struct`s becomes large. All non-member functions are protected by the `gocart` namepsace in order to eliminate the risk of external name collisions.

The bridged **EAM** code is organized into subdirectories with the naming pattern: `*bridgename*-interface` and correspond to individual libraries specified in the `CMakeList.txt` file, which are linked with the primary `gocart` library during compilation. Each bridge directory contains at a minimum:
  - A `*bridgename*_c2f.F90` file containing the Fortran-side interface exposing the desired **EAM** subprograms, parameters, and variables to external *C/C++* code.
  - A `*bridgename*.hpp` file containing the *C++* side interface declaring the external Fortran function signatures and optionally additional interface functions to improve code clarity in the primary **GOCART** process interface implementation.
The `setsox-interface` and `gocart-interface` additionally contain minimally refactored versions of their respective module files in order to trim or replace dependencies that are infeasible to initialize in **EAMxx**.
The `eam-interface` is a slightly special case, as it attempts to bridge all **EAM** infrastructure that may be required directly or indirectly by multiple modules and subprograms requested for bridging by the collaborators. Each bridge is protected by an additional nested namespace, `eam_bridge`, `sslt_sections`, `mo_setsox`, and `mo_gocart`, shadowing module names where relevant, in order to provide greater clarity on the origin of a bridged function or variable.

The `stubs` directory, contains modified versions of external (mostly **EAM**) code with certain dependencies removed. These are intended to provide required functionality while removing usages of subroutines and data strucutres that are not feasible to bridge, such as the `physics_buffer`. Subprogram call signatures are largely unchanged but all exceptions are clearly commented in their stubbed source file.

## Project status and deliverables
As of 9/5/2025, functioning bridges exist for `sslt_sections` and `mo_setsox`. Additionally, bridges have been created to initialize required **EAM** data structures, allowing tracers to be declared from the output of the chemistry preprocessor (`pp_chemuci_linozv3_mam5` in this instance) to provide greater flexibility for modal aeorosol models, as well as MMR<->VMR conversions. The provided files should be heavily commented to provide clarity on the program logic as well as the intended purpose of each variable. A summary of the deliverables and their status is provided in the following table:
| **Item**               | **Status**            | **Notes**                                                        |
| :--------------------- | :-------------------- | :--------------------------------------------------------------- |
| Surface emissions      | *Partially Complete*  | Needs remapping, proper time interpolation, and unit conversion. |
| Elevated emissions     | *Incomplete*          | Similar challenges to Surface Emissions                          |
| `sslt_sections` bridge | *Complete*            | Should be refactored to use 10m windspeed provided by MCT        |
| `mo_setsox` bridge     | *Complete*            | A C++ translation was also created, but incorrectly assumed *bulk aero* rather than *modal aero* models and must be revised. |
| `mo_gocart` bridge     | *Incomplete*          | Dry deposition requires file I/O and remapping<sup>1</sup>. Both wet and dry deposition may indirectly require **EAM** modal aeorol data to be initialized or processes to run<sup>2</sup>. |

1. If the old **EAM** remapping and file I/O is to be bridged, the user will also need to provide the grid as well in the old NetCDF file format (*i.e.*  *lat*,*lon* coordinates rather than *col* coordinates). This also may provide some challenges to subroutines expecting certain data to be structured with *lat* and *lon* indices rather than the flattened *col* indices. Alternatively, this I/O and mapping could be done on the **EAMxx** side using *scorpio*, which would require a C++ translation of some (possibly just one) of the `mo_gocart` subprograms and some modification to some `mo_gocart` subroutines to allow this data to be passed from C++.
2. Some data required by **EAM** processes such as indices of arrays, names, and classes of modal aerosols are set in `modal_aero_initialize_data`, which in turn call the initializers of the `modal_aero_amicphys`, `modal_aero_calcsize`, `modal_aero_coag`, `modal_aero_deposition`, `modal_aero_gasaerexch`, `modal_aero_newnuc`, `modal_aero_rename`, and `modal_aero_convproc` modules. These are presently supressed, but if data from these modules is indirectly required by `mo_gocart` for wet and dry deposition, these would need to be stubbed or the data set manually. Because these modules make heavy use of the physics buffer, it could be a considerable task creating equivalent stubs, likely on the order of a few weeks of effort. It's not clear yet if these are necessary so this represents a larger uncertainty on a hypothetical timeline for the GOCART bridge project.

# Additional Notes for Development
The following sections provide more specific notes and explanations that may be useful to developers seeking to complete or extend this bridge in the future.

## EAM Aerosol Models and GOCART-2G
**EAM** appears to provide several Aerosol models located in subdirectories of `${EAM_ROOT}/src/chemistry`, which provide information on the tracked chemical species, reactions, and supporting data. These models provide a version of the `chem_mods` module, which provides metadata, such as number of gas-phase species (`gas_pcnst`), used to size arrays among other things, and `mo_sim_dat.F90`, which provides data such as species names, invariant species, advective mass, *etc*. These models also provide implementations of various chemical processes in modules and subroutines not mentioned in the scope of the bridging process. Presumably, this means that the `mo_gocart` module offers either an alternative to these processes, or they are assumed to be accounted in some other **EAMxx** processes. It is similarly possible, that our collaborators have some misconceptions about how physical and chemical processes are handled in **EAMxx**, in which case, the procedures may also need to be bridged in order for the full atmospheric chemistry model to take effect.

The aerosol models provide names of the chemical species accounted for, which are marked either as *Explicit* or *Implicit*. The sum total of these adds to the total number of gas-phase species (`gas_pcnst`). While unconfirmed, this tagging may indicate that some species are advected, while others are computed locally according to some quasi-equilibrium assumption or prescribed as *Invariants*. If this is true, then the *Implicit* species may need to be registered as *Fields* rather than *tracers* in EAMxx.

The provided example case utilizing the implementation of **GOCART-2G** in **EAM** also uses the 5-mode aerosol model *pp_chemuci_linozv3_mam5_vbs*. This prescribes `gas_pcnst=73` gas-phase species consisting of 36 *Explicit* and 37 *Implicit* species. The present **EAMxx** bridge that is the `GOCART` `AtmosphereProcess` assumes this aerosol model, though many design decisions were made in hope of supporting alternate modal aerosol models as specified in the **EAM** source code. In some cases this may be as simple as supplying an alternative `chem_mods.F90` and `mo_sim_dat.F90` files in the `CMakeLists.txt`, though this has not been tested. This would in principle, open the possibility of users specifying the underlying aerosol model from some input YAML file prior to compiling **E3SM**.


## Aqueous Chemistry Bridge
The requested cloud aqueous chemistry component is provided in the `mo_setsox` module. The original (found in `${EAM_ROOT}/src/chemistry/aerosol/mo_setsox.F90`) has been refactored to remove usage of the `cam_history` and `phys_control`, which are used to write simulation data and are not essential to the function of the subprogram. Due to the change in logging infrastructure between **EAM** and **EAMxx**, a direct bridge of these modules may not be feasible. Creating stub for these modules as an interface for the new logging system may be possible, but is outside the scope of the present project.

Based on context and structure, the `invariants` dummy variable array in the `SETSOX` subroutine likely refers to species in the aerosol model marked as *Invariant*. This array is populated indirectly through the physics buffer so it is challenging to trace the path of this data, though it appears to have identical indexing as the tracer array, so it may be acceptable to pass the tracer as the dummy argument for `invariants` as well. Because the *pp_chemuci_linozv3_mam5_vbs* aerosol model does not appear to specify any invariants relevant to the `SETSOX` procedure, the stubbed version of `invariants` has been removed for clarity.



## mo_gocart module procedures
The following are a list of procedures in `mo_gocart`. Broadly speaking, these have been marked with `public` access rights, but are not necessarily used by other subprograms or at all in some cases. The list has been divided into **External Procedures**, which ***ARE*** called by subprograms external to `mo_gocart`, **Internal Procedures**, which are called only by other procedures in `mo_gocart`, and **Unused Procedures**, which ***ARE NOT*** presently called anywhere in **EAM**. Some of the internal procedures are only called by unused procedures and some either shadow or are perhaps copies of procedures defined in other modules.

### External Procedures:
- `gc_dvel_inti_fromlnd`: part of interface `gc_drydep_inti`
- `gc_drydep_fromlnd`: part of interface `gc_drydep`
- `gc_dvel_inti_table`: part of interface `gc_drydep_inti`
- `gc_drydep_table`: part of interface `gc_drydep`
- `interpdvel`: called in `mo_gocart` and in `mo_drydep`
- `intp2d`: called in `mo_gocart` and in `mo_drydep`
- `get_landuse_and_soilw_from_file`: called in `mo_gocart` and in `mo_drydep`
- `interp_map`: called in `mo_gocart` and in `mo_drydep`
- `gc_wetdep_init`: called in `mo_chemini`
- `gc_wetdep_inputs_set`: called in `mo_gas_phase_chemdr`
- `gc_wetdep`: called in `mo_gas_phase_chemdr`
- `aerosol_depvel_compute`: called in `modal_aero/aero_model.F90`
- `NIthermo`: called in `mo_gas_phase_chemdr`
### Internal Procedures:
- `gc_dvel_inti_fromlnd`: only called by `gc_dvel_xactive`
- `gc_drydep_xactive`: only called by `gc_drydp_fromlnd`
- `soilw_inti`: only called by `gc_dvel_inti_xactive`
- `chk_soilw`: called by `set_soilw`. Shares call signature with `chk_soilw` from `mo_drydep`
- `gc_has_drydep`: called by `gc_dvel_inti_table` and `gc_dvel_inti_xactive`
- `gc_clddiag`: called by `gc_wetdep_inputs_set`
- `flux_precnum_vs_flux_prec_mpln`: Shares call signature with function in `wetdep` module
- `faer_resusp_vs_fprec_evap_mpln`: Shares call signature with function in `wetdep` module
- `fprecn_resusp_vs_fprec_evap_mpln`: Shares call signature with function in `wetdep` module
- `flux_precnum_vs_flux_prec_mp`: Shares call signature with function in `wetdep` module
- `flux_precnum_vs_flux_prec_ln`: Shares call signature with function in `wetdep` module
- `faer_resusp_vs_fprec_evap_mp`: Shares call signature with function in `wetdep` module
- `faer_resusp_vs_fprec_evap_ln`: Shares call signature with function in `wetdep` module
- `fprecn_resusp_vs_fprec_evap_mp`: Shares call signature with function in `wetdep` module
- `fprecn_resusp_vs_fprec_evap_ln`: Shares call signature with function in `wetdep` module
- `RPMARES`
- `AWATER`
- `nh3_POLY4`
- `CUBIC`
- `ACTCOF`
- `HNO3_reaction_rate`
- `sktrs_hno3`
- `sktrs_sslt`


### Unused Procedures

- `gc_dvel_inti_xactive`
- `set_soilw`: This shares a call signature with `set_soilw` in `mo_drydep`, which is called in `mo_gas_phase_chemdr`.
- `gc_wetdep_inputs_unset`referenced in a `use` statement in `mo_gas_phase_chemdr` but not called
- `SSLT_reaction_rate`
- `apportion_reaction_rate`




## Refactor Notes:
- **Explicit-shape arrays in Fortran:** Many subroutines in **EAM** utilize explicitly-shaped local arrays (e.g. `real(r8) :: rainmr(pcols,pver)`) which often are based on the local grid subset size `pcols` and `pver`, which may be set as compile-time Fortran `parameters` in **EAM** or set at runtime. Prior to what I previously believed, even if a non-parameter integer is made available through a module (e.g. `ppgrid`), a local array utilizing this for size information is treated as a *automatic array*, which is automatically allocated and deallocated upon entering and leaving the subprogram. It's unclear if this is part of the language standard or if it is compiler-specific, but it is presently relied upon in **EAM**. So long as grid information is initialized, wrapped Fortran procedures will continue to use this design pattern.

  *Note:* the compiler should have discretion to allocate automatic arrays on the stack or on the heap. While not expected, some compilers may still stack-allocate larger arrays potentially leading to stack overflow issues. If this occurs, the procedures may need to be refactored as dynamic arrays (e.g. `allocate(rainmr(pcols,pver))` ... `deallocate(rainmr)`) or `pcols` and `pver` must somehow be supplied at compile-time as in **EAM**.

- In **EAM** `pcols` is the maximum number of columns per chunk, whereas `ncol` is the number of columns actually used in a given chunk. Thus far, there does not appear to be a compelling reason to assign different values to these variables, however if some of the wrapped **EAM** subprograms have MPI transactions that rely on these, it could cause problems. As a precaution, upon intialization, the largest `ncol` value is taken as `pcols`.




## May require additional investigation:
- **seq_drydep_mod module:** There are multiple potential sources of the `seq_drydep_mod` module required by `mo_gocart` located in:
    - `${E3SM_ROOT}/driver-mct/shr/seq_drydep_mod.F90`
    - `${E3SM_ROOT}/components/elm/src/utils/seq_drydep_mod_elm.F90`
    - `${E3SM_ROOT}/driver-moab/shr/seq_drydep_mod.F90`
    - `${E3SM_ROOT}/share/nuopc/seq_drydep_mod.F90`
  
  The `${E3SM_ROOT}/share/nuopc/seq_drydep_mod.F90` seems to be the most likely source as it is external to any specific **E3SM** or **EAM**/**EAMxx** component or process, however there are some deprecated procedures/data from the ESMF timing module that do not appear to be set in this version of **E3SM**. Therefore, a stub of this file is used in the GOCART process that removes these deprecated procedures. See `stubs/seq_drydep_mod.F90` for specifics.

- **shr_utils:** The `shr_utils` library referenced from the `CMakeLists.txt` file recreated the `${E3SM_ROOT}/share/util` dependencies. Some stubs of the shared utility files have been generated (and may be regenerated) from shell scripts which replace some of the template metaprogramming that would occur during a normal `cmake`. These utilities should be available in a more standardized form in the `csm_share` library and by proxy through the `scream_share` library. Each library that uses these utilities ***MUST*** link `csm_share` or `scream_share` directly before linking with the final `gocart` library. Ideally the `shr_utils` library and source stubs should not be used, but have been left as a backup in case necessary utilities are not compiled for some reason.

- The provided `atm_in` and `Filepath` files indicate the *rrtmg* radiation scheme is being used rather than the *rrtmgp* scheme. Therefore the `${EAM_ROOT}/src/physics/rrtmg/radconstants.F90` is used in the `eam_share` library rather than the version in the `${EAM_ROOT}/src/physics/rrtmgp` directory.

- At present, some preprocessor directives such as defining `MODAL_AERO_5MODE` are hard-coded in `CMakeLists.txt` to match settings seen in the provided example case. Ideally this should probably be settable by the user without altering the source code.



