#!/bin/bash

#==================================================================================================
# For posterity:
#
# There is no 'shr_assert_mod.F90' provided in the E3SM 'share' directory, rather a template file
# 'shr_assert_mod.F90.in' exists that CMake uses for metaprogramming. This shell script attempts
# to replicate the metaprogramming and generate a *.F90 file that can be used for dependencies 
# without building and linking the E3SM 'share' libraries (which is probably possible, but I don't
# know how to do)
#
# - Daniel Fust (Aug 2025)
#==================================================================================================

FILE="shr_assert_mod.F90"
TYPES=("real" "double" "int" "long")
VTYPES=("real(r4)" "real(r8)" "integer(i4)" "integer(i8)")
DIMS=(0 1 2 3 4 5 6 7)
DIMSTR=("" 
        "(:)"
        "(:,:)"
        "(:,:,:)"
        "(:,:,:,:)"
        "(:,:,:,:,:)"
        "(:,:,:,:,:,:)"
        "(:,:,:,:,:,:,:)")

MASKDIM=(""
         "size(mask,1)"
         "size(mask,1), size(mask,2)" 
         "size(mask,1), size(mask,2), size(mask,3)"
         "size(mask,1), size(mask,2), size(mask,3), size(mask,4)"
         "size(mask,1), size(mask,2), size(mask,3), size(mask,4), size(mask,5)"
         "size(mask,1), size(mask,2), size(mask,3), size(mask,4), size(mask,5), size(mask,6)"
         "size(mask,1), size(mask,2), size(mask,3), size(mask,4), size(mask,5), size(mask,6), size(mask,7)")

# It's not clear precisely what {ITYPE} should be, though it appears to
# specify 4-byte vs 8-byte.
#ITYPE="shr_kind_i4"
ITYPE="shr_kind_i8"


#-- Make or erase file --
> $FILE

cat << EOF>> $FILE
#include "preprocessor_defines.h"

module shr_assert_mod

! Assert subroutines for common debugging operations.

use shr_kind_mod, only: &
     r4 => shr_kind_r4, &
     r8 => shr_kind_r8, &
     i4 => shr_kind_i4, &
     i8 => shr_kind_i8

use shr_sys_mod, only: &
     shr_sys_abort

use shr_log_mod, only: &
     shr_log_Unit

use shr_infnan_mod, only: shr_infnan_isnan

use shr_strconvert_mod, only: toString

implicit none
private
save

! Assert that a logical is true.
public :: shr_assert
public :: shr_assert_all
public :: shr_assert_any

! Assert that a numerical value satisfies certain constraints.
public :: shr_assert_in_domain

interface shr_assert_all
   module procedure shr_assert
   ! DIMS 1,2,3,4,5,6,7
EOF
for (( j=1; j<${#DIMS[@]}; j++ )); do   
   echo "   module procedure shr_assert_all_${DIMS[j]}d" >> $FILE
done
cat << EOF>> $FILE 
end interface

interface shr_assert_any
   module procedure shr_assert
   ! DIMS 1,2,3,4,5,6,7
EOF
for (( j=1; j<${#DIMS[@]}; j++ )); do   
   echo "   module procedure shr_assert_any_${DIMS[j]}d" >> $FILE
done
cat << EOF>> $FILE   
end interface

interface shr_assert_in_domain
   ! TYPE double,real,int,long
   ! DIMS 0,1,2,3,4,5,6,7
EOF
for (( i=0; i<${#TYPES[@]}; i++ )); do
    for (( j=0; j<${#DIMS[@]}; j++ )); do
   echo "   module procedure shr_assert_in_domain_${DIMS[j]}d_${TYPES[i]}" >> $FILE
    done
done
cat << EOF>> $FILE
end interface

! Private utilities.

interface print_bad_loc
   ! TYPE double,real,int,long
   ! DIMS 0,1,2,3,4,5,6,7
EOF
for (( i=0; i<${#TYPES[@]}; i++ )); do
    for (( j=0; j<${#DIMS[@]}; j++ )); do   
   echo "   module procedure print_bad_loc_${DIMS[j]}d_${TYPES[i]}" >> $FILE
    done
done
cat << EOF>> $FILE   
end interface

interface find_first_loc
   ! DIMS 0,1,2,3,4,5,6,7
EOF
for (( j=0; j<${#DIMS[@]}; j++ )); do   
   echo "   module procedure find_first_loc_${DIMS[j]}d" >> $FILE
done
cat << EOF>> $FILE
end interface

interface within_tolerance
   ! TYPE double,real,int,long
EOF
for (( i=0; i<${#TYPES[@]}; i++ )); do   
   echo "   module procedure within_tolerance_${TYPES[i]}" >> $FILE
done 
cat << EOF>> $FILE   
end interface

contains

subroutine shr_assert(var, msg, file, line)

  ! Logical being asserted
  logical, intent(in) :: var
  ! Optional error message if assert fails
  character(len=*), intent(in), optional :: msg
  ! Optional file and line of the caller, written out if given
  ! (line is ignored if file is absent)
  character(len=*), intent(in), optional :: file
  integer         , intent(in), optional :: line

  character(len=:), allocatable :: full_msg

  if (.not. var) then
     full_msg = 'ERROR'
     if (present(file)) then
        full_msg = full_msg // ' in ' // trim(file)
        if (present(line)) then
           full_msg = full_msg // ' at line ' // toString(line)
        end if
     end if
     if (present(msg)) then
        full_msg = full_msg // ': ' // msg
     end if
     call shr_sys_abort(full_msg)
  end if

end subroutine shr_assert

! DIMS 1,2,3,4,5,6,7
EOF
for (( j=1; j<${#DIMS[@]}; j++ )); do
cat << EOF>> $FILE
subroutine shr_assert_all_${DIMS[j]}d(var, msg, file, line)

  ! Logical being asserted
  logical, intent(in) :: var${DIMSTR[j]}
  ! Optional error message if assert fails
  character(len=*), intent(in), optional :: msg
  ! Optional file and line of the caller, written out if given
  ! (line is ignored if file is absent)
  character(len=*), intent(in), optional :: file
  integer         , intent(in), optional :: line

  call shr_assert(all(var), msg=msg, file=file, line=line)

end subroutine shr_assert_all_${DIMS[j]}d

EOF
done
cat << EOF>> $FILE

! DIMS 1,2,3,4,5,6,7
EOF
for (( j=1; j<${#DIMS[@]}; j++ )); do
cat << EOF>> $FILE
subroutine shr_assert_any_${DIMS[j]}d(var, msg, file, line)

  ! Logical being asserted
  logical, intent(in) :: var${DIMSTR[j]}
  ! Optional error message if assert fails
  character(len=*), intent(in), optional :: msg
  ! Optional file and line of the caller, written out if given
  ! (line is ignored if file is absent)
  character(len=*), intent(in), optional :: file
  integer         , intent(in), optional :: line

  call shr_assert(any(var), msg=msg, file=file, line=line)

end subroutine shr_assert_any_${DIMS[j]}d

EOF
done
cat << EOF>> $FILE
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

! TYPE double,real,int,long
! DIMS 0,1,2,3,4,5,6,7
EOF
for (( i=0; i<${#TYPES[@]}; i++ )); do
    for (( j=0; j<${#DIMS[@]}; j++ )); do
cat << EOF>> $FILE
subroutine shr_assert_in_domain_${DIMS[j]}d_${TYPES[i]}(var, varname, msg, &
     is_nan, lt, gt, le, ge, eq, ne, abs_tol)

!-----------------------------
! BEGIN defining local macros
!-----------------------------

! Flag for floating point types.

EOF
# macro replaced with sketch bash logic (DFust)
#if (${ITYPE} == TYPEREAL) || (${ITYPE} == TYPEDOUBLE)
#define TYPEFP
#else
#undef TYPEFP
#endif
cat << EOF >> $FILE

! "Generalized" macro functions allow transformational intrinsic functions
! to handle both scalars and arrays.

#if (${DIMS[j]} != 0)
! When given an array, use the intrinsics.
#define GEN_SIZE(x) size(x)
#define GEN_ALL(x) all(x)
#else

! Scalar extensions:
!   GEN_SIZE always returns 1 for a scalar.
!   GEN_ALL (logical reduction) is a no-op for a scalar.
!   GEN_[MAX,MIN]LOC should return a 1D, size 0 (empty), integer array.
#define GEN_SIZE(x) 1
#define GEN_ALL(x) x

#endif

!-----------------------------
! END macro section
!-----------------------------

  ! Variable being checked.
  ${VTYPES[i]}, intent(in) :: var${DIMSTR[j]}
  ! Variable name to be used in error messages.
  character(len=*), intent(in), optional :: varname
  ! Optional error message if assert fails.
  character(len=*), intent(in), optional :: msg
  ! Assert that the variable is not (or is) NaN.
  logical, intent(in), optional :: is_nan
  ! Limits for (in)equalities.
  ${VTYPES[i]}, intent(in), optional :: lt
  ${VTYPES[i]}, intent(in), optional :: gt
  ${VTYPES[i]}, intent(in), optional :: le
  ${VTYPES[i]}, intent(in), optional :: ge
  ${VTYPES[i]}, intent(in), optional :: eq
  ${VTYPES[i]}, intent(in), optional :: ne
  ${VTYPES[i]}, intent(in), optional :: abs_tol

  ! Note that the following array is size 0 for scalars.
  integer :: loc_vec(${DIMS[j]})

  logical :: is_nan_passed
  logical :: lt_passed
  logical :: gt_passed
  logical :: le_passed
  logical :: ge_passed
  logical :: eq_passed
  logical :: ne_passed

  ${VTYPES[i]} :: abs_tol_loc

  ! Handling of abs_tol makes a couple of fairly safe assumptions.
  !  1. It is not the most negative integer.
  !  2. It is finite (not a floating point infinity or NaN).
  if (present(abs_tol)) then
     abs_tol_loc = abs(abs_tol)
  else
     abs_tol_loc = 0_i4
  end if

  is_nan_passed = .true.
  lt_passed = .true.
  gt_passed = .true.
  le_passed = .true.
  ge_passed = .true.
  eq_passed = .true.
  ne_passed = .true.

  ! Do one pass just to find out if we can return with no problem.


EOF
#ifdef TYPEFP
# compiler macro replaced with sketchy bash logic (DFust)
if [[ ${TYPES[i]} == "real(r4)" ]] || [[ ${TYPES[i]} == "real(r8)" ]]; then
cat << EOF >> $FILE
  ! Only floating-point values can actually be Inf/NaN.
  if (present(is_nan)) &
     is_nan_passed = GEN_ALL(shr_infnan_isnan(var) .eqv. is_nan)
EOF
#else
else
cat << EOF >> $FILE
  if (present(is_nan)) &
     is_nan_passed = .not. is_nan .or. GEN_SIZE(var) == 0
EOF
#endif
fi
cat << EOF >> $FILE

  if (present(lt)) &
     lt_passed = GEN_ALL(var < lt)

  if (present(gt)) &
     gt_passed = GEN_ALL(var > gt)

  if (present(le)) &
     le_passed = GEN_ALL(var <= le)

  if (present(ge)) &
     ge_passed = GEN_ALL(var >= ge)

  if (present(eq)) then
     eq_passed = GEN_ALL(within_tolerance(eq, var, abs_tol_loc))
  end if

  if (present(ne)) then
     ne_passed = GEN_ALL(.not. within_tolerance(ne, var, abs_tol_loc))
  end if

  if ( is_nan_passed .and. &
       lt_passed .and. &
       gt_passed .and. &
       le_passed .and. &
       ge_passed .and. &
       eq_passed .and. &
       ne_passed) &
       return

  ! If we got here, assert will fail, so find out where so that we
  ! can try to print something useful.

  if (.not. is_nan_passed) then
EOF
#ifdef TYPEFP
# macro replaced with sketchy bash logic (DFust)
if [[ ${TYPES[i]} == "real(r4)" ]] || [[ ${TYPES[i]} == "real(r8)" ]]; then
cat << EOF >> $FILE
     loc_vec = find_first_loc(shr_infnan_isnan(var) .neqv. is_nan)
     call print_bad_loc(var, loc_vec, varname)
     if (is_nan) then
        write(shr_log_Unit,*) "Expected value to be NaN."
     else
        write(shr_log_Unit,*) "Expected value to be a number."
     end if
EOF
#else
else
cat << EOF >> $FILE
     loc_vec = spread(1,1,${DIMS[j]})
     call print_bad_loc(var, loc_vec, varname)
     if (is_nan) then
        write(shr_log_Unit,*) &
             "Asserted NaN, but the variable is not floating-point!"
     end if
EOF
#endif
fi
cat << EOF >> $FILE
  end if

  if (.not. lt_passed) then
     loc_vec = find_first_loc(var >= lt)
     call print_bad_loc(var, loc_vec, varname)
     write(shr_log_Unit,*) "Expected value to be less than ",lt
  end if

  if (.not. gt_passed) then
     loc_vec = find_first_loc(var <= gt)
     call print_bad_loc(var, loc_vec, varname)
     write(shr_log_Unit,*) "Expected value to be greater than ",gt
  end if

  if (.not. le_passed) then
     loc_vec = find_first_loc(var > le)
     call print_bad_loc(var, loc_vec, varname)
     write(shr_log_Unit,*) "Expected value to be less than or &
          &equal to ",le
  end if

  if (.not. ge_passed) then
     loc_vec = find_first_loc(var < ge)
     call print_bad_loc(var, loc_vec, varname)
     write(shr_log_Unit,*) "Expected value to be greater than or &
          &equal to ",ge
  end if

  if (.not. eq_passed) then
     loc_vec = find_first_loc(.not. within_tolerance(eq, var, abs_tol_loc))
     call print_bad_loc(var, loc_vec, varname)
     write(shr_log_Unit,*) "Expected value to be equal to ",eq
     if (abs_tol_loc > 0) &
          write(shr_log_Unit,*) "Asserted with tolerance ", abs_tol_loc
  end if

  if (.not. ne_passed) then
     loc_vec = find_first_loc(within_tolerance(ne, var, abs_tol_loc))
     call print_bad_loc(var, loc_vec, varname)
     write(shr_log_Unit,*) "Expected value to never be equal to ",ne
     if (abs_tol_loc > 0) &
          write(shr_log_Unit,*) "Asserted with tolerance ", abs_tol_loc
  end if

  call shr_sys_abort(msg)

EOF
# macros replaced with sketch bash logic (DFust)
#! Undefine local macros.
#undef TYPEFP
#undef GEN_SIZE
#undef GEN_ALL
cat << EOF >> $FILE

end subroutine shr_assert_in_domain_${DIMS[j]}d_${TYPES[i]}
EOF
    done 
done
cat << EOF>> $FILE

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

! TYPE double,real,int,long
! DIMS 0,1,2,3,4,5,6,7
EOF
for (( i=0; i<${#TYPES[@]}; i++ )); do
    for (( j=0; j<${#DIMS[@]}; j++ )); do
cat << EOF>> $FILE

subroutine print_bad_loc_${DIMS[j]}d_${TYPES[i]}(var, loc_vec, varname)
  ! Print information about a bad location in an variable.
  ! For scalars, just print value.

  ${VTYPES[i]}, intent(in) :: var${DIMSTR[j]}
  integer, intent(in) :: loc_vec(${DIMS[j]})

  character(len=*), intent(in), optional :: varname

  character(len=:), allocatable :: varname_to_write

  if (present(varname)) then
     allocate(varname_to_write, source=varname)
  else
     allocate(varname_to_write, source="input variable")
  end if

  write(shr_log_Unit,*) &
       "ERROR: shr_assert_in_domain: ",trim(varname_to_write), &
       " has invalid value ", &
#if (${DIMS[j]} != 0)
EOF
        echo -n "       var(loc_vec(1)" >> $FILE

#       {REPEAT:loc_vec(#)}), &
        for (( k=2; k<=${DIMS[j]}; k++ )); do
            echo -n ",loc_vec($k)" >> $FILE
        done
        echo -n ")," >> $FILE
cat << EOF >> $FILE
       " at location: ",loc_vec
#else
       var

  ! Kill compiler spam for unused loc_vec.
  if (.false.) write(*,*) loc_vec
#endif

end subroutine print_bad_loc_${DIMS[j]}d_${TYPES[i]}

EOF
    done
done
cat << EOF>> $FILE 

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

! DIMS 0,1,2,3,4,5,6,7
EOF
for (( j=0; j<${#DIMS[@]}; j++ )); do
cat << EOF>> $FILE
pure function find_first_loc_${DIMS[j]}d(mask) result (loc_vec)
  ! Inefficient but simple subroutine for finding the location of
  ! the first .true. value in an array.
  ! If no true values, returns first value.

  logical, intent(in) :: mask${DIMSTR[j]}
  integer :: loc_vec(${DIMS[j]})

#if (${DIMS[j]} != 0)
  integer :: flags(${MASKDIM[j]})

  where (mask)
     flags = 1
  elsewhere
     flags = 0
  end where

  loc_vec = maxloc(flags)
#else

! Remove compiler warnings (statement will be optimized out).

#if (! defined CPRPGI && ! defined CPRCRAY)
  if (.false. .and. mask) loc_vec = loc_vec
#endif

#endif

end function find_first_loc_${DIMS[j]}d
EOF
done
cat << EOF>> $FILE


! TYPE double,real,int,long
EOF
for (( i=0; i<${#TYPES[@]}; i++ )); do
cat << EOF>> $FILE
elemental function within_tolerance_${TYPES[i]}(expected, actual, tolerance) &
     result(is_in_tol)
  ! Precondition: tolerance must be >= 0.
  ${VTYPES[i]}, intent(in) :: expected
  ${VTYPES[i]}, intent(in) :: actual
  ${VTYPES[i]}, intent(in) :: tolerance
  logical :: is_in_tol

  ! The following conditionals are to ensure that we don't overflow.

  ! This takes care of two identical infinities.
  if (actual == expected) then
     is_in_tol = .true.
  else if (actual > expected) then
     if (expected >= 0) then
        is_in_tol = (actual - expected) <= tolerance
     else
        is_in_tol = actual <= (expected + tolerance)
     end if
  else
     if (expected < 0) then
        is_in_tol = (expected - actual) <= tolerance
     else
        is_in_tol = actual >= (expected - tolerance)
     end if
  end if

end function within_tolerance_${TYPES[i]}
EOF
done
cat << EOF>> $FILE

end module shr_assert_mod
EOF