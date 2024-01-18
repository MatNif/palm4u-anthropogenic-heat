#!/usr/bin/env bash

# kpp4palm - script for creating gasphase module

#------------------------------------------------------------------------------#
# This file is part of the PALM model system.
#
# PALM is free software: you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
#
# PALM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PALM. If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2017-2023  Klaus Ketelsen and MPI-CH (April 2007)
# Copyright 2017-2023  Karlsruhe Institute of Technology
# Copyright 2017-2023  Leibniz Universitaet Hannover
#------------------------------------------------------------------------------#
# Nov. 2016: Initial Version of KPP chemistry convertor adapted for PALM 
# by Klaus Ketelsen
#
# This code is a modified version of KP4 (Jöckel, P., Kerkweg, A., Pozzer, A., 
# Sander, R., Tost, H., Riede, H., Baumgaertner, A., Gromov, S., and Kern, B.,
# 2010: Development cycle 2 of the Modular Earth Submodel System (MESSy2), 
# Geosci. Model Dev., 3, 717-752, https://doi.org/10.5194/gmd-3-717-2010).
# KP4 is part of the Modular Earth Submodel System (MESSy), which is is
# available under the  GNU General Public License (GPL).
#
#------------------------------------------------------------------------------#
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
   DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
   SOURCE="$(readlink "$SOURCE")"
   [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
SCRIPT_LOCATION="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

project_root_dir=$(readlink -f "${SCRIPT_LOCATION}/..")
project_src_dir=$(readlink -f "${project_root_dir}/src/")
project_kpp_dir=$(readlink -f "${project_root_dir}/kpp/")

program_name=$(basename -s ".sh" ${BASH_SOURCE[0]})

cd ${project_root_dir}

set -eu

########################### User SetUp ####################################

export KPP_HOME=${project_kpp_dir}
export KPP=${project_kpp_dir}/bin/kpp

########################## End User Setup ################################

WORK=tmp_${program_name}

# Default

MECH=phstatp  # use photo-stationary mechanism as default
OUTDIR=${project_root_dir}/build
OUTDIR_SET="NO"
OUTFILE=chem_gasphase_mod
PREFIX=chem_gasphase_mod
MODE="scalar"
VLEN=1
KEEP="NO"
UPDT="NO"
DE_INDEX=2
DE_INDEX_FAST="YES"

export KPP_SOLVER=Rosenbrock


show_usage() {
   echo "Usage: $0 [-h] -m <mechanism_name> -o <path> [-i <method_number>] [-v] [-l <vector_length>] [-k] [-u] [-f] [-p <name>] [-s <name>]"
}

show_help() {
   show_usage
   echo "      -h                     show this help message"
   echo "      -m <mechanism_name>    set chemical mechanism (phstat, phstatp, smog, simple,cbm4, etc, see directory mechanisms)"
   echo "      -o <path>              set output directory for the generated code"
   echo "      -i <method_number>     set deindexing method (0,1,2,3; 0=no deindexing)"
   echo "      -v                     switch to enable vector mode (not working completely yet)"
   echo "      -l <vector_length>     set vector length (not working completely yet)"
   echo "      -k                     switch on to keep working directory"
   echo "      -u                     switch to update the f90 code in the def_MECH directory (existing files will be overwritten)"
   echo "      -f                     switch to enable fast deindexing"
   echo "      -p <name>              set name prefix for output file (default: chem_gasphase_mod) should not be changed"
   echo "      -s <name>              set name of solver (only Rosebrock solvers work for vector mode)"
}

# get Command line option

while  getopts :hm:o:i:fkup:s:vl:w:  c     # get options
do case $c in
      m)   MECH=$OPTARG;;            # mechanism

      o)   OUTDIR=$(readlink -f ${OPTARG})  # Output directory of Generated Code
           OUTDIR_SET="YES";;

      i)   DE_INDEX=$OPTARG;;        # if set, deindexing

      f)   DE_INDEX_FAST="YES";;     # if set, fast deindexing

      k)   KEEP="YES";;              # keep Working directory

      p)   PREFIX=$OPTARG;;          # Name Prefix (chem_gasphase_mod, do not change)

      s)   KPP_SOLVER=$OPTARG;;      # Chosen solver (only Rosebrock solvers work for vector mode)

      u)   UPDT="YES";;              # update mechanisms/def_$MECH/chem_gasphase_mod.f90

      v)   MODE="vector";;           # Set to vector Mode

      l)   VLEN=$OPTARG;;            # Set vector length

      w)   WORK=$OPTARG;;            # Working directory

      h)
         show_help
         exit 0
         ;;
      *)
         show_usage
         exit 1
         ;;

   esac
done

if [ "${OUTDIR_SET}" != "YES" ]; then
   echo "ERROR: output directory for the generated code was not set. Please use option \"-o <path>\""
   exit 1
fi

echo MECHANISM = $MECH
echo DE_INDEX = $DE_INDEX
echo KEEP = $KEEP
echo UPDT = $UPDT
echo MODE = $MODE
echo VLEN = $VLEN

DEF_PREFIX=${PREFIX}.kpp
DEFDIR=${project_root_dir}/mechanisms/def_$MECH
echo DEFDIR = $DEFDIR

# Create or clean working directory

cd ${project_root_dir}
WORK=$(readlink -f "$WORK")
mkdir -p $WORK
rm -rf $WORK/*
cd $WORK

# kpp dependend, may be changed

KPP_FILE_LIST="Initialize Integrator LinearAlgebra Jacobian Function Rates Util"
if [[ $MODE = "vector" ]]
 then
KPP_FILE_LIST="$KPP_FILE_LIST kp4_compress_subroutines"
fi

KPP_SUBROUTINE_LIST="Initialize"
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST INTEGRATE"
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST Fun"
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST KppSolve"
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST KppDecomp"
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST Jac_SP"
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST Update_RCONST"
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST initialize_kpp_ctrl"
if [[ $MODE != "vector" ]]
then
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST error_output"
fi

# parse all function names from UserRateLaws.f90 and add them to KPP_SUBROUTINE_LIST

KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST $(sed -Enr 's/^[ \t]*[^\!].+FUNCTION *(\w*) *\(.+\).*$/\1/igp' ${DEFDIR}/UserRateLaws.f90)"

# if [[ $MODE = "vector" && $KPP_SOLVER = "ROS2" ]]
# then
#   cp ${project_root_dir}/templates/${KPP_SOLVER}_vec.f90 ${KPP_SOLVER}.f90    # get vector Solver
# else
# #  KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST FunTemplate JacTemplate Update_SUN "
#   KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
#   if [[ $MODE = "vector" ]]
#   then
#     cp ${project_root_dir}/templates/${KPP_SOLVER}_vec.f90 ${KPP_SOLVER}.f90  # get vector Solver
#   else
#     KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST Rosenbrock  FunTemplate JacTemplate Update_SUN"
#   fi
# fi
 if [[ $MODE = "vector" ]]
 then
   # get vector Solver 
   cp ${project_root_dir}/templates/${KPP_SOLVER}_vec.f90 ${KPP_SOLVER}.f90
   cp ${project_root_dir}/templates/kp4_compress_header ${PREFIX}_kp4_compress_header.f90
   cp ${project_root_dir}/templates/kp4_compress_subroutines ${PREFIX}_kp4_compress_subroutines.f90
fi

# Interface ignore list
KPP_INTERFACE_IGNORE=" "

echo " "
echo KPP_SOLVER $KPP_SOLVER
echo " "

case $KPP_SOLVER in
    ROS2) ;;

    Rosenbrock)   
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WAXPY"
      if [[ $MODE != "vector" ]]
      then
         KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WSCAL Rosenbrock  FunTemplate JacTemplate"
        KPP_INTERFACE_IGNORE="WAXPY"

      else 
        KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST FunTemplate JacTemplate"
        KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST kco_initialize kco_compress kco_finalize"
      fi;;

    rosenbrock_mz)
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST Rosenbrock  FunTemplate JacTemplate Update_SUN";;

    rosenbrock)
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST Rosenbrock  FunTemplate JacTemplate";;

    kpp_lsode)
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST KppLsode DLSODE JAC_CHEM FUN_CHEM"
      KPP_INTERFACE_IGNORE="$KPP_INTERFACE_IGNORE JAC_CHEM KppDecomp KppSolve";;

    kpp_radau5)
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY FUN_CHEM JAC_CHEM SET2ZERO"
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST RADAU5 Update_SUN"
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST KppSolveCmplx KppDecompCmplx";;

    kpp_sdirk)
       KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
       KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST SDIRK JAC_CHEM SET2ZERO FUN_CHEM"
       KPP_INTERFACE_IGNORE="$KPP_INTERFACE_IGNORE Set2zero SET2ZERO FUN_CHEM";;

    kpp_seulex)
       KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
       KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST ATMSEULEX"
       KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST SEULEX_ErrorMsg SEULEX_Integrator FUN_CHEM JAC_CHEM SEUL"
       KPP_INTERFACE_IGNORE="$KPP_INTERFACE_IGNORE SEULEX_Integrator SDIRK FUN_CHEM SEUL";;

   \?)  print "SORRY ONLY ROSENBROCK METHODS WORK AT THE MOMENT:" $KPP_SOLVER
        exit 1;;
esac
#mz-ak-20070509+

KPP_INCLUDE_LIST="Parameters Global JacobianSP Monitor"
if [[ $MODE = "vector" ]]
 then
KPP_INCLUDE_LIST="$KPP_INCLUDE_LIST kp4_compress_header"
fi

#Get definition Files

cp $DEFDIR/*.eqn         .
cp $DEFDIR/*.spc         .
cp $DEFDIR/${PREFIX}.kpp     .

# make sure kpp is using the correct UserRateLaws.f90 file

cp -p ${DEFDIR}/UserRateLaws.f90 ${project_root_dir}/kpp/util/

# Global variable are defined here 
# This has the advantage that it is not necessary to include these variables in all .kpp definition files

cat  >> ${PREFIX}.kpp  <<  EOF
#INLINE F90_GLOBAL
! QVAP - Water vapor
  REAL(kind=dp) :: QVAP
! FAKT - Conversion factor
  REAL(kind=dp) :: FAKT

! CS_MECH for check of mechanism name with namelist
  CHARACTER(LEN=30) :: CS_MECH
#ENDINLINE
EOF

# Store mechanism name in file mech_list
cat  >> mech_list  <<  EOF
!   Mechanism: $MECH
!
EOF

# Store mechanism name for cs_mech
cat  >> set_cm  <<  EOF

! Set cs_mech for check with mechanism name from namelist
    cs_mech = '$MECH'
EOF

# Run kpp

$KPP $DEF_PREFIX

# Get templates for C++ program

cp ${project_root_dir}/templates/module_header* .           # Use fixed Module_header
cp ${project_root_dir}/templates/initialize_kpp_ctrl_template.f90 .  # CTRL kpp time stepping

# file with subroutine list for c++ program create_kpp_module

for i in $KPP_FILE_LIST
do
  echo ${PREFIX}_${i} >> file_list
done
echo initialize_kpp_ctrl_template >> file_list

# file with subroutine list for c++ program create_kpp_module

for i in $KPP_SUBROUTINE_LIST
do
  echo $i >> subroutine_list
done

# file with include list for c++ program create_kpp_module

for i in $KPP_INCLUDE_LIST
do
  echo ${PREFIX}_${i} >> include_list
done

touch interface_ignore_list
for i in $KPP_INTERFACE_IGNORE
do
  echo $i >> interface_ignore_list
done

echo start ${program_name}.exe with arguments
echo $PREFIX $MODE $VLEN $DE_INDEX $DE_INDEX_FAST

${project_root_dir}/bin/${program_name}.exe $PREFIX $MODE $VLEN $DE_INDEX $DE_INDEX_FAST

# Add dummy statements in order to prevent warnings due to unused variables
#
sed -i -e '/cfactor =/a !  ' kk_kpp.f90
sed -i -e '/cfactor =/a BLANKS  IF ( lu_crow(1) == 1  .OR.  lu_icol(1) == 1  .OR.  lu_irow(1) == 1 )  CONTINUE ' kk_kpp.f90
sed -i -e '/cfactor =/a BLANKS  IF ( time >= -1.0_dp )  CONTINUE' kk_kpp.f90
sed -i -e '/cfactor =/a ! Following lines are just to avoid compiler message about unused variables' kk_kpp.f90
sed -i -e '/cfactor =/a !  ' kk_kpp.f90

if [[ $MODE = "vector" ]]
then
sed -i -e '/! Computation of equation rates/i ! The following lines are just to avoid compiler message about unused variables' kk_kpp.f90
sed -i -e '/! Computation of equation rates/i ! (some of the are only required if there only passive tracers)' kk_kpp.f90
sed -i -e '/! Computation of equation rates/i BLANKS  IF ( f(vl,nfix) >= 0.0_dp )  CONTINUE' kk_kpp.f90
sed -i -e '/! Computation of equation rates/i BLANKS  IF ( v(vl,nvar) >= 0.0_dp )  CONTINUE' kk_kpp.f90
sed -i -e '/! Computation of equation rates/i BLANKS  IF ( rct(vl,nreact) >= 0.0_dp )  CONTINUE' kk_kpp.f90
if ! grep -q ' a(1:' kk_kpp.f90 ; then
sed -i -e '/! Computation of equation rates/i BLANKS  IF ( a(vl,nreact) >= 0.0_dp )  CONTINUE' kk_kpp.f90
fi
sed -i -e '/! Computation of equation rates/i !  ' kk_kpp.f90
sed -i -e '/!   i = 1/i ! Following line is just to avoid compiler message about unused variables' kk_kpp.f90
sed -i -e '/!   i = 1/i ! if there are only passive tracers' kk_kpp.f90
sed -i -e '/!   i = 1/i BLANKS   IF ( jvs(vl,1) >= 0.0_dp )  CONTINUE' kk_kpp.f90
sed -i -e '/!   i = 1/i !  ' kk_kpp.f90

else

sed -i -e '/! Computation of equation rates/i ! The following lines are just to avoid compiler message about unused variables' kk_kpp.f90
sed -i -e '/! Computation of equation rates/i ! (some of the are only required if there only passive tracers)' kk_kpp.f90
sed -i -e '/! Computation of equation rates/i BLANKS  IF ( f(nfix) >= 0.0_dp )  CONTINUE' kk_kpp.f90
sed -i -e '/! Computation of equation rates/i BLANKS  IF ( v(nvar) >= 0.0_dp )  CONTINUE' kk_kpp.f90
sed -i -e '/! Computation of equation rates/i BLANKS  IF ( rct(nreact) >= 0.0_dp )  CONTINUE' kk_kpp.f90
if ! grep -q ' a(1)' kk_kpp.f90 ; then
sed -i -e '/! Computation of equation rates/i BLANKS  IF ( a(nreact) >= 0.0_dp )  CONTINUE' kk_kpp.f90
fi
sed -i -e '/! Computation of equation rates/i !  ' kk_kpp.f90
sed -i -e '/!   i = 1/i ! Following line is just to avoid compiler message about unused variables' kk_kpp.f90
sed -i -e '/!   i = 1/i ! if there are only passive tracers' kk_kpp.f90
sed -i -e '/!   i = 1/i BLANKS   IF ( jvs(1) >= 0.0_dp )  CONTINUE' kk_kpp.f90
sed -i -e '/!   i = 1/i !  ' kk_kpp.f90

sed -i -e '/IF ( alpha .eq. zero ) RETURN/i !  ' kk_kpp.f90
sed -i -e '/IF ( alpha .eq. zero ) RETURN/i ! Following line is just to avoid compiler message about unused variables' kk_kpp.f90
sed -i -e '/IF ( alpha .eq. zero ) RETURN/i BLANKS     IF ( incx == 0  .OR.  incy == 0 )  CONTINUE' kk_kpp.f90

fi

if [[ $MODE = "vector" ]]
then
sed -i -e '/REAL(kind=dp) :: b/a BLANKS  IF ( v(vl,nvar) >= 0.0_dp )  CONTINUE' kk_kpp.f90
sed -i -e '/REAL(kind=dp) :: b/a BLANKS  IF ( f(vl,nfix) >= 0.0_dp )  CONTINUE' kk_kpp.f90
sed -i -e '/REAL(kind=dp) :: b/a ! The following lines are just to avoid compiler message about unused variables' kk_kpp.f90
sed -i -e '/REAL(kind=dp):: b/a !' kk_kpp.f90
else
sed -i -e '/REAL(kind=dp):: b/a BLANKS  IF ( v(nvar) >= 0.0_dp )  CONTINUE' kk_kpp.f90
sed -i -e '/REAL(kind=dp):: b/a BLANKS  IF ( f(nfix) >= 0.0_dp )  CONTINUE' kk_kpp.f90
sed -i -e '/REAL(kind=dp):: b/a ! The following lines are just to avoid compiler message about unused variables' kk_kpp.f90
sed -i -e '/REAL(kind=dp):: b/a !' kk_kpp.f90
fi

sed -i -e '/one=1.0_dp/a BLANKS  IF ( incx == 0 )  CONTINUE' kk_kpp.f90
sed -i -e '/one=1.0_dp/a ! Following line is just to avoid compiler message about unused variables' kk_kpp.f90
sed -i -e '/one=1.0_dp/a !  ' kk_kpp.f90

sed -i -e '/IF ( sum(alpha(1:vl)) .eq. zero ) RETURN/i !  ' kk_kpp.f90
sed -i -e '/IF ( sum(alpha(1:vl)) .eq. zero ) RETURN/i ! The following line is just to avoid compiler message about unused variables' kk_kpp.f90
sed -i -e '/IF ( sum(alpha(1:vl)) .eq. zero ) RETURN/i BLANKS     IF ( incx == 0  .OR.  incy == 0 )  CONTINUE' kk_kpp.f90

sed -i -e '/INTENT(INOUT):: b(n)/a BLANKS IF ( pivot(1) == 0 )  CONTINUE' kk_kpp.f90
sed -i -e '/INTENT(INOUT):: b(n)/a ! Following line is just to avoid compiler message about unused variables' kk_kpp.f90
sed -i -e '/INTENT(INOUT):: b(n)/a !  ' kk_kpp.f90

sed -i -e '/singular = .FALSE./a BLANKS  IF ( direction==0 )  CONTINUE' kk_kpp.f90
sed -i -e '/singular = .FALSE./a ! The following line is just to avoid compiler message about unused variables' kk_kpp.f90
sed -i -e '/singular = .FALSE./a !  ' kk_kpp.f90


sed -i -e '1,$s/BLANKS /  /  ' kk_kpp.f90

sed -i -e '/! Time/ d ' kk_kpp.f90
sed -i -e '/! Working directory/ d ' kk_kpp.f90

mkdir -p ${OUTDIR}
if [[ -e $OUTDIR/${OUTFILE}.f90 ]] 
then 
 mv $OUTDIR/${OUTFILE}.f90 $OUTDIR/${OUTFILE}.f90.sav
fi
cp -p kk_kpp.f90    $OUTDIR/${OUTFILE}.f90
echo " "
echo "Write kpp module -- > " $OUTDIR/${OUTFILE}.f90

if [[ $UPDT = "YES" ]]
then
cp -p kk_kpp.f90    $DEFDIR/${OUTFILE}.f90
echo " "
echo "Write kpp module -- > " $DEFDIR/${OUTFILE}.f90
fi

if [[ $KEEP = "NO" ]]
then
  rm -rf $WORK
fi
exit

