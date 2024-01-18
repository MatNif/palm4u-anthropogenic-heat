# This is the CMake configuration file for INIFOR

# on macOS, use GNU readlink instead of the system BSD version, can be installed
# via Homebrew: 'brew install coreutils'
set(READLINK_CMD "readlink")
if(APPLE)
   set(READLINK_CMD "greadlink")
endif()

cmake_minimum_required (VERSION 2.8)

# set program name
set(INIFOR_PROGRAM_NAME "INIFOR")

project(${INIFOR_PROGRAM_NAME} NONE)
exec_program( hostname OUTPUT_VARIABLE INIFOR_HOSTNAME)

# .palm.config locations
set(makefile_in ${INIFOR_MAKEFILE_TEMPLATE})
set(makefile ${CMAKE_INSTALL_PREFIX}/Makefile)

# include local cmake files
list(APPEND CMAKE_MODULE_PATH ${INIFOR_CMAKE_PATH})

# enabling Fortran language support
enable_language(Fortran) # required to compile the main model and all utilities

# check for netCDF
#set(NETCDF_C "YES")
set(NETCDF_FORTRAN "YES")
find_package (NetCDF REQUIRED)
find_program(NETCDF_FORTRAN_COMPILER_FULL NAMES ${NETCDF_FORTRAN_COMPILER})
execute_process(COMMAND ${READLINK_CMD} -e ${NETCDF_FORTRAN_COMPILER_FULL} OUTPUT_VARIABLE NETCDF_FORTRAN_COMPILER_RESOLVED)
execute_process(COMMAND ${READLINK_CMD} -e ${CMAKE_Fortran_COMPILER} OUTPUT_VARIABLE FORTRAN_COMPILER_RESOLVED)
if(NOT ${FORTRAN_COMPILER_RESOLVED} STREQUAL ${NETCDF_FORTRAN_COMPILER_RESOLVED} )
   message(WARNING "Fortran compiler \"${FORTRAN_COMPILER_RESOLVED}\" does not match netCDF Fortran compiler \"${NETCDF_FORTRAN_COMPILER_RESOLVED}\".")
else()
   message(STATUS "Fortran compiler matches netCDF Fortran compiler.")
endif()

# compiler flag management
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
   set(INIFOR_COMPILER_OPTIONS "-em -K trap=fp -O3 -hnoomp -hnoacc -hfp3 -hdynamic")
   set(INIFOR_LINKER_OPTIONS "${INIFOR_COMPILER_OPTIONS} -dynamic")
   set(INIFOR_CPP_FLAGS "-eZ")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
   set(INIFOR_COMPILER_OPTIONS "-O3 -fp-model source -fno-alias -fpe0 -ftz -no-prec-div -no-prec-sqrt -ip -nbs -diag-disable 8290,8291")
   set(INIFOR_LINKER_OPTIONS "${INIFOR_COMPILER_OPTIONS}")
   set(INIFOR_CPP_FLAGS "-cpp")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
   set(INIFOR_COMPILER_OPTIONS "-O3 -ffree-line-length-none -ffpe-trap=invalid,zero,overflow")
   set(INIFOR_LINKER_OPTIONS "${INIFOR_COMPILER_OPTIONS}")
   set(INIFOR_CPP_FLAGS "-cpp")
endif()

if(NETCDF_FOUND)
   list(APPEND INIFOR_CPP_OPTIONS __netcdf)
   list(APPEND INIFOR_CPP_OPTIONS __netcdf4)
   set(INIFOR_COMPILER_OPTIONS "${INIFOR_COMPILER_OPTIONS} -I ${NETCDF_INCLUDES}")
   set(INIFOR_LINKER_OPTIONS "${INIFOR_LINKER_OPTIONS} ${NETCDF_LIBRARIES}")
endif(NETCDF_FOUND)

string(REPLACE ";" " -D" INIFOR_CPP_OPTIONS_STR_1 "${INIFOR_CPP_OPTIONS}")
set(INIFOR_COMPILER_OPTIONS "${INIFOR_COMPILER_OPTIONS} ${INIFOR_CPP_FLAGS} -D${INIFOR_CPP_OPTIONS_STR_1}")

# configuring the Makefile
configure_file(${makefile_in} ${makefile} @ONLY)
