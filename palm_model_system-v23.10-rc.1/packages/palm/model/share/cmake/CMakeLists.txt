# This is the CMake configuration file for PALM

cmake_minimum_required (VERSION 2.8)

# set program name
set(PALM_PROGRAM_NAME "PALM")

project(${PALM_PROGRAM_NAME} NONE)
exec_program( hostname OUTPUT_VARIABLE PALM_HOSTNAME)

# .palm.config locations
set(config_in ${PALM_CONFIG_TEMPLATE})
set(config ${CMAKE_INSTALL_PREFIX}/.palm.config.${PALM_CONFIGURATION_NAME})

# palmtest.yml locations
set(palmtest_in ${PALMTEST_YAML_TEMPLATE})
set(palmtest ${CMAKE_INSTALL_PREFIX}/palmtest.yml)

# build_rrtmg locations
set(build_rrtmg_in ${CMAKE_CURRENT_SOURCE_DIR}/rrtmg_build_options.in)
set(build_rrtmg ${CMAKE_CURRENT_SOURCE_DIR}/../rrtmg/.rrtmg_build_options)

# build_tenstream locations
set(build_tenstream_in ${CMAKE_CURRENT_SOURCE_DIR}/tenstream_build_options.in)
set(build_tenstream ${CMAKE_CURRENT_SOURCE_DIR}/../tenstream/.tenstream_build_options)

# include local cmake files
list(APPEND CMAKE_MODULE_PATH ${PALM_CMAKE_PATH})

# enabling C language support
enable_language(C) # required to compile the FASTv8 Server

# enabling Fortran language support
enable_language(Fortran) # required to compile the main model and all utilities

# check for Fortran MPI support
find_package(MPI REQUIRED)

# check for netCDF
#set(NETCDF_C "YES")
set(NETCDF_FORTRAN "YES")
find_package (NetCDF REQUIRED)
find_program(NETCDF_FORTRAN_COMPILER_FULL NAMES ${NETCDF_FORTRAN_COMPILER})
execute_process(COMMAND readlink -e ${NETCDF_FORTRAN_COMPILER_FULL} OUTPUT_VARIABLE NETCDF_FORTRAN_COMPILER_RESOLVED)
execute_process(COMMAND readlink -e ${CMAKE_Fortran_COMPILER} OUTPUT_VARIABLE FORTRAN_COMPILER_RESOLVED)
if(NOT ${FORTRAN_COMPILER_RESOLVED} STREQUAL ${NETCDF_FORTRAN_COMPILER_RESOLVED} )
   message(WARNING "Fortran compiler \"${FORTRAN_COMPILER_RESOLVED}\" does not match netCDF Fortran compiler \"${NETCDF_FORTRAN_COMPILER_RESOLVED}\".")
else()
   message(STATUS "Fortran compiler matches netCDF Fortran compiler.")
endif()

# check for fftw
find_package(FFTW)
if(FFTW_FOUND)
   message(STATUS "${PALM_PROGRAM_NAME} is using an external fftw library.")
else(FFTW_FOUND)
   message(STATUS "${PALM_PROGRAM_NAME} is using the buildin fft algorithm.")
endif(FFTW_FOUND)

# check for rrtmg
find_package(RRTMG)
if(RRTMG_FOUND)
   message(STATUS "${PALM_PROGRAM_NAME} is using the rrmtg library.")
else(RRTMG_FOUND)
   message(STATUS "${PALM_PROGRAM_NAME} is not using the rrtmg library.")
endif(RRTMG_FOUND)

# check for TenStream
find_package(TenStream)
if(TENSTREAM_FOUND)
  message(STATUS "${PALM_PROGRAM_NAME} is using the TenStream library.")
else(TENSTREAM_FOUND)
  message(STATUS "${PALM_PROGRAM_NAME} is not using the TenStream library.")
endif(TENSTREAM_FOUND)

# compiler flag management
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
   set(PALM_COMPILER_OPTIONS "-em -K trap=fp -O3 -hnoomp -hnoacc -hfp3 -hdynamic")
   set(PALM_LINKER_OPTIONS "${PALM_COMPILER_OPTIONS} -dynamic")
   set(PALM_CPP_FLAGS "-eZ")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
   set(PALM_COMPILER_OPTIONS "-O3 -fp-model source -fno-alias -fpe0 -ftz -no-prec-div -no-prec-sqrt -ip -nbs -diag-disable 8290,8291")
   set(PALM_LINKER_OPTIONS "${PALM_COMPILER_OPTIONS}")
   set(PALM_CPP_FLAGS "-cpp")
   list(APPEND PALM_CPP_OPTIONS __intel_compiler)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
   set(PALM_COMPILER_OPTIONS "-Ofast -ffree-line-length-none")
   set(PALM_LINKER_OPTIONS "${PALM_COMPILER_OPTIONS}")
   set(PALM_CPP_FLAGS "-cpp")
   list(APPEND PALM_CPP_OPTIONS __gfortran)
endif()
if(${CMAKE_C_COMPILER_ID} STREQUAL "Cray")
   set(PALM_C_COMPILER_OPTIONS "-O3")
elseif(${CMAKE_C_COMPILER_ID} STREQUAL "Intel")
   set(PALM_C_COMPILER_OPTIONS "-O3")
elseif(${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
   set(PALM_C_COMPILER_OPTIONS "-Ofast -Wno-format -Wno-unused-result")
endif()

# workaround
list(APPEND PALM_CPP_OPTIONS __parallel)

# adding cpp options
list(APPEND PALM_CPP_OPTIONS MPI_REAL=MPI_DOUBLE_PRECISION)
list(APPEND PALM_CPP_OPTIONS MPI_2REAL=MPI_2DOUBLE_PRECISION)
if(NETCDF_FOUND)
   list(APPEND PALM_CPP_OPTIONS __netcdf)
   set(PALM_COMPILER_OPTIONS "${PALM_COMPILER_OPTIONS} -I ${NETCDF_INCLUDES}")
   set(PALM_LINKER_OPTIONS "${PALM_LINKER_OPTIONS} ${NETCDF_LIBRARIES}")
   set(PALMTEST_NETCDF_INCLUDES "netcdf: \"-I ${NETCDF_INCLUDES}\"")
   set(PALMTEST_NETCDF_LIBRARIES "netcdf: \"${NETCDF_LIBRARIES}\"")
endif(NETCDF_FOUND)

set(RRTMG_COMPILER_OPTIONS "${PALM_COMPILER_OPTIONS}")
set(RRTMG_LINKER_OPTIONS "${PALM_LINKER_OPTIONS}")

if(FFTW_FOUND)
   list(APPEND PALM_CPP_OPTIONS __fftw)
   set(PALM_COMPILER_OPTIONS "${PALM_COMPILER_OPTIONS} -I ${FFTW_INCLUDES}")
   set(PALM_LINKER_OPTIONS "${PALM_LINKER_OPTIONS} ${FFTW_LIBRARIES}")
   set(PALMTEST_FFTW_INCLUDES "fftw: \"-I ${FFTW_INCLUDES}\"")
   set(PALMTEST_FFTW_LIBRARIES "fftw: \"${FFTW_LIBRARIES}\"")
endif(FFTW_FOUND)

if(RRTMG_FOUND)
   list(APPEND PALM_CPP_OPTIONS __rrtmg)
   set(PALM_COMPILER_OPTIONS "${PALM_COMPILER_OPTIONS} -I ${RRTMG_INCLUDES}")
   set(PALM_LINKER_OPTIONS "${PALM_LINKER_OPTIONS} ${RRTMG_LIBRARIES}")
   set(PALMTEST_RRTMG_INCLUDES "rrtmg: \"-I ${RRTMG_INCLUDES}\"")
   set(PALMTEST_RRTMG_LIBRARIES "rrtmg: \"${RRTMG_LIBRARIES}\"")
endif(RRTMG_FOUND)

if(TENSTREAM_FOUND)
   list(APPEND PALM_CPP_OPTIONS __tenstream)
   string(REPLACE ";" " " PALM_TENSTREAM_CFLAGS "${TENSTREAM_STATIC_CFLAGS}")
   string(REPLACE ";" " " PALM_TENSTREAM_LDFLAGS "${TENSTREAM_STATIC_LDFLAGS}")
   set(PALM_COMPILER_OPTIONS "${PALM_COMPILER_OPTIONS} ${PALM_TENSTREAM_CFLAGS}")
   set(PALM_LINKER_OPTIONS "${PALM_LINKER_OPTIONS} ${PALM_TENSTREAM_LDFLAGS}")
   set(PALMTEST_TENSTREAM_INCLUDES "tenstream: \"${TENSTREAM_STATIC_CFLAGS}\"")
   set(PALMTEST_TENSTREAM_LIBRARIES "tenstream: \"${TENSTREAM_STATIC_LDFLAGS}\"")
endif(TENSTREAM_FOUND)

string(REPLACE ";" " -D" PALM_CPP_OPTIONS_STR_1 "${PALM_CPP_OPTIONS}")
set(PALM_CPP_OPTIONS_STR "${PALM_CPP_FLAGS} -D${PALM_CPP_OPTIONS_STR_1}")

set(PALM_EXECUTE_COMMAND "mpirun -n {{mpi_tasks}}")

# configuring the .palm.config
configure_file(${config_in} ${config} @ONLY)

# configuring the palmtest.yml
configure_file(${palmtest_in} ${palmtest} @ONLY)

# export configuring for rrtmg
configure_file(${build_rrtmg_in} ${build_rrtmg} @ONLY)

# export configuring for TenStream
configure_file(${build_tenstream_in} ${build_tenstream} @ONLY)
