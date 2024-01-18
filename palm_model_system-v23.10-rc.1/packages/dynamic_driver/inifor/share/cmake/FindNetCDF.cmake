# - Find NetCDF
# Find the native NetCDF includes and library
#
#  NETCDF_INCLUDES    - where to find netcdf.h, etc
#  NETCDF_LIBRARIES   - Link these libraries when using NetCDF
#  NETCDF_FOUND       - True if NetCDF found including required interfaces (see below)
#
# Your package can require certain interfaces to be FOUND by setting these
#
#  NETCDF_C           - require the C interface and link the C library
#  NETCDF_CXX         - require the C++ interface and link the C++ library
#  NETCDF_FORTRAN     - require the Fortran interface and link the Fortran library
#
# The following are not for general use and are included in
# NETCDF_LIBRARIES if the corresponding option above is set.
#
#  NETCDF_LIBRARIES_C        - Just the C interface
#  NETCDF_LIBRARIES_CXX      - C++ interface, if available
#  NETCDF_LIBRARIES_FORTRAN  - Fortran 90 interface, if available
#
# Normal usage would be:
#  set (NETCDF_FORTRAN "YES")
#  find_package (NetCDF REQUIRED)
#  target_link_libraries (uses_f90_interface ${NETCDF_LIBRARIES})

if (NETCDF_INCLUDES AND NETCDF_LIBRARIES)
  # Already in cache, be silent
  set (NETCDF_FIND_QUIETLY TRUE)
endif (NETCDF_INCLUDES AND NETCDF_LIBRARIES)


macro(NETCDF_CONFIG flag output)
    if(NETCDF_CONFIG_EXECUTABLE)
        exec_program( ${NETCDF_CONFIG_EXECUTABLE} ARGS ${flag}
            OUTPUT_VARIABLE ${output} RETURN_VALUE return_value)
        if(NOT ${return_value} EQUAL 0 )
            message( STATUS "Unable to determine ${flag} from ${NETCDF_CONFIG_EXECUTABLE}." )
        endif()
    endif(NETCDF_CONFIG_EXECUTABLE)
endmacro()

if(NETCDF_C_ROOT)
   list(APPEND NETCDF_HINTS "${NETCDF_C_ROOT}")
else()
   list(APPEND NETCDF_HINTS "$ENV{NETCDF_ROOT}")
endif()

if(NETCDF_FORTRAN_ROOT)
   list(APPEND NETCDF_HINTS "${NETCDF_FORTRAN_ROOT}")
else()
   list(APPEND NETCDF_HINTS "$ENV{NETCDF_ROOT}")
endif()

if(NETCDF_C_ROOT)
   find_program(NETCDF_C_CONFIG_EXECUTABLE NAMES nc-config
      HINTS ${NETCDF_HINTS} PATH_SUFFIXES bin Bin NO_DEFAULT_PATH
      DOC "NETCDF CONFIG PROGRAM. Used to detect NETCDF compile flags." )
else()
   find_program(NETCDF_C_CONFIG_EXECUTABLE NAMES nc-config
      HINTS ${NETCDF_HINTS} PATH_SUFFIXES bin Bin
      DOC "NETCDF CONFIG PROGRAM. Used to detect NETCDF compile flags." )
endif()

set(NETCDF_CONFIG_EXECUTABLE ${NETCDF_C_CONFIG_EXECUTABLE})
   if(NETCDF_C_CONFIG_EXECUTABLE)
      NETCDF_CONFIG(--cc NETCDF_C_COMPILER_C)
      NETCDF_CONFIG(--fc NETCDF_C_COMPILER_FORTRAN)
      NETCDF_CONFIG(--prefix NETCDF_C_ROOT)
      NETCDF_CONFIG(--includedir NETCDF_C_INCLUDE)
      NETCDF_CONFIG(--version NETCDF_C_VERSION)
      #NETCDF_CONFIG(--has-c++ NETCDF_C_CXX)
      #NETCDF_CONFIG(--has-f77 NETCDF_C_F77)
      NETCDF_CONFIG(--has-f90 NETCDF_C_F90)
      #NETCDF_CONFIG(--has-dap NETCDF_C_DAP)
      #NETCDF_CONFIG(--has-nc2 NETCDF_C_NC2)
      #NETCDF_CONFIG(--has-nc4 NETCDF_C_NC4)
      #NETCDF_CONFIG(--has-hdf4 NETCDF_C_HDF4)
      #NETCDF_CONFIG(--has-hdf5 NETCDF_C_HDF5)
      #NETCDF_CONFIG(--has-pnetcdf NETCDF_C_PARALLEL)
      list(APPEND NETCDF_INCLUDE_HINTS "${NETCDF_C_INCLUDE}")
      list(APPEND NETCDF_HINTS "${NETCDF_C_ROOT}")
      message(STATUS "Found ${NETCDF_C_VERSION} compiled with ${NETCDF_C_COMPILER_C}")
   else(NETCDF_C_CONFIG_EXECUTABLE)
      message(STATUS "nc-config not found")
   endif(NETCDF_C_CONFIG_EXECUTABLE)

if(NETCDF_C_ROOT AND NETCDF_FORTRAN_ROOT)
   find_program(NETCDF_FORTRAN_CONFIG_EXECUTABLE NAMES nf-config
       HINTS ${NETCDF_HINTS} PATH_SUFFIXES bin Bin NO_DEFAULT_PATH
       DOC "NETCDF CONFIG PROGRAM. Used to detect NETCDF compile flags." )
else()
   find_program(NETCDF_FORTRAN_CONFIG_EXECUTABLE NAMES nf-config
       HINTS ${NETCDF_HINTS} PATH_SUFFIXES bin Bin
       DOC "NETCDF CONFIG PROGRAM. Used to detect NETCDF compile flags." )
endif()

set(NETCDF_CONFIG_EXECUTABLE ${NETCDF_FORTRAN_CONFIG_EXECUTABLE})
   if(NETCDF_FORTRAN_CONFIG_EXECUTABLE)
      NETCDF_CONFIG(--cc NETCDF_FORTRAN_COMPILER_C)
      NETCDF_CONFIG(--fc NETCDF_FORTRAN_COMPILER_FORTRAN)
      NETCDF_CONFIG(--prefix NETCDF_FORTRAN_ROOT)
      NETCDF_CONFIG(--includedir NETCDF_FORTRAN_INCLUDE)
      NETCDF_CONFIG(--version NETCDF_FORTRAN_VERSION)
      #NETCDF_CONFIG(--has-c++ NETCDF_FORTRAN_CXX)
      #NETCDF_CONFIG(--has-f77 NETCDF_FORTRAN_F77)
      NETCDF_CONFIG(--has-f90 NETCDF_FORTRAN_F90)
      #NETCDF_CONFIG(--has-dap NETCDF_FORTRAN_DAP)
      #NETCDF_CONFIG(--has-nc2 NETCDF_FORTRAN_NC2)
      #NETCDF_CONFIG(--has-nc4 NETCDF_FORTRAN_NC4)
      #NETCDF_CONFIG(--has-hdf4 NETCDF_FORTRAN_HDF4)
      #NETCDF_CONFIG(--has-hdf5 NETCDF_FORTRAN_HDF5)
      #NETCDF_CONFIG(--has-pnetcdf NETCDF_FORTRAN_PARALLEL)
      list(APPEND NETCDF_INCLUDE_HINTS "${NETCDF_FORTRAN_INCLUDE}")
      list(APPEND NETCDF_HINTS "${NETCDF_FORTRAN_ROOT}")
      message(STATUS "Found ${NETCDF_FORTRAN_VERSION} compiled with ${NETCDF_FORTRAN_COMPILER_FORTRAN}")
   else(NETCDF_FORTRAN_CONFIG_EXECUTABLE)
      #message(STATUS "nf-config not found")
      set(NETCDF_FORTRAN_COMPILER_C ${NETCDF_C_COMPILER_C})
      set(NETCDF_FORTRAN_COMPILER_FORTRAN ${NETCDF_C_COMPILER_FORTRAN})
      set(NETCDF_FORTRAN_ROOT ${NETCDF_C_ROOT})
      set(NETCDF_FORTRAN_INCLUDE ${NETCDF_C_INCLUDE})
      set(NETCDF_FORTRAN_VERSION ${NETCDF_C_VERSION})
      #set(NETCDF_FORTRAN_CXX ${NETCDF_C_CXX})
      #set(NETCDF_FORTRAN_F77 ${NETCDF_C_F77})
      set(NETCDF_FORTRAN_F90 ${NETCDF_C_F90})
      #set(NETCDF_FORTRAN_DAP ${NETCDF_C_DAP})
      #set(NETCDF_FORTRAN_NC2 ${NETCDF_C_NC2})
      #set(NETCDF_FORTRAN_NC4 ${NETCDF_C_NC4})
      #set(NETCDF_FORTRAN_HDF4 ${NETCDF_C_HDF4})
      #set(NETCDF_FORTRAN_HDF5 ${NETCDF_C_HDF5})
      #set(NETCDF_FORTRAN_PARALLEL ${NETCDF_C_PARALLEL})
      if(NETCDF_FORTRAN_F90)
         message(STATUS "Found ${NETCDF_FORTRAN_VERSION} compiled with ${NETCDF_FORTRAN_COMPILER_FORTRAN}")
      else(NETCDF_FORTRAN_F90)
         message(STATUS "nc-config found no netCDF Fortran libraries")
      endif(NETCDF_FORTRAN_F90)
   endif(NETCDF_FORTRAN_CONFIG_EXECUTABLE)

# find netcdf c
if(NOT NETCDF_C_INCLUDE)
   find_path(NETCDF_C_INCLUDE netcdf.h HINTS ${NETCDF_HINTS} PATH_SUFFIXES include Include)
endif()
find_library(NETCDF_C_LIB netcdf HINTS ${NETCDF_HINTS} PATH_SUFFIXES lib lib64)
#message(STATUS "NETCDF_C_INCLUDE so far: ${NETCDF_C_INCLUDE}")
#message(STATUS "NETCDF_C_LIB so far: ${NETCDF_C_LIB}")

# find netcdf fortran
if(NOT NETCDF_FORTRAN_INCLUDE)
   find_path(NETCDF_FORTRAN_INCLUDE netcdf.mod HINTS ${NETCDF_HINTS} PATH_SUFFIXES include Include)
endif()
find_library(NETCDF_FORTRAN_LIB netcdff HINTS ${NETCDF_HINTS} PATH_SUFFIXES lib lib64)
if(NOT NETCDF_FORTRAN_LIB)
   find_library(NETCDF_FORTRAN_LIB netcdf HINTS ${NETCDF_HINTS} PATH_SUFFIXES lib lib64)
endif()
#message(STATUS "NETCDF_FORTRAN_INCLUDE so far: ${NETCDF_FORTRAN_INCLUDE}")
#message(STATUS "NETCDF_FORTRAN_LIB so far: ${NETCDF_FORTRAN_LIB}")

if ((NOT NETCDF_C_LIB) OR (NOT NETCDF_C_INCLUDE))
   message(STATUS "Trying to find NetCDF using LD_LIBRARY_PATH (we're desperate)...")
   file(TO_CMAKE_PATH "$ENV{LD_LIBRARY_PATH}" LD_LIBRARY_PATH)
   find_library(NETCDF_C_LIB NAMES netcdf HINTS ${LD_LIBRARY_PATH})

   if (NETCDF_C_LIB)
      get_filename_component(NETCDF_LIB_DIR ${NETCDF_C_LIB} PATH)
      string(REGEX REPLACE "/lib/?$" "/include" NETCDF_H_HINT ${NETCDF_LIB_DIR})
      find_path (NETCDF_C_INCLUDE netcdf.h HINTS ${NETCDF_H_HINT} DOC "Path to netcdf.h")
      message(STATUS "found netcdf.h in: ${NETCDF_C_INCLUDE}")
      list(APPEND NETCDF_INCLUDE_HINTS "${NETCDF_C_INCLUDE}")
   endif()
endif()

get_filename_component (NETCDF_C_LIB_DIR "${NETCDF_C_LIB}" PATH)
get_filename_component (NETCDF_FORTRAN_LIB_DIR "${NETCDF_FORTRAN_LIB}" PATH)
list(APPEND NETCDF_LIB_HINTS "${NETCDF_C_LIB_DIR}")
list(APPEND NETCDF_LIB_HINTS "${NETCDF_FORTRAN_LIB_DIR}")

#message(STATUS "All include Hints: ${NETCDF_INCLUDE_HINTS}")
#message(STATUS "All lib Hints: ${NETCDF_LIB_HINTS}")

macro(NetCDF_add_interface lang)
   if(NETCDF_${lang})
      if(NETCDF_${lang}_INCLUDE AND NETCDF_${lang}_LIB)
         list(INSERT NetCDF_includes 0 ${NETCDF_${lang}_INCLUDE})
         list(INSERT NetCDF_libs 0 ${NETCDF_${lang}_LIB}) # prepend so that -lnetcdf is last
      else()
         set(NetCDF_has_interfaces "NO")
         message(STATUS "Failed to find NetCDF interface for ${lang}")
      endif()
   endif(NETCDF_${lang})
endmacro(NetCDF_add_interface)

set(NetCDF_has_interfaces "YES") # will be set to NO if we're missing any interfaces
NetCDF_add_interface(C)
NetCDF_add_interface(CXX)
NetCDF_add_interface(FORTRAN)

# macro (NetCDF_check_interface lang header libs)
#    if (NETCDF_${lang})
#       find_path (NETCDF_INCLUDES_${lang} NAMES ${header} HINTS ${NETCDF_HINTS} PATH_SUFFIXES include Include NO_DEFAULT_PATH)
#       find_library (NETCDF_LIBRARIES_${lang} NAMES ${libs} HINTS ${NETCDF_HINTS} PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
#       mark_as_advanced (NETCDF_INCLUDES_${lang} NETCDF_LIBRARIES_${lang})
#       if (NETCDF_INCLUDES_${lang} AND NETCDF_LIBRARIES_${lang})
#          list (INSERT NetCDF_libs 0 ${NETCDF_LIBRARIES_${lang}}) # prepend so that -lnetcdf is last
#       else (NETCDF_INCLUDES_${lang} AND NETCDF_LIBRARIES_${lang})
#          set (NetCDF_has_interfaces "NO")
#          message (STATUS "Failed to find NetCDF interface for ${lang}")
#       endif (NETCDF_INCLUDES_${lang} AND NETCDF_LIBRARIES_${lang})
#    endif (NETCDF_${lang})
# endmacro (NetCDF_check_interface)
#
# set (NetCDF_has_interfaces "YES") # will be set to NO if we're missing any interfaces
# NetCDF_check_interface (C netcdf.h netcdf)
# NetCDF_check_interface (CXX netcdfcpp.h netcdf_c++)
# NetCDF_check_interface (FORTRAN netcdf.mod  netcdff)

set (NETCDF_C_COMPILER "${NETCDF_C_COMPILER_C}" CACHE STRING "The C compiler used to build netCDF")
set (NETCDF_FORTRAN_COMPILER "${NETCDF_FORTRAN_COMPILER_FORTRAN}" CACHE STRING "The Fortran compiler used to build netCDF")
set (NETCDF_INCLUDES "${NetCDF_includes}" CACHE STRING "All NetCDF includes required for interface level")
set (NETCDF_LIBRARIES "${NetCDF_libs}" CACHE STRING "All NetCDF libraries required for interface level")

# handle the QUIETLY and REQUIRED arguments and set NETCDF_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (NetCDF DEFAULT_MSG NETCDF_LIBRARIES NETCDF_INCLUDES NetCDF_has_interfaces)

mark_as_advanced (NETCDF_LIBRARIES NETCDF_INCLUDES)

#message(STATUS "netCDF library: ${NETCDF_LIBRARIES}")
#message(STATUS "netCDF include: ${NETCDF_INCLUDES}")

