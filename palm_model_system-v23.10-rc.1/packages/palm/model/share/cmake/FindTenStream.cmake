# - Find TenStream
# Find the TenStream includes and library
# We use the pkgconfig file emitted from TenStream.
# First add the path to TenStream to the CMAKE_PREFIX_PATH (new since cmake 3.1)
# or prepend to the PKG_CONFIG_PATH environment variable.
#
# This exports a couple of variables with the TENSTREAM prefix.
# To show all of them uncomment the below debug lines.
#
# Additionally add the TenStream path to the binary RPATH to make sure we can find
# all shared libs from TenStream/PETSc when running PALM
#
# The ones we use are:
#
#  TENSTREAM_STATIC_CFLAGS  - compile flags used to build TenStream
#  TENSTREAM_STATIC_LDFLAGS - linker flags used to build TenStream
#  TENSTREAM_FOUND          - True if found.

find_path (TENSTREAM_PKGCONFIG
  TenStream.pc
  HINTS ${TENSTREAM_HINTS}
  PATH_SUFFIXES pkgconfig)
list(PREPEND CMAKE_PREFIX_PATH ${TENSTREAM_PKGCONFIG})
set(ENV{PKG_CONFIG_PATH} "${TENSTREAM_PKGCONFIG}:$ENV{PKG_CONFIG_PATH}")

include(FindPkgConfig)
pkg_search_module(TENSTREAM TenStream)

if(TENSTREAM_FOUND)
  # list_all_cmake_vars for debugging purpose
  # get_cmake_property(_variableNames VARIABLES)
  # list (SORT _variableNames)
  # foreach (_variableName ${_variableNames})
  #   message(STATUS "${_variableName}=${${_variableName}}")
  # endforeach()

  # TODO: check if this way of providing the rpath is portable
  list(PREPEND TENSTREAM_STATIC_LDFLAGS "-Wl,-rpath,${TENSTREAM_LIBDIR}")

  message(STATUS "TENSTREAM_STATIC_CFLAGS ${TENSTREAM_STATIC_CFLAGS}")
  message(STATUS "TENSTREAM_STATIC_LDFLAGS ${TENSTREAM_STATIC_LDFLAGS} ")
endif(TENSTREAM_FOUND)
