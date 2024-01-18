# - Find FFTW
# Find the native FFTW includes and library
#
#  FFTW_INCLUDES    - where to find fftw3.h
#  FFTW_LIBRARIES   - List of libraries when using FFTW.
#  FFTW_FOUND       - True if FFTW found.

if (FFTW_INCLUDES)
  # Already in cache, be silent
  set (FFTW_FIND_QUIETLY TRUE)
endif (FFTW_INCLUDES)

#exec_program(sed ARGS -i '/%fftw_inc.*/d' ${config})
#list(APPEND NETCDF_HINTS "$ENV{LD_LIBRARY_PATH}")

find_path (FFTW_INCLUDES fftw3.f03 HINTS ${FFTW_HINTS} ENV LD_LIBRARY_PATH  PATH_SUFFIXES include Include)
find_library (FFTW_LIBRARIES NAMES fftw3 HINTS ${FFTW_HINTS} ENV LD_LIBRARY_PATH PATH_SUFFIXES lib lib64)

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDES)

mark_as_advanced (FFTW_LIBRARIES FFTW_INCLUDES)

