# - Find RRTMG
# Find the native RRTMG includes and library
#
#  RRTMG_INCLUDES    - where to find rrtmg module files
#  RRTMG_LIBRARIES   - List of libraries when using RRTMG.
#  RRTMG_FOUND       - True if RRTMG found.

if (RRTMG_INCLUDES)
  # Already in cache, be silent
  set (RRTMG_FIND_QUIETLY TRUE)
endif (RRTMG_INCLUDES)

find_path (RRTMG_INCLUDES parkind.mod HINTS ${RRTMG_HINTS} ENV LD_LIBRARY_PATH  PATH_SUFFIXES include Include)
find_library (RRTMG_LIBRARIES NAMES rrtmg HINTS ${RRTMG_HINTS} ENV LD_LIBRARY_PATH PATH_SUFFIXES lib lib64)

# handle the QUIETLY and REQUIRED arguments and set RRTMG_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (RRTMG DEFAULT_MSG RRTMG_LIBRARIES RRTMG_INCLUDES)

mark_as_advanced (RRTMG_LIBRARIES RRTMG_INCLUDES)
