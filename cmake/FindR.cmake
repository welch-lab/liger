# FindR.cmake from MLPACK
# Make sure find package macros are included
set(TEMP_CMAKE_FIND_APPBUNDLE ${CMAKE_FIND_APPBUNDLE})
set(CMAKE_FIND_APPBUNDLE "NEVER")

# Find R.
find_program(R_EXECUTABLE R DOC "R executable.")

if(R_EXECUTABLE)
    # Get the location of R.
    execute_process(
        WORKING_DIRECTORY .
        COMMAND ${R_EXECUTABLE} RHOME
        OUTPUT_VARIABLE R_BASE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # Get the R version.
    execute_process(
        COMMAND ${R_EXECUTABLE} --version
        OUTPUT_VARIABLE R_VERSION_STRING
        RESULT_VARIABLE RESULT
    )

    if(RESULT EQUAL 0)
        string(REGEX REPLACE ".*([0-9]+\\.[0-9]+\\.[0-9]+).*" "\\1"
            R_VERSION_STRING ${R_VERSION_STRING})
    endif()

    set(R_HOME ${R_BASE_DIR} CACHE PATH "R home directory obtained from R RHOME")
    mark_as_advanced(R_HOME)
endif()

# Find the Rscript program.
find_program(RSCRIPT_EXECUTABLE Rscript DOC "Rscript executable.")

set(CMAKE_FIND_APPBUNDLE ${TEMP_CMAKE_FIND_APPBUNDLE})

# Search for non-standard R.h include path if header missing
CHECK_INCLUDE_FILE("R.h" _haveR_h)
if(_haveR_h)
    set(R_INCLUDE_DIR ${_haveR_H})
else()
execute_process(COMMAND ${RSCRIPT_EXECUTABLE} "-e" "R.home('include')"
    RESULT_VARIABLE _haveR_h
    OUTPUT_VARIABLE _R_INCLUDE_location
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

    # Some cleanup in location of R.
    string(REGEX MATCHALL "\".*\"" _R_INCLUDE_location "${_R_INCLUDE_location}")
    string(REGEX REPLACE "\"" "" _R_INCLUDE_location "${_R_INCLUDE_location}")
    set(R_INCLUDE_DIR ${_R_INCLUDE_location})
endif()

mark_as_advanced(RSCRIPT_EXECUTABLE R_EXECUTABLE)
set(_REQUIRED_R_VARIABLES R_EXECUTABLE RSCRIPT_EXECUTABLE R_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    R
    REQUIRED_VARS ${_REQUIRED_R_VARIABLES}
    VERSION_VAR R_VERSION_STRING
    FAIL_MESSAGE "R not found"
)
