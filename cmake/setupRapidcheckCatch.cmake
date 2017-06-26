# setup the rapidcheck and catch libraries in
# external/rapidcheck for build.

#
# Check that recursive submodule checkout has been performed
#
set(rapidcheck_DIR "${PROJECT_SOURCE_DIR}/modules/rapidcheck")
if (NOT EXISTS "${rapidcheck_DIR}/CMakeLists.txt")
	message(FATAL_ERROR "Could not find rapidcheck submodule in expected directory \
${rapidcheck_DIR}. Try \"git submodule update --init --recursive\"")
endif()

set(catch_INCLUDE_DIR "${rapidcheck_DIR}/ext/catch/include")
if (NOT EXISTS "${catch_INCLUDE_DIR}/catch.hpp")
	message(FATAL_ERROR "Could not find catch header at expected place \
${catch_INCLUDE_DIR}. Try \"git submodule update --init --recursive\"")
endif()

#
# Configure rapidcheck
#
set(rapidcheck_TARGET rapidcheck)

# Set option such that rapidcheck tests are built.
set(RC_ENABLE_TESTS ON CACHE BOOL "Build RapidCheck tests")

# Change compiler flags (CMAKE_CXX_FLAGS) to make fresh build config
set(CMAKE_CXX_FLAGS_STORED_TMP ${CMAKE_CXX_FLAGS})
set(CMAKE_CXX_FLAGS "")

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU"
		AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "6.0")
	# TODO Problems in rapidcheck's code
	#      Fix them some day
	enable_if_compiles(CMAKE_CXX_FLAGS "-Wno-error=shift-negative-value")
	enable_if_compiles(CMAKE_CXX_FLAGS_STORED_TMP "-Wno-error=shift-negative-value")
	enable_if_compiles(CMAKE_CXX_FLAGS "-Wno-error=misleading-indentation")
	enable_if_compiles(CMAKE_CXX_FLAGS_STORED_TMP "-Wno-error=misleading-indentation")
endif()

stdlib_cxx_flags(${DRB_CXX_STANDARD_LIBRARY} CMAKE_CXX_FLAGS dummy)
unset(dummy)

# Add the rapidcheck subdirectory and configure its built.
message(STATUS "Configuring rapidcheck in dir ${rapidcheck_DIR}")
add_subdirectory(${rapidcheck_DIR})

# undo the changes to the compiler flags:
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_STORED_TMP})
unset(CMAKE_CXX_FLAGS_STORED_TMP)

#
# Configure catch
#
set(catch_TARGET catch_hdr)
add_library(${catch_TARGET} INTERFACE IMPORTED)
message(STATUS "Using catch from dir ${catch_INCLUDE_DIR}")
include_directories(${catch_INCLUDE_DIR})
