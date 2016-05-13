# Setup build type for examples
if (CMAKE_BUILD_TYPE MATCHES "Release")
	set(EXAMPLES_BUILD_TYPE "RELEASE")
else()
	set(EXAMPLES_BUILD_TYPE "DEBUG")
endif()

# Include subdirectories:
add_subdirectory(hard_coded_hf)
