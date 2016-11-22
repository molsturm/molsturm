#
# macro to setup an example target
#
function(setup_example_target TARGET)
	# usually the examples are built and linked against RELEASE
	# if this version is available
	# Supply PREFER_DEBUG as ARGV1 if you want to have the opposite
	# preference

	# Select build type for example target:
	if ("${ARGV1}" STREQUAL "PREFER_DEBUG" AND CMAKE_BUILD_TYPE MATCHES "Debug")
		set(EXAMPLE_BUILD_TYPE "DEBUG")
	elseif(CMAKE_BUILD_TYPE MATCHES "Release")
		set(EXAMPLE_BUILD_TYPE "RELEASE")
	else()
		set(EXAMPLE_BUILD_TYPE "DEBUG")
	endif()

	set_target_properties(
		${TARGET}
		PROPERTIES
		#${CMAKE_EXE_LINKER_FLAGS} are added
		LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS_${EXAMPLE_BUILD_TYPE}}"
		LINKER_LANGUAGE "CXX"
		#${CMAKE_CXX_FLAGS} are added automaticaly
		COMPILE_FLAGS "${CMAKE_CXX_FLAGS_${EXAMPLE_BUILD_TYPE}}"
		COMPILE_DEFINITIONS "${EXAMPLE_BUILD_TYPE}"
	)

	# setup definitions.
	target_compile_definitions(${TARGET} PRIVATE "${MOLSTURM_DEFINITIONS};${MOLSTURM_DEFINITIONS_${EXAMPLES_BUILD_TYPE}}")

	# link it with the appropriate molsturm library target
	target_link_libraries(${TARGET}
		${molsturm_${EXAMPLE_BUILD_TYPE}_TARGET}
		${MOLSTURM_DEPENDENCIES}
		${MOLSTURM_DEPENDENCIES_${EXAMPLE_BUILD_TYPE}}
	)
endfunction()

# Include subdirectories:
add_subdirectory(hf)
