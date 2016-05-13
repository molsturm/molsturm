#
# The sources for this example executable
#
set(EXAMPLE_HARD_CODED_HF_SOURCES
	main.cc)

# Compile and link:
add_executable(hard_coded_hf ${EXAMPLE_HARD_CODED_HF_SOURCES})

set_target_properties(
	hard_coded_hf
	PROPERTIES
	LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS_${EXAMPLES_BUILD_TYPE}}" #${CMAKE_EXE_LINKER_FLAGS} are added
	LINKER_LANGUAGE "CXX"
	COMPILE_FLAGS "${CMAKE_CXX_FLAGS_${EXAMPLES_BUILD_TYPE}}" #${CMAKE_CXX_FLAGS} are added automaticaly
	COMPILE_DEFINITIONS "${EXAMPLES_BUILD_TYPE}"
)

# setup definitions.
target_compile_definitions(hard_coded_hf PRIVATE "${MOLSTURM_DEFINITIONS};${MOLSTURM_DEFINITIONS_${EXAMPLES_BUILD_TYPE}}")

# link it with the appropriate linalgwrap library target
target_link_libraries(hard_coded_hf
	${molsturm_${EXAMPLES_BUILD_TYPE}_TARGET}
	${MOLSTURM_DEPENDENCIES}
	${MOLSTURM_DEPENDENCIES_${EXAMPLES_BUILD_TYPE}}
)
