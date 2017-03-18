#
# The sources for the test executable
#
set(MOLSTURM_TESTS_SOURCES
	FunctionalityTest.cc
	HfEnergyTest.cc
	main.cc
)

drb_add_executable(molsturm_tests
	DBGSUFFIX  "_DEBUG"
	RELSUFFIX  "_RELEASE" 
	FILES ${MOLSTURM_TESTS_SOURCES}
)

drb_target_compile_definitions(ALL     molsturm_tests 
	PRIVATE ${MOLSTURM_DEFINITIONS})
drb_target_compile_definitions(DEBUG   molsturm_tests 
	PRIVATE ${MOLSTURM_DEFINITIONS_DEBUG})
drb_target_compile_definitions(RELEASE molsturm_tests 
	PRIVATE ${MOLSTURM_DEFINITIONS_RELEASE})

drb_target_link_libraries(ALL     molsturm_tests 
	${MOLSTURM_DEPENDENCIES} ${MOLSTURM_DEPENDENCIES_TEST})
drb_target_link_libraries(DEBUG   molsturm_tests 
	${molsturm_DEBUG_TARGET} ${MOLSTURM_DEPENDENCIES_DEBUG})
drb_target_link_libraries(RELEASE molsturm_tests 
	${molsturm_RELEASE_TARGET} ${MOLSTURM_DEPENDENCIES_RELEASE})

# Add the test executables
FOREACH(target ${molsturm_tests_TARGETS})
	add_test(NAME ${target} COMMAND ${target})
ENDFOREACH()

