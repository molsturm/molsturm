# sets these things
#
# 	MOLSTURM_DEPENDENCIES			everyone needs these libraries
# 	MOLSTURM_DEPENDENCIES_DEBUG		debug mode needs these extras
# 	MOLSTURM_DEPENDENCIES_RELEASE		release mode needs these extras
# 	MOLSTURM_DEPENDENCIES_TEST		tests need these extra libraries
#
#       MOLSTURM_DEFINITIONS			definitions for all compilation
#       MOLSTURM_DEFINITIONS_DEBUG		definitions for debug mode
#       MOLSTURM_DEFINITIONS_RELEASE		definitions for release mode
#       

####################
#-- Empty it all --#
####################
set(MOLSTURM_DEPENDENCIES "")
set(MOLSTURM_DEPENDENCIES_DEBUG "")
set(MOLSTURM_DEPENDENCIES_RELEASE "")
set(MOLSTURM_DEPENDENCIES_TEST "")
set(MOLSTURM_DEFINITIONS "")
set(MOLSTURM_DEFINITIONS_DEBUG "")
set(MOLSTURM_DEFINITIONS_RELEASE "")

##############
#-- Macros --#
##############

macro(add_submodule_dependency LIBRARY VERSION)
	set(${LIBRARY}_DIR "${PROJECT_SOURCE_DIR}/modules/${LIBRARY}")
	if(NOT EXISTS "${${LIBRARY}_DIR}/CMakeLists.txt")
		message(FATAL_ERROR "Did not find expected submodule ${LIBRARY} in ${${LIBRARY}_DIR}. \
Try \"git submodule update --init --recursive\"")
	endif()

	# Extract version from CMakeLists.txt:
	file(STRINGS "${${LIBRARY}_DIR}/CMakeLists.txt" VERSION_RAW
		REGEX "${LIBRARY} VERSION [0-9.]+"
		LIMIT_COUNT 1)
	string(REGEX MATCH "[0-9.]+" ${LIBRARY}_VERSION "${VERSION_RAW}")

	# Compare against what is needed
	if("${${LIBRARY}_VERSION}" VERSION_LESS "${VERSION}")
		message(FATAL_ERROR "Inconsistency in the repo: Version ${VERSION} of ${LIBRARY} \
was requested, but only version ${${LIBRARY}_VERSION} was found. Maybe a \
\"git submodule update --init --recursive\" helps?")
	endif()

	# Set the project up:
	add_subdirectory("${${LIBRARY}_DIR}")
	include_directories("${${LIBRARY}_DIR}/src")

	# Add dependencies:
	foreach(build ${DRB_BUILD_TYPES})
		set(MOLSTURM_DEPENDENCIES_${build} ${MOLSTURM_DEPENDENCIES_${build}} ${${LIBRARY}_${build}_TARGET})
	endforeach()
endmacro(add_submodule_dependency)

############################
#-- rapidcheck and catch --#
############################
if (MOLSTURM_ENABLE_TESTS)
	# We need rapidcheck and catch for the tests:
	include(cmake/setupRapidcheckCatch.cmake)
	set(MOLSTURM_DEPENDENCIES_TEST ${MOLSTURM_DEPENDENCIES_TEST} ${rapidcheck_TARGET} ${catch_TARGET})
endif()

#######################
#-- Other libraries --#
#######################

add_submodule_dependency(krims 0.0.0)
add_submodule_dependency(linalgwrap 0.2.0)
add_submodule_dependency(sturmint 0.0.0)
add_submodule_dependency(gint 0.0.0)
add_submodule_dependency(gscf 0.0.0)
