# sets these things
#
#       MOLSTURM_DEPENDENCIES			everyone needs these libraries
#       MOLSTURM_DEPENDENCIES_DEBUG		debug mode needs these extras
#       MOLSTURM_DEPENDENCIES_RELEASE		release mode needs these extras
#       MOLSTURM_DEPENDENCIES_TEST		tests need these extra libraries
#

####################
#-- Empty it all --#
####################
set(MOLSTURM_DEPENDENCIES "")
set(MOLSTURM_DEPENDENCIES_DEBUG "")
set(MOLSTURM_DEPENDENCIES_RELEASE "")
set(MOLSTURM_DEPENDENCIES_TEST "")

##############
#-- Macros --#
##############

macro(add_submodule_dependency LIBRARY VERSION)
	set(${LIBRARY}_SOURCE_DIR "${PROJECT_SOURCE_DIR}/modules/${LIBRARY}")
	set(${LIBRARY}_BINARY_DIR "${PROJECT_BINARY_DIR}/modules/${LIBRARY}")
	if(NOT EXISTS "${${LIBRARY}_SOURCE_DIR}/CMakeLists.txt")
		message(FATAL_ERROR "Did not find expected submodule ${LIBRARY} in ${${LIBRARY}_DIR}. \
Try \"git submodule update --init --recursive\"")
	endif()

	# Extract version from CMakeLists.txt:
	file(STRINGS "${${LIBRARY}_SOURCE_DIR}/CMakeLists.txt" VERSION_RAW
		REGEX "${LIBRARY} VERSION [0-9.]+"
		LIMIT_COUNT 1)
	string(REGEX MATCH "[0-9.]+" ${LIBRARY}_VERSION "${VERSION_RAW}")

	# Compare against what is needed
	if("${${LIBRARY}_VERSION}" VERSION_LESS "${VERSION}")
		message(FATAL_ERROR "Inconsistency in the repo: Version ${VERSION} of ${LIBRARY} \
was requested, but only version ${${LIBRARY}_VERSION} was found. Maybe a \
\"git submodule update --init ${${LIBRARY}_SOURCE_DIR}\" helps?")
	endif()

	# Set the project up:
	add_subdirectory("${${LIBRARY}_SOURCE_DIR}")
	include_directories("${${LIBRARY}_SOURCE_DIR}/src")  # For sources
	include_directories("${${LIBRARY}_BINARY_DIR}/src")  # For generated files

	# Add dependencies:
	foreach(build ${DRB_BUILD_TYPES})
		set(MOLSTURM_DEPENDENCIES_${build} ${MOLSTURM_DEPENDENCIES_${build}} ${${LIBRARY}_${build}_TARGET})
	endforeach()
endmacro(add_submodule_dependency)

###############
#--  Types  --#
###############
# Determine scalar types we want to use here:
include(ScalarTypes)
setup_scalar_types()
set(MOLSTURM_DEPENDENCIES_ ${MOLSTURM_DEPENDENCIES} ${SCALAR_TYPES_LIBRARIES})

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

add_submodule_dependency(krims 0.1.0)
add_submodule_dependency(lazyten 0.3.0)
add_submodule_dependency(gint 0.0.0)
add_submodule_dependency(gscf 0.1.0)
