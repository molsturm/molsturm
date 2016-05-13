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

# macro to find a library from the system
# and set it up.
macro(add_required_upstream_library LIBRARY VERSION)
	find_package(${LIBRARY} ${VERSION} REQUIRED CONFIG)

	# Set the apropriate target names:
	set(${LIBRARY}_DEBUG_TARGET   "${LIBRARY}.g" 
		CACHE INTERNAL "Target name of debug version of ${LIBRARY}")
	set(${LIBRARY}_RELEASE_TARGET "${LIBRARY}"
		CACHE INTERNAL "Target name of release version of ${LIBRARY}")

	message(STATUS "Found ${LIBRARY} config at ${${LIBRARY}_CONFIG}")

	# Check that all required targets are available.
	foreach(build ${DRB_BUILD_TYPES})
		if(NOT TARGET "${${LIBRARY}_${build}_TARGET}")
			message(FATAL_ERROR "We could not find a ${build} version of linalwrap at this location. \
Either disable building a ${build} version of ${CMAKE_PROJECT_NAME} or else \
rebuild linalgwrap with a ${build} version as well.")
		endif()

		# Add dependencies to appropriate versions of molsturm
		set(MOLSTURM_DEPENDENCIES_${build} ${MOLSTURM_DEPENDENCIES_${build}} ${${LIBRARY}_${build}_TARGET})
	endforeach()
endmacro(add_required_upstream_library)

# macro to setup a library from the git submodules in modules/
# TODO Version is silently ignored.
macro(add_required_submodule_library LIBRARY VERSION)
	if(NOT EXISTS "${molsturm_SOURCE_DIR}/${MOLSTURM_SUBMODULES_DIR}/${LIBRARY}/CMakeLists.txt")
		message(FATAL_ERROR "Did not find expected submodule in ${molsturm_SOURCE_DIR}/${MOLSTURM_SUBMODULES_DIR}/${LIBRARY}/. \
Try \"git submodule update --init --recursive\"")
	endif()

	message(STATUS "Found ${LIBRARY} as a submodule in ${molsturm_SOURCE_DIR}/modules/${LIBRARY}/.")

	# Set the project up:
	add_subdirectory("${MOLSTURM_SUBMODULES_DIR}/${LIBRARY}")

	# Now add dependencies:
	foreach(build ${DRB_BUILD_TYPES})
		set(MOLSTURM_DEPENDENCIES_${build} ${MOLSTURM_DEPENDENCIES_${build}} ${${LIBRARY}_${build}_TARGET})
	endforeach()
endmacro(add_required_submodule_library)


###################
#--  libraries  --#
###################

#
# common
#
if(MOLSTURM_WITH_SYSTEM_COMMON)
	add_required_upstream_library(common 0.0.0)
else()
	add_required_submodule_library(common 0.0.0)
endif()

#
# linalgwrap
#
if(MOLSTURM_WITH_SYSTEM_LINALGWRAP)
	add_required_upstream_library(linalgwrap 0.1.0)
else()
	add_required_submodule_library(linalgwrap 0.1.0)
endif()

#
# sturmint
#
if(MOLSTURM_WITH_SYSTEM_STURMINT)
	add_required_upstream_library(sturmint 0.0.0)
else()
	add_required_submodule_library(sturmint 0.0.0)
endif()

#
# gint
#
if(MOLSTURM_WITH_SYSTEM_GINT)
	add_required_upstream_library(gint 0.0.0)
else()
	add_required_submodule_library(gint 0.0.0)
endif()

#
# gscf
#
if(MOLSTURM_WITH_SYSTEM_GSCF)
	add_required_upstream_library(gscf 0.0.0)
else()
	add_required_submodule_library(gscf 0.0.0)
endif()

#
# catch and rapidcheck
#
if(MOLSTURM_ENABLE_TESTS)
	# We need to have a way for the common library to be forced to expose
	# rapidcheck and catch in this way.

	if (NOT TARGET common_catch)
		message(FATAL_ERROR "No target common_catch defined.")
	endif()

	if (NOT TARGET common_rapidcheck)
		message(FATAL_ERROR "No target common_rapidcheck defined.")
	endif()

	# Add them to dependencies:
	set(MOLSTURM_DEPENDENCIES_TEST ${MOLSTURM_DEPENDENCIES_TEST} common_catch common_rapidcheck)
endif()
