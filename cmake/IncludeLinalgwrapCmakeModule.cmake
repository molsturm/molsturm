# Try to find a cmake module which usually gets shipped with linalgwrap.

function(include_linalgwrap_cmake_module MODULE)
	# First try to load it plainly as a module:
	include(${MODULE} OPTIONAL RESULT_VARIABLE RES)

	if ("${RES}" STREQUAL "NOTFOUND")
		# We could not "just" find it.
		# Now try some hints:
		set(ModuleHints 
			# if the user specified the location of the module explictly
			"$ENV{${MODULE}_DIR}/${MODULE}.cmake"
			# if the user provided a hint for linalgwrap:
			"$ENV{linalgwrap_DIR}/share/cmake/modules/${MODULE}.cmake"
		)

		if(NOT MOLSTURM_WITH_SYSTEM_LINALGWRAP)
			# add linalgwrap in modules directory.
			set(ModuleHints 
				${ModuleHints}
				"${CMAKE_SOURCE_DIR}/modules/linalgwrap/cmake/modules/${MODULE}.cmake"
			)
		endif()

		foreach(hfile ${ModuleHints})
			include(${hfile} OPTIONAL RESULT_VARIABLE RES)
			if (NOT "${RES}" STREQUAL "NOTFOUND")
				message(STATUS "Found ${MODULE} file at ${RES}")
				return()
			endif()
		endforeach()

		if (MOLSTURM_WITH_SYSTEM_LINALGWRAP)
			message(FATAL_ERROR "Could not find the ${MODULE} module.
Try specifying the enviroment variable ${MODULE}_DIR for a hint toward the directory \
containing the ${MODULE} module. 
Alternatively you can also specify in the env. variable linalgwrap_DIR the installation prefix\
which was used when installing the linalgwrap binaries (which contains this module).")
		else()
			message(FATAL_ERROR "Could not find the ${MODULE} module.
Did you checkout the git submodules in the modules directory? \
Try \"git submodule update --init --recursive\" to achieve this.")
		endif()
	else()
		message(STATUS "Using system-provided ${MODULE} file.")
	endif()
endfunction(include_linalgwrap_cmake_module)
