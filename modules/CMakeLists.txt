# krims
if (NOT MOLSTURM_WITH_SYSTEM_KRIMS)
	# We need to build krims lib from subdirectory
	add_subdirectory(krims)
endif()

# linalgwrap
if (NOT MOLSTURM_WITH_SYSTEM_LINALGWRAP)
	# We need to build linalgwrap lib from subdirectory
	add_subdirectory(linalgwrap)
endif()

# sturmint
if (NOT MOLSTURM_WITH_SYSTEM_STURMINT)
	# We need to build sturmint lib from subdirectory
	add_subdirectory(sturmint)
endif()

# gint
if (NOT MOLSTURM_WITH_SYSTEM_GINT)
	# We need to build gint lib from subdirectory
	add_subdirectory(gint)
endif()

# gscf
if (NOT MOLSTURM_WITH_SYSTEM_GSCF)
	# We need to build gscf lib from subdirectory
	add_subdirectory(gscf)
endif()
