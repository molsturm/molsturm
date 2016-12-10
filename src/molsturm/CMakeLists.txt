#
# set the source files.
#
set(MOLSTURM_SOURCES
	scf_guess.cc
	version.cc
	IopScfKeys.cc
)

# Dump the current version as an include file.
include(WriteVersionHeader)
write_version_header("${CMAKE_BINARY_DIR}/src/molsturm/version_defs.hh")

# Add the libraries for Debug and Release builds:
drb_add_library(molsturm
	DBGSUFFIX  ".g"
	FILES ${MOLSTURM_SOURCES}
)

drb_set_target_properties(
	ALL molsturm
	PROPERTIES
	VERSION "${PROJECT_VERSION}"
)

drb_target_compile_definitions(ALL     molsturm PUBLIC ${MOLSTURM_DEFINITIONS})
drb_target_compile_definitions(DEBUG   molsturm PUBLIC ${MOLSTURM_DEFINITIONS_DEBUG})
drb_target_compile_definitions(RELEASE molsturm PUBLIC ${MOLSTURM_DEFINITIONS_RELEASE})

drb_target_link_libraries(ALL     molsturm ${MOLSTURM_DEPENDENCIES})
drb_target_link_libraries(DEBUG   molsturm ${MOLSTURM_DEPENDENCIES_DEBUG})
drb_target_link_libraries(RELEASE molsturm ${MOLSTURM_DEPENDENCIES_RELEASE})

# set the installation property for these targets:
install(TARGETS ${molsturm_TARGETS}
	EXPORT MolsturmTargets
	DESTINATION lib
	COMPONENT library
)

# install the includes:
install(DIRECTORY .
	DESTINATION "include/molsturm"
	COMPONENT devel
	FILES_MATCHING PATTERN "*.hh"
)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/version_defs.hh"
	DESTINATION "include/molsturm"
	COMPONENT devel
)

# Export the target specifications for molsturm.
export(EXPORT MolsturmTargets
	FILE "${molsturm_BINARY_DIR}/molsturmTargets.cmake"
	NAMESPACE Upstream::
)

# and install them appropriately:
install(EXPORT MolsturmTargets
	FILE "molsturmTargets.cmake"
	NAMESPACE Upstream::
	DESTINATION "${PackageModuleLocation}/molsturm"
)

