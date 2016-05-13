# Installs the cmake apckage information this project provides
#
# Requires the variable PackageModuleLocation to be set.

# Write a basic version file for molsturm
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
	"${molsturm_BINARY_DIR}/molsturmConfigVersion.cmake"
	COMPATIBILITY AnyNewerVersion
)

# Adjust a configure file
configure_file(cmake/molsturmConfig.cmake.in
	"${molsturm_BINARY_DIR}/molsturmConfig.cmake"
	COPYONLY
)

# Set an export location:
install(FILES
	"${molsturm_BINARY_DIR}/molsturmConfig.cmake"
	"${molsturm_BINARY_DIR}/molsturmConfigVersion.cmake"
	DESTINATION "${PackageModuleLocation}/molsturm"
	COMPONENT devel
)

