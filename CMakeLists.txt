cmake_minimum_required(VERSION 3.0.0)

# The DebugReleaseBuild module provides a nice infrastructure to 
# simultaneously build a Debug and a Release version of the same program
# or library, exactly what we need. Try to load it:
include("${CMAKE_SOURCE_DIR}/cmake/IncludeLinalgwrapCmakeModule.cmake")
include_linalgwrap_cmake_module(DebugReleaseBuild)

# initialise DebugReleaseBuild:
drb_init()

##########################################################################
# Global options

# Path where cmake wants the modules and config packages:
set(PackageModuleLocation share/cmake)

# setup global options (needs to be done before call to project())
option(MOLSTURM_ENABLE_TESTS "Build molsturm tests" ON)
option(MOLSTURM_ENABLE_EXAMPLES "Build molsturm examples" ON)

option(MOLSTURM_WITH_SYSTEM_LINALGWRAP
	"Compile and install the linalgwrap library
from the git submodule (OFF) or use system-provided version (ON)" OFF)
option(MOLSTURM_WITH_SYSTEM_KRIMS
	"Compile and install the krims library from the git
submodule (OFF) or use system-provided (ON)" OFF)
option(MOLSTURM_WITH_SYSTEM_GINT
	"Compile and install the gint library from the git
submodule (OFF) or use system-provided (ON)" OFF)
option(MOLSTURM_WITH_SYSTEM_GSCF
	"Compile and install the gscf library from the git
submodule (OFF) or use system-provided (ON)" OFF)
option(MOLSTURM_WITH_SYSTEM_STURMINT
	"Compile and install the sturmint library from the git
submodule (OFF) or use system-provided (ON)" OFF)

##########################################################################
# Project setup

# enable the new version feature introduced in cmake 3.0.0
cmake_policy(SET CMP0048 NEW)
project(molsturm VERSION 0.0.0)

# Global include directories for this project:
include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_BINARY_DIR}/src")

# enable testing of this project
enable_testing()

##########################################################################
# Setup hard and optional dependencies and find components

# Determine and setup compiler flags
drb_setup_compiler_flags()

# Set modules directory:
set(MOLSTURM_SUBMODULES_DIR "modules")

# setup hard dependencies
include(cmake/setup_dependencies.cmake)

# setup optional dependencies and features.
include(cmake/setup_optional.cmake)

##########################################################################
# Pass on to subdirectories:

# The src directory does most of the work: Here the library is build and 
# installation is configured.
add_subdirectory(src)

# Add subdirectories for the test.
if(MOLSTURM_ENABLE_TESTS) 
	add_subdirectory(tests)
endif()

# Add subdirectories for the test.
if(MOLSTURM_ENABLE_EXAMPLES) 
	add_subdirectory(examples)
endif()

##########################################################################
# Finally install package:
include(cmake/install_package.cmake)
