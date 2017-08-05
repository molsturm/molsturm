## ---------------------------------------------------------------------
##
## Copyright (C) 2017 by the molsturm authors
##
## This file is part of molsturm.
##
## molsturm is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## molsturm is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with molsturm. If not, see <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------

cmake_minimum_required(VERSION 3.0.0)
cmake_policy(VERSION 3.0.0)

# Try to locate module location of krims submodule:
if (NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/modules/krims/CMakeLists.txt")
	message(FATAL_ERROR "Could not find krims submodule. \
Try \"git submodule update --init --recursive\"")
endif()
list(APPEND CMAKE_MODULE_PATH  "${CMAKE_CURRENT_LIST_DIR}/modules/krims/cmake/modules")

# Load the DebugReleaseBuild module for a simultaneous Debug and Release
# and release build and initialise it.
include(DebugReleaseBuild)
drb_init()

# Prefer shared to static libs (cmake does the reverse by default.
# We to set this option *before* calling the first 'project')
option(BUILD_SHARED_LIBS "Build shared libraries instead of static ones" ON)

##########################################################################
# Project setup

project(molsturm VERSION 0.0.2)
message("\n# Configuring ${PROJECT_NAME}\n#")

# Global include directories for this project:
include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_BINARY_DIR}/src")

# enable testing of this project
enable_testing()

##########################################################################
# Global options

# Path where cmake wants the modules and config packages:
set(PackageModuleLocation share/cmake)

# setup global options
include(GlobalOption)
global_option(ENABLE_DOCUMENTATION "Build documentation"                       OFF )
global_option(ENABLE_EXAMPLES      "Build example executables"                 ON  )
global_option(ENABLE_TESTS         "Build unit test executables"               ON  )
global_option(INTERFACE_PYTHON     "Build and install python interface module" ON  )

##########################################################################
# Setup hard and optional dependencies and find components

# Determine and setup compiler flags, enforcing C++11
drb_setup_compiler_flags(11)

# Export compile commands if cmake supports it.
if (CMAKE_VERSION VERSION_GREATER 3.5.0)
	set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

# If a user wants a python module, setup the interface helper class
# This might alter the compiler flags in order for the compiled
# code to be compatible with the interfacing, hence this call is
# needed right here.
if (MOLSTURM_INTERFACE_PYTHON)
	include(CppInterface)
endif()

# setup hard dependencies
include(cmake/setup_dependencies.cmake)

# setup optional dependencies and features.
include(cmake/setup_optional.cmake)

##########################################################################
# Pass on to subdirectories:

# The src directory does most of the work: Here the library is build and
# installation is configured.
add_subdirectory(src)

# Add subdirectories for the testing.
if(MOLSTURM_ENABLE_TESTS)
	add_subdirectory(tests)
endif()

# Add subdirectories for the examples.
if(MOLSTURM_ENABLE_EXAMPLES)
	add_subdirectory(examples)
endif()

if(MOLSTURM_ENABLE_DOCUMENTATION)
	add_subdirectory(doc)
endif()

##########################################################################
# Finally install package:
include(cmake/install_package.cmake)

# and print a closing message
message("#\n# Configuring ${PROJECT_NAME} done\n")
