## ---------------------------------------------------------------------
##
## Copyright (C) 2016 by the molsturm authors
##
## This file is part of molsturm.
##
## molsturm is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## molsturm is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
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

##########################################################################
# Global options

# Path where cmake wants the modules and config packages:
set(PackageModuleLocation share/cmake)

# setup global options (needs to be done before call to project())
option(MOLSTURM_ENABLE_TESTS "Build molsturm tests" ON)
option(MOLSTURM_ENABLE_EXAMPLES "Build molsturm examples" ON)

##########################################################################
# Project setup

project(molsturm VERSION 0.0.0)
message("\n# Configuring ${PROJECT_NAME}\n#")

# Global include directories for this project:
include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_BINARY_DIR}/src")

# enable testing of this project
enable_testing()

##########################################################################
# Setup hard and optional dependencies and find components

# Determine and setup compiler flags, enforcing C++11
drb_setup_compiler_flags(11)

# Export compile commands if cmake supports it.
if (CMAKE_VERSION VERSION_GREATER 3.5.0)
	set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
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

##########################################################################
# Finally install package:
include(cmake/install_package.cmake)

# and print a closing message
message("#\n# Configuring ${PROJECT_NAME} done\n")
