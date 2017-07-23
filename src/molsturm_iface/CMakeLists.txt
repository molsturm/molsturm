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

set(SRC
	available_features.cc
	export_hf_results.cc
	hartree_fock.cc
	Parameters.cc
	parse_parameters.cc
)

#
# Parse build type
#
set(MOLSTURM_IFACE_BUILD_TYPE "Release" CACHE STRING
"The build type to use for the molsturm python bindings. \
Valid are Debug and Release.")

set(module_BUILD_TYPE "DEBUG")
if (MOLSTURM_IFACE_BUILD_TYPE STREQUAL "Release")
	set(module_BUILD_TYPE "RELEASE")
elseif(MOLSTURM_IFACE_BUILD_TYPE  STREQUAL "Debug")
	set(module_BUILD_TYPE "DEBUG")
else()
	message(FATAL_ERROR "The value of MOLSTURM_IFACE_BUILD_TYPE needs to be either \"Debug\" or \"Release\".")
endif()

if (NOT MOLSTURM_IFACE_BUILD_TYPE STREQUAL GINT_INTERFACE_BUILD_TYPE)
	message(FATAL_ERROR "MOLSTURM_IFACE_BUILD_TYPE and GINT_INTERFACE_BUILD_TYPE need to agree.")
endif()

#
# Find swig and python
#
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

set(Python_ADDITIONAL_VERSIONS 3)
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})

#
# Setup cpp flags and includes for swig
# (Note, that swig makes heavy use of old-style casts)
#
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${module_BUILD_TYPE}}")
enable_if_compiles(CMAKE_CXX_FLAGS -Wno-missing-field-initializers)
enable_if_compiles(CMAKE_CXX_FLAGS -Wno-pedantic)
enable_if_compiles(CMAKE_CXX_FLAGS -Wno-old-style-cast)
enable_if_compiles(CMAKE_CXX_FLAGS -Wno-sign-compare)
enable_if_compiles(CMAKE_CXX_FLAGS -Wno-maybe-uninitialized)
enable_if_compiles(CMAKE_CXX_FLAGS -Wno-deprecated)
enable_if_compiles(CMAKE_CXX_FLAGS -Wno-error)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#
# molsturm_iface swig module
#
set(CMAKE_SWIG_FLAGS "")
set_source_files_properties(molsturm_iface.i PROPERTIES CPLUSPLUS ON)
swig_add_module(molsturm_iface python molsturm_iface.i ${SRC})

swig_link_libraries(molsturm_iface
	${PYTHON_LIBRARIES}
	${molsturm_${module_BUILD_TYPE}_TARGET}
)
target_compile_definitions(${SWIG_MODULE_molsturm_iface_REAL_NAME} PRIVATE
        ${MOLSTURM_DEFINITIONS} ${MOLSTURM_DEFINITIONS_${molsturm_BUILD_TYPE}}
	${module_BUILD_TYPE}=${module_BUILD_TYPE}
)

