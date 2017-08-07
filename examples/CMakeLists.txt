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

function(setup_example_target TARGET)
	# usually the examples are built and linked against RELEASE
	# if this version is available
	# Supply PREFER_DEBUG as ARGV1 if you want to have the opposite
	# preference
	#
	# Use this function like
	#    add_executable(test test.cc main.cc)
	#    setup_example_target(test)
	#

	# Select build type for example target:
	if ("${ARGV1}" STREQUAL "PREFER_DEBUG" AND CMAKE_BUILD_TYPE MATCHES "Debug")
		set(EXAMPLE_BUILD_TYPE "DEBUG")
	elseif(CMAKE_BUILD_TYPE MATCHES "Release")
		set(EXAMPLE_BUILD_TYPE "RELEASE")
	else()
		set(EXAMPLE_BUILD_TYPE "DEBUG")
	endif()

	set_target_properties(
		${TARGET}
		PROPERTIES
		#${CMAKE_EXE_LINKER_FLAGS} are added
		LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS_${EXAMPLE_BUILD_TYPE}}"
		LINKER_LANGUAGE "CXX"
		#${CMAKE_CXX_FLAGS} are added automaticaly
		COMPILE_FLAGS "${CMAKE_CXX_FLAGS_${EXAMPLE_BUILD_TYPE}}"
		COMPILE_DEFINITIONS "${EXAMPLE_BUILD_TYPE}"
	)

	# link it with the appropriate molsturm library target
	target_link_libraries(${TARGET} ${molsturm_${EXAMPLE_BUILD_TYPE}_TARGET})
endfunction()

# Include subdirectories:
add_subdirectory(hf)
