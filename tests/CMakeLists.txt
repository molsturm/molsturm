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

#
# The sources for the C++ test executable
#
set(MOLSTURM_TESTS_SOURCES
	FunctionalityTest.cc
	HfEnergyTest.cc
	main.cc
)

drb_add_executable(molsturm_tests
	DBGSUFFIX  "_DEBUG"
	RELSUFFIX  "_RELEASE"
	FILES ${MOLSTURM_TESTS_SOURCES}
)

drb_target_compile_definitions(ALL     molsturm_tests
	PRIVATE ${MOLSTURM_DEFINITIONS})
drb_target_compile_definitions(DEBUG   molsturm_tests
	PRIVATE ${MOLSTURM_DEFINITIONS_DEBUG})
drb_target_compile_definitions(RELEASE molsturm_tests
	PRIVATE ${MOLSTURM_DEFINITIONS_RELEASE})

drb_target_link_libraries(ALL     molsturm_tests
	${MOLSTURM_DEPENDENCIES} ${MOLSTURM_DEPENDENCIES_TEST})
drb_target_link_libraries(DEBUG   molsturm_tests
	${molsturm_DEBUG_TARGET} ${MOLSTURM_DEPENDENCIES_DEBUG})
drb_target_link_libraries(RELEASE molsturm_tests
	${molsturm_RELEASE_TARGET} ${MOLSTURM_DEPENDENCIES_RELEASE})

FOREACH(target ${molsturm_tests_TARGETS})
	add_test(NAME ${target} COMMAND ${target})
	set_tests_properties(${target} PROPERTIES TIMEOUT 300)
ENDFOREACH()

#
# Add the python test script
#
set(Python_ADDITIONAL_VERSIONS 3)
find_package(PythonInterp REQUIRED)
add_test(NAME molsturm_iface_tests
	COMMAND ${PYTHON_EXECUTABLE} -m unittest discover
-s "${CMAKE_CURRENT_LIST_DIR}/molsturm_iface_tests"
)
set_tests_properties(molsturm_iface_tests
	PROPERTIES
	ENVIRONMENT "PYTHONPATH=${molsturm_SOURCE_DIR}/src/molsturm_iface:\
${molsturm_BINARY_DIR}/src/molsturm_iface"
	TIMEOUT 300
)
