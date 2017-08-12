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
	export_hf_results.cc
	hartree_fock.cc
	Parameters.cc
	parse_parameters.cc
)

# Setup molsturm_iface module
set(MODULE_NAME molsturm_iface)
cpp_interface_python3(${MODULE_NAME} LINK_TARGETS molsturm
	SOURCES molsturm_iface.i ${SRC}
)

#
# Installing
#
# The site packages folder of the python version we use
find_package(PythonInterp REQUIRED)
set(PYSITE_PACKAGES "lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages")

# Install the python module files (except _iface, which we override)
install(DIRECTORY molsturm
	COMPONENT python-module
	DESTINATION ${PYSITE_PACKAGES}
	PATTERN molsturm/_iface.py EXCLUDE
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_NAME}.py
	COMPONENT python-module
	RENAME _iface.py
	DESTINATION ${PYSITE_PACKAGES}/molsturm
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/_${MODULE_NAME}.so
	COMPONENT python-module
	DESTINATION ${PYSITE_PACKAGES}/molsturm
)
