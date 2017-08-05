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

include(Documentation)

# Exclude patterns for the actual Doxyfile generated with configure_file
set(DOXY_EXCLUDE_PATTERNS "")
set(DOXY_LATEX "OFF")
if (DOCUMENTATION_DOXYGEN_LATEX)
	set(DOXY_LATEX "ON")
endif()

configure_file("Doxyfile.in" "Doxyfile" @ONLY)
setup_doxygen(doc-molsturm-doxygen "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile"
	      "${CMAKE_CURRENT_BINARY_DIR}/codedoc"
	      "share/doc/molsturm")

add_custom_target(doc-molsturm ALL)
add_dependencies(doc-molsturm doc-molsturm-doxygen)
