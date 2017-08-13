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

set(MOLSTURM_SOURCES
	scf_guess.cc
	Version.cc
	IopScfKeys.cc
	read_xyz.cc
	MolecularSystem.cc
	GuessAlgorithms/guess_hcore.cc
	GuessAlgorithms/GuessAlgorithmsKeysBase.cc
)

# Write the config file config.hh (into the binary dir)
configure_file("config.hh.in" "config.hh")

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
install(FILES "${molsturm_BINARY_DIR}/src/molsturm/config.hh"
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
	DESTINATION "share/cmake/molsturm"
)

