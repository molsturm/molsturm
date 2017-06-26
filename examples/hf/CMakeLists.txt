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

set(SRC main.cc parse_args.cc)

# TODO Only compile release version of hf:
# add_executable(hf parse_args.cc ${SRC})
# setup_example_target(hf)

# XXX: For now compile hf with DRB
set(TARGET hf)
drb_add_executable(${TARGET} 
	DBGSUFFIX  "_DEBUG"
	FILES ${SRC}
)

drb_target_compile_definitions(ALL     ${TARGET}  PRIVATE ${MOLSTURM_DEFINITIONS})
drb_target_compile_definitions(DEBUG   ${TARGET}  PRIVATE ${MOLSTURM_DEFINITIONS_DEBUG})
drb_target_compile_definitions(RELEASE ${TARGET}  PRIVATE ${MOLSTURM_DEFINITIONS_RELEASE})

drb_target_link_libraries(ALL     ${TARGET}  ${MOLSTURM_DEPENDENCIES})
drb_target_link_libraries(DEBUG   ${TARGET}  ${molsturm_DEBUG_TARGET} ${MOLSTURM_DEPENDENCIES_DEBUG})
drb_target_link_libraries(RELEASE ${TARGET}  ${molsturm_RELEASE_TARGET} ${MOLSTURM_DEPENDENCIES_RELEASE})
