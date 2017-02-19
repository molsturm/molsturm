#!/bin/bash
## ---------------------------------------------------------------------
##
## Copyright (C) 2016 by the gscf authors
##
## This file is part of molsturm.
##
## gscf is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## gscf is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with gscf. If not, see <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------

. update_from_sisters.lib.sh || exit 1


update_file "linalgwrap" "templates/cc.template" "keep_header" || exit 1
update_file "linalgwrap" "templates/cmake.template" "keep_header" || exit 1
update_file "linalgwrap" "templates/hh.template" "keep_header" || exit 1
update_file "linalgwrap" "templates/README.md" || exit 1

update_file "krims" ".clang-format" || exit 1
update_file "linalgwrap" "update_from_sisters.lib.sh" || exit 1
