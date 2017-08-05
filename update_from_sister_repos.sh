#!/bin/bash
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

. update_from_sisters.lib.sh || exit 1

update_file "gint" ".travis/update_cached_libint.sh" || exit 1

update_file "linalgwrap" "templates/py.template" "keep_header" || exit 1
update_file "linalgwrap" "templates/cc.template" "keep_header" || exit 1
update_file "linalgwrap" "templates/cmake.template" "keep_header" || exit 1
update_file "linalgwrap" "templates/hh.template" "keep_header" || exit 1
update_file "linalgwrap" "templates/README.md" || exit 1

update_file "krims" ".clang-format" || exit 1
update_file "krims" ".clang-tidy" || exit 1
update_file "linalgwrap" "update_from_sisters.lib.sh" || exit 1
