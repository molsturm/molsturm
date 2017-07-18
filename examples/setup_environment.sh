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

# The relative path to the build path
BUILD_PATH="../build"

# ------------------------------------------------

for dir in $PWD/../src/molsturm_iface $PWD/$BUILD_PATH/src/molsturm_iface; do
	dir=$(readlink -f $dir)
	if [ ! -d "$dir" ]; then
		echo "Not a valid path:  $dir" >&2
		echo "Probably the examples won't work" >&2
	fi
	export PYTHONPATH="$PYTHONPATH:$dir"
done

# ------------------------------------------------

for file in $PWD/../../adcc/build/src/adcc_iface/adcc_iface.py $PWD/../../adcc/src/adcc_iface/adcc_iface.i; do
	if [ -f "$file" ]; then
		dir=$(readlink -f $(dirname $file))
		export PYTHONPATH="$PYTHONPATH:$dir"
	else
		break
	fi
done

for file in $PWD/../../pyscf/__init__.py; do
	if [ -f "$file" ]; then
		dir=$(readlink -f $(dirname $file))
		export PYTHONPATH="$PYTHONPATH:$dir"
	else
		break
	fi
done

export PYTHONPATH="/usr/lib/python3/dist-packages:$PYTHONPATH"

python3 <<- EOF
from molsturm.posthf import available_methods
print("Available post-HF methods:")
for m in available_methods:
  print("  ",m)
EOF
