#!/bin/sh
## ---------------------------------------------------------------------
##
## Copyright (C) 2017 by the gint authors
##
## This file is part of gint.
##
## gint is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## gint is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with gint. If not, see <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------

# Replaces the libint cache folder by the version build
# in this built.

# ----------------------------------------------------------------------
. ${TRAVIS_BUILD_DIR}/.travis/common.sh

# The subdirectories of the LIBINT_INSTALL_DIR directory into which
# the readily built libint files are placed.
SUBDIRS="include lib share"

# ----------------------------------------------------------------------

# Go to the top build directory
cd "${TRAVIS_BUILD_DIR}"

# Check all folders exist
[ ! -d "$LIBINT_INSTALL_DIR" ] && exit 0
for subdir in $SUBDIRS; do
	[ ! -d "$LIBINT_INSTALL_DIR/$subdir" ] && exit 0
done

echo "Updating libint cache at $LIBINT_CACHE_DIR"
rm -rf "$LIBINT_CACHE_DIR" || exit 1
mkdir "$LIBINT_CACHE_DIR" || exit 1

for subdir in $SUBDIRS; do
	[ ! -d "$LIBINT_INSTALL_DIR/$subdir" ] && continue
	cp -a $LIBINT_INSTALL_DIR/$subdir $LIBINT_CACHE_DIR || exit 1
done

exit 0
