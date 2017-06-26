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

## DO NOT EDIT
## This file is automatically generated from a file in the repository "linalgwrap".
## Edit the original and call the script "update_from_sister_repos.sh" instead.

# A script file which contains a script licence header
SCRIPT_FILE="CMakeLists.txt"

if [ ! -r "$SCRIPT_FILE" ]; then
	echo "SCRIPT_FILE ($SCRIPT_FILE) not readable. Are you running" \
		"the script from the top level of the repo?" >&2
	return 1
fi

# extract licence header from SCRIPT_FILE
script_header() {
	# The header is delimited by "-------" characters
	awk '
		$1 == "##" && $2 ~ /^-+$/ {
			# invert printing flag and print line
			pr = !pr
			print
			next
		}

		# print the line
		pr
	' "$SCRIPT_FILE"
}

# Replace the licence header on the stream
replace_header() {
	local REPO="$1"   # The original repo the file is from
	awk -v "newheader=`script_header`" -v "repo=$REPO" '
		!headerprinted && noprint{
			# print header exactly once if inside
			# the noprint part of the file (where the
			# header to replace is located)
			headerprinted=1
			print newheader
			print "\n## DO NOT EDIT"
			print "## This file is automatically generated from " \
				"a file in the repository \"" repo "\"."
			print "## Edit the original and call the script " \
				"\"update_from_sister_repos.sh\" instead."
		}

		$1 == "##" && $2 ~ /^-+$/ {
			# Invert printing flag => no printing
			# if inside header
			noprint = !noprint
			next
		}

		# print the line if required
		!noprint
	'
}

# Try to find another repository
find_repo() {
	local REPO="$1"
	if [ -e "../$REPO/.git" ]; then
		echo "../$REPO"
		return 0
	fi
	if [ -e "external/$REPO/.git" ]; then
		echo "external/$REPO"
		return 0
	fi
	if [ -e "modules/$REPO/.git" ]; then
		echo "modules/$REPO"
		return 0
	fi
	return 1
}

# Update a file if necessary and replace its licence header
update_file() {
	local ORIGREPO=$1       # repo to get the file from
	local FILE=$2           # File to update (same location in both repos)
	local REPLACEHEADER=$3  # Replace the licence header or not

	# Try to find the other repository to get the file from
	local REPODIR
	if ! REPODIR=`find_repo "$ORIGREPO"`; then
		echo "Could not find repo $ORIGREPO" >&2
		return 1
	fi

	# Path to get the file from:
	local FROMPATH="$REPODIR/$FILE"
	if [ ! -f "$FROMPATH" ]; then
		echo "FILE file ($FILE) does not exist in source repository $REPODIR." >&2
		echo "Are you running the script from the top level directory of this repo?" >&2
		return 1
	fi

	# relative dir in which file will be placed
	local TODIR=`dirname "$FILE"`
	if [ ! -d "$TODIR" ]; then
		echo "TODIR ($TODIR) is not a valid directory." >&2
		echo "Are you running the script from the top level directory of this repo?" >&2
		return 1
	fi

	if [ ! -f "$FILE" -o "$FROMPATH" -nt "$FILE" ]; then
		echo "Updating $FILE from $FROMPATH"
		if [ "$REPLACEHEADER" == "keep_header" ]; then
			cp "$FROMPATH" "$FILE" || return 1
		else
			< "$FROMPATH" replace_header "$ORIGREPO"> "$FILE" || return 1
		fi
		touch --reference="$FROMPATH" "$FILE"
		chmod --reference="$FROMPATH" "$FILE"
	fi
	return 0
}
