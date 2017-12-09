#!/bin/bash

if [ ! -d modules ] || [ ! -f CMakeLists.txt ]; then
	echo "Run from top dir of molsturm." >&2
fi

SKIP_DIRS=(build rapidcheck external deprecated)
FINDARGS=$(for d in ${SKIP_DIRS[@]}; do echo -n " ! -ipath \"**/$d\" ! -ipath \"**/$d/**\""; done)
DIRS=$(eval "find -type d $FINDARGS \( -iname src -o -iname doc -o -iname tests -o -iname examples -o -iname .travis \)")
cloc $DIRS
