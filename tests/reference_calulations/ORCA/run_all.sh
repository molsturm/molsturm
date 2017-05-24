#!/bin/sh

for i in *.in; do
	OUT=$(basename "$i" .in).out
	[ -f "$OUT" ] && continue
	orca "$i" > "$OUT" 2>&1
done
rm -f *.gbw *.prop *_property.txt *.ges
