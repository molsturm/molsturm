# The relative path to the build path
BUILD_PATH="../build"

# ------------------------------------------------

for dir in $PWD/../src/molsturm_iface $PWD/$BUILD_PATH/src/molsturm_iface; do
	dir=$(realpath $dir)
	if [ ! -d "$dir" ]; then
		echo "Not a valid path:  $dir" >&2
		echo "Probably the examples won't work" >&2
	fi
	export PYTHONPATH="$PYTHONPATH:$dir"
done
