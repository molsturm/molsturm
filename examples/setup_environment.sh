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

python3 <<- EOF
from molsturm.posthf import available_methods
print("Available post-HF methods:")
for m in available_methods:
  print("  ",m)
EOF
