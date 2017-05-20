# Try to find adc:
FOUND=1
for file in $PWD/../../../adcc/build/src/pyadc/pyadc_iface.py $PWD/../../../adcc/src/pyadc/pyadc.py; do
	if [ -f "$file" ]; then
		dir=$(realpath $(dirname $file))
		export PYTHONPATH="$PYTHONPATH:$dir"
	else
		FOUND=0
		break
	fi
done
[ $FOUND == 1 ] && echo "Found pyadc"

return 0
