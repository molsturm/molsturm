#!/bin/bash

KEXP=4.6120
Z=4
ELEC=$((Z/2))
ERROR="1e-2"

#

NMIN=3
NMAX=8

LMIN=2
LMAX=4

BINARY="../build/examples/hf/hf"

# ----------------

for n in $(seq $NMIN $NMAX); do
	lm=$LMAX
	((n<LMAX)) && lm=$((n-1))

	for l in $(seq $LMIN $LMAX); do
		echo
		echo "######################"
		echo "#--  n=$n   l=$l   --#"
		echo "######################"
		$BINARY --Z_charge $Z --alpha $ELEC --beta $ELEC \
			--n_max $n --l_max $l  --basis_type cs_naive \
			--kexp $KEXP --error "$ERROR" \
			--n_eigenpairs 1000 \
			| tee ${Z}_n${n}_l${l}_k${KEXP}.out
	done
done
