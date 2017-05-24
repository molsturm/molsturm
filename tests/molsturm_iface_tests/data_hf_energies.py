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
## vi: tabstop=2 shiftwidth=2 softtabstop=2 expandtab

# Default parameters, which apply for all runs
params_all = {
  "basis_type":     "gaussian/libint",
  "n_eigenpairs":   100,
  "error":          1e-9,
  "eigensolver":    "lapack",
  "guess_method":   "hcore",
  "guess_esolver":  "lapack",
  "export_repulsion_integrals_max_norbs": 40,
}

# ---------------------------------
# The molecular systems we look at below
system_he = {
  "atom_numbers":  [ 2 ],
  "coords":        [ [0,0,0] ],
}

system_h2o = {
  "atoms":  [ "O", "H", "H" ],
  "coords":
  [
    [ 0,0,0 ],
    [ 0, 0, 1.795239827225189 ],
    [ 1.693194615993441, 0, -0.599043184453037 ],
  ],
}

system_ethane = {
  "atoms":  2* [ "C" ] + 6* [ "H" ],
  "coords":
  [
    [0.11773122601810967, 2.706579152570035e-07, -0.36282169701109324],
    [0.8516848710413312, -3.5297816510288056e-07, 2.420736070352015],
    [2.9018095114074667, -2.3439803917193614e-06, 2.656181084227907],
    [0.11360588178037678, -1.6647921043201028, 3.3913802657097327],
    [0.1136091723264787, 1.6647931911034923, 3.39137973979706],
    [-1.9323875392990797, 3.854368570695653e-06, -0.5983078420216195],
    [0.8558352286352563, 1.66479071436795, -1.3334507521041237],
    [0.8558296748964446, -1.664793229221158, -1.3334498295870683],
  ],
}

# ---------------------------------
# The computational cases covered

cases = [
  {
    "description": "He_aug-cc-pvdz",
    "system":   system_he,
    "params": {
      "basis_set": "aug-cc-pvdz",
      "max_iter":  15,
    },
    "reference": {  # ORCA
      "energy_total":  -2.855704667706,
      "energy_1e":     -3.87716137,
      "energy_2e":      1.02145670,
      "energy_nucrep":  0,
      "energy_mos":
                  2* [ -0.917124,  0.174366,  0.530376,
                       0.530376,  0.530376,  1.713453,
                       3.024883,  3.024883,  3.024883, ],
    },
  },
# ---------------------------------------------------------------
  {
    "description": "H2O_sto-3g",
    "system":   system_h2o,
    "params": {
      "basis_set": "sto-3g",
      "max_iter":  15,
    },
    "reference": {  # ORCA
      "energy_total": -74.959319286910,
      "energy_1e":    -122.50621280,
      "energy_2e":      38.29541424,
      "energy_nucrep":  9.25147927,
      "energy_mos":
                2* [ -20.233397, -1.265715, -0.629267, -0.441668,
                     -0.387645,  0.602839,  0.765918 ],
    },
  },
# ---------------------------------------------------------------
  {
    "description": "ethane_3-21g.out",
    "system":   system_ethane,
    "params": {
      "basis_set": "3-21g",
      "max_iter":  20,
      "error":     1e-7,
    },
    "reference": {  # ORCA
      "energy_total":  -78.793532514352,
      "energy_1e":     -188.47879917,
      "energy_2e":     67.39701724,
      "energy_nucrep": 42.28824941,
      "energy_mos":
                  2*[ -11.153926, -11.153668,  -1.021042,
                       -0.842799,  -0.598213,  -0.598213,
                       -0.510102,  -0.482209,  -0.482209,
                        0.278467,   0.335282,   0.335282,
                        0.340851,   0.370603,   0.370604,
                        0.432701,   0.909415,   0.909416,
                        0.966257,   1.084323,   1.092219,
                        1.092221,   1.287797,   1.360790,
                        1.360791,   1.363882,   1.363882,
                        1.443343,   1.935394,   2.297531, ],
    },
  },
]
