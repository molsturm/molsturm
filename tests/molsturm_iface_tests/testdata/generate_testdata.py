#!/usr/bin/env python3
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

import molsturm
import os
import yaml

def build_input_params():
  """Gather a dictionary of all input parameters"""
  dir_of_this_script = os.path.dirname( os.path.abspath( __file__ ) )

  # Build the list of input data
  inputs=dict()
  for f in os.listdir(dir_of_this_script):
    fn,ext = os.path.splitext(f)

    # Only consider input files:
    if ext != ".yaml": continue
    if fn[-3:] != ".in": continue

    key = fn[:-3]
    with open(f,"r") as stream:
      inputs[key] = yaml.safe_load(stream)
  return inputs

# Cache of calculations we already did
calculation_hf_cache=dict()

def run_hf_calculation(name,params):
  if not name in calculation_hf_cache:
    calculation_hf_cache[name] = molsturm.hartree_fock(**params)
  return calculation_hf_cache[name]

# --------------------------------------------------------------------

def job_dump_yaml(name,params):
  """Run a full calculation and dump the result as a yaml file"""
  output = name+".yaml"
  if not os.path.exists(output):
    res = run_hf_calculation(name,params)
    molsturm.dump_yaml(res,output)

# --------------------------------------------------------------------

if __name__ == "__main__":
  inputs = build_input_params()
  for name in inputs:
    params = inputs[name]

    # The include list of jobs which should be performed
    # on this input:
    include = params["include"]
    del params["include"]

    for job in include:
      try:
        locals()["job_" + job](name,params)
      except KeyError:
        raise SystemExit("Unknown job name '" + job + "' in input file '" +
                         name + ".in.yaml'")

