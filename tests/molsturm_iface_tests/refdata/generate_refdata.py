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

import itertools
import molsturm
import numpy as np
import os
import re
import subprocess
import tempfile
import yaml

def dir_of_this_script():
  return os.path.dirname( os.path.abspath( __file__ ) )

FLPAT = "[-.0-9eE]+"  # Pattern for floating point numbers
ORCADIR = os.path.join(dir_of_this_script(), "orca")

def build_input_params():
  """Gather a dictionary of all input parameters"""
  # Build the list of input data
  inputs=dict()
  for f in os.listdir(dir_of_this_script()):
    fn,ext = os.path.splitext(f)

    # Only consider input files:
    if ext != ".yaml": continue
    if fn[-3:] != ".in": continue

    key = fn[:-3]
    with open(f,"r") as stream:
      inputs[key] = yaml.safe_load(stream)
  return inputs


def build_orca_input(params, simple_keywords, input_file):
  """
  params: The params array
  simple_keywords: Extra keywords to put on the simple input line
  input_file: The file to write
  """
  os.makedirs(os.path.dirname(input_file), exist_ok=True)

  for key in [ "atoms", "basis_set", "coords", "multiplicity" ]:
    if not key in params:
      raise ValueError("Missing required key '" + key + "' in input parameters.")
  if len(params["atoms"]) != len(params["coords"]):
    raise ValueError("The length of the atoms array and of the coords array "+
                     "does not agree.")

  simple_keywords = ["bohrs", "extremescf", "nofrozencore", "nopop", "noprintmos",
                     "largeprint"] + simple_keywords
  simple_keywords.append(params["basis_set"])

  with open(input_file, "w") as f:
    f.write("! " + " ".join(simple_keywords)+"\n")
    f.write("* xyz " + str(params.get("charge", 0)) + " " +
            str(params["multiplicity"]) + "\n")

    for atom, coords in itertools.zip_longest(params["atoms"], params["coords"]):
      f.write(str(atom))
      for c in coords:
        f.write("  " + str(c))
      f.write("\n")

    f.write("*\n")


def run_orca(input_file):
  """Run an orca calculation and return the output string"""

  basename = os.path.basename(input_file)
  basenoext = os.path.splitext(basename)[0]
  orca_out_file = os.path.join(ORCADIR, basenoext + ".out")

  if os.path.exists(orca_out_file):
    # Calculation has been done before, just return result
    with open(orca_out_file, "r") as out:
      return out.read()
  else:
    # Copy file to calculation folder:
    if not os.path.isfile(input_file):
      raise IOError("File does not exist: " + input_file)

    # Run calculation and dump
    res = subprocess.check_output(["orca", basename], cwd=ORCADIR)
    with open(orca_out_file, "wb") as out:
      out.write(res)

    # Cleanup
    for ext in [ ".gbw", ".prop", "_property.txt" ]:
      tmpfile = os.path.join(ORCADIR, basenoext + ext)
      if os.path.isfile(tmpfile):
        os.remove(tmpfile)

    return res.decode()

def output_find_patterns(output, patterns):
  """
  patterns: Dict of keys -> regexes to search for
  """
  result = dict()
  for patdict in patterns:
    key = patdict["key"]
    match = re.search(patdict["regex"], output)
    if patdict["default"] is None and not match:
      raise Exception("Could not find value for " + key + " in orca output. "+
                      "No default is given either.")
    elif not match:
      result[key] = patdict["default"]
    else:
      result[key] = patdict["convert"](match.group("value"))
  return result

def orca_extract_orben(output):
  # Parse MOs: Find beginning and end of the MO block
  mo_begin = re.search(".*----\nORBITAL ENERGIES\n-----.*", output)
  mo_end = output[mo_begin.end():].find("--")
  mo_text = output[mo_begin.end():] [:mo_end]

  up = mo_text.find("SPIN UP ORBITALS")
  down = mo_text.find("SPIN DOWN ORBITALS")
  if up == -1 and down == -1:
    retricted = True
    # Restricted => Use the full table twice
    loops=[ (0,None), (0,None) ]
  elif up < down and up >= 0:
    restricted = False
    # Restricted => Use the alpha table for the alphas
    # and the beta table for the betas
    loops=[ (up, down), (down, None) ]
  else:
    raise ValueError("output has unrecognisied MO energy table format")

  mos=[]
  re_mo_energy = re.compile("\s*[0-9]+\s*" + FLPAT + "\s*(?P<value>" + FLPAT + ")")
  for begin, end in loops:
    for line in mo_text[begin:end].split("\n"):
      match = re_mo_energy.match(line)
      if match:
        mos.append(float(match.group("value")))
  return np.array(mos)

# --------------------------------------------------------------------

def job_orca_hf(name, params):
  """Run an orca HF calculation and extract reference values from it."""

  orca_in_file = os.path.join(ORCADIR, name + ".orca_hf.in")
  build_orca_input(params, [ "hf" ], orca_in_file)
  res = run_orca(orca_in_file)

  def parse_restricted(value):
    return False if value == "UHF" else True

  orca_patterns = [
    {
      "key":      "energy_ground_state",
      "regex":    "FINAL SINGLE POINT ENERGY\s*(?P<value>" + FLPAT + ")",
      "convert":  float,
      "default":  None,
    }, {
      "key":      "energy_1e",
      "regex":    "One Electron Energy\s*:\s*(?P<value>" + FLPAT + ") Eh",
      "convert":  float,
      "default":  None,
    }, {
      "key":      "energy_2e",
      "regex":    "Two Electron Energy\s*:\s*(?P<value>" + FLPAT + ") Eh",
      "convert":  float,
      "default":  None,
    }, {
      "key":      "energy_kinetic",
      "regex":    "Kinetic Energy\s*:\s*(?P<value>" + FLPAT + ") Eh",
      "convert":  float,
      "default":  None,
    }, {
      "key":      "energy_potential",
      "regex":    "Potential Energy\s*:\s*(?P<value>" + FLPAT + ") Eh",
      "convert":  float,
      "default":  None,
    }, {
      "key":      "energy_nuclear_repulsion",
      "regex":    "Nuclear Repulsion\s*:\s*(?P<value>" + FLPAT + ") Eh",
      "convert":  float,
      "default":  None,
    }, {
      "key":      "restricted",
      "regex":    "Hartree-Fock type\s*HFTyp\s*\.*\s*(?P<value>\w+)",
      "convert":  lambda value: False if value == "UHF" else True,
      "default":  None,
    }, {
      "key":      "spin_squared",
      "regex":    "Expectation value of <S\*\*2>\s*:\s*(?P<value>" + FLPAT + ")",
      "convert":  float,
      "default":  0,
    },
  ]
  reference = output_find_patterns(res, orca_patterns)
  print("TODO: It would be nice if the MO coefficients would be extracted from ORCA, too.")

  # Extra parsing
  reference["orben_f"] = orca_extract_orben(res)

  with open(os.path.join(dir_of_this_script(), name + ".hf.yaml"), "w") as f:
    # TODO It would be really nice if this was an actual hfres yaml
    #      data file and not just a faked one as it is now.
    #      This would require extracting a couple of more keys from
    #      the ORCA output.
    f.write("# Data from ORCA calculation " + orca_in_file + "\n")
    molsturm.dump_yaml(reference, f)


def job_orca_mp2(name, params):
  orca_in_file = os.path.join(ORCADIR, name + ".orca_mp2.in")
  build_orca_input(params, [ "mp2" ], orca_in_file)
  res = run_orca(orca_in_file)

  orca_patterns = [
    {
      "key":      "energy_ground_state",
      "regex":    "FINAL SINGLE POINT ENERGY\s*(?P<value>" + FLPAT + ")",
      "convert":  float,
      "default":  None,
    }, {
      "key":      "energy_mp2",
      "regex":    "MP2 CORRELATION ENERGY\s*:\s*(?P<value>" + FLPAT + ") Eh",
      "convert":  float,
      "default":  None,
    }
  ]
  reference = output_find_patterns(res, orca_patterns)
  reference[molsturm.INPUT_PARAMETER_KEY] = {}

  with open(os.path.join(dir_of_this_script(), name + ".mp2.yaml"), "w") as f:
    f.write("# Data from ORCA calculation " + orca_in_file + "\n")
    yaml.safe_dump(reference, f)

def job_orca_fci(name, params):
  orca_in_file = os.path.join(ORCADIR, name + ".orca_fci.in")
  build_orca_input(params, [ "fci" ], orca_in_file)
  res = run_orca(orca_in_file)

  orca_patterns = [
    {
      "key":      "energy_root0",
      "regex":    "FINAL SINGLE POINT ENERGY Root 0\s*=\s*(?P<value>" + FLPAT + ")",
      "convert":  float,
      "default":  None,
    }
  ]
  parsed = output_find_patterns(res, orca_patterns)

  reference = {
    molsturm.INPUT_PARAMETER_KEY: { "n_roots": 1 },
    "states": [ {
      "energy":        parsed["energy_root0"],
      "civector":      None,  # TODO No clue how to get that
      "spin_squared":  None,  # TODO No clue how to get that
      "multiplicity":  None,  # TODO No clue how to get that
    } ],
  }
  with open(os.path.join(dir_of_this_script(), name + ".fci.yaml"), "w") as f:
    f.write("# Data from ORCA calculation " + orca_in_file + "\n")
    yaml.safe_dump(reference, f)


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
        locals()["job_" + job](name, params)
      except KeyError:
        raise SystemExit("Unknown job name '" + job + "' in input file '" +
                         name + ".in.yaml'")

