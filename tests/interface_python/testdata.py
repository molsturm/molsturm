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
import molsturm.yaml_utils
import os
import re
import yaml

def dir_of_this_script():
  return os.path.dirname( os.path.abspath( __file__ ) )

"""Directory with the reference cases"""
REFERENCE_DIR = os.path.join(dir_of_this_script(), "refdata")

"""Directory with the test cases"""
TESTDATA_DIR = os.path.join(dir_of_this_script(), "testdata")


def __parse_test_case(infile):
  # The possible posthf cases:
  posthf_cases = [ "mp2", "fci" ]

  # Testcase directory is composed of the files
  #   case.in.yaml  -> molsturm hartree_fock input and test metadata
  #   case.hf.yaml  -> hartree fock reference result
  #   case.mp2.yaml -> mp2 reference result
  #   case.fci.yaml -> fci reference result
  # So this next variable gives the part excluding
  # the varying extension:
  basepath = infile[:-len(".in.yaml")]
  ret = { "testing": { "name": os.path.basename(basepath) }, }

  # Load "hf" subtree from hf file
  if os.path.exists(basepath + ".hf.yaml"):
    ret["hf"] = molsturm.load_yaml(basepath + ".hf.yaml")
  elif os.path.exists(basepath + ".hf.hdf5"):
    ret["hf"] = molsturm.load_hdf5(basepath + ".hf.hdf5")
  else:
    raise ValueError("No hf results found for " + infile)

  # Load "params" and "testing" subtree from in file
  with open(basepath + ".in.yaml", "r") as f:
    ret["params"] = yaml.safe_load(f)

  # Remove the "include" subbranch since it is only
  # needed for test data generation
  del ret["params"]["include"]

  # Bring the existing testing metadata one level up
  if "testing" in ret["params"]:
    ret["testing"].update(ret["params"]["testing"])
    del ret["params"]["testing"]

  # Read posthf subcases:
  molsturm.yaml_utils.install_constructors()
  for sub in posthf_cases:
    yamlfile = basepath + "." + sub + ".yaml"
    if os.path.isfile(yamlfile):
      with open(yamlfile) as f:
        ret[sub] = yaml.safe_load(f)

  # Store if we have any posthf case:
  ret["testing"]["any_posthf"] = any( sub in ret for sub in posthf_cases )

  return ret

__reference_cases_cache = None
def reference_cases():
  global __reference_cases_cache

  """Return the list of test cases against reference data"""
  if __reference_cases_cache is None:
    __reference_cases_cache = [
      __parse_test_case(os.path.join(REFERENCE_DIR,fname))
      for fname in os.listdir(REFERENCE_DIR)
      if fname.endswith(".in.yaml")
    ]
  return __reference_cases_cache


__test_cases_cache = None
def test_cases():
  """Return the list of test cases against previously computed data"""
  global __test_cases_cache

  if __test_cases_cache is None:
    __test_cases_cache = [
      __parse_test_case(os.path.join(TESTDATA_DIR,fname))
      for fname in os.listdir(TESTDATA_DIR)
      if fname.endswith(".in.yaml")
    ]
  return __test_cases_cache


def test_cases_by_pred(pred):
  """Return those test cases matching a given predicate"""
  return [ case for case in test_cases() if pred(case) ]


def test_cases_by_name(nameRegex):
  """Return the test case which matches the given name"""

  if isinstance(nameRegex, str):
    nameRegex = re.compile(nameRegex)
  return test_cases_by_pred(lambda x : nameRegex.match(x["testing"]["name"]))
