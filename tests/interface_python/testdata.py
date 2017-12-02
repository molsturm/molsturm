#!/usr/bin/env python3
## vi: tabstop=4 shiftwidth=4 softtabstop=4 expandtab
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

import molsturm
import molsturm.yaml_utils
import os
import re
import yaml
import directories


def __parse_test_case(subdir, infile):
    """
    subdir:  The subdirectory in which the infile was found
             (either "testdata" or "refdata")
    infile:  The input filename (path relative to the directory above)
    """
    # The possible posthf cases:
    posthf_cases = ["mp2", "fci"]

    # A test case is split up over these files:
    #   case.in.yaml  -> molsturm hartree_fock input and test metadata
    #   case.hf.yaml  -> hartree fock reference result
    #   case.hf.hdf5  -> hartree fock reference result
    #   case.mp2.yaml -> mp2 reference result
    #   case.fci.yaml -> fci reference result
    # Since some of these files might get rather large,
    # they might not exist inside the actual repository,
    # but in fact reside on a webserver from which they are downloaded
    # into the binary directory by CMake.
    # In either case the part first part of the file name gives
    # the name of the test case:
    name = infile[:-len(".in.yaml")]

    basepath = os.path.join(subdir, name)
    ret = {"testing": {"name": name}, }

    # Find and load "hf" subtree:
    for ext in ["yaml", "hdf5"]:
        filename = basepath + ".hf." + ext
        if directories.exists(filename):
            ret["hf"] = molsturm.load_state(directories.find(filename))
            break
    if "hf" not in ret:
        raise ValueError("No hf results found for " + infile)

    # Load "input_parameters" and "testing" subtree from in file
    infile_full = directories.find(basepath + ".in.yaml")
    with open(infile_full, "r") as f:
        params = yaml.safe_load(f)
    ret["input_parameters"] = params["input_parameters"]

    # Forward testing metadata and ScfParameters
    if "testing" in params:
        ret["testing"].update(params["testing"])

    # Read posthf subcases:
    molsturm.yaml_utils.install_constructors()
    for sub in posthf_cases:
        yamlfile = basepath + "." + sub + ".yaml"
        if directories.exists(yamlfile):
            with open(directories.find(yamlfile)) as f:
                ret[sub] = yaml.safe_load(f)

    return ret


# Caches
__reference_cases_cache = None
__test_cases_cache = None


def reference_cases():
    global __reference_cases_cache

    """Return the list of test cases against reference data"""
    if __reference_cases_cache is None:
        testnames = [os.path.basename(f)
                     for f in directories.iglob(os.path.join("refdata", "*.in.yaml"))]
        __reference_cases_cache = [
            __parse_test_case("refdata", name) for name in sorted(set(testnames))
        ]
    return __reference_cases_cache


def test_cases():
    """Return the list of test cases against previously computed data"""
    global __test_cases_cache

    if __test_cases_cache is None:
        testnames = [os.path.basename(f)
                     for f in directories.iglob(os.path.join("testdata", "*.in.yaml"))]
        __test_cases_cache = [
            __parse_test_case("testdata", name) for name in sorted(set(testnames))
        ]
    return __test_cases_cache


class predicates:
    @staticmethod
    def is_expensive():
        return lambda test: test["testing"]["expensive"]

    @staticmethod
    def is_not_expensive():
        return lambda test: not test["testing"]["expensive"]

    @staticmethod
    def all():
        return lambda test: True

    @staticmethod
    def by_name(nameRegex):
        if isinstance(nameRegex, str):
            nameRegex = re.compile(nameRegex)
        return lambda test: nameRegex.match(test["testing"]["name"])

    @staticmethod
    def tests_method(method):
        """Should a particular method be tested for this test case"""
        return lambda test: method in test


def test_cases_by_pred(pred):
    """Return those test cases matching a given predicate"""
    return [case for case in test_cases() if pred(case)]


def test_cases_by_name(nameRegex):
    """Return the test case which matches the given name"""
    return test_cases_by_pred(predicates.by_name(nameRegex))
