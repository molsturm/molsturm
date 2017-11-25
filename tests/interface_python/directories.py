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
import os
import itertools
import glob


def dir_of_this_script():
    return os.path.dirname(os.path.abspath(__file__))


def find_tests_source_dir():
    """Find the 'tests' directory in the source tree of the repository"""
    # If appropriate environment variable is set by ctest, use it
    if "TESTS_SOURCE_DIR" in os.environ:
        trial_dir = os.environ["TESTS_SOURCE_DIR"]
    else:
        # Try to guess from the location of this file
        trial_dir = os.path.join(dir_of_this_script(), "..")

    # Check that we are correct
    for checkfile in ["CMakeLists.txt", "DownloadTestData.cmake"]:
        if not os.path.isfile(os.path.join(trial_dir, checkfile)):
            raise RuntimeError("Test source dir not found, since checkfile '" +
                               checkfile + "' not at the expected location '" +
                               trial_dir + "'")
    return os.path.abspath(trial_dir)


def find_tests_binary_dir():
    """Find the 'tests' directory in the binary tree of the repository"""
    # If appropriate environment variable is set by ctest, use it
    if "TESTS_BINARY_DIR" in os.environ:
        trial_dir = os.environ["TESTS_BINARY_DIR"]
    else:
        # Try to guess from the location of this file
        trial_dir = os.path.join(dir_of_this_script(), "..", "..", "build", "tests")

    for checkfile in ["CTestTestfile.cmake", "molsturm-download-testdata_config.cmake"]:
        if not os.path.isfile(os.path.join(trial_dir, checkfile)):
            raise RuntimeError("Test binary dir not found, since checkfile '" +
                               checkfile + "' not at the expected location '" +
                               trial_dir + "'")
    return os.path.abspath(trial_dir)


"""Directory where the test source files are located (absolute path)"""
TESTS_SOURCE_DIR = find_tests_source_dir()

"""Directory where the test binary files are located (absolute path)"""
TESTS_BINARY_DIR = find_tests_binary_dir()

"""Directory where the python test source files are located (absolute path)"""
PYTHON_TESTS_SOURCE_DIR = dir_of_this_script()

"""Directory where the python test binary files are located (absolute path)"""
PYTHON_TESTS_BINARY_DIR = os.path.join(TESTS_BINARY_DIR,
                                       os.path.relpath(dir_of_this_script(),
                                                       TESTS_SOURCE_DIR))

"""Helper dict to transform the keywords to the actual paths."""
__keyword_to_path = {"source": PYTHON_TESTS_SOURCE_DIR, "binary": PYTHON_TESTS_BINARY_DIR}


def iglob(globexpr, trees=["source", "binary"]):
    """Glob inside the directory trees associated to the python tests
    for a particular file. By default first the source and then the binary
    tree are globbed. This behaviour can be changed by the trees
    variable, which may take 'source' and 'binary' in any arbitrary order"""
    try:
        iglobs = [glob.iglob(os.path.join(__keyword_to_path[t], globexpr))
                  for t in trees]
    except KeyError as e:
        raise ValueError("Unknown tree specifier: " + str(e))
    return itertools.chain(*iglobs)


def find(path, trees=["source", "binary"]):
    """
    Try to find exactly the path specified below the source
    or the binary tree of the tests.

    If the file cannot be found a FileNotFoundError is raised.
    """
    for t in trees:
        if os.path.exists(os.path.join(__keyword_to_path[t], path)):
            return os.path.join(__keyword_to_path[t], path)
    raise FileNotFoundError(path + " does not exist inside the "
                            "source or binary tree")


def exists(path, trees=["source", "binary"]):
    """
    Try to find exactly the path specified below the source
    or the binary tree of the tests. If it does not exist,
    return False, else return True.
    """
    for t in trees:
        if os.path.exists(os.path.join(__keyword_to_path[t], path)):
            return True
    return False
