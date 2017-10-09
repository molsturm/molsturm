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

from ._constants import HFRES_ARRAY_KEYS, INPUT_PARAMETER_KEY
from datetime import datetime
from distutils.version import StrictVersion
from ._hdf5 import h5py
from ._iface import Version
from . import _hdf5 as hdf5
from io import IOBase
import numpy as np
import yaml


def metadata_common():
    """
    Return metadata about the file format and the
    version of molsturm.
    """
    return {
        "molsturm_version": Version.as_string(),
        "date": datetime.now().isoformat(),
    }


#
# Yaml
#
def dump_yaml(hfres, stream):
    """Take a HartreeFock result and dump the data
       in yaml format.

       It can later be picked up using load_yaml.

       Note that files are stored fully in plain text,
       so they can get quite large if you choose to
       dump the two electron integral tensor as well.

       This is mainly intended for debugging and
       to 'look' at the data at hand. For archiving
       and transferring data dump_hdf5 should be preferred.

       stream:   Path or file stream to write the data to
    """
    if isinstance(stream, str):
        with open(stream, "w") as f:
            return dump_yaml(hfres, f)
    elif not isinstance(stream, IOBase):
        raise TypeError("stream parameter needs to be a string or a stream object")

    # Make a copy and populate with metadata:
    res = dict(hfres)
    res["meta"] = metadata_common()
    res["meta"]["format_type"] = "yaml"
    res["meta"]["format_version"] = "0.1.0"

    # Convert numpy arrays to plain list of lists
    for k in HFRES_ARRAY_KEYS:
        if k in res:
            res[k] = res[k].tolist()

    yaml.safe_dump(res, stream)


def load_yaml(stream):
    """Read an input file or stream in the molsturm yaml format
       and return the parsed hfres dictionary.

       Note that ``dump_yaml(dict,file); dict = load_yaml(file)``
       is an identity.

       stream:   Path or file stream

       throws: ValueError if the file is not in the expected format.
    """
    # Minimal and maximal version this parser understands
    parser_min_version = "0.1.0"
    parser_max_version = parser_min_version

    if isinstance(stream, str):
        with open(stream, "r") as f:
            return load_yaml(f)
    elif not isinstance(stream, IOBase):
        raise TypeError("stream parameter needs to be a string or a stream object")

    res = yaml.load(stream)

    # Check metadata:
    try:
        meta_type = res["meta"]["format_type"]
        meta_version = res["meta"]["format_version"]
    except KeyError as e:
        raise ValueError("Yaml stream metadata seems to be missing or corrupted. "
                         "The key " + str(e.args[0]) + " is missing.")

    if meta_type != "yaml":
        raise ValueError("Yaml stream format is not 'yaml', but '" +
                         str(meta_type) + "'.")

    if StrictVersion(meta_version) < StrictVersion(parser_min_version):
        raise ValueError("Parser not compatible to versions below " +
                         parser_min_version + ". Encountered version " +
                         meta_version + ".")
    elif StrictVersion(meta_version) > StrictVersion(parser_max_version):
        raise ValueError("Parser not compatible to versions beyond " +
                         parser_max_version + ". Encountered version " +
                         meta_version + ".")

    # Remove metadata block
    del res["meta"]

    # Convert plain list of lists to numpy arrays
    for k in HFRES_ARRAY_KEYS:
        if k in res:
            res[k] = np.array(res[k])
    return res


def metadata_yaml(stream):
    """Extract meta data block of a molsturm yaml file and return as a dictionary
       Returns None no meta data block found
    """
    if isinstance(stream, str):
        with open(stream, "r") as f:
            return metadata_yaml(f)
    elif not isinstance(stream, IOBase):
        raise TypeError("stream parameter needs to be a string or a stream object")

    res = yaml.load(stream)
    try:
        return res["meta"]
    except KeyError:
        return None


#
# HDF5
#
def dump_hdf5(hfres, path):
    """Take a HartreeFock result and dump the data in hdf5 format.
       Most of the data will be plain text, but the large numpy
       arrays are stored in binary form.

       It can later be picked up using load_hdf5.

       path: File path (Note: *not* a stream)
    """
    if not isinstance(path, str):
        raise TypeError("path needs to be a valid path in form of a string")

    with h5py.File(path, 'w') as h5f:
        # Split into the dict corresponding to array values (which are gzipped)
        # and the rest
        array_dict = {k: hfres[k] for k in hfres if k in HFRES_ARRAY_KEYS}
        rest_dict = {k: hfres[k] for k in hfres if k not in HFRES_ARRAY_KEYS}

        hdf5.emplace_dict(array_dict, h5f, compression="gzip")
        hdf5.emplace_dict(rest_dict,  h5f)

        # Write meta data as attributes to root group
        meta = metadata_common()
        meta["format_type"] = "hdf5"
        meta["format_version"] = "0.1.0"

        for attr in meta:
            h5f.attrs[attr] = meta[attr]


def load_hdf5(path, meta_check=True):
    """Read an input file or stream in the molsturm hdf5 format
       and return the parsed hfres dictionary.

       Note that ``dump_hdf5(dict,file); dict = load_hdf5(file)``
       is an identity.

       path:         File path (Note: *not* a stream)
       meta_check:   If False skips the check of the metadata
    """
    # Minimal and maximal version this parser understands
    parser_min_version = "0.1.0"
    parser_max_version = parser_min_version

    def check_metadata(h5f):
        try:
            meta_type = h5f.attrs["format_type"]
            meta_version = h5f.attrs["format_version"]
        except KeyError:
            raise ValueError("Could not find required metadata from HDF5 file for "
                             "metadata check. You can disable checking the metadata "
                             "using 'meta_check=False', but strange things might "
                             "happen ...")

        if meta_type != "hdf5":
            raise ValueError("Hdf5 file format is not 'hdf5', but '" +
                             str(meta_type) + "'.")

        if StrictVersion(meta_version) < StrictVersion(parser_min_version):
            raise ValueError("Parser not compatible to versions below " +
                             parser_min_version + ". Encountered version " +
                             meta_version + ".")
        elif StrictVersion(meta_version) > StrictVersion(parser_max_version):
            raise ValueError("Parser not compatible to versions beyond " +
                             parser_max_version + ". Encountered version " +
                             meta_version + ".")

    with h5py.File(path, "r") as h5f:
        if meta_check:
            check_metadata(h5f)
        return hdf5.extract_group(h5f)


def metadata_hdf5(path):
    """Extract meta data block of an hdf5 file and return as a dictionary
       Returns None no meta data block found
    """
    with h5py.File(path, "r") as h5f:
        return {at: h5f.attrs[at] for at in h5f.attrs}
