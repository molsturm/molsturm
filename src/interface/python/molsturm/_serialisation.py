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

from ._constants import HFRES_ARRAY_KEYS
from datetime import datetime
from distutils.version import StrictVersion
from .State import State
from ._hdf5 import h5py
from ._iface import Version
from . import _hdf5 as hdf5
from io import IOBase
import numpy as np
import os
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
def _dump_yaml(hfres, stream):
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
            return _dump_yaml(hfres, f)
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


def _load_yaml(stream, meta_check=True):
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
            return _load_yaml(f)
    elif not isinstance(stream, IOBase):
        raise TypeError("stream parameter needs to be a string or a stream object")

    res = yaml.load(stream)

    # Check metadata:
    if meta_check:
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
    res.pop("meta", None)

    # Convert plain list of lists to numpy arrays
    for k in HFRES_ARRAY_KEYS:
        if k in res:
            res[k] = np.array(res[k])
    return State(res)


def _metadata_yaml(stream):
    """Extract meta data block of a molsturm yaml file and return as a dictionary
       Returns None no meta data block found
    """
    if isinstance(stream, str):
        with open(stream, "r") as f:
            return _metadata_yaml(f)
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
def _dump_hdf5(hfres, path):
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


def _load_hdf5(path, meta_check=True):
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
        return State(hdf5.extract_group(h5f))


def _metadata_hdf5(path):
    """Extract meta data block of an hdf5 file and return as a dictionary
       Returns None no meta data block found
    """
    with h5py.File(path, "r") as h5f:
        return {at: h5f.attrs[at] for at in h5f.attrs}


#
# Common interface
#
__dump_funcs = {
    "hdf5": _dump_hdf5,
    "h5": _dump_hdf5,
    "hdf": _dump_hdf5,
    "yaml": _dump_yaml,
    "yml": _dump_yaml,
}

__load_funcs = {
    "hdf5": _load_hdf5,
    "h5": _load_hdf5,
    "hdf": _load_hdf5,
    "yaml": _load_yaml,
    "yml": _load_yaml,
}

__metadata_funcs = {
    "hdf5": _metadata_hdf5,
    "h5": _metadata_hdf5,
    "hdf": _metadata_hdf5,
    "yaml": _metadata_yaml,
    "yml": _metadata_yaml,
}


def dump_state(state, path, type="auto"):
    """
    state:   State object resulting from an SCF calculation
    path:    Path to dump the state
    type:    Type to use for dumping. If "auto" it will be determined
             from the file extension.
    """
    _, ext = os.path.splitext(path)
    ext = ext[1:]  # remove leading .

    if ext not in __dump_funcs:
        raise ValueError("Unrecognised dump file extension: " + ext + ".")
    else:
        return __dump_funcs[ext](state, path)


def load_state(path, type="auto", meta_check=True):
    """
    path:         Path to the file from which to load the state
    type:         Type to use for loading. If "auto" it will be determined
                  from the file extension.
    meta_check:   If False skips the check of the metadata

    Note that ``dump_state(state, file); state = load_state(file)``
    is an identity.

    Returns the loaded state object
    """
    _, ext = os.path.splitext(path)
    ext = ext[1:]  # remove leading .

    if ext not in __load_funcs:
        raise ValueError("Unrecognised dump file extension: " + ext + ".")
    else:
        return __load_funcs[ext](path, meta_check=meta_check)


def load_metadata(path, type="auto"):
    """Extract meta data block of a file produced with `dump_state` and
       return as a dictionary.
       Returns None no meta data block found
    """
    _, ext = os.path.splitext(path)
    ext = ext[1:]  # remove leading .

    if ext not in __load_funcs:
        raise ValueError("Unrecognised dump file extension: " + ext + ".")
    else:
        return __metadata_funcs[ext](path)
