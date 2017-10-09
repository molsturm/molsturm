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

from . import _iface as iface
from .ParameterMap import ParameterMap
import collections
import numpy as np


class ParamSpecial:
    """Named tuple to contain value and type of an parameter,
       which is treated with a special update function on the C++ side."""
    def __init__(self, value, type):
        self.value = value
        self.type = type


def __to_iface_parameters(params, interface_type):
    """Make an interface parametermap from a ParameterMap parameter tree."""

    if not isinstance(params, ParameterMap):
        raise TypeError("params needs to be a ParameterMap object.")

    # Map which changes the keys string from the one in the tree
    # to the one actually exported to the C++ side
    # (this is for legacy and compatiblity reasons)
    map_key_remap = {
        "scf/diis_size": "scf/diis_n_prev_steps",
        "discretisation/k_exp": "discretisation/k_exponent"
    }

    # Construct the actual object, which we will return
    # The type of the object constructed is determined by the string
    # interface_type
    ret = getattr(iface, interface_type)()

    # Map which maps the actual type of a data value to the precise function name
    # inside the ret object to call and the type which is expected by the SWIG
    # interface code.
    ParamUpdate = collections.namedtuple("ParamUpdate", ["function", "typecast"])
    map_update = {
        int: ParamUpdate("update_int", int),
        np.uint64: ParamUpdate("update_size_t", int),
        float: ParamUpdate("update_scalar", float),
        str: ParamUpdate("update_string", str),
        bool: ParamUpdate("update_bool", bool),
    }

    # This call is special since it needs two inputs from the tree
    if "system/atom_numbers" in params or "system/coords" in params:
        ret.update_structure("system/structure",
                             params["system/atom_numbers"].value,
                             params["system/coords"].value)

    # Deal with the rest:
    for key_in in params:
        key = map_key_remap.get(key_in, key_in)

        if isinstance(params[key_in], ParamSpecial):
            # Deal with special nodes
            value = params[key_in].value
            typestr = params[key_in].type
            if typestr not in ["structure", "ignore"]:
                # Setting the value has not been already done above
                getattr(ret, "update_" + typestr)(key, value)
        else:
            value = params[key_in]
            try:
                updatefct, typecast = map_update[type(value)]
            except KeyError as e:
                raise ValueError("Did not understand type " + str(type(value)) +
                                 " of key " + key_in)
            getattr(ret, updatefct)(key, typecast(value))
    return ret
