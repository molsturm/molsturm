#!/usr/bin/env python3
## vi: tabstop=4 shiftwidth=4 softtabstop=4 expandtab
## ---------------------------------------------------------------------
##
## Copyright (C) 2018 by the molsturm authors
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

from ._iface_conversion import __to_iface_parameters
from . import _iface as iface
from .ScfParameters import ScfParameters
import numpy as np


def __get_integral(params, kind):
    """
    Obtain a one-electron integral by calling the molsturm interface
    """
    if isinstance(params, ScfParameters):
        # Make a shallow copy before we modify the parameters
        params = params.copy()
    else:
        raise TypeError("params needs to be an ScfParameters object.")

    # Check and normalise parameters:
    try:
        params.normalise()
    except (ValueError, KeyError, TypeError) as e:
        raise ValueError("ScfParameters params not valid: " + str(e))

    iface_params = __to_iface_parameters(params, "ScfParameters")

    n_bas, _, _ = params.scf_sizes
    if kind != "electron_repulsion_bbbb":
        output = np.zeros((n_bas, n_bas))
    else:
        output = np.zeros((n_bas, n_bas, n_bas, n_bas))

    # Call appropriate function to get integral values
    getattr(iface, kind)(iface_params, output)
    return output  # Return result


def overlap_bb(params):
    """Get overlap matrix in terms of basis functions"""
    return __get_integral(params, "overlap_bb")


def kinetic_bb(params):
    """Get kinetic matrix in terms of basis functions"""
    return __get_integral(params, "kinetic_bb")


def nuclear_attraction_bb(params):
    """Get nuclear attraction matrix in terms of basis functions"""
    return __get_integral(params, "nuclear_attraction_bb")


def electron_repulsion_bbbb(params):
    """Get electron repulsion integrals in terms of basis functions"""
    return __get_integral(params, "electron_repulsion_bbbb")
