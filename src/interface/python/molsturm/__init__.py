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

from gint import construct_basis, available_basis_types
from ._scf import self_consistent_field
from ._scf_hliface import hartree_fock
from ._scf_hliface import compute_derived_hartree_fock_energies
from ._scf_hliface import compute_exchange_ff, compute_coulomb_ff
from .ScfParameters import ScfParameters
from ._print import *
from ._serialisation import dump_state, load_state, load_metadata
from ._iface import Version
from ._constants import INPUT_PARAMETER_KEY
from .System import System


# Mappings to previous names and deprecation
def MolecularSystem(*args, **kwargs):
    from warnings import warn
    warn("MolecularSystem is deprecated. Use System instead. "
         "This alias will be removed after the next release.")
    return System(*args, **kwargs)


def load_hdf5(*args, **kwargs):
    from warnings import warn
    warn("Use load_state function instead. "
         "This function will be removed after the next release.",
         DeprecationWarning)
    return load_state(*args, type="hdf5", **kwargs)


def load_yaml(*args, **kwargs):
    from warnings import warn
    warn("Use load_state function instead. "
         "This function will be removed after the next release.",
         DeprecationWarning)
    return load_state(*args, type="yaml", **kwargs)


def dump_hdf5(*args, **kwargs):
    from warnings import warn
    warn("Use dump_state function instead. "
         "This function will be removed after the next release.",
         DeprecationWarning)
    return dump_state(*args, type="hdf5", **kwargs)


def dump_yaml(*args, **kwargs):
    from warnings import warn
    warn("Use dump_state function instead. "
         "This function will be removed after the next release.",
         DeprecationWarning)
    return dump_state(*args, type="yaml", **kwargs)
