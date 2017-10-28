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

from ._iface_conversion import __to_iface_parameters
from . import _iface as iface
from . import yaml_utils
from .ScfParameters import ScfParameters
from .State import State
import numpy as np

# Old stuff ... will probably go some day
from ._constants import HFRES_ARRAY_KEYS, INPUT_PARAMETER_KEY


def self_consistent_field(params):
    """
    Run an SCF procedure from a ScfParameters object,
    which should contain all the parameters needed for the SCF.

    The routine will check the params object for validity
    before starting the SCF.

    For a more high-level entry point to start a calculation see
    also the function hartree_fock
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

    scf_sizes = params.scf_sizes
    if params["guess/method"] == "external":
        orben = params["guess/orben_f"].value
        orbcoeff = params["guess/orbcoeff_bf"].value
    else:
        orben = np.empty((scf_sizes.n_spin, scf_sizes.n_fock))
        orbcoeff = np.empty((scf_sizes.n_spin, scf_sizes.n_bas,
                             scf_sizes.n_fock))
    solution_view = iface.ScfSolutionView(orben, orbcoeff)

    # Convert to the interface parameters and call the C++ side
    # The convention used on the python side and the C++ side regarding
    # scf/n_eigenpairs differ. Python assums its always divisible by 2,
    # i.e. refers to the actual number of spin orbital eigenpairs computed.
    # C++ refers to the implicit number computed. Since for restricted
    # we only use the alpha block this value is hence half. This is
    # expressed here explicitly by using the scf_sizes object.
    iface_params = __to_iface_parameters(params, "ScfParameters")
    iface_params.update_size_t("scf/n_eigenpairs", scf_sizes.n_fock * scf_sizes.n_spin)
    res = iface.self_consistent_field(iface_params, solution_view)

    #
    # Output
    #
    # TODO This part has not been altered since the move to the new interface
    #      on the input side. Rework it.
    #
    assert res.n_bas == scf_sizes.n_bas
    res_keys = [k for k in dir(iface.ScfResults) if k[0] != "_"]
    shape_lookup = {"f": res.n_orbs_alpha + res.n_orbs_beta,
                    "b": res.n_bas}

    # Return orbital energies as (n_spin, n_fock)
    # Return orbcoeff as (n_spin, n_bas, n_fock)

    out = {k: getattr(res, k) for k in res_keys if k not in HFRES_ARRAY_KEYS}
    for k in HFRES_ARRAY_KEYS:
        if k in ["fock_bb", "overlap_bb"]:
            continue

        # Build the shape to cast the numpy arrays into from the
        # suffixes (e.g. _ffff, _bf) and the shape lookup object
        # we created above
        target_shape = tuple(shape_lookup[c] for c in k[k.rfind("_") + 1:])
        ary = np.array(getattr(res, k))
        if ary.size != 0:
            # If the size is 0, then the data has not been computed,
            # so we can ignore it
            out[k] = ary.reshape(target_shape)

    f_bb = np.array(getattr(res, "fock_bb"))
    if f_bb.size != 0:
        out["fock_bb"] = f_bb.reshape(2 * res.n_bas, 2 * res.n_bas)

    s_bb = np.array(getattr(res, "overlap_bb"))
    if s_bb.size != 0:
        out["overlap_bb"] = s_bb.reshape(2 * res.n_bas, 2 * res.n_bas)

    # Forward input parameters to output
    # TODO Later store as a ScfParameters object in the state.
    out[INPUT_PARAMETER_KEY] = yaml_utils.strip_special(params.to_dict(),
                                                        convert_np_arrays=True,
                                                        convert_np_scalars=True)

    return State(out)
