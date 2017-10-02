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

from .State import State
from .ScfParameters import ScfParameters
from .MolecularSystem import MolecularSystem
import numpy as np
from . import _iface as iface
from ._iface_conversion import __to_iface_parameters
from .scf_guess import extrapolate_from_previous
import inspect
from gint.util import split_basis_type

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
        pass
    elif isinstance(params, dict):
        return self_consistent_field(ScfParameters.from_dict(params))
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
    #
    # TODO Better return a dict-like class instead of a dict. That way
    #      we can use the class more easily and distinguish between results
    #      at different levels better.
    res_keys = [k for k in dir(iface.ScfResults) if k[0] != "_"]
    shape_lookup = {"f": res.n_orbs_alpha + res.n_orbs_beta,
                    "b": res.n_bas}

    # Return orbital energies as (n_spin, n_fock)
    # Return orbcoeff as (n_spin, n_bas, n_fock)

    out = {k: getattr(res, k) for k in res_keys if k not in HFRES_ARRAY_KEYS}
    for k in HFRES_ARRAY_KEYS:
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

    # Forward input parameters to output
    out[INPUT_PARAMETER_KEY] = params.to_dict()

    return State(out)


def hartree_fock(system, basis=None, basis_type=None,
                 conv_tol=None, max_iter=None, n_eigenpairs=None,
                 restricted=None,
                 guess=None, guess_esolver=None,
                 eigensolver=None, diis_size=None,
                 print_iterations=None,
                 **kwargs):
    """
    Run a Hartree-Fock calculation with molsturm.

    system    The molecular system to model
    basis    A valid basis object. If None the basis will be constructed
             on the fly from teh basis_type and the kwargs.
    basis_type    String describing the type of basis function to use.
    params        Path to a yaml file containing the parameter hierachy
                  or a ParameterMap object containing the parameter hierachy.

    Examples:
        hartree_fock(system=("Be",), basis_type="gaussian",
                     basis_set_name="sto-3g")

    For a more low-level entry point starting directly from a set of parameters
    see the function self_consistent_field
    """
    #
    # Input normalisation and building of parameter tree
    #
    if isinstance(system, str):
        # Only a single atom
        system = MolecularSystem(system)
    elif not isinstance(system, MolecularSystem):
        raise TypeError("The first argument needs to be a MolecularSystem object or a "
                        "tuple to setup a MolecularSystem object on the fly.")

    if basis is None:
        if basis_type is None:
            raise ValueError("Either the basis or the basis_type needs to be given.")
        Basis, backend = split_basis_type(basis_type)

        # Get the name of the parameters which are accepted by the constructor
        # of the Basis object
        init_params = inspect.signature(Basis.__init__).parameters.keys()

        bas_kwargs = dict()
        for p in kwargs:
            if p in init_params:
                bas_kwargs[p] = kwargs[p]
        for key in bas_kwargs:
            del kwargs[key]

        try:
            basis = Basis(system, **bas_kwargs, backend=backend)
        except (TypeError, ValueError) as e:
            raise ValueError("Invalid kwarg for basis construction: " + str(e))
    else:
        if basis_type is not None:
            raise ValueError("Only one of basis or basis_type may be given.")

    # TODO Check for wrongful or unused kwargs.

    # Build a parameter tree from the commandline arguments
    param_tree = ScfParameters()
    param_tree.system = system
    param_tree.basis = basis

    # Add the scf parameters
    for frm, to, typ in [
        ("conv_tol", "scf/conv_tol", float),
        ("max_iter", "scf/max_iter", np.uint64),
        ("n_eigenpairs", "scf/n_eigenpairs", np.uint64),
        ("diis_size", "scf/diis_size", np.uint64),
        ("print_iterations", "scf/print_iterations", bool),
        ("eigensolver", "scf/eigensolver/method", str),
        ("restricted", "scf/restricted", bool),
        ("guess_esolver", "guess/eigensolver/method", str),
    ]:
        val = locals()[frm]
        if val is not None:
            param_tree[to] = typ(val)

    # Set all other key-value pairs given on the commandline
    for k, v in kwargs.items():
        param_tree[k] = v

    # Build guess parameters
    if isinstance(guess, str):
        if guess == "external":
            raise ValueError("External guess is only available via the "
                             "ScfParameters object.")
        else:
            param_tree["guess/method"] = str(guess)
    elif isinstance(guess, State):
        # Normalise parameter tree such that the extrapolate_from_previous
        # has access to a proper set of parameters.
        param_tree.normalise()

        # Extrapolate and add the external guess:
        orben_f, orbcoeff_bf = extrapolate_from_previous(guess, param_tree)
        param_tree.set_guess_external(orben_f, orbcoeff_bf)
    elif guess is None:
        pass
    else:
        raise TypeError("guess can only be a string or a State object.")

    return self_consistent_field(param_tree)


def compute_derived_hartree_fock_energies(hfres):
    """Compute various derived hartree-fock energy terms."""
    # TODO It would be better to have this in a hfres class,
    #      which is returned by the hartree_fock function
    res = dict()

    # Prefix all energy keys use:
    prefix = "energy_"

    # Classify the different keys:
    zeroElectron = ["nuclear_repulsion"]    # No electrons involved
    twoElectron = ["coulomb", "exchange"]   # 2 electron terms

    # Keys with special treatment
    special = zeroElectron + twoElectron + ["ground_state"]
    oneElectron = sorted([k[len(prefix):]
                          for k in hfres
                          if k.startswith(prefix) and not k[len(prefix):] in special
                          ])

    # All energy terms:
    energies = zeroElectron + oneElectron + twoElectron

    # Store individual terms in returned dictionary
    res["terms"] = {ene: hfres[prefix + ene] for ene in energies}

    # Derived energies:
    res[prefix + "ground_state"] = hfres[prefix + "ground_state"]
    res[prefix + "1e"] = sum([hfres[prefix + ene] for ene in oneElectron])
    res[prefix + "2e"] = sum([hfres[prefix + ene] for ene in twoElectron])
    res[prefix + "electronic"] = res[prefix + "1e"] + res[prefix + "2e"]
    res[prefix + "nuclear"] = hfres[prefix + "nuclear_repulsion"]
    res[prefix + "potential"] = sum([hfres[prefix + ene] for ene in energies
                                     if ene not in ["kinetic"]])
    res[prefix + "kinetic"] = hfres[prefix + "kinetic"]
    res["virial_ratio"] = - res[prefix + "potential"] / res[prefix + "kinetic"]

    return res


def compute_coulomb_ff(hfres):
    """Compute the coulomb matrix in MO space"""
    noa = hfres["n_orbs_alpha"]
    na = hfres["n_alpha"]
    nb = hfres["n_beta"]
    jirep = hfres["eri_ffff"]

    return np.trace(jirep[:na, :na, :, :], axis1=0, axis2=1) + \
        np.trace(jirep[noa:noa + nb, noa:noa + nb, :, :], axis1=0, axis2=1)


def compute_exchange_ff(hfres):
    """Compute the exchange matrix in MO space"""
    noa = hfres["n_orbs_alpha"]
    na = hfres["n_alpha"]
    nb = hfres["n_beta"]
    jirep = hfres["eri_ffff"]

    return np.trace(jirep[:na, :, :, :na], axis1=0, axis2=3) + \
        np.trace(jirep[noa:noa + nb, :, :, noa:noa + nb], axis1=0, axis2=3)
