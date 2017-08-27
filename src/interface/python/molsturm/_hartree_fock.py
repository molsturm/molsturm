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

from gint.util import split_basis_type
from . import _iface as iface
from .MolecularSystem import MolecularSystem
from .State import State
from .ParameterMap import ParameterMap
from ._iface_conversion import ParamSpecial, __to_iface_parameters
import gint.gaussian
import gint.sturmian.atomic
import inspect
import numpy as np

from .scf_guess import extrapolate_from_previous
from ._constants import HFRES_ARRAY_KEYS, INPUT_PARAMETER_KEY


def __to_double_vector(val):
    ret = iface.DoubleVector()
    for v in val.flatten():
        ret.push_back(float(v))
    return ret


def __setup_user_provided_guess(kwargs, inputargs):
    """
    This is triggered if guess_external_orben_f or guess_external_orbcoeff_bf
    is provided by the user. This may only be done if "guess" is also set
    to external. Alternatively we set guess to external and check that both
    values are provided.
    """
    if "guess" not in kwargs or kwargs["guess"] != "external":
        raise ValueError("The parameters 'guess_external_orben_f' or "
                         "'guess_external_orbcoeff_bf' may only be provided iff 'guess'"
                         " is 'external'")

    if "guess_external_orben_f" not in kwargs \
       or "guess_external_orbcoeff_bf" not in kwargs:
        raise ValueError("If the guess parameter is set to 'external', the parameters "
                         "'guess_external_orben_f' (for the guess orbital energies) and "
                         "'guess_external_orbcoeff_bf' (and the guess orthonormal SCF "
                         "coefficients) are both required.")

    inputargs["guess_external_orben_f"] = kwargs["guess_external_orben_f"]
    inputargs["guess_external_orbcoeff_bf"] = kwargs["guess_external_orbcoeff_bf"]
    inputargs["guess"] = "external"


def __setup_guess(kwargs, inputargs):
    """
    This is triggered if guess is provided by the user.
    """
    guess_key = "guess"
    guess = kwargs[guess_key]

    if isinstance(guess, str):
        if guess == "external":
            __setup_user_provided_guess(kwargs, inputargs)
        else:
            inputargs[guess_key] = guess
        return
    elif not isinstance(guess, dict):
        raise TypeError("The value provided to guess must be a previous hartree_fock "
                        "result or a string describing a valid guess method.")

    inputargs["guess_external_orben_f"], inputargs["guess_external_orbcoeff_bf"] = \
        extrapolate_from_previous(guess, **kwargs)
    inputargs["guess"] = "external"


# Map from input kwargs to inputargs list
__kwargs_parse_map = {
    "guess":                       __setup_guess,
    "guess_external_orben_f":      __setup_user_provided_guess,
    "guess_external_orbcoeff_bf":  __setup_user_provided_guess,
}

# Special input parameters for which the above conversion functions need
# to be used before assignment
__params_transform_map = {
    "guess_external_orben_f":      __to_double_vector,
    "guess_external_orbcoeff_bf":  __to_double_vector,
}


# TODO The n_bas_tmp parameter should go.
#      It is not the way this should be done.
#
#      Note, that we cannot just put it into the tree, since
#      the user (who might specify the tree via a yaml file)
#      certainly does not know n_bas by himself.
def self_consistent_field(param_tree, n_bas_tmp):
    """
    Run an SCF procedure from a param_tree object, which should
    contain all the parameters needed for the SCF.

    The parameter tree should be ParameterMap object and it will
    be checked for basic consistency before running the scf
    """
    # TODO We need to set defaults here as well in case we get
    #      an incomplete parameter tree, where some parameters
    #      are missing. Or shall we rely on the C++ side. Probably
    #      not because there could be different values compared
    #      to the ones we choos in the hartree_fock python function.

    # TODO To determine n_bas from the param_tree we probably need to
    #      construct a basis object here first and call its size function.

    #
    # Checking the param_tree object
    #
    if param_tree["scf/kind"] == "restricted" and \
       param_tree["system/n_alpha"] != param_tree["system/n_beta"]:
        raise ValueError("Currently restricted is only possible for closed-shell systems")

    if param_tree["scf/kind"] == "restricted-open":
        raise ValueError("restricted-open is currently not implemented.")

    if param_tree["scf/n_eigenpairs"] % 2 != 0:
        raise ValueError("The n_eigenpairs parameter applies to the accumulated number "
                         "of eigenpairs in the SCF calculations, i.e. the number of "
                         "alpha plus the number of beta orbitals. This is the done even "
                         "for restricted calculations. For now we further require this "
                         "number to be even number.")

    # Number of spin components
    n_spin = 1 if param_tree["scf/kind"] == "restricted" else 2
    n_bas = n_bas_tmp
    n_fock = min(n_bas, param_tree["scf/n_eigenpairs"] // 2)
    param_tree["scf/n_eigenpairs"] = np.uint64(n_fock * n_spin)

    # Check enough eigenpairs are requested.
    if n_fock < max(param_tree["system/n_alpha"], param_tree["system/n_beta"]):
        raise ValueError("Cannot treat a system with " +
                         str(param_tree["system/n_alpha"]) + " alpha and " +
                         str(param_tree["system/n_beta"]) +
                         " beta electrons with computing only " + str(n_fock) +
                         " eigenpairs. Either choose a larger basis or a larger "
                         "value for n_eigenpairs.")

    #
    # Run scf
    #
    orben = np.empty((n_spin, n_fock))
    orbcoeff = np.empty((n_spin, n_bas, n_fock))
    solution_view = iface.ScfSolutionView(orben, orbcoeff)

    iface_params = __to_iface_parameters(param_tree, "ScfParameters")
    res = iface.self_consistent_field(iface_params, solution_view)

    #
    # Output
    #
    # TODO This part has not been altered since the move to the new interface
    #      on the input side. Rework it.
    #
    assert res.n_bas == n_bas
    #
    # TODO Better return a dict-like class instead of a dict. That way
    #      we can use the class more easily and distinguish between results
    #      at different levels better.
    res_keys = [k for k in dir(iface.ScfResults) if k[0] != "_"]
    shape_lookup = {"f": res.n_orbs_alpha + res.n_orbs_beta,
                    "b": res.n_bas}

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

    # TODO not yet possible
    print("WARNING: Exporting the input parameters is not re-implemented at the moment!")
    # Forward input parameters to output
    # out[INPUT_PARAMETER_KEY] = inputargs
    # TODO sericalise the dict_tree to something yaml and hdf5 can shamelessly write

    return State(out)


def __scf_tree_from_args(molecular_system, basis, conv_tol, max_iter, n_eigenpairs,
                         restricted, guess, guess_esolver, eigensolver, diis_size,
                         print_iterations):
    """
    Read the input arguments of hartree_fock and build a param_tree object from it.
    """

    param_tree = ParameterMap()

    #
    # System parameters
    #
    param_tree["system/n_alpha"] = np.uint64(molecular_system.n_alpha)
    param_tree["system/n_beta"] = np.uint64(molecular_system.n_beta)
    param_tree["system/coords"] = ParamSpecial(molecular_system.coords, "structure")
    param_tree["system/atom_numbers"] = \
        ParamSpecial(molecular_system.atom_numbers, "structure")

    #
    # Integral parameters (parts hackish)
    #
    param_tree["integrals/basis_type"] = str(basis.basis_type)
    if isinstance(basis, gint.gaussian.Basis):
        # TODO Instead set basis functions here directly and omit passing
        #      the basis set name here and allow to pass the full description
        #      down to gint
        param_tree["integrals/basis_set"] = str(basis.basis_set_name)

        # TODO This is still some sort of legacy stuff we kind of need
        #      to do at the moment unfortunately. One should remove that soon.
        param_tree["integrals/orbital_type"] = ParamSpecial("real_molecular",
                                                            type="orbital_type")
    elif isinstance(basis, gint.sturmian.atomic.Basis):
        if (molecular_system.n_atoms > 1):
            raise ValueError("Invalid basis: Atomic Coulomb-Strumians can only be used "
                             "on atoms and not on molecules.")

        param_tree["integrals/k_exp"] = float(basis.k_exp)
        param_tree["integrals/nlm_basis"] = ParamSpecial(basis.functions, "nlm_basis")

        # TODO This is still some sort of legacy stuff we kind of need
        #      to do at the moment unfortunately. One should remove that soon.
        param_tree["integrals/orbital_type"] = ParamSpecial("complex_atomic",
                                                            type="orbital_type")

        # TODO This is only temporary and until the gint layer has fully moved
        #      to using nlm_basis.
        print("WARNING: Right now we assume all sturmian basis sets to be dense")
        n_max, l_max, m_max = np.max(basis.functions, axis=0)
        param_tree["integrals/n_max"] = int(n_max)
        param_tree["integrals/l_max"] = int(l_max)
        param_tree["integrals/m_max"] = int(m_max)
    else:
        raise TypeError("basis has an unrecognised type.")

    #
    # Build guess parameters
    #
    param_tree["guess/eigensolver/method"] = str(guess_esolver)
    if isinstance(guess, str):
        if guess == "external":
            # TODO We would need to get the guess data into the
            #      coefficients from which we start the calculation somehow
            raise NotImplementedError("external guess not yet implemented.")
        else:
            param_tree["guess/method"] = str(guess)
    else:
        raise NotImplementedError("guess from previous not yet implemented.")
    print("WARNING: External guess of guess from previous is not yet re-implemented.")

    #
    # Build scf parameters
    #
    param_tree["scf/max_error_norm"] = float(conv_tol)
    param_tree["scf/max_1e_energy_change"] = float(conv_tol * 100.)
    param_tree["scf/max_tot_energy_change"] = float(conv_tol / 4.)
    param_tree["scf/max_iter"] = np.uint64(max_iter)
    param_tree["scf/n_eigenpairs"] = np.uint64(n_eigenpairs)
    param_tree["scf/diis_size"] = np.uint64(diis_size)
    param_tree["scf/print_iterations"] = bool(print_iterations)

    # Set the precise scf kind done.
    param_tree["scf/kind"] = "restricted-open" if restricted else "unrestricted"
    if restricted and molecular_system.is_closed_shell:
        param_tree["scf/kind"] = "restricted"

    return param_tree


def hartree_fock(molecular_system=None, basis=None, basis_type=None,
                 conv_tol=5e-7, max_iter=25, n_eigenpairs=10000,
                 restricted=None,
                 guess="hcore", guess_esolver="auto",
                 eigensolver="auto", diis_size=4,
                 print_iterations=False,
                 params=None,
                 **kwargs):
    """
    Run a Hartree-Fock calculation with molsturm.

    molecular_system    The molecular system to model
    basis    A valid basis object. If None the basis will be constructed
             on the fly from teh basis_type and the kwargs.
    basis_type    String describing the type of basis function to use.
    params        Path to a yaml file containing the parameter hierachy
                  or a ParameterMap object containing the parameter hierachy.

    Examples:
        hartree_fock(molecular_system=("Be",), basis_type="gaussian",
                     basis_set_name="sto-3g")
    """

    if params is not None:
        if isinstance(params, ParameterMap):
            # If the argument of hartree_fock is a ParameterMap,
            # short-circuit all of the normalisation and run the inner function.
            return self_consistent_field(params)
        elif isinstance(params, str):
            # TODO Read parameters from a yaml file and run using those
            raise NotImplementedError("Cannot read parameters from yaml file yet.")
        else:
            raise TypeError("params needs to be a ParameterMap object.")

    #
    # Input normalisation and building of parameter tree
    #
    if isinstance(molecular_system, str):
        # Only a single atom
        molecular_system = MolecularSystem(molecular_system)
    elif not isinstance(molecular_system, MolecularSystem):
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

        try:
            basis = Basis(molecular_system, **bas_kwargs, backend=backend)
        except (TypeError, ValueError) as e:
            raise ValueError("Invalid kwarg for basis construction: " + str(e))
    else:
        if basis_type is not None:
            raise ValueError("Only one of basis or basis_type may be given.")

    if restricted is None:
        restricted = molecular_system.is_closed_shell

    # TODO Check for wrongful or unused kwargs.

    # Build a parameter tree from the commandline arguments
    param_tree = __scf_tree_from_args(molecular_system, basis, conv_tol, max_iter,
                                      n_eigenpairs, restricted, guess, guess_esolver,
                                      eigensolver, diis_size, print_iterations)
    return self_consistent_field(param_tree, basis.size)


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
