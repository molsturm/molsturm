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
import numpy as np
from gint.util import basis_class_from_name
import collections
import gint.gaussian
import gint.sturmian.atomic
import inspect

from .scf_guess import extrapolate_from_previous
from .MolecularSystem import MolecularSystem
from ._constants import HFRES_ARRAY_KEYS, INPUT_PARAMETER_KEY
from collections import Iterable
from .sturmian import CoulombSturmianBasis

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


class soarray(np.ndarray):
    """Spin-Orbital array"""

    # Subclassing numpy arrays is a little different than usual.
    # See https://docs.scipy.org/doc/numpy/user/basics.subclassing.html
    # for a guideline.

    def __new__(cls, n_orbs, n_electrons, inner):
        # Create numpy array instance
        order = "C" if not np.isfortran(inner) else "F"
        obj = super(soarray, cls).__new__(cls, shape=inner.shape, dtype=inner.dtype,
                                          buffer=inner, strides=inner.strides, order=order)

        if isinstance(n_orbs, int):
            obj.n_orbs = (n_orbs // 2, n_orbs - n_orbs // 2)
        elif isinstance(n_orbs, tuple):
            obj.n_orbs = n_orbs
        else:
            raise TypeError("n_orbs needs to be a tuple or an int.")

        if isinstance(n_electrons, int):
            obj.n_electrons = (n_electrons // 2, n_electrons - n_electrons // 2)
        elif isinstance(n_electrons, tuple):
            obj.n_electrons = n_electrons
        else:
            raise TypeError("n_electrons needs to be a tuple or an int.")

        if obj.n_orbs[0] < obj.n_electrons[0] or obj.n_orbs[1] < obj.n_electrons[1]:
            raise ValueError("Need less electrons than orbitals.")

        if any(sh != sum(obj.n_orbs) for sh in inner.shape):
            raise ValueError("Shape of the inner tensor needs to be equal to the number "
                             "of orbitals.")
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return

        if not hasattr(obj, "n_electrons"):
            raise ValueError("Sorry, cannot construct soarray from a plain numpy ndarray")
        else:
            self.n_electrons = getattr(obj, "n_electrons")
            self.n_orbs = getattr(obj, "n_orbs")

    def block(self, block):
        if len(block) != self.ndim:
            raise ValueError("Tensor dimension and number of letters in "
                             "the block string needs to agree")

        n_alpha = self.n_electrons[0]
        n_beta = self.n_electrons[1]
        n_orbs_alpha = self.n_orbs[0]
        n_orbs_beta = self.n_orbs[1]

        occa = slice(0, n_alpha)
        occb = slice(n_orbs_alpha, n_orbs_alpha + n_beta)
        virta = slice(n_alpha, n_orbs_alpha)
        virtb = slice(n_orbs_alpha + n_beta, n_orbs_alpha + n_orbs_beta)

        ranges = {
            "o": (occa, occb),
            "v": (virta, virtb),
        }

        if any(a not in ranges for a in block):
            raise ValueError("Block string may only contain letters " +
                             ",".join(ranges.keys()))

        ret = self
        for i, b in enumerate(block):
            # Build numpy index:
            idxa = i * (slice(None), ) + (ranges[b][0], Ellipsis)
            idxb = i * (slice(None), ) + (ranges[b][1], Ellipsis)
            ret = np.concatenate((ret[idxa], ret[idxb]), axis=i)
        return ret


# TODO Quick and dirty wrapper to mimic the new syntax
class TmpState(dict):
    def __init__(self, other):
        for key in other:
            self.__setitem__(key, other[key])

    def eri_block(self, block="ffff"):
        """Get a block of the electron repulsion tensor
           in shell pair or chemists notation.
           No antisymmetrisation done yet.
        """
        if block != "ffff":
            raise NotImplementedError("Only ffff can be obtained at the moment.")
        return self.make_soarray(self.__getitem__("eri_ffff"))

    @property
    def eri(self):
        """
        Return the full eri tensor.
        """
        return self.eri_block()

    @property
    def fock(self):
        """
        Return the full fock matrix.
        """
        return self.make_soarray(self.__getitem__("fock_ff"))

    @property
    def n_orbs(self):
        return (self.__getitem__("n_orbs_alpha"), self.__getitem__("n_orbs_beta"))

    @property
    def n_electrons(self):
        return (self.__getitem__("n_alpha"), self.__getitem__("n_beta"))

    def make_soarray(self, inner):
        a = soarray(self.n_orbs, self.n_electrons, inner)
        return a


class ParameterMap(dict):
    # TODO Three like structure which can be built from
    #      the input parameters to hartree_fock
    #
    # Probably this guy should be directly linked to a GenMap
    # object back on the C++ side or made such that a GenMap
    # in C++ can be easily generated from this.
    pass

    def update(key, value, typestr=None):
        if type(value) == bool:
            typestr = "bool"
        elif type(value) == int:
            raise ValueError("For int type string is needed")
        elif type(value) == str:
            typestr = "string"


def hartree_fock(molecular_system, basis=None, basis_type=None,
                 conv_tol=5e-7, max_iter=25, n_eigenpairs=10000,
                 restricted=None,
                 guess="hcore", guess_esolver="auto",
                 eigensolver="auto", diis_size=4,
                 print_iterations=False, **kwargs):
    """
    Run a Hartree-Fock calculation with molsturm.

    molecular_system    The molecular system to model
    basis    A valid basis object. If None the basis will be constructed
             on the fly from teh basis_type and the kwargs.
    basis_type    String describing the type of basis function to use.


    Examples:
        hartree_fock(molecular_system=("Be",),    )

    # TODO old docs:
    The list of valid input
    parameters can be retrieved by the means of the list "hartree_fock_keys".

    A couple of selected kwargs:
      - basis_type        The type of the basis used for the calculation
      - coords            List of iterables of size 3: Coordinats of the atoms
      - atoms             List of the atom symbols (in the same order as coords)
      - guess             The guess method to use
    """

    #
    # Input normalisation
    #
    if isinstance(molecular_system, tuple):
        # Construct from the provided tuple:
        molecular_system = MolecularSystem(*molecular_system)
    if isinstance(molecular_system, str):
        # Only a single atom
        molecular_system = MolecularSystem(molecular_system)
    elif not isinstance(molecular_system, MolecularSystem):
        raise TypeError("The first argument needs to be a MolecularSystem object or a "
                        "tuple to setup a MolecularSystem object on the fly.")

    if basis is None:
        if basis_type is None:
            raise ValueError("Either the basis or the basis_type needs to be given.")
        Basis = basis_class_from_name(basis_type)

        # Get the name of the parameters which are accepted by the constructor
        # of the Basis object
        init_params = inspect.signature(Basis.__init__).parameters.keys()

        bas_kwargs = dict()
        for p in kwargs:
            if p in init_params:
                bas_kwargs[p] = kwargs[p]

        try:
            basis = Basis(molecular_system, **bas_kwargs)
        except (TypeError, ValueError) as e:
            raise ValueError("Invalid kwarg for basis construction: " + str(e))
    else:
        if basis_type is not None:
            raise ValueError("Only one of basis or basis_type may be given.")

    if restricted is None:
        restricted = molecular_system.is_closed_shell

    # TODO Check for wrongful or unused kwargs.

    #
    # Move system parameters into dict tree
    #
    ScfParam = collections.namedtuple("ScfParam", ["value", "type"])
    param_tree = dict()
    param_tree["restricted"] = ScfParam(restricted, "bool")

    #
    # System parameters
    #
    param_tree["system/n_alpha"] = ScfParam(molecular_system.n_alpha, "size_t")
    param_tree["system/n_beta"] = ScfParam(molecular_system.n_beta, "size_t")
    param_tree["system/coords"] = ScfParam(molecular_system.coords, "structure")
    param_tree["system/atom_numbers"] = \
        ScfParam(molecular_system.atom_numbers, "structure")

    #
    # Integral parameters (parts hackish)
    #
    param_tree["integrals/basis_type"] = ScfParam(basis.basis_type, "string")
    if isinstance(basis, gint.gaussian.Basis):
        # TODO Instead set basis functions here directly and omit passing
        #      the basis set name here and allow to pass the full description
        #      down to gint
        param_tree["integrals/basis_set"] = ScfParam(basis.basis_set_name, "string")

        # TODO This is still some sort of legacy stuff we kind of need
        #      to do at the moment unfortunately. One should remove that soon.
        param_tree["integrals/orbital_type"] = ScfParam("real_molecular", "orbital_type")
    elif isinstance(basis, gint.sturmian.atomic.Basis):
        if (molecular_system.n_atoms > 1):
            raise ValueError("Invalid basis: Atomic Coulomb-Strumians can only be used "
                             "on atoms and not on molecules.")

        param_tree["integrals/k_exponent"] = ScfParam(basis.k_exp, "scalar")
        param_tree["integrals/nlm_basis"] = ScfParam(basis.functions, "nlm_basis")

        # TODO This is still some sort of legacy stuff we kind of need
        #      to do at the moment unfortunately. One should remove that soon.
        param_tree["integrals/orbital_type"] = ScfParam("complex_atomic", "orbital_type")

        # TODO This is only temporary and until the gint layer has fully moved
        #      to using nlm_basis.
        print("WARNING: Right now we assume all sturmian basis sets to be dense")
        n_max, l_max, m_max = np.max(basis.functions, axis=0)
        param_tree["integrals/n_max"] = ScfParam(n_max, "int")
        param_tree["integrals/l_max"] = ScfParam(l_max, "int")
        param_tree["integrals/m_max"] = ScfParam(m_max, "int")
    else:
        raise TypeError("basis has an unrecognised type.")

    #
    # Build guess parameters
    #
    param_tree["guess/eigensolver/method"] = ScfParam(guess_esolver, "string")
    if isinstance(guess, str):
        if guess == "external":
            # TODO We would need to get the guess data into the
            #      coefficients from which we start the calculation somehow
            raise NotImplementedError("external guess not yet implemented.")
        else:
            param_tree["guess/method"] = ScfParam(guess, "string")
    else:
        raise NotImplementedError("guess from previous not yet implemented.")

    #
    # Build scf parameters
    #
    param_tree["scf/max_error_norm"] = ScfParam(conv_tol, "scalar")
    param_tree["scf/max_1e_energy_change"] = ScfParam(conv_tol * 100., "scalar")
    param_tree["scf/max_tot_energy_change"] = ScfParam(conv_tol / 4., "scalar")
    param_tree["scf/max_iter"] = ScfParam(max_iter, "size_t")
    param_tree["scf/n_eigenpairs"] = ScfParam(n_eigenpairs, "size_t")
    param_tree["scf/diis_size"] = ScfParam(diis_size, "size_t")
    param_tree["scf/print_iterations"] = ScfParam(print_iterations, "bool")

    # TODO import param_tree from a yaml file
    #      This will set n_bas, n_fock and n_spin from the cached input

    #
    # Checking and normalisation
    #
    if param_tree["restricted"].value and \
       param_tree["system/n_alpha"].value != param_tree["system/n_beta"].value:
        raise ValueError("Currently restricted is only possible for closed-shell systems")

    if param_tree["scf/n_eigenpairs"].value % 2 != 0:
        raise ValueError("The n_eigenpairs parameter applies to the accumulated number "
                         "of eigenpairs in the SCF calculations, i.e. the number of "
                         "alpha plus the number of beta orbitals. This is the done even "
                         "for restricted calculations. For now we further require this "
                         "number to be even number.")

    n_spin = 1 if param_tree["restricted"].value else 2  # Number of spin components
    n_bas = basis.size
    n_fock = min(n_bas, param_tree["scf/n_eigenpairs"].value // 2)
    param_tree["scf/n_eigenpairs"] = ScfParam(n_fock * n_spin, "size_t")

    # Check enough eigenpairs are requested.
    if n_fock < max(param_tree["system/n_alpha"].value, param_tree["system/n_beta"].value):
        raise ValueError("Cannot treat a system with " + str(param_tree["system/n_alpha"].value) + " alpha and " + str(param_tree["system/n_beta"].value) + " beta electrons with computing only " + str(n_fock) + " eigenpairs. Either choose a larger basis or a larger value for n_eigenpairs.")

    #
    # Set the keys inside the ScfParameters
    # TODO This section is really bäh
    #
    # Map which changes the keys from the tree to the keys exported to the C++ side
    key_remap = {"scf/diis_size": "scf/diis_n_prev_steps"}

    scfparams = iface.ScfParameters()
    scfparams.update_structure("system/structure",
                               param_tree["system/atom_numbers"].value,
                               param_tree["system/coords"].value)

    for key in param_tree:
        if param_tree[key].type in ["structure"]:
            pass  # Already done
        else:
            # TODO This is an ugly hack
            # Should be done transparently in some param_tree manager class
            # (see ParameterMap above)
            if param_tree[key].type == "size_t":
                print("WARNIG  using size_t hacke on " + key)
                value = np.asscalar(np.array([param_tree[key].value], dtype=np.uint64))
            elif param_tree[key].type == "bool":
                print("WARNIG  using bool hacke on " + key)
                value = bool(param_tree[key].value)
            elif param_tree[key].type == "int":
                print("WARNIG  using int hacke on " + key)
                value = int(param_tree[key].value)
            else:
                value = param_tree[key].value

            # If the key is in the key_remap it needs to be remapped
            if key in key_remap:
                key_out = key_remap[key]
            else:
                key_out = key

            # Use the type to construct the version of the update
            # function to call and call it with the value we have stored
            try:
                getattr(scfparams, "update_" + param_tree[key].type)(key_out, value)
            except:
                print(key, value, type(value), param_tree[key].type)
                raise

    #
    # Run scf
    #
    orben = np.empty((n_spin, n_fock))
    orbcoeff = np.empty((n_spin, n_bas, n_fock))
    solution_view = iface.ScfSolutionView(orben, orbcoeff)
    scf_kind = iface.RHF if param_tree["restricted"] else iface.UHF
    res = iface.self_consistent_field(scf_kind, scfparams, solution_view)

    #
    # Output
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

    return TmpState(out)

#   #
#   # Input
#   #
#   inputargs = dict()
#   if molecular_system:
#     if not isinstance(molecular_system, MolecularSystem):
#       raise TypeError("molecular_system needs to be of type MolecularSystem")
#     inputargs.update(molecular_system.as_hartree_fock_parameters())
# 
#   if basis:
#     if not isinstance(basis, CoulombSturmianBasis):
#       raise TypeError("basis needs to be of type CoulombSturmianBasis")
#     inputargs.update(basis.as_hartree_fock_parameters())
# 
#   # Keys which need to be parsed *after* all other ones have been.
#   delayed_keys = [ "guess", "guess_external_orben_f", "guess_external_orbcoeff_bf" ]
# 
#   for key in kwargs:
#     if not key in hartree_fock_keys:
#       raise ValueError("Keyword " + key + " is unknown to hartree_fock")
#     if key in delayed_keys:
#       continue
# 
#     # Copy key and (possibly transformed) value:
#     if key in __kwargs_parse_map:
#       __kwargs_parse_map[key](kwargs,inputargs)
#     else:
#       inputargs[key] = kwargs[key]
# 
#   # Deal with the delayed keys:
#   for key in delayed_keys:
#     if key in kwargs:
#       __kwargs_parse_map[key](kwargs,inputargs)
# 
#   # Setup parameters:
#   params = iface.Parameters()
#   for key in inputargs:
#     if key in __params_transform_map:
#       setattr(params, key, __params_transform_map[key](inputargs[key]))
#     else:
#       setattr(params, key, inputargs[key])
# 
#   if "restricted" in kwargs:
#     # Make a note that the user specified the restricted keyword
#     params.internal_restricted_set_by_user = True
# 
#   #
#   res = iface.hartree_fock(params)
#   #


def compute_derived_hartree_fock_energies(hfres):
  """Compute various derived hartree-fock energy terms."""
  # TODO It would be better to have this in a hfres class,
  #      which is returned by the hartree_fock function
  res=dict()

  # Prefix all energy keys use:
  prefix = "energy_"

  # Classify the different keys:
  zeroElectron = [ "nuclear_repulsion" ]    # No electrons involved
  twoElectron = [ "coulomb", "exchange" ]   # 2 electron terms

  # Keys with special treatment
  special = zeroElectron + twoElectron + [ "ground_state" ]
  oneElectron = sorted([ k[len(prefix):] for k in hfres
                         if k.startswith(prefix) and \
                           not k[len(prefix):] in special
                       ])

  # All energy terms:
  energies = zeroElectron + oneElectron + twoElectron

  # Store individual terms in returned dictionary
  res["terms"] = { ene : hfres[prefix+ene] for ene in energies }

  # Derived energies:
  res[prefix + "ground_state"] = hfres[prefix + "ground_state"]
  res[prefix + "1e"]           = sum([ hfres[prefix+ene] for ene in oneElectron ])
  res[prefix + "2e"]           = sum([ hfres[prefix+ene] for ene in twoElectron ])
  res[prefix + "electronic"]   = res[prefix + "1e"] + res[prefix + "2e"]
  res[prefix + "nuclear"]      = hfres[prefix + "nuclear_repulsion"]
  res[prefix + "potential"]    = sum([ hfres[prefix+ene] for ene in energies
                                     if not ene in [ "kinetic" ] ])
  res[prefix + "kinetic"]      = hfres[prefix + "kinetic"]
  res["virial_ratio"]          = - res[prefix + "potential"] / res[prefix + "kinetic"]

  return res

def compute_coulomb_ff(hfres):
  """Compute the coulomb matrix in MO space"""
  noa   = hfres["n_orbs_alpha"]
  na    = hfres["n_alpha"]
  nb    = hfres["n_beta"]
  jirep = hfres["eri_ffff"]

  return np.trace(jirep[   :na,          :na,       :, : ], axis1=0, axis2=1) + \
         np.trace(jirep[noa:noa + nb, noa:noa + nb, :, : ], axis1=0, axis2=1)


def compute_exchange_ff(hfres):
  """Compute the exchange matrix in MO space"""
  noa   = hfres["n_orbs_alpha"]
  na    = hfres["n_alpha"]
  nb    = hfres["n_beta"]
  jirep = hfres["eri_ffff"]

  return np.trace(jirep[   :na      , :, :,    :na      ], axis1=0, axis2=3) + \
         np.trace(jirep[noa:noa + nb, :, :, noa:noa + nb], axis1=0, axis2=3)

