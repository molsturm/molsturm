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
from ._constants import HFRES_ARRAY_KEYS, INPUT_PARAMETER_KEY
from collections import Iterable
from .scf_guess import extrapolate_from_previous
from .MolecularSystem import MolecularSystem
from .sturmian import CoulombSturmianBasis
from gint.util import basis_class_from_name

def __to_double_vector(val):
  ret = iface.DoubleVector()
  for v in val.flatten():
    ret.push_back(float(v))
  return ret


def __to_coords(arr):
  if not isinstance(arr,Iterable):
    raise TypeError("Argument passed to \"coords\" needs to be iterable")

  ret = iface.DoubleVector()
  def parse_coord(c):
    if len(c) != 3:
      raise ValueError("All items of the coordinates array need have exactly 3 items.")
    for i in range(3):
      try:
        ret.push_back(float(c[i]))
      except ValueError:
        raise ValueError("All elements of coord vectors need to be castible to float")

  if isinstance(arr[0], Iterable):
    for c in arr:
      parse_coord(c)
  else:
    parse_coord(arr)
  return ret


def __to_nlm(arr):
  ret = iface.IntVector()
  for c in arr:
    if len(c) != 3:
      raise ValueError("Nlm array needs to be of shape n times 3")
    for i in range(3):
      ret.push_back(c[i])
  return ret


def __to_atom_numbers(li):
  ret = iface.IntVector()
  if isinstance(li,int):
    ret.push_back(int(li))
  elif isinstance(li,Iterable):
    for n in li:
      ret.push_back(int(n))
  else:
    raise ValueError("atom numbers needs be an int or a list of ints")
  return ret


def __to_atoms(li):
  ret = iface.StringVector()
  if isinstance(li,str):
    ret.push_back(str(li))
  elif isinstance(li,Iterable):
    for n in li:
      ret.push_back(str(n))
  else:
    raise ValueError("atoms needs to be a string or a list of strings")
  return ret


def __setup_user_provided_guess(kwargs, inputargs):
  """
  This is triggered if guess_external_orben_f or guess_external_orbcoeff_bf
  is provided by the user. This may only be done if "guess" is also set
  to external. Alternatively we set guess to external and check that both
  values are provided.
  """
  if not "guess" in kwargs or kwargs["guess"] != "external":
    raise ValueError("The parameters 'guess_external_orben_f' or "
                     "'guess_external_orbcoeff_bf' may only be provided iff 'guess'"
                     " is 'external'")

  if not "guess_external_orben_f" in kwargs \
     or not "guess_external_orbcoeff_bf" in kwargs:
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
    raise TypeError("The value provided to guess must be a previous hartree_fock result "+
                    " or a string describing a valid guess method.")

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
  "coords":                      __to_coords,
  "atom_numbers":                __to_atom_numbers,
  "atoms":                       __to_atoms,
  "nlm_basis":                   __to_nlm,
  "guess_external_orben_f":      __to_double_vector,
  "guess_external_orbcoeff_bf":  __to_double_vector,
}


# Extra keys witch which we deal with on the python level
HF_EXTRA_KEYS = [ ]


# These keys exist in the Parameters struct, but should not be exposed
# via the hartree_fock interface
HF_PARAMS_EXCLUDE_KEYS = [ "all" ] + [ k for k in dir(iface.Parameters)
                                       if k[0] == "_" or k.startswith("internal_") ]

"""The list of keys understood by the hartree_fock function"""
hartree_fock_keys = HF_EXTRA_KEYS + [ k for k in dir(iface.Parameters)
                                      if not k in HF_PARAMS_EXCLUDE_KEYS ]

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
        return self.__getitem__("eri_ffff")

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
        return self.__getitem__("fock_ff")


def hartree_fock(molecular_system, basis=None, basis_type=None,
                 conv_tol=5e-7, max_iter=25, n_eigenpairs=10000,
                 restricted=None,
                 guess="hcore", guess_esolver="auto",
                 eigensolver="auto",
                 print_iterations=False, **kwargs):
    """
    Run a Hartree-Fock calculation with molsturm.

    basis    A valid basis object. If None the basis will be constructed
             on the fly from teh basis_type and the kwargs.



    Examples:
        hartree_fock(molecular_system=("Be"),    )


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
    elif not isinstance(molecular_system, MolecularSystem):
        raise TypeError("The first argument needs to be a MolecularSystem object or a "
                        "tuple to setup a MolecularSystem object on the fly.")

    if basis is None:
        if basis_type is None:
            raise ValueError("Either the basis or the basis_type needs to be given.")
        Basis = basis_class_from_name(basis_type)
        bas_args = []
        bas_kwargs = {}

        # TODO use inspect to check the interface of the Basis
        #      object and use the kwargs to construct it.
        raise NotImplementedError("Basis construction on the fly is not yet implemented.")

        basis = Basis(molecular_system, *bas_args, **bas_kwargs)

    if restricted is None:
        restricted = molecular_system.is_closed_shell
    if restricted and not molecular_system.is_closed_shell:
        raise ValueError("Currenlty restricted is only possible for closed-shell systems")

    #
    # Input parameters
    #
    params = iface.ScfParameters()
    params.set_molecular_system(molecular_system.atom_numbers, molecular_system.coords,
                                molecular_system.n_alpha, molecular_system.n_beta)

    # The following functions set three kinds of parameters:
    # a) params.INTEGRAL
    #       Integral parameters. These are interpreted by gint::IntegralLookup and
    #       basically determine everything about the integral backend to use and
    #       the parameters it needs (molecular structure and so on.)
    #
    # b) params.GUESS
    #       Guess parameters, these are interpreted by the molsturm::scf_guess function
    #
    # c) params.SCF
    #       SCF parameters, these are read by molsturm::run_scf and from there
    #       subsequently passed onto the SCF algorithms inside gscf.
    SCF = params.SCF


    #
    # Integral parameters (parts hackish)
    #
    INTEGRAL = params.INTEGRAL
    params.set_param_string(INTEGRAL, "basis_type", basis.basis_type)
    if isinstance(basis, gint.gaussian.Basis):
        # TODO Instead set basis functions here directly and omit passing
        #      the basis set name here and allow to pass the full description
        #      down to gint
        params.set_param_string(INTEGRAL, "basis_set", basis.basis_set_name)

        # TODO This is still some sort of legacy stuff we kind of need
        #      to do at the moment unfortunately. One should remove that soon.
        params.set_integral_param_orbital_type("real_molecular")
    elif isinstance(basis, gint.sturmian.atomic.Basis):
        if (molecular_system.n_atoms > 1):
            raise ValueError("Invalid basis: Atomic Coulomb-Strumians can only be used "
                             "on atoms and not on molecules.")

        params.set_param_double(INTEGRAL, "k_exponent", basis.k_exp)
        params.set_integral_param_nlm_basis(np.array(basis.functions))

        # TODO Do this only if we are sure ... or better leave it off entirely
        #      It is some kind of legacy anyway.
        print("WARNING: Right now we assume all sturmian basis sets to be dense")
        params.set_param_int(INTEGRAL, "n_max", basis.n_max)
        params.set_param_int(INTEGRAL, "l_max", basis.l_max)
        params.set_param_int(INTEGRAL, "m_max", basis.m_max)

        # TODO This is still some sort of legacy stuff we kind of need
        #      to do at the moment unfortunately. One should remove that soon.
        params.set_integral_param_orbital_type("complex_atomic")
    else:
        raise TypeError("basis has an unrecognised type.")

    #
    # Build guess parameters
    #
    GUESS = params.GUESS

    # TODO The way to set the guess eigensolver parameters from C++ is:
    # guess_params.update(ScfGuessKeys::eigensolver_params,
    #                     {{EigensystemSolverKeys::method, params.guess_esolver}});
    # This is not yet implemented
    print("WARNING: Setting the guess eigensolver is not yet implemented.")

    if isinstance(guess, str):
        if guess == "external":
            # TODO We would need to get the guess data into the
            #      coefficients from which we start the calculation somehow
            raise NotImplementedError("external guess not yet implemented.")
        else:
            params.set_param_string(GUESS, "method", guess)
    else:
        raise NotImplementedError("guess from previous not yet implemented.")

    #
    # Build scf parameters
    #

    # TODO see build_scf_params in parse_parameters.cc




    #
    # Run scf
    #
    scf_kind = iface.RHF if restricted else iface.UHF
    results = iface.self_consistent_field(kind, params, solution_view)

    print("WARNING: Exporting the results or the input parameters is not implemented at the moment!")

  

    # TODO old stuff follows




  #
  # Input
  #
  inputargs = dict()
  if molecular_system:
    if not isinstance(molecular_system, MolecularSystem):
      raise TypeError("molecular_system needs to be of type MolecularSystem")
    inputargs.update(molecular_system.as_hartree_fock_parameters())

  if basis:
    if not isinstance(basis, CoulombSturmianBasis):
      raise TypeError("basis needs to be of type CoulombSturmianBasis")
    inputargs.update(basis.as_hartree_fock_parameters())

  # Keys which need to be parsed *after* all other ones have been.
  delayed_keys = [ "guess", "guess_external_orben_f", "guess_external_orbcoeff_bf" ]

  for key in kwargs:
    if not key in hartree_fock_keys:
      raise ValueError("Keyword " + key + " is unknown to hartree_fock")
    if key in delayed_keys:
      continue

    # Copy key and (possibly transformed) value:
    if key in __kwargs_parse_map:
      __kwargs_parse_map[key](kwargs,inputargs)
    else:
      inputargs[key] = kwargs[key]

  # Deal with the delayed keys:
  for key in delayed_keys:
    if key in kwargs:
      __kwargs_parse_map[key](kwargs,inputargs)

  # Setup parameters:
  params = iface.Parameters()
  for key in inputargs:
    if key in __params_transform_map:
      setattr(params, key, __params_transform_map[key](inputargs[key]))
    else:
      setattr(params, key, inputargs[key])

  if "restricted" in kwargs:
    # Make a note that the user specified the restricted keyword
    params.internal_restricted_set_by_user = True

  #
  res = iface.hartree_fock(params)
  #

  #
  # Output
  #
  # TODO Better return a dict-like class instead of a dict. That way
  #      we can use the class more easily and distinguish between results
  #      at different levels better.
  res_keys = [ k for k in dir(iface.HfResults) if k[0] != "_" ]
  shape_lookup = { "f": res.n_orbs_alpha + res.n_orbs_beta,
                   "b": res.n_bas }

  out = { k :getattr(res,k) for k in res_keys if not k in HFRES_ARRAY_KEYS }
  for k in HFRES_ARRAY_KEYS:
    # Build the shape to cast the numpy arrays into from the
    # suffixes (e.g. _ffff, _bf) and the shape lookup object
    # we created above
    target_shape = tuple( shape_lookup[c] for c in k[k.rfind("_")+1:] )
    ary = np.array(getattr(res,k))
    if ary.size != 0:
      # If the size is 0, then the data has not been computed,
      # so we can ignore it
      out[k] = ary.reshape(target_shape)

  # Forward input parameters to output
  out[INPUT_PARAMETER_KEY] = inputargs

  return TmpState(out)


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

