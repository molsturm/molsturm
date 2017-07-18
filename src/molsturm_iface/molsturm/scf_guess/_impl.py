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

from .._basis import is_sturmian, is_gaussian
from .._constants import INPUT_PARAMETER_KEY
from .. import sturmian
import numpy as np

def __crop_orben_orbcoeff(restricted, n_alpha, n_beta, orben, orbcoeff):
  """
  Take an orben and an orbcoeff array and crop it in a way to serve as
  a guess for a system with n_alpha and n_beta electrons which is modelled
  in restricted or unrestricted manor
  """

  if restricted:
    return orben[:n_alpha], orbcoeff[:,:n_alpha]
  else:
    n = max(n_alpha, n_beta)
    n_orbs_alpha = orben.shape[0]//2
    if 2*n_orbs_alpha != orben.shape[0]:
      raise ValueError("orben does not have an even number of entries.")
    return np.concatenate( (orben[:n], orben[n_orbs_alpha:n_orbs_alpha+n]) ), \
        np.concatenate( (orbcoeff[:,:n], orbcoeff[:, n_orbs_alpha:n_orbs_alpha+n]),
                       axis=1)


def __extrapolate_from_previous_gaussian(old_hfres, kwargs):
  old_kwargs  = old_hfres[INPUT_PARAMETER_KEY]
  old_orben    = old_hfres["orben_f"]
  old_orbcoeff = old_hfres["orbcoeff_bf"]

  # TODO Here we assume the n_bas and n_alpha/n_beta values will
  #      be the same such that we can plainly copy the old results
  #      over
  n_elec = max(old_hfres["n_alpha"], old_hfres["n_beta"])
  return __crop_orben_orbcoeff(old_hfres["restricted"], old_hfres["n_alpha"],
                               old_hfres["n_beta"], old_orben, old_orbcoeff)

def __extrapolate_from_previous_sturmian(old_hfres, kwargs):
  old_kwargs  = old_hfres[INPUT_PARAMETER_KEY]
  old_orben    = old_hfres["orben_f"]
  old_orbcoeff = old_hfres["orbcoeff_bf"]

  # TODO
  # This is really hackish, since we reimplement the decisions
  # which go on inside the C++ code right here!
  #
  # The way to get around this is to only transfer nlm_basis to the C++
  # and do all the parsing on the python side, see also Parameters.hh
  if "nlm_basis" in old_kwargs or "nlm_basis" in kwargs:
    raise NotImplementedError("Cannot deal with sturmian bases, "
                              "which contain 'nlm_basis' at the moment")

  def get_nlm_basis(kwargs):
    # TODO Here we assume nlm order!
    order="nlm"

    n = kwargs["n_max"]
    l = kwargs["l_max"] if "l_max" in kwargs else None
    m = kwargs["m_max"] if "m_max" in kwargs else None
    return sturmian.CoulombSturmianBasis(kwargs["k_exp"], n, l, m, order=order)

  old_nlm_basis = get_nlm_basis(old_kwargs)
  nlm_basis     = get_nlm_basis(kwargs)

  assert len(old_nlm_basis) == old_hfres["n_bas"], \
         "Expected number of basis functions differ: n_bas == " + \
         str(old_hfres["n_bas"]) + " vs. len(old_nlm_basis) == " + \
         str(len(old_nlm_basis))

  if old_nlm_basis == nlm_basis:
    return __crop_orben_orbcoeff(old_hfres["restricted"], old_hfres["n_alpha"],
                                 old_hfres["n_beta"], old_orben, old_orbcoeff)
  else:
    proj = old_nlm_basis.obtain_projection_to(nlm_basis)
    return __crop_orben_orbcoeff(old_hfres["restricted"], old_hfres["n_alpha"],
                                 old_hfres["n_beta"], old_orben,
                                 np.matmul(proj, old_orbcoeff))

def extrapolate_from_previous(old_hfres, kwargs):
  """
  Extrapolate the old SCF results onto the new parameters to build a
  guess for the new SCF procedure
  """
  old_kwargs  = old_hfres[INPUT_PARAMETER_KEY]

  def check_agreement(key):
    both_have_key = key in old_kwargs and key in kwargs
    both_without_key = not key in old_kwargs and not key in kwargs
    if (both_without_key): return

    if not both_have_key or old_kwargs[key] != kwargs[key]:
      raise ValueError("Cannot extrapolate an SCF guess if the old and new "
                       "value for '" + key + "' differ.")

  check_agreement("basis_type")
  if is_gaussian(**kwargs):
    check_agreement("basis_set")
    return __extrapolate_from_previous_gaussian(old_hfres, kwargs)
  elif is_sturmian(**kwargs):
    return __extrapolate_from_previous_sturmian(old_hfres, kwargs)
  else:
    raise ValueError("Did not understand basis_type: '"+basis_type+"'.")

