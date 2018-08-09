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

from ._adcc import run_adcc_adcman, adcc_found
from ._pyscf import pyscf_found, run_fci_pyscf
from .._constants import INPUT_PARAMETER_KEY
from .._mp2 import mp2 as run_mp2

def __build_available_methods():
  ret=[ "mp2" ]
  if adcc_found:
    ret += [ "mp3", "adc0", "adc1", "adc2", "adc2x", "adc3" ]
  if pyscf_found:
    ret += [ "fci" ]
  return ret
available_methods = __build_available_methods()

def __assert_available(method):
  """Check if the method is available and bail out if not"""
  if not method in available_methods:
    raise RuntimeError("molsturm.posthf did not find a way to do an " + method +
                       " calculation. Check that this method is available.")

def __forward_to_adcc(hfres, method, **params):
  ret = run_adcc_adcman(hfres, method=method, **params)
  del ret[INPUT_PARAMETER_KEY]["method"]
  return ret


def mp2(hfres, **params):
  """
  Take hf results and extra parameters as an input
  and run an Møller Plesset 2nd order calculation.

  Return the resulting dictionary of computed data.
  """
  __assert_available("mp2")
  e_mp2, t2 = run_mp2(hfres, **params)
  return {
    "energy_mp2":           e_mp2,
    "energy_ground_state":  hfres["energy_ground_state"] + e_mp2,
    "t2_ovov":              t2,
    INPUT_PARAMETER_KEY:    params,
  }

def mp3(hfres, **params):
  """
  Take hf results and extra parameters as an input
  and run an Møller Plesset 3rd order calculation.

  Return the resulting dictionary of computed data.
  """
  __assert_available("mp3")
  return __forward_to_adcc(hfres, "mp3", **params)

def adc0(hfres, **params):
  """
  Take hf results and extra parameters as an input
  and run an Algebraic Diagrammatic Construction
  0th order calculation.

  Return the resulting dictionary of computed data.
  """
  __assert_available("adc0")
  return __forward_to_adcc(hfres, "adc0", **params)

def adc1(hfres, **params):
  """
  Take hf results and extra parameters as an input
  and run an Algebraic Diagrammatic Construction
  1st order calculation.

  Return the resulting dictionary of computed data.
  """
  __assert_available("adc1")
  return __forward_to_adcc(hfres, "adc1", **params)

def adc2(hfres, **params):
  """
  Take hf results and extra parameters as an input
  and run an Algebraic Diagrammatic Construction
  strict 2nd order calculation.

  Return the resulting dictionary of computed data.
  """
  __assert_available("adc2")
  return __forward_to_adcc(hfres, "adc2s", **params)

def adc2x(hfres, **params):
  """
  Take hf results and extra parameters as an input
  and run an Algebraic Diagrammatic Construction
  extended 2nd order calculation.

  Return the resulting dictionary of computed data.
  """
  __assert_available("adc2x")
  return __forward_to_adcc(hfres, "adc2x", **params)

def adc3(hfres, **params):
  """
  Take hf results and extra parameters as an input
  and run an Algebraic Diagrammatic Construction
  3rd order calculation.

  Return the resulting dictionary of computed data.
  """
  __assert_available("adc3")
  return __forward_to_adcc(hfres, "adc3", **params)

def fci(hfres, **kwargs):
  """
  Take the hf results and extra parameters and
  run a Full-CI calculation. Return the dictionary
  of computed data.

  A few kwargs:
    verbosity     Control the verbosity of the output
    n_roots       Number of FCI roots to compute
    conv_tol      Tolerance fo consiter root to be converged
  """
  __assert_available("fci")
  return run_fci_pyscf(hfres, **kwargs)
