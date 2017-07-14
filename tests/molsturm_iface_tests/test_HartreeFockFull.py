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

from NumCompTestCase import NumCompTestCase
import molsturm
import data_cs_be as data
import numpy as np
import unittest

@unittest.skipUnless(data.params["basis_type"] in molsturm.available_basis_types,
                     "Required basis type is not available")
class TestHartreeFockFull(NumCompTestCase):
  """This test should ensure, that we get exactly the same data
     through the python interface as we get without it by directly
     invoking the library via a C++ program.
  """
  @classmethod
  def setUpClass(cls):
    cls._hf_result = molsturm.hartree_fock(**data.params)

  def test_energies(self):
    for ene in data.ref_energies:
      self.assertAlmostEqual(self._hf_result[ene], data.ref_energies[ene],
                             tol=data.params["conv_tol"],prefix=ene+": ")

  def test_scf_convergence(self):
    self.assertGreaterEqual(data.ref_n_iter,self._hf_result["n_iter"])
    self.assertLessEqual(self._hf_result["final_error_norm"],data.params["conv_tol"])

    for key in data.ref_convergence_result:
      self.assertEqual(self._hf_result[key], data.ref_convergence_result[key])

  def test_orbital_energies(self):
    self.assertArrayAlmostEqual(self._hf_result["orben_f"],
                                data.ref_orbital_energies,
                                tol=data.params["conv_tol"],
                                prefix="Orbital energies: ")

  def test_coefficients(self):
    self.assertArrayAlmostEqual(self._hf_result["orbcoeff_bf"],
                                data.ref_coefficients,
                                tol=data.params["conv_tol"],
                                prefix="Coefficients: ")

  def test_fock(self):
    n_oa = self._hf_result["n_orbs_alpha"]
    slicemap = { "a": slice(None,n_oa), "b": slice(n_oa,None), }
    sizemap = { "a": n_oa, "b": self._hf_result["n_orbs_beta"], }

    for i in [ "a", "b" ]:
      for j in [ "a", "b" ]:
        fij = self._hf_result["fock_ff"][slicemap[i],slicemap[j]]
        try:
          ref_fij = data.ref_fock[i+j]
        except KeyError as e:
          ref_fij = np.zeros((sizemap[i],sizemap[j]))
        self.assertArrayAlmostEqual(fij, ref_fij, tol=data.params["conv_tol"],
                                    prefix="Fock "+i+"-"+j+": ")

  def test_repulsion_integrals(self):
    n_oa = self._hf_result["n_orbs_alpha"]
    slicemap = { "a": slice(None,n_oa), "b": slice(n_oa,None), }
    sizemap = { "a": n_oa, "b": self._hf_result["n_orbs_beta"], }

    for i in [ "a", "b" ]:
      for j in [ "a", "b" ]:
        for k in [ "a", "b" ]:
          for l in [ "a", "b" ]:
            Jijkl = self._hf_result["eri_ffff"][slicemap[i],slicemap[j],
                                                slicemap[k], slicemap[l]]
            try:
              ref_Jijkl = data.ref_repulsion_integrals[i+j+k+l]
            except KeyError as e:
              ref_Jijkl = np.zeros( (sizemap[i], sizemap[j], sizemap[k], sizemap[l]) )

            self.assertArrayAlmostEqual(Jijkl, ref_Jijkl, tol=data.params["conv_tol"],
                                        prefix="Repulsion tensor "+i+"-"+j+"-"+k+\
                                        "-"+l+": ")
