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

from NumCompTestCase import NumCompTestCase
import molsturm
import data_cs_be as data
import numpy as np
import unittest


basis_type = data.input_parameters["integrals"]["basis_type"]


@unittest.skipUnless(basis_type in molsturm.available_basis_types,
                     "Required basis type " + basis_type +
                     " is not available")
class TestHartreeFockFull(NumCompTestCase):
    """This test should ensure, that we get exactly the same data
       through the python interface as we get without it by directly
       invoking the library via a C++ program.
    """
    @classmethod
    def setUpClass(cls):
        scfparams = molsturm.ScfParameters.from_dict(data.input_parameters)
        cls._hf_result = molsturm.self_consistent_field(scfparams)

    def test_energies(self):
        conv_tol = data.input_parameters["scf"]["conv_tol"]
        for ene in data.ref_energies:
            self.assertAlmostEqual(self._hf_result[ene], data.ref_energies[ene],
                                   tol=conv_tol, prefix=ene + ": ")

    def test_scf_convergence(self):
        conv_tol = data.input_parameters["scf"]["conv_tol"]
        self.assertGreaterEqual(data.ref_n_iter, self._hf_result["n_iter"])
        self.assertLessEqual(self._hf_result["final_error_norm"], conv_tol)

        for key in data.ref_convergence_result:
            self.assertEqual(self._hf_result[key], data.ref_convergence_result[key])

    def test_orbital_energies(self):
        conv_tol = data.input_parameters["scf"]["conv_tol"]
        self.assertArrayAlmostEqual(self._hf_result["orben_f"],
                                    data.ref_orbital_energies,
                                    tol=conv_tol,
                                    prefix="Orbital energies: ")

    def build_coefficient_rotation_matrix(self):
        """
        Return the rotation matrix, which rotates the
        coefficients from the obtained coefficients to the reference
        coefficients. The shape of the returned array is
        (n_spin, n_fock, n_fock)
        """
        n_oa = self._hf_result["n_orbs_alpha"]
        n_ba = self._hf_result["n_bas"]
        slicemap = {"a": slice(None, n_oa), "b": slice(n_oa, None), }
        slicemap_b = {"a": slice(None, n_ba), "b": slice(n_ba, None), }

        assert n_oa == self._hf_result["n_orbs_beta"]

        overlap_bb = self._hf_result["overlap_bb"]
        orbcoeff_bf = self._hf_result["orbcoeff_bf"]
        refcoeff_bf = data.ref_coefficients

        ret = np.empty((2, n_oa, n_oa))
        for i, s in enumerate(["a", "b"]):
            c_bf = orbcoeff_bf[:, slicemap[s]]
            cref_bf = refcoeff_bf[:, slicemap[s]]
            s_bb = overlap_bb[slicemap_b[s], slicemap_b[s]]
            ret[i, :, :] = np.dot(cref_bf.transpose(), np.dot(s_bb, c_bf))
        return ret

    def test_coefficients(self):
        """
        Test that reference coefficients and the coefficients
        actually obtained are related by unitary rotation.
        """
        conv_tol = data.input_parameters["scf"]["conv_tol"]
        ui_sff = self.build_coefficient_rotation_matrix()

        for i in range(2):
            s = ["a", "b"][i]
            ui_ff = ui_sff[i, :, :]
            self.assertArrayAlmostEqual(np.dot(ui_ff.transpose(), ui_ff),
                                        np.eye(*ui_ff.shape),
                                        tol=conv_tol,
                                        prefix=s + " coefficient rotation matrix: ")

    def test_fock(self):
        conv_tol = data.input_parameters["scf"]["conv_tol"]
        n_oa = self._hf_result["n_orbs_alpha"]
        slicemap = {"a": slice(None, n_oa), "b": slice(n_oa, None), }
        sizemap = {"a": n_oa, "b": self._hf_result["n_orbs_beta"], }

        for i in ["a", "b"]:
            for j in ["a", "b"]:
                fij = self._hf_result["fock_ff"][slicemap[i], slicemap[j]]
                try:
                    ref_fij = data.ref_fock[i + j]
                except KeyError as e:
                    ref_fij = np.zeros((sizemap[i], sizemap[j]))
                self.assertArrayAlmostEqual(fij, ref_fij, tol=conv_tol,
                                            prefix="Fock " + i + "-" + j + ": ")

    def test_repulsion_integrals(self):
        conv_tol = data.input_parameters["scf"]["conv_tol"]
        n_oa = self._hf_result["n_orbs_alpha"]
        slicemap = {"a": slice(None, n_oa), "b": slice(n_oa, None), }
        sizemap = {"a": n_oa, "b": self._hf_result["n_orbs_beta"], }

        u_sff = self.build_coefficient_rotation_matrix()
        umap = {"a": u_sff[0, :, :], "b": u_sff[1, :, :], }

        need_unitary_trans = True
        if np.allclose(umap["a"], np.eye(*umap["a"].shape)) and \
           np.allclose(umap["b"], np.eye(*umap["b"].shape)):
            need_unitary_trans = False

        for i in ["a", "b"]:
            for j in ["a", "b"]:
                for k in ["a", "b"]:
                    for l in ["a", "b"]:
                        Jijkl = self._hf_result["eri_ffff"][slicemap[i], slicemap[j],
                                                            slicemap[k], slicemap[l]]
                        try:
                            ref_Jijkl = data.ref_repulsion_integrals[i + j + k + l]
                        except KeyError as e:
                            ref_Jijkl = np.zeros((sizemap[i], sizemap[j],
                                                  sizemap[k], sizemap[l]))

                        if need_unitary_trans:
                            # Transform from the reference coefficients to the expected
                            # result for the actual coefficients obtained in this
                            # calculation
                            # Note, that the actual coeffiecients might be a unitary
                            # rotation of the reference by numerical noise
                            ref_Jijkl = np.einsum("ijkl,ia,jb,kc,ld->abcd", ref_Jijkl,
                                                  umap[i], umap[j], umap[k], umap[l])

                        self.assertArrayAlmostEqual(Jijkl, ref_Jijkl, tol=conv_tol,
                                                    prefix="Repulsion tensor " + i + "-" +
                                                    j + "-" + k + "-" + l + ": ")
