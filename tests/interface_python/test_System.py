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

import unittest
import molsturm
import gint


class TestSystem(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        water_atoms = ["H", "O", "H"]
        water_coords = [[0, 0, 1], [0, 0, 0], [1, 0, 0]]
        no_atoms = ["N", "O"]
        no_coords = [[0, 0, 0], [1, 0, 0]]

        cls.empty = molsturm.System()
        cls.water = molsturm.System(water_atoms, water_coords)
        cls.water_triplet = molsturm.System(water_atoms, water_coords, (6, 4))
        cls.no = molsturm.System(no_atoms, no_coords)
        cls.oxide = molsturm.System(["O"], [[0, 0, 0]], (5, 5))

    def test_0_constructor(self):
        # Test what the setup constructed
        self.assertEqual(len(self.empty.atom_numbers), 0)
        self.assertEqual(len(self.empty.coords), 0)
        self.assertEqual(self.empty.n_alpha, 0)
        self.assertEqual(self.empty.n_beta, 0)

        self.assertListEqual(self.water.atom_numbers.tolist(), [1, 8, 1])
        self.assertListEqual(self.water.coords[0].tolist(), [0, 0, 1])
        self.assertListEqual(self.water.coords[1].tolist(), [0, 0, 0])
        self.assertListEqual(self.water.coords[2].tolist(), [1, 0, 0])
        self.assertEqual(self.water.n_alpha, 5)
        self.assertEqual(self.water.n_beta, 5)

        self.assertListEqual(self.water_triplet.atom_numbers.tolist(), [1, 8, 1])
        self.assertListEqual(self.water_triplet.coords[0].tolist(), [0, 0, 1])
        self.assertListEqual(self.water_triplet.coords[1].tolist(), [0, 0, 0])
        self.assertListEqual(self.water_triplet.coords[2].tolist(), [1, 0, 0])
        self.assertEqual(self.water_triplet.n_alpha, 6)
        self.assertEqual(self.water_triplet.n_beta, 4)

        self.assertListEqual(self.no.atom_numbers.tolist(), [7, 8])
        self.assertListEqual(self.no.coords[0].tolist(), [0, 0, 0])
        self.assertListEqual(self.no.coords[1].tolist(), [1, 0, 0])
        self.assertEqual(self.no.n_alpha, 8)
        self.assertEqual(self.no.n_beta, 7)

        self.assertListEqual(self.oxide.atom_numbers.tolist(), [8])
        self.assertListEqual(self.oxide.coords[0].tolist(), [0, 0, 0])
        self.assertEqual(self.oxide.n_alpha, 5)
        self.assertEqual(self.oxide.n_beta, 5)

        # Construction from atom names and atom numbers
        wat = molsturm.System(["hydrogen", "O", 1], [[0, 0, 1], [0, 0, 0], [1, 0, 0]])
        self.assertListEqual(wat.atom_numbers.tolist(), [1, 8, 1])
        self.assertListEqual(wat.coords[0].tolist(), [0, 0, 1])
        self.assertListEqual(wat.coords[1].tolist(), [0, 0, 0])
        self.assertListEqual(wat.coords[2].tolist(), [1, 0, 0])
        self.assertEqual(wat.n_alpha, 5)
        self.assertEqual(wat.n_beta, 5)

        # Try to get a triplet oxygen wrongly
        with self.assertRaises(ValueError):
            molsturm.System(["O"], [0, 0, 0], (3, 5))

        # Missing the coord is ok if only a single atom
        ox = molsturm.System(["O"])
        self.assertListEqual(ox.coords[0].tolist(), [0, 0, 0])
        self.assertListEqual(self.oxide.atom_numbers.tolist(), [8])

        # But not if two or more
        with self.assertRaises(ValueError):
            molsturm.System(["O", "O"])

        # May miss the liss for one atom:
        ox = molsturm.System("O")
        self.assertListEqual(ox.coords[0].tolist(), [0, 0, 0])
        self.assertListEqual(self.oxide.atom_numbers.tolist(), [8])

    def test_1_atoms(self):
        water_list = [gint.element.by_symbol("H"), gint.element.by_symbol("O"),
                      gint.element.by_symbol("H")]

        self.assertListEqual(self.empty.atoms, [])
        self.assertListEqual(self.water.atoms, water_list)
        self.assertListEqual(self.water_triplet.atoms, water_list)
        self.assertListEqual(self.no.atoms, [gint.element.by_symbol("N"),
                                             gint.element.by_symbol("O")])
        self.assertListEqual(self.oxide.atoms, [gint.element.by_symbol("O")])

    def test_1_n_atoms(self):
        self.assertEqual(self.empty.n_atoms, 0)
        self.assertEqual(self.water.n_atoms, 3)
        self.assertEqual(self.water_triplet.n_atoms, 3)
        self.assertEqual(self.no.n_atoms, 2)
        self.assertEqual(self.oxide.n_atoms, 1)

    def test_1_empty(self):
        self.assertTrue(self.empty.empty)
        self.assertFalse(self.water.empty)
        self.assertFalse(self.water_triplet.empty)
        self.assertFalse(self.no.empty)
        self.assertFalse(self.oxide.empty)

    def test_1_charge(self):
        self.assertEqual(self.empty.charge, 0)
        self.assertEqual(self.water.charge, 0)
        self.assertEqual(self.water_triplet.charge, 0)
        self.assertEqual(self.no.charge, 0)
        self.assertEqual(self.oxide.charge, -2)

        self.assertEqual(self.empty.total_charge, 0)
        self.assertEqual(self.water.total_charge, 0)
        self.assertEqual(self.water_triplet.total_charge, 0)
        self.assertEqual(self.no.total_charge, 0)
        self.assertEqual(self.oxide.total_charge, -2)

    def test_1_n_electrons(self):
        self.assertEqual(self.empty.n_electrons, 0)
        self.assertEqual(self.water.n_electrons, 10)
        self.assertEqual(self.water_triplet.n_electrons, 10)
        self.assertEqual(self.no.n_electrons, 15)
        self.assertEqual(self.oxide.n_electrons, 10)

    def test_1_is_closed_shell(self):
        self.assertTrue(self.empty.is_closed_shell)
        self.assertTrue(self.water.is_closed_shell)
        self.assertFalse(self.water_triplet.is_closed_shell)
        self.assertFalse(self.no.is_closed_shell)
        self.assertTrue(self.oxide.is_closed_shell)

    def test_1_multiplicity(self):
        self.assertEqual(self.empty.multiplicity, 1)
        self.assertEqual(self.water.multiplicity, 1)
        self.assertEqual(self.water_triplet.multiplicity, 3)
        self.assertEqual(self.no.multiplicity, 2)
        self.assertEqual(self.oxide.multiplicity, 1)

    def test_2_adjust_electrons(self):
        oxy = molsturm.System("O")
        self.assertEqual(oxy.n_alpha, 4)
        self.assertEqual(oxy.n_beta, 4)
        self.assertEqual(oxy.multiplicity, 1)
        self.assertEqual(oxy.charge, 0)

        # Adjust multiplicity to make triplet
        oxy.adjust_electrons(multiplicity=3)
        self.assertEqual(oxy.n_alpha, 5)
        self.assertEqual(oxy.n_beta, 3)
        self.assertEqual(oxy.multiplicity, 3)
        self.assertEqual(oxy.charge, 0)

        # Adjust electrons to make doublet O^-
        oxy.adjust_electrons(charge=-1)
        self.assertEqual(oxy.n_alpha, 5)
        self.assertEqual(oxy.n_beta, 4)
        self.assertEqual(oxy.multiplicity, 2)
        self.assertEqual(oxy.charge, -1)

        # Adjust electrons to make quartet O^-
        oxy.adjust_electrons(multiplicity=4)
        self.assertEqual(oxy.n_alpha, 6)
        self.assertEqual(oxy.n_beta, 3)
        self.assertEqual(oxy.multiplicity, 4)
        self.assertEqual(oxy.charge, -1)

        # And go back to triplet oxygen
        oxy.adjust_electrons(charge=0, multiplicity=3)
        self.assertEqual(oxy.n_alpha, 5)
        self.assertEqual(oxy.n_beta, 3)
        self.assertEqual(oxy.multiplicity, 3)
        self.assertEqual(oxy.charge, 0)

        # Strip off all electrons:
        oxy.adjust_electrons(charge=8)
        self.assertEqual(oxy.n_alpha, 0)
        self.assertEqual(oxy.n_beta, 0)

        # Too large charge
        with self.assertRaises(ValueError):
            oxy.adjust_electrons(charge=9)

        # Invalid multiplicity number
        with self.assertRaises(ValueError):
            oxy.adjust_electrons(charge=0, multiplicity=2)
        with self.assertRaises(ValueError):
            oxy.adjust_electrons(charge=1, multiplicity=3)

        # All electrons in alpha:
        oxy.adjust_electrons(charge=0, multiplicity=9)
        self.assertEqual(oxy.n_alpha, 8)
        self.assertEqual(oxy.n_beta, 0)

        # Request larger multiplicity
        with self.assertRaises(ValueError):
            oxy.adjust_electrons(charge=0, multiplicity=11)

        # Change from triplet oxygen to triplet O^{2+}
        oxy.adjust_electrons(charge=0, multiplicity=3)
        oxy.adjust_electrons(charge=2)
        self.assertEqual(oxy.n_alpha, 4)
        self.assertEqual(oxy.n_beta, 2)
        self.assertEqual(oxy.multiplicity, 3)
        self.assertEqual(oxy.charge, 2)

    def test_2_charge_setter(self):
        oxy = molsturm.System("O")
        self.assertEqual(oxy.charge, 0)
        self.assertEqual(oxy.multiplicity, 1)

        oxy.charge = -2
        self.assertEqual(oxy.charge, -2)
        self.assertEqual(oxy.multiplicity, 1)

        oxy = molsturm.System("O", electrons=(5, 3))
        self.assertEqual(oxy.charge, 0)
        self.assertEqual(oxy.multiplicity, 3)

        oxy.charge = -2
        self.assertEqual(oxy.charge, -2)
        self.assertEqual(oxy.multiplicity, 3)

        oxy.charge = -1
        self.assertEqual(oxy.charge, -1)
        self.assertEqual(oxy.multiplicity, 2)

    def test_2_multiplicity_setter(self):
        oxy = molsturm.System("O")
        self.assertEqual(oxy.charge, 0)
        self.assertEqual(oxy.multiplicity, 1)

        oxy.multiplicity = 3
        self.assertEqual(oxy.charge, 0)
        self.assertEqual(oxy.multiplicity, 3)

        with self.assertRaises(ValueError):
            oxy.multiplicity = 2

    def test_2_n_electrons_setter(self):
        oxy = molsturm.System("O")
        self.assertEqual(oxy.charge, 0)
        self.assertEqual(oxy.multiplicity, 1)

        oxy.n_electrons += 2
        self.assertEqual(oxy.charge, -2)
        self.assertEqual(oxy.multiplicity, 1)

        oxy = molsturm.System("O", electrons=(5, 3))
        self.assertEqual(oxy.charge, 0)
        self.assertEqual(oxy.multiplicity, 3)

        oxy.n_electrons += 2
        self.assertEqual(oxy.charge, -2)
        self.assertEqual(oxy.multiplicity, 3)

        oxy.n_electrons = 9
        self.assertEqual(oxy.charge, -1)
        self.assertEqual(oxy.multiplicity, 2)
