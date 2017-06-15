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

# This files include constants which are needed all over the code

"""The keys corresponding to numpy arrays in the HfResults class"""
HFRES_ARRAY_KEYS = [ "orben_f", "eri_ffff", "fock_ff", "orbcoeff_fb", "hcore_ff" ]

"""Keys which could be numpy arrays in the input Parameters to hartree_fock"""
INPUT_PARAMETER_ARRAY_KEYS = [ "atom_numbers", "atoms", "coords", "nlm_basis" ]

"""The key used to indicate the input parameters in the dict returned by hartree_fock"""
HFRES_INPUT_PARAMETER_KEY= "input_parameters"
