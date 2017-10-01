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
from .ParameterMap import ParameterMap
from .MolecularSystem import MolecularSystem
import gint.gaussian
import gint.sturmian.atomic
from ._iface_conversion import ParamSpecial
import numpy as np
import warnings
import collections


def __assert_c_contiguous_numpy(arr, shape, dtype=None, prefix="numpy array"):
    if not isinstance(arr, np.ndarray):
        raise TypeError(prefix + " is not a numpy array, but a " + str(type(arr)))
    if arr.shape != shape:
        raise ValueError(prefix + " has not the expected shape " + str(shape) +
                         ", but " + str(arr.shape))
    if dtype:
        if arr.dtype != dtype:
            raise TypeError(prefix + " has not the expected dtype " + str(dtype) +
                            ", but " + str(arr.dtype))

    if not arr.flags["C_CONTIGUOUS"]:
        raise ValueError(prefix + " needs to be a C-contiguous numpy array.")


"""
A tuple with the following named components:
    n_bas     Number of basis functions
    n_spin    Number of spin components (1 for restricted, else 2)
    n_fock    Number of spatial orbitals the computation contains
              (equals n_eigenpairs // 2)
"""
ScfSizes = collections.namedtuple("ScfSizes", ["n_bas", "n_spin", "n_fock"])


class ScfParameters(ParameterMap):
    __size_t_keys = ["scf/max_iter", "scf/diis_size"]
    __special_keys = {
        # key : special type
        "guess/orben_f": "ignore",
        "guess/orbcoeff_bf": "ignore",
        "integrals/orbital_type": "orbital_type",
        "system/coords": "structure",
        "system/atom_numbers": "structure",
        "integrals/nlm_basis": "nlm_basis",
    }

    def __from_dict_inner(self, d, prefix):
        for k, v in d.items():
            fullkey = prefix + "/" + k

            if isinstance(v, dict):
                self.__from_dict_inner(v, fullkey)
            elif fullkey in self._special_keys:
                self[fullkey] = ParamSpecial(v, type=self.__special_keys[fullkey])
            elif isinstance(v, int):
                if k.startswith("n_") or fullkey in self.__size_t_keys:
                    self[fullkey] = np.uint64(v)
                else:
                    self[fullkey] = v
            elif isinstance(v, (str, bool, float)):
                self[fullkey] = v

    @classmethod
    def from_dict(cls, d):
        """
        Read the dict tree and construct
        a ParameterMap from it, then fill with
        default values, then check for validity.
        """
        params = cls()
        params.__from_dict_inner(d, "")
        params.fill_defaults()
        params.assert_valid()
        return params

    def __assert_system(self):
        """
        Check wether the system subtree is valid. Throws if not.
        """

        # TODO This could be made much simpler if the MolecularSystem
        #      class had a from_params function which constructed
        #      (and checked) the class from such a subtree.

        if "system" not in self:
            raise KeyError("No system subtree found in ScfParameters.")

        system = self["system"]
        if len(system["coords"].value) != len(system["atom_numbers"].value):
            raise ValueError("Length of the system/coords and length of "
                             "system/atom_numbers needs to agree.")
        n_atom = len(system["coords"].value)
        __assert_c_contiguous_numpy(system["coords"].value, (n_atom, 3),
                                    "system/coords")
        __assert_c_contiguous_numpy(system["atom_numbers"].value, (n_atom, 3),
                                    "system/atom_numbers")

        if "n_alpha" not in system or "n_beta" not in system:
            raise KeyError("system/n_alpha and system/n_beta need to be present.")
        if type(system["n_alpha"]) != np.uint64:
            raise TypeError("system/n_alpha needs to be pof type np.uint64")
        if type(system["n_beta"]) != np.uint64:
            raise TypeError("system/n_beta needs to be pof type np.uint64")

    def __normalise_integrals(self, system):
        # TODO This could be made much simpler if the Basis classes
        #      had a from_params function which constructed
        #      (and checked) the class from such a subtree.

        if "integrals" not in self:
            raise KeyError("No integrals subtree found in ScfParameters.")

        integrals = self["integrals"]
        if "basis_type" not in self:
            raise KeyError("No key basis_type found in integrals subtree.")

        if integrals["basis_type"] not in gint.available_basis_types:
            raise ValueError("Basis type " + integrals["basis_type"] +
                             " is not available.")

        # Get the backend string and the basis class type
        Basis, backend = split_basis_type(integrals["basis_type"])
        if isinstance(basis, gint.gaussian.Basis):
            if "basis_set_name" not in integrals:
                raise KeyError("Required key basis_set_name not found in integrals "
                               "subtree.")
            if type(integrals["basis_set_name"]) != str:
                raise TypeError("integrals/basis_set_name needs to be a str.")

            # TODO This is still some sort of legacy stuff we kind of need
            #      to do at the moment unfortunately. One should remove that soon.
            integrals.setdefault("orbital_type", ParamSpecial("real_molecular",
                                                              type="orbital_type"))
        elif isinstance(basis, gint.sturmian.atomic.Basis):
            if (system.n_atoms > 1):
                raise ValueError("Invalid basis: Atomic Coulomb-Strumians can only be"
                                 " used on atoms and not on molecules.")

            # TODO This is still some sort of legacy stuff we kind of need
            #      to do at the moment unfortunately. One should remove that soon.
            integrals.setdefault("orbital_type", ParamSpecial("complex_atomic",
                                                              type="orbital_type"))
            if "k_exp" not in integrals:
                raise KeyError("Required key k_exp not found in integrals subtree.")
            if type(integrals["k_exp"]) != float:
                raise TypeError("integrals/k_exp needs to be a float")

            # TODO This is only temporary and until the gint layer has fully moved
            #      to using nlm_basis.
            if "n_max" in integrals:
                integrals.setdefault("l_max") = integrals["n_max"] - 1
                integrals.setdefault("m_max") = integrals["l_max"]

                if "nlm_basis" in integrals:
                    warnings.warn("Overriding integrals/nlm_basis in ScfParameters")
                tmpbas = Basis(system, 1, integrals["n_max"], integrals["l_max"],
                               integrals["m_max"])
                integrals["nlm_basis"] = ParamSpecial(tmpbas.functions, type="nlm_basis")
            else:
                if "nlm_basis" not in integrals:
                    raise KeyError("Required key nlm_basis not found in integrals "
                                   "subtree.")

                n_max, l_max, m_max = np.max(integrals["nlm_basis"].value, axis=0)
                integrals["n_max"] = n_max
                integrals["l_max"] = l_max
                integrals["m_max"] = m_max

            if not isinstance(integrals["nlm_basis"].value, np.ndarray):
                raise TypeError("integrals/nlm_basis needs to be a numpy array")

            n_bas = len(integrals["nlm_basis"].value)
            __assert_c_contiguous_numpy(integrals["nlm_basis"].value, (n_bas, 3))
        else:
            raise AssertionError("Unrecognised Basis type")

    def __normalise_scf(self, system, basis):
        """
        Normalise and check the scf parameters.

        Return an ScfSizes object (since the __normalise_guess function needs this)
        """
        if "scf" not in self:
            raise KeyError("No scf subtree found in ScfParameters.")
        scf = self["scf"]

        # Convergence tolerance parameters:
        if "conv_tol" in scf:
            if "max_error_norm" in scf:
                warnings.warn("Overriding scf/nlm_basis in ScfParameters, since scf/conv_tol is given.")
            scf["max_error_norm"] = scf["conv_tol"]
            del scf["conv_tol"]
        scf.setdefault("max_error_norm", 5e-7)
        scf.setdefault("max_1e_energy_change", float(scf["max_error_norm"] / 100.))
        scf.setdefault("max_tot_energy_change", float(scf["max_error_norm"] / 4.))

        # Iteration control parameters:
        scf.setdefault("max_iter", np.uint64(25))
        scf.setdefault("diis_size", np.uint64(4))
        scf.setdefault("print_iterations", False)
        scf.setdefault("eigensolver/method", "auto")

        # Scf kind
        if "restricted" in scf:
            if "scf/kind" in scf:
                warnings.warn("Overriding scf/kind in ScfParameters, since "
                              "scf/restricted is given.")
            scf["kind"] = "restricted-open" if scf["restricted"] else "unrestricted"
            if system.is_closed_shell:
                param_tree["scf/kind"] = "restricted"
        scf.setdefault("scf/kind",
                       "restricted" if system.is_closed_shell else "unrestricted")

        if scf["kind"] == "restricted-open":
            raise NotImplementedError("restricted-open is currently not implemented.")
        n_spin = 1 if self["scf/kind"] == "restricted" else 2

        # Number of eigenpairs to obtain:
        # Our convention here is that we count both alpha and beta eigenpairs,
        # i.e. if this is 4 than for restricted calculations only 2 eigenpairs are found
        # and they are used for both alpha and beta.
        scf.setdefault("n_eigenpairs", 2 * basis.size)
        scf["n_eigenpairs"] = np.uint64(min(2 * basis.size, scf["n_eigenpairs"]))

        if scf["n_eigenpairs"] % 2 != 0:
            raise ValueError("The scf/n_eigenpairs parameter applies to the accumulated "
                             "number of eigenpairs in the SCF calculations, i.e. the "
                             "number of alpha plus the number of beta orbitals. This is "
                             "the done even for restricted calculations. For now we "
                             "further require this number to be even number.")

        # Check enough eigenpairs are requested.
        n_fock = scf["n_eigenpairs"] // 2
        if n_fock < max(system.n_alpha, system.n_beta):
            raise ValueError("Cannot treat a system with " +
                             str(system.n_alpha) + " alpha and " + str(system.n_beta) +
                             " beta electrons with computing only " + 
                             str(scf["n_eigenpairs"]) +
                             " eigenpairs. You need to use a large enough basis and "
                             "request enough eigenpairs (scf/n_eigenpairs parameter).")

        return ScfSizes(n_bas=basis.size, n_fock=n_fock, n_spin=n_spin)

    def __normalise_guess(self, system, basis, scf_sizes):
        self.setdefault("guess/method", "hcore")
        self.setdefault("guess/eigensolver/method", "auto")
        guess = self["guess"]

        if guess["method"] == "external":
            if "orben_f" not in guess:
                raise KeyError("For external guesses the key "
                               "guess/orben_f is required")
            if "orbcoeff_bf" not in guess:
                raise KeyError("For external guesses the key "
                               "guess/orbcoeff_bf is required")

            orben = param_tree["guess/orben_f"].value
            orbcoeff = param_tree["guess/orbcoeff_bf"].value

            __assert_c_contiguous_numpy(orben, (scf_sizes.n_spin, scf_sizes.n_fock), "guess/orben_f")
            __assert_c_contiguous_numpy(orbceoff, (scf_sizes.n_spin, scf_sizes.n_bas, scf_sizes.n_fock),
                                        "guess/orbcoeff_bf")

    def normalise(self):
        """
        Look at the current set of values in the ScfParameters object
        and try to normalise them, i.e. fill in default values where appropriate
        and replace meta-parameters (which set a couple of different parameters at once)
        by their actual primitive values.

        After the function call the parameters are either in a valid state
        or the function will throw KeyError, ValueError or TypeError exceptions
        to indicate that this is not the case.
        """
        self.__assert_system()
        system = self.system

        self.__normalise_integrals(system)
        basis = self.basis

        sizes = self.__normalise_scf(system, basis)
        self.__normalise_guess(system, basis, sizes)

    def set_guess_external(self, orben_f, orbcoeff_bf):
        """
        Convenience function to set an external guess by supplying
        a guess for the orbital energies and the orbital coefficients.

        This function only works proprly after the relevant integrals
        and scf parameters have been set up and will throw errors
        if this has not been done yet.

        The guess is expected in the format
           orben_f             (n_spin, n_fock)
           guess/orbcoeff_bf   (n_spin, n_bas, n_fock)
        i.e. as blocks for the individual spin components. For restricted
        we expect n_spin == 1 and for unrestricted n_spin == 2
        """
        param_tree["guess/method"] = "external"
        orben_f = np.ascontiguousarray(orben_f)
        orbcoeff_bf = np.ascontiguousarray(orbcoeff_bf)
        param_tree["guess/orben_f"] = ParamSpecial(orben_f, type="ignore")
        param_tree["guess/orbcoeff_bf"] = ParamSpecial(orbcoeff_bf, type="ignore")

        # Check and normalise what we have:
        self.normalise()

    @property
    def scf_sizes(self):
        """
        Return a tuple with the following named components:
            n_bas     Number of basis functions
            n_spin    Number of spin components (1 for restricted, else 2)
            n_fock    Number of spatial orbitals the computation contains
                      (equals n_eigenpairs // 2)
        """
        assert self["scf/kind"] in ["restricted", "unrestricted"]
        assert self["scf/n_eigenpairs"] % 2 == 0
        return Sizes(self.basis.size, 1 if self["scf/kind"] == "restricted" else 2,
                     self["scf/n_eigenpairs"] // 2)


    @property
    def system(self):
        """
        Return the MolecularSystem which is represented by the
        internal parameters.
        """
        self.__assert_tree_system()
        return MolecularSystem(
            atoms=self["system/atom_numbers"].value,
            coords=self["system/coords"].value,
            electrons=(self["system/n_alpha"], self["system/n_beta"])
        )

    @system.setter
    def system(self, system):
        """Set the molecular system parameters"""
        self["system/n_alpha"] = np.uint64(system.n_alpha)
        self["system/n_beta"] = np.uint64(system.n_beta)
        self["system/coords"] = ParamSpecial(system.coords, "structure")
        self["system/atom_numbers"] = ParamSpecial(system.atom_numbers, "structure")

    @property
    def basis(self):
        """
        Return the basis object represented by the internal parameters.
        """
        self.__normalise_integrals(self.system)

        # Get the backend string and the basis class type
        integrals = self["integrals"]
        Basis, backend = split_basis_type(integrals["basis_type"])

        # TODO This explicit branching is a bit ugly
        if Basis == gint.gaussian.Basis:
            return gint.gaussian.Basis(self.system(), integrals["basis_set_name"],
                                       backend=backend)
        elif Basis == gint.sturmian.atomic.Basis:
            warnings.warn("ScfParameters.basis assumes that the atomic Coulomb-Sturmian"
                         " basis is dense and in nlm order.")
            n_max, l_max, m_max = np.max(integrals["nlm_basis"].value, axis=0)

            return gint.sturmian.atomic.Basis(self.system(), integrals["k_exp"],
                                              integrals["n_max"], integrals["l_max"],
                                              integrals["m_max"])
        else:
            raise AssertionError("Unrecognised Basis type")

    @basis.setter
    def basis(self, basis):
        """
        Set the basis into the parameters.
        """
        self["integrals/basis_type"] = str(basis.basis_type)
        if isinstance(basis, gint.gaussian.Basis):
            # TODO Instead set basis functions here directly and omit passing
            #      the basis set name here and allow to pass the full description
            #      down to gint
            self["integrals/basis_set_name"] = str(basis.basis_set_name)

            # TODO This is still some sort of legacy stuff we kind of need
            #      to do at the moment unfortunately. One should remove that soon.
            self["integrals/orbital_type"] = ParamSpecial("real_molecular",
                                                          type="orbital_type")
        elif isinstance(basis, gint.sturmian.atomic.Basis):
            if (self.system.n_atoms > 1):
                raise ValueError("Invalid basis: Atomic Coulomb-Strumians can only be"
                                 " used on atoms and not on molecules.")

            # TODO This is still some sort of legacy stuff we kind of need
            #      to do at the moment unfortunately. One should remove that soon.
            self["integrals/orbital_type"] = ParamSpecial("complex_atomic",
                                                          type="orbital_type")

            self["integrals/k_exp"] = float(basis.k_exp)
            self["integrals/nlm_basis"] = ParamSpecial(basis.functions, "nlm_basis")

            warnings.warn("Right now we assume all sturmian basis sets to be dense.")
            n_max, l_max, m_max = np.max(basis.functions, axis=0)
            self["integrals/n_max"] = int(n_max)
            self["integrals/l_max"] = int(l_max)
            self["integrals/m_max"] = int(m_max)
        else:
            raise TypeError("basis has an unrecognised type.")
