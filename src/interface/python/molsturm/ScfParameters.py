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


import gint.element
from gint import split_basis_type
from .ParameterMap import ParameterMap
from .System import System
import gint.gaussian
import gint.sturmian.atomic
from ._iface_conversion import ParamSpecial
import numpy as np
import warnings
import collections


"""
A tuple with the following named components:
    n_bas     Number of basis functions
    n_spin    Number of spin components (1 for restricted, else 2)
    n_fock    Number of spatial orbitals the computation contains
              (equals n_eigenpairs // 2)
"""
ScfSizes = collections.namedtuple("ScfSizes", ["n_bas", "n_spin", "n_fock"])


class ScfParameters(ParameterMap):

    """The keys which do not start with "n_" and
       which still need to be of type size_t
    """
    __size_t_keys = ["scf/max_iter", "scf/diis_size", "scf/diis_startup_iter"]

    """The keys which need a ParamSpecial wrapper"""
    __special_keys = {
        # key : special type
        "guess/orben_f": "ignore",
        "guess/orbcoeff_bf": "ignore",
        "discretisation/orbital_type": "orbital_type",
        "system/coords": "structure",
        "system/atom_numbers": "structure",
        "system/atoms": "ignore",
        "discretisation/nlm_basis": "nlm_basis",
    }

    def __from_dict_inner(self, d, prefix):
        inttypes = (int,) + tuple(np.sctypes["int"] + np.sctypes["uint"])
        scalartypes = (bool, float)

        for k, v in d.items():
            fullkey = prefix + "/" + k if len(prefix) > 0 else k

            if isinstance(v, dict):
                self.__from_dict_inner(v, fullkey)
            elif isinstance(v, ParamSpecial):
                self[fullkey] = v
            elif fullkey in self.__special_keys:
                self[fullkey] = ParamSpecial(v, type=self.__special_keys[fullkey])
            elif isinstance(v, inttypes):
                if k.startswith("n_") or fullkey in self.__size_t_keys:
                    self[fullkey] = np.uint64(v)
                else:
                    self[fullkey] = int(v)
            elif isinstance(v, (str,) + scalartypes):
                self[fullkey] = v
            else:
                raise TypeError("Unknown type " + str(type(v)) + " for key " + fullkey)

    @classmethod
    def from_dict(cls, d):
        """
        Read the dict tree and construct a ParameterMap from it.

        The input parameters are not checked for valididy.
        In order to do so, call the function "normalise"
        after constructing the object.
        """
        params = ScfParameters()
        params.__from_dict_inner(d, "")
        return params

    def __copy__(self):
        """
        Make an exact shallow copy of the ScfParameters
        """
        return ScfParameters.from_dict(self)

    def copy(self):
        return self.__copy__()

    @classmethod
    def from_args(cls, system, basis, conv_tol=None, max_iter=None, n_eigenpairs=None,
                  restricted=None, eigensolver=None, print_iterations=None):
        """
        Construct an ScfParameter set from some arguments.
        The resulting structure will not be normalised and not checked for validity.
        You can enforce this by calling normalise() on the returned object.

        system        The molecular system to model
        basis         A valid basis object. If None the basis will be constructed
                      on the fly from teh basis_type and the kwargs.
        basis_type    String describing the type of basis function to use.
        conv_tol      SCF convergence tolerance
        max_iter      Maximum number of SCF iterations
        n_eigenpairs  Number of orbitals to compute
        eigensolver   SCF eigensolver to use
        print_iterations  Shall some diagnostics about the SCF iterations be printed
        restricted    Shall an restricted or an unrestricted SCF be run
        """

        if isinstance(system, str):
            system = System(system)
        elif not isinstance(system, System):
            raise TypeError("The first argument needs to be a System object or "
                            "a string of exactly one atom.")

        # Build a parameter tree from the commandline arguments
        params = ScfParameters()
        params.system = system
        params.basis = basis

        # Add the scf parameters
        for frm, to, typ in [
            ("conv_tol", "scf/conv_tol", float),
            ("max_iter", "scf/max_iter", np.uint64),
            ("n_eigenpairs", "scf/n_eigenpairs", np.uint64),
            ("print_iterations", "scf/print_iterations", bool),
            ("eigensolver", "scf/eigensolver/method", str),
            ("restricted", "scf/restricted", bool),
        ]:
            val = locals()[frm]
            if val is not None:
                params[to] = typ(val)

        return params

    def __normalise_numpy_array(self, key, shape, dtype=None):
        try:
            # The .value is around, since all keys going through
            # this point to ParameterSpecial objects.
            arr = np.array(self[key].value, copy=False, dtype=dtype)
        except TypeError as e:
            raise TypeError(key + " is not a numpy array and casting returned an "
                            "error: " + str(e))
        except ValueError as e:
            raise ValueError(key + " cannot be uses as a numpy array: " + str(e))

        if arr.shape != shape:
            raise ValueError(key + " has not the expected shape " + str(shape) +
                             ", but " + str(arr.shape))

        if dtype:
            if arr.dtype != dtype:
                raise TypeError(key + " has not the expected dtype " + str(dtype) +
                                ", but " + str(arr.dtype))

        if not arr.flags["C_CONTIGUOUS"]:
            raise ValueError(key + " needs to be a C-contiguous numpy array.")

        self[key].value = arr

    def __normalise_system(self):
        """
        Check wether the system subtree is valid. Throws if not.

        Transforms these keys:
            system/atoms -> system/atom_numbers
            system/charge -> system/n_alpha, system/n_beta
            system/multipliciy -> system/n_alpha, system/n_beta
        """
        if "system" not in self:
            raise KeyError("No system subtree found in ScfParameters.")

        system = self["system"]

        # TODO This could be made much simpler if the System
        #      class had a from_params function which constructed
        #      (and checked) the class from such a subtree.

        if "atoms" in system:
            if "atom_numbers" in system:
                warnings.warn("Overriding system/atom_numbers in ScfParameters since "
                              "system/atoms is present.")
            atnums = np.array([gint.element.by_anything(at).atom_number
                               for at in system["atoms"].value], dtype=int)
            system["atom_numbers"] = ParamSpecial(atnums, "structure")
            del system["atoms"]

        if len(system["coords"].value) != len(system["atom_numbers"].value):
            raise ValueError("Length of the system/coords and length of "
                             "system/atom_numbers needs to agree.")

        n_atom = len(system["coords"].value)
        self.__normalise_numpy_array("system/coords", (n_atom, 3), dtype=float)
        self.__normalise_numpy_array("system/atom_numbers", (n_atom,), dtype=int)

        if "multiplicity" in system or "charge" in system or \
           ("n_alpha" not in system and "n_beta" not in system):
            sysobj = System(atoms=system["atom_numbers"].value,
                            coords=system["coords"].value)

            try:
                sysobj.adjust_electrons(charge=system.pop("charge", None),
                                        multiplicity=system.pop("multiplicity", None))
            except ValueError as e:
                raise ValueError("system/multiplicity or system/charge erroneous: " +
                                 str(e))

            if "n_alpha" in system or "n_beta" in system:
                warnings.warn("Overriding system/n_alpha and system/n_beta in "
                              "ScfParameters, since system/multiplicity or "
                              "system/charge is present.")

            system["n_alpha"] = np.uint64(sysobj.n_alpha)
            system["n_beta"] = np.uint64(sysobj.n_beta)

        if "n_alpha" not in system or "n_beta" not in system:
            raise KeyError("system/n_alpha and system/n_beta or alternatively "
                           "system/multiplicity need to be present.")
        if type(system["n_alpha"]) != np.uint64:
            raise TypeError("system/n_alpha needs to be pof type np.uint64")
        if type(system["n_beta"]) != np.uint64:
            raise TypeError("system/n_beta needs to be pof type np.uint64")

    def __normalise_discretisation(self, system):
        # TODO This could be made much simpler if the Basis classes
        #      had a from_params function which constructed
        #      (and checked) the class from such a subtree.

        if "discretisation" not in self:
            raise KeyError("No discretisation subtree found in ScfParameters.")

        discretisation = self["discretisation"]
        if "basis_type" not in discretisation:
            raise KeyError("No key basis_type found in discretisation subtree.")

        if discretisation["basis_type"] not in gint.available_basis_types:
            raise ValueError("Basis type " + discretisation["basis_type"] +
                             " is not available.")

        # Get the backend string and the basis class type
        Basis, backend = split_basis_type(discretisation["basis_type"])
        if Basis == gint.gaussian.Basis:
            if "basis_set_name" not in discretisation:
                raise KeyError("Required key basis_set_name not found in discretisation "
                               "subtree.")
            if type(discretisation["basis_set_name"]) != str:
                raise TypeError("discretisation/basis_set_name needs to be a str.")

            # TODO This is still some sort of legacy stuff we kind of need
            #      to do at the moment unfortunately. One should remove that soon.
            discretisation.setdefault("orbital_type",
                                      ParamSpecial("real_molecular", type="orbital_type"))
        elif Basis == gint.sturmian.atomic.Basis:
            if (system.n_atoms > 1):
                raise ValueError("Invalid basis: Atomic Coulomb-Strumians can only be"
                                 " used on atoms and not on molecules.")

            # TODO This is still some sort of legacy stuff we kind of need
            #      to do at the moment unfortunately. One should remove that soon.
            discretisation.setdefault("orbital_type",
                                      ParamSpecial("complex_atomic", type="orbital_type"))
            if "k_exp" not in discretisation:
                raise KeyError("Required key k_exp not found in discretisation subtree.")
            if type(discretisation["k_exp"]) != float:
                raise TypeError("discretisation/k_exp needs to be a float")

            # TODO This is only temporary and until the gint layer has fully moved
            #      to using nlm_basis.
            if "n_max" in discretisation:
                discretisation["n_max"] = int(discretisation["n_max"])
                discretisation.setdefault("l_max", discretisation["n_max"] - 1)
                discretisation.setdefault("m_max", discretisation["l_max"])

                if "nlm_basis" in discretisation:
                    warnings.warn("Overriding discretisation/nlm_basis in ScfParameters")
                tmpbas = Basis(system, k_exp=1, n_max=discretisation["n_max"],
                               l_max=discretisation["l_max"],
                               m_max=discretisation["m_max"])
                discretisation["nlm_basis"] = ParamSpecial(tmpbas.functions,
                                                           type="nlm_basis")
            else:
                if "nlm_basis" not in discretisation:
                    raise KeyError("Required key discretisation/nlm_basis not found in "
                                   "ScfParameters.")

                n_max, l_max, m_max = np.max(discretisation["nlm_basis"].value, axis=0)
                discretisation["n_max"] = n_max
                discretisation["l_max"] = l_max
                discretisation["m_max"] = m_max

            if not isinstance(discretisation["nlm_basis"].value, np.ndarray):
                raise TypeError("discretisation/nlm_basis needs to be a numpy array")

            n_bas = len(discretisation["nlm_basis"].value)
            self.__normalise_numpy_array("discretisation/nlm_basis",
                                         (n_bas, 3), dtype=int)
        else:
            raise AssertionError("Unrecognised Basis type")

    def __normalise_scf(self, system, basis):
        """
        Normalise and check the scf parameters.

        Return an ScfSizes object (since the __normalise_guess function needs this)

        Transforms these keys:
            scf/conv_tol -> scf/max_error_norm, scf/max_1e_energy_change,
                            scf/max_tot_energy_change
            scf/restricted -> scf/kind
        """
        self.setdefault("scf/eigensolver/method", "auto")
        scf = self["scf"]

        # Convergence tolerance parameters:
        if "conv_tol" in scf:
            if "max_error_norm" in scf:
                warnings.warn("Overriding scf/nlm_basis in ScfParameters, "
                              "since scf/conv_tol is given.")
            scf["max_error_norm"] = float(scf["conv_tol"])
            del scf["conv_tol"]
        scf.setdefault("max_error_norm", 5e-7)
        scf.setdefault("max_1e_energy_change", float(scf["max_error_norm"] * 100.))
        scf.setdefault("max_tot_energy_change", float(scf["max_error_norm"] / 4.))

        # Iteration control parameters:
        scf.setdefault("diis_size", np.uint64(4))
        scf.setdefault("max_iter", np.uint64(25))
        scf.setdefault("diis_startup_iter", np.uint64(12))
        scf.setdefault("diis_startup_error_norm", 0.25)
        scf.setdefault("diis_shutdown_error_norm", 5e-7)

        # I have no clue why, but sometimes the type information is lost here.
        scf.setdefault("print_iterations", False)
        scf["print_iterations"] = bool(scf["print_iterations"])

        # Scf kind
        if "restricted" in scf:
            if "kind" in scf:
                warnings.warn("Overriding scf/kind in ScfParameters, since "
                              "scf/restricted is given.")
            scf["kind"] = "restricted-open" if scf["restricted"] else "unrestricted"
            if system.is_closed_shell and scf["restricted"]:
                scf["kind"] = "restricted"
            del scf["restricted"]
        scf.setdefault("kind", "restricted" if system.is_closed_shell else "unrestricted")

        if scf["kind"] == "restricted-open":
            raise NotImplementedError("restricted-open is currently not implemented.")
        n_spin = 1 if scf["kind"] == "restricted" else 2

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
        n_fock = int(scf["n_eigenpairs"] // 2)
        if n_fock < max(system.n_alpha, system.n_beta):
            raise ValueError("Cannot treat a system with " +
                             str(system.n_alpha) + " alpha and " + str(system.n_beta) +
                             " beta electrons with computing only " +
                             str(scf["n_eigenpairs"]) +
                             " eigenpairs. You need to use a large enough basis and "
                             "request enough eigenpairs (scf/n_eigenpairs parameter).")

        return ScfSizes(n_bas=basis.size, n_fock=n_fock, n_spin=n_spin)

    def __normalise_guess(self, system, basis, scf_sizes):
        n_bas, n_spin, n_fock = scf_sizes

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

            self.__normalise_numpy_array("guess/orben_f", (n_spin, n_fock), dtype=float)
            self.__normalise_numpy_array("guess/orbcoeff_bf", (n_spin, n_bas, n_fock),
                                         dtype=float)
        elif guess["method"] == "hcore":
            guess.setdefault("nuclear_attraction_factor", 1.0)

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
        self.__normalise_system()
        system = self.system

        self.__normalise_discretisation(system)
        basis = self.basis

        sizes = self.__normalise_scf(system, basis)
        self.__normalise_guess(system, basis, sizes)

    def set_guess_external(self, orben_f, orbcoeff_bf):
        """
        Convenience function to set an external guess by supplying
        a guess for the orbital energies and the orbital coefficients.

        No checking wether the guess is consistent with the other
        internal parameters is performed. This can be done via a call
        to normalise() after setting up the guess like this.

        The guess is expected in the format
           orben_f             (n_spin, n_fock)
           guess/orbcoeff_bf   (n_spin, n_bas, n_fock)
        i.e. as blocks for the individual spin components. For restricted
        we expect n_spin == 1 and for unrestricted n_spin == 2
        """
        # Clear guess parameters first:
        self.pop("guess", None)

        self["guess/method"] = "external"
        orben_f = np.ascontiguousarray(orben_f)
        orbcoeff_bf = np.ascontiguousarray(orbcoeff_bf)
        self["guess/orben_f"] = ParamSpecial(orben_f, type="ignore")
        self["guess/orbcoeff_bf"] = ParamSpecial(orbcoeff_bf, type="ignore")

    def clear_guess(self):
        """
        Clear any guess parameters, which are currently set.

        This is deprecated. Plainly use

        >>> del params["guess"]

        instead.
        """
        from warnings import warn
        warn("Simply use 'del params[\"guess\"]' instead of this function. "
             "This function will be removed after the next release.",
             DeprecationWarning)
        self.pop("guess", None)

    @property
    def scf_sizes(self):
        """
        Return a tuple with the following named components:
            n_bas     Number of basis functions
            n_spin    Number of spin components (1 for restricted, else 2)
            n_fock    Number of spatial orbitals the computation contains
                      (equals n_eigenpairs // 2)

        This function only works properly if the system and scf parameters
        have been properly set up and throw KeyError or ValueError in case
        things are missing or invalid.
        """
        if self["scf/kind"] not in ["restricted", "unrestricted"]:
            raise ValueError("Invalid value for scf/kind: " + self["scf/kind"])
        if self["scf/n_eigenpairs"] % 2 != 0:
            raise ValueError("scf/n_eigenpairs needs to be an even number.")
        return ScfSizes(self.basis.size, 1 if self["scf/kind"] == "restricted" else 2,
                        int(self["scf/n_eigenpairs"] // 2))

    @property
    def system(self):
        """
        Normalise the system parameter subtree, then return the system
        object represented by the internal parameters.
        """
        self.__normalise_system()
        return System(
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
        Normalise the system and basis subtrees, then return the basis
        object represented by the internal parameters.
        """
        self.__normalise_discretisation(self.system)

        # Get the backend string and the basis class type
        discretisation = self["discretisation"]
        Basis, backend = split_basis_type(discretisation["basis_type"])

        # TODO This explicit branching is a bit ugly
        if Basis == gint.gaussian.Basis:
            return gint.gaussian.Basis(self.system, discretisation["basis_set_name"],
                                       backend=backend)
        elif Basis == gint.sturmian.atomic.Basis:
            warnings.warn("ScfParameters.basis assumes that the atomic Coulomb-Sturmian"
                          " basis is dense and in nlm order.")
            n_max, l_max, m_max = np.max(discretisation["nlm_basis"].value, axis=0)

            return gint.sturmian.atomic.Basis(self.system, discretisation["k_exp"],
                                              discretisation["n_max"],
                                              discretisation["l_max"],
                                              discretisation["m_max"], backend=backend)
        else:
            raise AssertionError("Unrecognised Basis type")

    @basis.setter
    def basis(self, basis):
        """
        Set the basis into the parameters.
        """
        self["discretisation/basis_type"] = str(basis.basis_type)
        if isinstance(basis, gint.gaussian.Basis):
            # TODO Instead set basis functions here directly and omit passing
            #      the basis set name here and allow to pass the full description
            #      down to gint
            self["discretisation/basis_set_name"] = str(basis.basis_set_name)

            # TODO This is still some sort of legacy stuff we kind of need
            #      to do at the moment unfortunately. One should remove that soon.
            self["discretisation/orbital_type"] = ParamSpecial("real_molecular",
                                                               type="orbital_type")
        elif isinstance(basis, gint.sturmian.atomic.Basis):
            if (self.system.n_atoms > 1):
                raise ValueError("Invalid basis: Atomic Coulomb-Strumians can only be"
                                 " used on atoms and not on molecules.")

            # TODO This is still some sort of legacy stuff we kind of need
            #      to do at the moment unfortunately. One should remove that soon.
            self["discretisation/orbital_type"] = ParamSpecial("complex_atomic",
                                                               type="orbital_type")

            self["discretisation/k_exp"] = float(basis.k_exp)
            self["discretisation/nlm_basis"] = ParamSpecial(basis.functions, "nlm_basis")

            warnings.warn("Right now we assume all sturmian basis sets to be dense.")
            n_max, l_max, m_max = np.max(basis.functions, axis=0)
            self["discretisation/n_max"] = int(n_max)
            self["discretisation/l_max"] = int(l_max)
            self["discretisation/m_max"] = int(m_max)
        else:
            raise TypeError("basis has an unrecognised type.")

    def __strip_special(self, dicti):
        for k, v in dicti.items():
            if isinstance(v, dict):
                dicti[k] = self.__strip_special(v)
            if isinstance(v, ParamSpecial):
                dicti[k] = v.value
        return dicti

    def to_dict(self, strip_special=True):
        """
        Return a dict of dicts which represents the same data

        strip_special  Should the ParamSpecial instances be stripped off
                       when returning the values (default: True)
        """
        if strip_special:
            return self.__strip_special(super().to_dict())
        else:
            return super().to_dict()
