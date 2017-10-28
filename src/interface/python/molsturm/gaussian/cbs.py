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

from scipy.optimize import curve_fit
import gint
import molsturm
import numpy as np
import re


def fit_cbs_curve(orders, quantities):
    """
    Performs the CBS fitting explained in

    Jensen, F. "Estimating the HF libit from finite basis set calculations"
    Theor. Chem. Acc (2005), 113: 267-273

    i.e. it tries to fit the function CBS + A * exp(-B * sqrt(order)).

    in order to estimate a quantity at the HF limit.
    The zeta orders of the basis should be given in the numpy array
    orders and the quantities for each order should be specified in
    quantities.

    Typical use cases are the energies of a HF calculation.

    Returns the popt and pcov results from the fit, i.e.
    popt, pcov = fit_cbs(orders, quantities) followed by a
    popt[0] returns the fitted CBS limit.
    """

    if len(orders) != len(quantities):
        raise ValueError("Length of the orders and the quantities array "
                         "need to be the same.")
    if len(orders) < 3:
        raise ValueError("Need at least three points for the CBS fit.")

    # Minimum of the quantities
    shift = np.min(quantities)

    def functional(x, E, A, B):
        return shift + E + A * np.exp(-B * np.sqrt(x))

    try:
        popt, pcov = curve_fit(functional, orders, quantities)
    except RuntimeError as e:
        raise RuntimeError("CBS curve fit failed, more info: " + str(e))

    # Correct shift:
    popt[0] += shift
    return popt, pcov


class CBSExtrapolationResult:
    """
    Class to capture the results of a extrapolate_cbs_limit call.
    The parameters correspond to the fitted curve
    CBS + A * exp(-B * sqrt(order))

    cbs   The extrapolated CBS limit
    A     The parameter A of the CBS curve
    B     The parameter B of the CBS curve
    orders  The orders used for the fit
    quantities   The quantities used to fit the cbs value
    states       The states of all calculations done.
                 Roughly speaking states[0] corresponds to orders[0]
    success   Has the fit been successful. If this is false
              the values of A, B, cbs are not sensible and
              extrapolation_function should not be used.
    """

    def __init__(self, orders, quantities, popt, pcov, states, success):
        if popt is not None and pcov is not None:
            self.cbs, self.A, self.B = popt
            self.cov = pcov
            self.success = success
        else:
            self.cbs = None
            self.A = None
            self.B = None
            self.cov = None
            self.success = False
        self.orders = orders
        self.quantities = quantities
        self.states = states

    def extrapolation_function(self, x):
        """
        Evaluate the extrapolated function for the CBS limit at a particular
        order x
        """
        if not self.success:
            raise ValueError("CBSExtrapolationResult does not contain a "
                             "successful cbsfit.")
        return self.cbs + self.A * np.exp(-self.B * np.sqrt(x))


def extrapolate_cbs_limit(scfparams, n_points=3, quantity="hfenergy"):
    """
    Perform an automatic CBS extrapolation for a particular quantity using
    the basis set family specified in the ScfParameters object as well
    as n_points points on the CBS surface, where the first point of the
    extrapolation is exactly the basis set specified in the scfparams
    object.

    Example: scfparams.basis is a cc-pVDZ basis. n_points=3.
    The quantity will be computed for cc-pVDZ, cc-pVTZ and cc-pVQZ
    and then extrapolated to the CBS.

    quantity may either be a string describing the quantity to compute
    or a function to take a state and return the quantity to extrapolate.

    If the extrapolation fails the function returns a CBSExtrapolationResult
    where the success flag is set to False. One may then retrieve
    the states of the performed calculations and the derived quantities
    for the fit, but the values of the CBS limit and the extrapolation_function
    are not reliable.
    """
    scfparams = scfparams.copy()
    origbasis = scfparams.basis
    if not isinstance(origbasis, gint.gaussian.Basis):
        raise ValueError("The basis of the scfparams needs to be a gaussian basis set.")

    if isinstance(quantity, str):
        if quantity == "hfenergy":
            def quantity(x):
                return x["energy_ground_state"]
        else:
            raise NotImplementedError("The quantity '" + quantity +
                                      "' is not supported at the moment. "
                                      "Currently we only support 'hfenergy'")

    #
    # Translate between zetaness as number and as letter
    #
    # Translation map for letter to zetaness of the basis set
    letter_to_zeta = {"D": 2, "T": 3, "Q": 4, "P": 5, "H": 6, }

    # And function to translate back in the conventional way:
    def zeta_to_string(zeta):
        if zeta <= 4:
            return [k for k, v in letter_to_zeta.items() if v == zeta][0]
        else:
            return str(zeta)
    #
    # Match for Dunning basis sets
    #
    dunning_re = re.compile(r"""
        ^(?P<pre>cc-pV)           # Prefix of the basis set
        (?P<zeta>[0-9]+|[DTQPH])  # Zetaness either as string or as letter
        (?P<post>Z)$              # Postfix
    """, flags=re.IGNORECASE | re.VERBOSE)

    is_dunning = dunning_re.match(origbasis.basis_set_name)
    if is_dunning:
        pre = is_dunning.group("pre")
        zetaness = letter_to_zeta[is_dunning.group("zeta").upper()]
        post = is_dunning.group("post")
    else:
        raise NotImplementedError("Currently automatic CBS limit extrapolation is only "
                                  "implemented for plain Dunning basis sets.")

    # Construct containers and basis sets
    orders = np.arange(zetaness, zetaness + n_points)

    basis_sets = []
    for i, zeta in enumerate(orders):
        basis_set_name = pre + zeta_to_string(zeta) + post

        try:
            basis = gint.gaussian.Basis(scfparams.system,
                                        basis_set_name=basis_set_name,
                                        backend=origbasis.backend)
        except RuntimeError as e:
            raise ValueError("Could not construct the basis set '" + basis_set_name +
                             "'. Check your input parameters, especially n_points. "
                             "Further info: " + str(e))
        basis_sets.append(basis)

    states = []
    quantities = []
    for i, basis in enumerate(basis_sets):
        # TODO One could very likely improve convergence if one took the
        #      just computed hf result of the previous step as a guess for
        #      the next. Right now we cannot do this, since the interpolation
        #      for this has not yet been implemented, but for the future ...

        scfparams.basis = basis
        res = molsturm.self_consistent_field(scfparams)
        states.append(res)
        quantities.append(quantity(res))
    quantities = np.array(quantities)

    try:
        popt, pcov = fit_cbs_curve(orders, quantities)
    except RuntimeError as e:
        return CBSExtrapolationResult(orders, quantities, None, None,
                                      states, success=False)

    return CBSExtrapolationResult(orders, quantities, popt, pcov, states, success=True)
