//
// Copyright (C) 2017 by the molsturm authors
//
// This file is part of molsturm.
//
// molsturm is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// molsturm is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with molsturm. If not, see <http://www.gnu.org/licenses/>.
//

#pragma once
#include <gint/Integral.hh>
#include <lazyten/Constants.hh>

namespace molsturm {

// TODO Better design:
//      - Every Fock matrix term does not only know how to apply itself, but also
//        how to compute its energy. A default implementation for this would be
//        tr(C^T Term C)

template <typename StoredMatrix>
struct IntegralTermContainer {
  typedef StoredMatrix stored_matrix_type;
  typedef typename stored_matrix_type::scalar_type scalar_type;

  //! The type to use for the integral terms.
  typedef typename gint::Integral<stored_matrix_type> int_term_type;

  /** The one electron terms */
  std::vector<int_term_type> integral_terms_1e;

  /** Set of coefficients for the one electron integral terms */
  std::vector<scalar_type> coefficients_1e;

  /** The coulomb term to use */
  int_term_type coulomb_term;

  /** The coefficient for the coulomb term */
  scalar_type coefficient_coulomb;

  /** The exchange term to use */
  int_term_type exchange_term;

  /** The coefficient for the HF exchange term */
  scalar_type coefficient_exchange;

  /** \brief Setup the class with the given integral terms.
   *  \note By default all coefficients, but the exchange coefficient, are set
   * to one. The exchange term coefficient is set to -1.
   *  This is done in order to make sure that by default this data structure is
   *  ready for a Hartree-Fock calculation.
   */
  IntegralTermContainer(std::vector<int_term_type> integral_terms_1e_,
                        int_term_type coulomb_term_, int_term_type exchange_term_)
        : integral_terms_1e{std::move(integral_terms_1e_)},
          coefficients_1e{std::vector<scalar_type>(integral_terms_1e.size(),
                                                   lazyten::Constants<scalar_type>::one)},
          coulomb_term{std::move(coulomb_term_)},
          coefficient_coulomb{lazyten::Constants<scalar_type>::one},
          exchange_term{std::move(exchange_term_)},
          coefficient_exchange{-lazyten::Constants<scalar_type>::one} {}
};

}  // namespace molsturm
