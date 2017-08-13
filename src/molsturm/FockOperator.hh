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
#include "RestrictedClosedIntegralOperator.hh"
#include "RestrictionType.hh"
#include "UnrestrictedIntegralOperator.hh"
#include <gint/IntegralLookup.hh>
#include <gint/IntegralType.hh>
#include <krims/TypeUtils/UsingLibrary.hh>

namespace molsturm {

namespace {
template <typename StoredMatrix>
IntegralTermContainer<StoredMatrix> build_hf_terms(
      const gint::IntegralLookup<StoredMatrix>& integrals,
      std::vector<gint::Integral<StoredMatrix>> extra_1e_terms) {
  using gint::IntegralTypeKeys;

  auto T_bb  = integrals.lookup_integral(IntegralTypeKeys::kinetic);
  auto V0_bb = integrals.lookup_integral(IntegralTypeKeys::nuclear_attraction);
  auto J_bb  = integrals.lookup_integral(IntegralTypeKeys::coulomb);
  auto K_bb  = integrals.lookup_integral(IntegralTypeKeys::exchange);

  std::vector<gint::Integral<StoredMatrix>> terms_1e = {
        {std::move(T_bb), std::move(V0_bb)}};
  for (auto&& term : extra_1e_terms) terms_1e.push_back(std::move(term));
  IntegralTermContainer<StoredMatrix> integral_container(
        std::move(terms_1e), std::move(J_bb), std::move(K_bb));

  return integral_container;
}
}  // namespace anonymous

/** Class representing a fock operator */
template <typename StoredMatrix, RestrictionType restriction>
class FockOperator final
      : public krims::conditional_t<restriction == RestrictionType::RestrictedClosed,
                                    RestrictedClosedIntegralOperator<StoredMatrix>,
                                    UnrestrictedIntegralOperator<StoredMatrix>> {
 public:
  static_assert(restriction != RestrictionType::RestrictedOpen,
                "RestrictedOpen is not implemented yet.");

  /** Access to the RestrictionType flag */
  static constexpr RestrictionType restriction_type = restriction;

  typedef krims::conditional_t<restriction == RestrictionType::RestrictedClosed,
                               RestrictedClosedIntegralOperator<StoredMatrix>,
                               UnrestrictedIntegralOperator<StoredMatrix>>
        base_type;
  typedef StoredMatrix stored_matrix_type;
  typedef gint::Integral<stored_matrix_type> int_term_type;
  typedef typename base_type::lazy_matrix_expression_ptr_type
        lazy_matrix_expression_ptr_type;

  FockOperator(const gint::IntegralLookup<StoredMatrix>& integrals,
               const MolecularSystem& system,
               std::vector<int_term_type> extra_1e_terms = {})
        : base_type{build_hf_terms(integrals, std::move(extra_1e_terms)), system} {}

  /** Clone the matrix */
  lazy_matrix_expression_ptr_type clone() const override {
    return lazy_matrix_expression_ptr_type(new FockOperator(*this));
  }
};

}  // namespace molsturm
