//
// Copyright (C) 2018 by the molsturm authors
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

#include "integrals.hh"
#include "ExcInvalidParameters.hh"
#include "config.hh"
#include <gint/IntegralLookup.hh>
#include <molsturm/MolecularSystem.hh>

namespace molsturm {
namespace iface {

gint::IntegralLookup<matrix_type> get_integrals(const ScfParameters& params, int n_bas1,
                                                int n_bas2) {
  // Add the structure for the integral library if not done yet:
  params.insert_default("discretisation/structure",
                        params.at_raw_value("system/structure"));
  gint::IntegralLookup<matrix_type> integrals(params.submap("discretisation"));
  const size_t n_bas = integrals.n_bas();

  assert_throw(n_bas1 == n_bas2,
               ExcInvalidParameters("Input array should have a quadratic shape."));
  assert_throw(
        n_bas == static_cast<size_t>(n_bas1),
        ExcInvalidParameters("The shape of the input array (= " + std::to_string(n_bas1) +
                             ") and the number of basis "
                             "functions returned by the integral library (= " +
                             std::to_string(n_bas) + ") does not agree."));

  return integrals;
}

void copy_matrix(const matrix_type& mat, double* out) {
  for (size_t i = 0; i < mat.n_rows(); ++i) {
    for (size_t j = 0; j < mat.n_rows(); ++j, ++out) {
      *out = mat(i, j);
    }
  }
}

void overlap_bb(const ScfParameters& params, double* out_bb, int n_bas1, int n_bas2) {
  auto integrals = get_integrals(params, n_bas1, n_bas2);
  auto sbb       = integrals.lookup_integral(gint::IntegralTypeKeys::overlap);
  copy_matrix(static_cast<matrix_type>(sbb), out_bb);
}

void kinetic_bb(const ScfParameters& params, double* out_bb, int n_bas1, int n_bas2) {
  auto integrals = get_integrals(params, n_bas1, n_bas2);
  auto tbb       = integrals.lookup_integral(gint::IntegralTypeKeys::kinetic);
  copy_matrix(static_cast<matrix_type>(tbb), out_bb);
}

void nuclear_attraction_bb(const ScfParameters& params, double* out_bb, int n_bas1,
                           int n_bas2) {
  auto integrals = get_integrals(params, n_bas1, n_bas2);
  auto vbb       = integrals.lookup_integral(gint::IntegralTypeKeys::nuclear_attraction);
  copy_matrix(static_cast<matrix_type>(vbb), out_bb);
}

}  // namespace iface
}  // namespace molsturm
