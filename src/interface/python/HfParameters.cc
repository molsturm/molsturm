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

#include "HfParameters.hh"
#include "ExcInvalidParameters.hh"
#include "gint/config.hh"
#include <gint/sturmian/atomic/NlmBasis.hh>
#include <lazyten/Base/Solvers/Eigensolution.hh>
#include <lazyten/SmallMatrix.hh>
#include <limits>
#include <molsturm/scf_guess.hh>

namespace molsturm {
namespace iface {

void HfParameters::set_molecular_system(long* atom_numbers, int n_atoms_an,
                                        double* coords, int n_atoms_c, int three_c,
                                        size_t n_alpha_, size_t n_beta_) {
  size_t n_atoms = static_cast<size_t>(n_atoms_an);
  assert_throw(n_atoms_an == n_atoms_c,
               krims::ExcSizeMismatch(n_atoms, static_cast<size_t>(n_atoms_c)));
  assert_throw(three_c == 3, krims::ExcSizeMismatch(static_cast<size_t>(three_c), 3ul));
  static_assert(std::is_same<gint::real_type, double>::value,
                "gint::real_type is expected to be double right now.");

  // Build the structure object from the plain arrays:
  structure.clear();
  structure.reserve(n_atoms);
  for (size_t i = 0; i < n_atoms; ++i) {
    std::array<double, 3> arr_coords{
          {coords[3 * i + 0], coords[3 * i + 1], coords[3 * i + 2]}};
    assert_throw(atom_numbers[i] > 0, krims::ExcTooLarge<long>(0, atom_numbers[i]));
    structure.emplace_back(atom_numbers[i], std::move(arr_coords));
  }

  n_alpha = n_alpha_;
  n_beta  = n_beta_;
}

void HfParameters::set_integral_param_nlm_basis(long* nlm, int n_nlm, int three_n) {
  using gint::sturmian::atomic::Nlm;
  using gint::sturmian::atomic::NlmBasis;

  assert_throw(three_n == 3,
               ExcInvalidParameters("The nlm array needs to be (n_bas, 3), i.e. the "
                                    "column dimension needs to be 3."));
  const int qnum_max = 1000;

  const size_t n_nlm_size = static_cast<size_t>(n_nlm);
  NlmBasis nlm_basis;
  nlm_basis.reserve(n_nlm_size);
  for (size_t i = 0; i < n_nlm_size; ++i) {
    Nlm tple{static_cast<int>(nlm[3 * i]), static_cast<int>(nlm[3 * i + 1]),
             static_cast<int>(nlm[3 * i + 2])};

    // Check values (mostly because we narrow long -> int)
    assert_throw(tple.n > 0 && tple.n < qnum_max,
                 ExcInvalidParameters("Values for n need to be in the range [1, " +
                                      std::to_string(qnum_max - 1) + "] and not " +
                                      std::to_string(tple.n) + "."));
    assert_throw(
          tple.l >= 0 && tple.l < tple.n,
          ExcInvalidParameters("Values for l need to be in the range [0, n-1] == [0, " +
                               std::to_string(tple.n - 1) + "] and not " +
                               std::to_string(tple.l) + "."));
    assert_throw(
          std::abs(tple.m) <= tple.l,
          ExcInvalidParameters("Values for m need to be in the range [-l, l] == [" +
                               std::to_string(-tple.l) + ", " + std::to_string(tple.l) +
                               "] and not " + std::to_string(tple.m) + "."));

    nlm_basis.push_back(std::move(tple));
  }
  integral_params.update("nlm_basis", std::move(nlm_basis));
}

void HfParameters::set_scf_param_external_guess(double* orbena_f, int n_fock_ae,
                                                double* orbenb_f, int n_fock_be,
                                                double* orbcoeffa_bf, int n_bas_ac,
                                                int n_fock_ac, double* orbcoeffb_bf,
                                                int n_bas_bc, int n_fock_bc) {
  assert_throw(
        n_fock_ae == n_fock_ac,
        ExcInvalidParameters("Size mismatch in the number of orbitals in orbena_f (" +
                             std::to_string(n_fock_ae) + ") and orbcoeffa_bf (" +
                             std::to_string(n_fock_ac) + ")."));
  assert_throw(
        n_fock_be == n_fock_bc,
        ExcInvalidParameters("Size mismatch in the number of orbitals in orbenb_f (" +
                             std::to_string(n_fock_be) + ") and orbcoeffb_bf (" +
                             std::to_string(n_fock_bc) + ")."));
  assert_throw(n_bas_ac == n_bas_bc,
               ExcInvalidParameters(
                     "Size mismatch in the number of basis functions in orbcoeffa_bf (" +
                     std::to_string(n_bas_ac) + ") and orbcoeffb_bf (" +
                     std::to_string(n_bas_bc) + ")."));

  const int n_fock_a = n_fock_ac;
  const int n_fock_b = n_fock_bc;
  const int n_bas    = n_bas_ac;

  // If we are unrestricted then we need to build the block-diagonal structure inside
  // the eigensolution, which we pass onto the guess_external method and hence to
  // the initial Fock matrix. We do this by shifting the beta block by n_fock_ac
  // and n_bas_ac and placing a block of size (n_fock_bc, n_bas_bc) there.
  // If (n_fock_bc, n_bas_bc) == (0, 0) i.e. if the beta block is absent,
  // then we are plainly RHF.
  lazyten::Eigensolution<double, lazyten::SmallVector<double>> esolution;

  auto& evalues = esolution.evalues();
  evalues.resize(n_fock_a + n_fock_b);
  auto it = std::copy(orbena_f, orbena_f + n_fock_a, evalues.begin());
  std::copy(orbenb_f, orbenb_f + n_fock_b, it);

  auto& evectors = esolution.evectors();
  evectors       = lazyten::MultiVector<lazyten::SmallVector<double>>(2 * n_bas,
                                                                n_fock_a + n_fock_b);
  for (int b = 0; b < n_bas; ++b) {
    // The alpha block
    for (int f = 0; f < n_fock_a; ++f) {
      evectors[f][b] = orbcoeffa_bf[b * n_fock_a + f];
    }

    // The beta block
    for (int f = 0; f < n_fock_b; ++f) {
      evectors[n_fock_a + f][n_bas + b] = orbcoeffb_bf[b * n_fock_b + f];
    }
  }

  guess_params.update(GuessExternalKeys::eigensolution, std::move(esolution));
}

template <typename T>
void HfParameters::set_param(HfParameters::ParameterKind kind, std::string key, T value) {
  switch (kind) {
    case HfParameters::ParameterKind::SCF:
      scf_params.update(key, value);
      break;
    case HfParameters::ParameterKind::GUESS:
      guess_params.update(key, value);
      break;
    case HfParameters::ParameterKind::INTEGRAL:
      integral_params.update(key, value);
      break;
    default:
      assert_throw(false, ExcInvalidParameters(
                                "Invalid value for HfParameters::ParameterKind: " +
                                std::to_string(kind)));
  }
}

#define INSTANTIATE_SETPARAM(TYPE)                                        \
  template void HfParameters::set_param(HfParameters::ParameterKind kind, \
                                        std::string key, TYPE value);
INSTANTIATE_SETPARAM(bool);
INSTANTIATE_SETPARAM(double);
INSTANTIATE_SETPARAM(int);
INSTANTIATE_SETPARAM(size_t);
INSTANTIATE_SETPARAM(std::string);

#undef INSTANTIATE_SET

}  // namespace iface
}  // namespace molsturm
