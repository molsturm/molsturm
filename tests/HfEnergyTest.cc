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

#include <catch.hpp>
#include <gint/IntegralLookup.hh>
#include <linalgwrap/TestingUtils.hh>
#include <molsturm/IopScf.hh>
#include <molsturm/RestrictedClosedIntegralOperator.hh>
#include <molsturm/scf_guess.hh>

namespace molsturm {
namespace tests {
using namespace linalgwrap;
using namespace krims;

//
// Types to use
//
typedef double scalar_type;
typedef SmallMatrix<scalar_type> stored_matrix_type;
typedef stored_matrix_type::vector_type vector_type;
typedef gint::Integral<stored_matrix_type> integral_type;

namespace hf_energy_test {

struct ReferenceData {
  size_t max_iter{0};
  double energy_total{0};
  double energy_1e{0};
  double energy_2e{0};
  std::vector<double> energies_mos{};
};

template <typename Fock = RestrictedClosedIntegralOperator<stored_matrix_type>>
void run_testscf(const MolecularSystem& sys, const GenMap& intparams,
                 const GenMap& scfparams, const ReferenceData& reference) {
  using gint::IntegralTypeKeys;

  // Types:
  const gint::OrbitalType otype = gint::OrbitalType::REAL_MOLECULAR;
  typedef gint::IntegralLookup<otype> int_lookup_type;

  // Update integral parameters:
  intparams.insert_default("structure", sys.structure);

  // Get the integrals:
  int_lookup_type integrals{intparams};
  integral_type S_bb = integrals.lookup_integral(IntegralTypeKeys::overlap);
  integral_type T_bb = integrals.lookup_integral(IntegralTypeKeys::kinetic);
  integral_type V0_bb = integrals.lookup_integral(IntegralTypeKeys::nuclear_attraction);
  integral_type J_bb = integrals.lookup_integral(IntegralTypeKeys::coulomb);
  integral_type K_bb = integrals.lookup_integral(IntegralTypeKeys::exchange);

  // The fock operator of this problem:
  Fock fock_init({{T_bb, V0_bb}, J_bb, K_bb}, sys);

  // Get a guess and run scf:
  auto guess = scf_guess(sys, fock_init, S_bb);
  auto res = run_scf(fock_init, S_bb, guess, scfparams);

  // Compare results:
  const auto& orben = res.orbital_energies();
  const auto& fock = res.problem_matrix();
  double tolerance = scfparams.at<double>(IopScfKeys::max_error_norm);

  CHECK(res.n_iter() < reference.max_iter);
  CHECK(fock.energy_total() == numcomp(reference.energy_total).tolerance(tolerance));

  // TODO Right now the reference value are unfortunately very crude
  //      Get better reference values!
  tolerance = std::max(tolerance, 5e-9);

  CHECK(fock.energy_1e_terms() == numcomp(reference.energy_1e).tolerance(tolerance));
  CHECK(fock.energy_2e_terms() == numcomp(reference.energy_2e).tolerance(tolerance));

  // TODO Right now the reference value are unfortunately very crude
  //      Get better reference values!
  tolerance = std::max(tolerance, 1e-7);
  CHECK(vector_type(orben) ==
        numcomp(vector_type(reference.energies_mos)).tolerance(10 * tolerance));
}

}  // namespace hf_energy_test

TEST_CASE("Test HF energies and MOs compared to ORCA", "[hf energies]") {
  using namespace hf_energy_test;

  //
  // Molecules we use
  //
  MolecularSystem he(
        {
              {"He", {{0, 0, 0}}},
        },
        /* charge */ 0);
  //
  MolecularSystem be(
        {
              {"Be", {{0, 0, 0}}},
        },
        /* charge */ 0);
  //
  MolecularSystem water(
        {
              {"O", {{0, 0, 0}}},                                   //
              {"H", {{0, 0, 1.795239827225189}}},                   //
              {"H", {{1.693194615993441, 0, -0.599043184453037}}},  //
        },
        /* charge */ 0);

//
// Common parameters
//
#ifdef DEBUG
  const double tolerance = 1e-7;  // Debug mode is slower
#else
  const double tolerance = 1e-9;  // Release mode
#endif  // DEBUG

  krims::GenMap intparams{{"basis_type", "gaussian/libint"}};
  krims::GenMap scfparams{{IopScfKeys::n_eigenpairs, size_t(10)},
                          {IopScfKeys::max_error_norm, tolerance}};

  //
  // =========================================================================
  //

  //
  // Water STO-3G
  //
  SECTION("H2o sto-3g") {
    ReferenceData d;
    d.max_iter = 15;
    d.energy_total = -74.959319286910;
    d.energy_1e = -122.50621280;
    d.energy_2e = 38.29541424;
    d.energies_mos = {-20.233397, -1.265715, -0.629267, -0.441668,
                      -0.387645,  0.602839,  0.765918};

    intparams.update("basis_set", "sto-3g");
    run_testscf(water, intparams, scfparams, d);
  }  // h2o sto-3g

  //
  // Helium aug-cc-pvdz
  //
  SECTION("he aug-cc-pvdz") {
    ReferenceData d;
    d.max_iter = 8;
    d.energy_total = -2.855704667706;
    d.energy_1e = -3.87716137;
    d.energy_2e = 1.02145670;
    d.energies_mos = {-0.917124, 0.174366, 0.530376, 0.530376, 0.530376,
                      1.713453,  3.024883, 3.024883, 3.024883};

    intparams.update("basis_set", "aug-cc-pvdz");
    scfparams.update(IopScfKeys::max_error_norm, tolerance);
    run_testscf(he, intparams, scfparams, d);
  }  // he aug-cc-pvdz

  //
  // Beryllium Def2-tzvp
  //
  SECTION("be def2-tzvp") {
    ReferenceData d;
    d.max_iter = 12;
    d.energy_total = -14.572579867386;
    d.energy_1e = -19.06163937;
    d.energy_2e = 4.48905951;
    d.energies_mos = {-4.732554, -0.309206, 0.058138, 0.058138, 0.058138,
                      0.087734,  0.349301,  0.349301, 0.349301, 0.492669};

    intparams.update("basis_set", "def2-tzvp");
    scfparams.update(IopScfKeys::max_error_norm, tolerance);
    run_testscf(be, intparams, scfparams, d);
  }  // be def2-tzvp

//
// =========================================================================
//

#ifndef DEBUG
  //
  // Water cc-pvdz
  //
  SECTION("H2o cc-pvdz") {
    // Decrease tolerance (harder problem)
    const double tolerance = 1e-8;
    scfparams.update(IopScfKeys::max_error_norm, tolerance);

    ReferenceData d;
    d.max_iter = 19;
    d.energy_total = -76.026348968453;
    d.energy_1e = -123.26657747;
    d.energy_2e = 37.98874924;
    d.energies_mos = {-20.546361, -1.336322, -0.712324, -0.557705, -0.491828,
                      0.187645,   0.256614,  0.812822,  0.837446,  1.160574};

    intparams.update({{"basis_set", "cc-pvdz"}});
    run_testscf(water, intparams, scfparams, d);
  }     // h2o cc-pvdz
#endif  // DEBUG

}  // hf energies

}  // namespace tests
}  // namespace molsturm
