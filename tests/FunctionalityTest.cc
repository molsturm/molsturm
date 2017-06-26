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
#include <gint/IntegralLookupKeys.hh>
#include <gint/OrbitalType.hh>
#include <linalgwrap/SmallMatrix.hh>
#include <linalgwrap/TestingUtils.hh>
#include <molsturm/IopScf.hh>
#include <molsturm/MolecularSystem.hh>
#include <molsturm/RestrictedClosedIntegralOperator.hh>

namespace molsturm {
namespace tests {
using namespace gscf;
using namespace gint;
using namespace linalgwrap;
using namespace krims;

TEST_CASE("HF functionality test", "[hf functionality]") {
  //
  // Parameters of the test problem
  //
  MolecularSystem sys(
        {
              {"Be", {{0, 0, 0}}},
        },
        /* charge */ 0);
  double k_exp = 1.;
  int n_max = 3;
  int l_max = 2;
  size_t n_eigenpairs = 4;
  const double tolerance = 1e-9;
  const std::string basis_type = "sturmian/atomic/cs_naive";
  const gint::OrbitalType otype = gint::OrbitalType::COMPLEX_ATOMIC;

  //
  // Types to use
  //
  // Types of scalar and matrix
  typedef double scalar_type;
  typedef SmallMatrix<scalar_type> stored_matrix_type;
  typedef stored_matrix_type::vector_type vector_type;

  // The lookup class type to get the actual integrals
  typedef gint::IntegralLookup<stored_matrix_type> int_lookup_type;

  // The type of the integral terms:
  typedef typename int_lookup_type::integral_type integral_type;

  //
  // Reference data
  //
  // The expected eigenvalues
  std::vector<scalar_type> eval_expected{
        -3.7496030113878,    -0.1755827076801148, -0.03942386983346891,
        0.08845420348863607, 0.08845420348863607, 0.2470756635901147,
        0.3644265733015022,  0.3644265733015022,  0.3663266541722434,
        0.3727104156955993,  0.3727104156955993,  0.3761194352770711,
        0.383552426594183,   0.383552426594183};

  // The expected eigenvectors
  MultiVector<vector_type> evec_expected{
        {-1.1884746467844811, 0., -0.24288050447076578, 0., 0., 0.16813161525502301, 0.,
         0., 0., 0., 0., 0.016396171817405908, 0., 0.},
        {-1.064786764522789, 0., 0.8777407505081162, 0., 0., -0.3081669311487548, 0., 0.,
         0., 0., 0., -0.028869768632584114, 0., 0.},
        {0., 0., 0., -3.69049544150639e-9, 0.8573394853277652, 0., 0., 0., 0., 0., 0., 0.,
         0.00002919059875836615, -0.6818863586007807},
        {0., 0.9857660367413854, 0., 0., 0., 0., 0., 0., 0.47777120131625944, 0., 0., 0.,
         0., 0.},
        {0., 0., 0., 0.8573394853277649, 3.69049544151069e-9, 0., 0., 0., 0., 0., 0., 0.,
         -0.6818863586007805, -0.000029190598758366127},
        {-0.5840485708369669, 0., 0.05174625401524502, 0., 0., -1.0729001918355632, 0.,
         0., 0., 0., 0., -0.07137766077631158, 0., 0.},
        {0., 0., 0., 1.1728582480320243e-9, -0.27246685510597846, 0., 0., 0., 0., 0., 0.,
         0., 0.000045420745388700296, -1.0610192320620837},
        {0., -0.033706141181753996, 0., 0., 0., 0., 0., 0., 1.0949264340797673, 0., 0.,
         0., 0., 0.},
        {0., 0., 0., -0.2724668551059781, -1.1728582480253274e-9, 0., 0., 0., 0., 0., 0.,
         0., -1.0610192320620837, -0.00004542074538870031},
        {0., 0., 0., 0., 0., 0., -1., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.000043369360575274524, 0.9999999990595501,
         0., 0., 0.},
        {-0.0019206466502236202, 0., -0.011672197660675484, 0., 0., 0.06685683586559842,
         0., 0., 0., 0., 0., -0.9976924548257627, 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., -0.9999999990595487, 0.00004336936057527441,
         0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.}};

  // The expected energies:
  scalar_type exp_energy_coulomb = 5.268082025825932;
  scalar_type exp_energy_exchange = -1.73626312131552;
  scalar_type exp_energy_kinetic = 5.720460090270233;
  scalar_type exp_energy_nucattr = -20.6344714235971;
  scalar_type exp_energy_1e_terms = exp_energy_nucattr + exp_energy_kinetic;
  scalar_type exp_energy_2e_terms = exp_energy_coulomb + exp_energy_exchange;
  scalar_type exp_energy_total = exp_energy_1e_terms + exp_energy_2e_terms;

  //
  // Setup integrals
  //
  krims::GenMap intparams{{IntegralLookupKeys::basis_type, basis_type},
                          {IntegralLookupKeys::orbital_type, otype},
                          {"k_exponent", k_exp},
                          {IntegralLookupKeys::structure, sys.structure},
                          {"n_max", n_max},
                          {"l_max", l_max},
                          {"m_max", l_max}};
  int_lookup_type integrals{intparams};

  // Get the integral as actual objects.
  integral_type S_bb = integrals.lookup_integral("overlap");
  integral_type T_bb = integrals.lookup_integral("kinetic");
  integral_type V0_bb = integrals.lookup_integral("nuclear_attraction");
  integral_type J_bb = integrals.lookup_integral("coulomb");
  integral_type K_bb = integrals.lookup_integral("exchange");

  // Combine 1e terms:
  std::vector<integral_type> terms_1e{T_bb, V0_bb};

  //
  // Problem setup and run.
  //
  auto guess_bf_ptr = std::make_shared<MultiVector<vector_type>>(S_bb.n_rows(), 4);
  (*guess_bf_ptr)[0](0) = -0.9238795325112872L;
  (*guess_bf_ptr)[0](1) = -1.306562964876376L;
  (*guess_bf_ptr)[0](5) = -0.9238795325112864L;
  (*guess_bf_ptr)[1](3) = (*guess_bf_ptr)[1](7) = -0.9192110607898044L;
  (*guess_bf_ptr)[2](2) = (*guess_bf_ptr)[2](6) = -0.9192110607898044L;
  (*guess_bf_ptr)[3](4) = (*guess_bf_ptr)[3](8) = -0.9192110607898044L;

  /*
   * TODO the guess builder is currently pretty bad and unreliable
   *      ... so short circuit it here.
   *
   * auto guess_bf_ptr = std::make_shared<linalgwrap::MultiVector<vector_type>>(
   *    loewdin_guess(S_bb, n_eigenpairs));
   */

  // The term container for the fock operator matrix
  IntegralTermContainer<stored_matrix_type> integral_container(std::move(terms_1e), J_bb,
                                                               K_bb);

  RestrictedClosedIntegralOperator<stored_matrix_type> fock_bb(integral_container, sys);
  fock_bb.update(guess_bf_ptr);

  const krims::GenMap params{{IopScfKeys::max_iter, 15ul},
                             {IopScfKeys::n_eigenpairs, n_eigenpairs},
                             {IopScfKeys::max_error_norm, tolerance}};

#ifdef DEBUG
  std::cout << "Running test SCF ... please wait." << std::endl;
#endif
  auto res = run_scf(fock_bb, S_bb, params);

  //
  // Check the results
  //
  CHECK(res.n_iter() <= 13);

  // Check the eigenvalues
  const auto& evalues = res.orbital_energies();
  eval_expected.resize(n_eigenpairs);
  CHECK(vector_type(evalues) ==
        numcomp(vector_type(eval_expected)).tolerance(1000. * tolerance));

  const auto& evectors = res.orbital_coeff();
  // TODO For comparing all of them one needs to take rotations
  //      inside degenerate subspaces into account
  for (size_t i = 0; i < sys.n_alpha; ++i) {
    linalgwrap::adjust_phase(evectors[i], evec_expected[i]);
    CHECK(evectors[i] == numcomp(evec_expected[i]).tolerance(100. * tolerance));
  }

  // Check the energies:
  double energytol = tolerance;
  const auto& fock = res.problem_matrix();
  auto energies = fock.energies();
  CHECK(energies[J_bb.id()] == numcomp(exp_energy_coulomb).tolerance(energytol));
  CHECK(energies[K_bb.id()] == numcomp(exp_energy_exchange).tolerance(energytol));
  CHECK(energies[V0_bb.id()] == numcomp(exp_energy_nucattr).tolerance(energytol));
  CHECK(energies[T_bb.id()] == numcomp(exp_energy_kinetic).tolerance(energytol));
  CHECK(fock.energy_1e_terms() == numcomp(exp_energy_1e_terms).tolerance(energytol));
  CHECK(fock.energy_2e_terms() == numcomp(exp_energy_2e_terms).tolerance(energytol));
  CHECK(fock.energy_total() == numcomp(exp_energy_total).tolerance(energytol));

}  // TEST_CASE

}  // namespace test
}  // namespace molsturm
