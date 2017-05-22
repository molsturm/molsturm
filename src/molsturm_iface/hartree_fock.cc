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

#include "hartree_fock.hh"
#include <gint/IntegralLookup.hh>
#include <gint/IntegralLookupKeys.hh>
#include <gint/IntegralType.hh>
#include <gint/OrbitalType.hh>
#include <gint/Structure.hh>
#include <gint/sturmian/atomic/NlmBasis.hh>
#include <iterator>
#include <linalgwrap/SmallMatrix.hh>
#include <molsturm/IopScf.hh>
#include <molsturm/scf_guess.hh>

namespace molsturm {
namespace iface {

// Type definitions (should move to a config file later)
typedef double scalar_type;
typedef linalgwrap::SmallMatrix<scalar_type> matrix_type;

MolecularSystem build_molecular_system(const Parameters& params) {
  assert_throw(params.coords.size() == params.atom_numbers.size(),
               ExcInvalidParameters("Sizes of the atomic numbers array and of the "
                                    "coordinates array need to agree."));
  assert_throw(params.coords.size() > 0,
               ExcInvalidParameters("Need at least one atom to do a calculation."));

  gint::Structure s;
  std::transform(params.coords.begin(), params.coords.end(), params.atom_numbers.begin(),
                 std::back_inserter(s),
                 [](const std::array<double, 3>& coord, unsigned int atom_number) {
                   return gint::Atom{atom_number, coord};
                 });

  if (params.multiplicity > 0) {
    return MolecularSystem(s, params.charge, params.multiplicity);
  } else {
    return MolecularSystem(s, params.charge);
  }
}

krims::GenMap build_int_params_sturmian(const Parameters& params,
                                        const MolecularSystem& system) {
  using gint::IntegralLookupKeys;

  krims::GenMap intparams;

  assert_implemented(system.structure.n_atoms() == 1);
  intparams.update({
        {IntegralLookupKeys::orbital_type, gint::OrbitalType::COMPLEX_ATOMIC},
        {"k_exponent", params.k_exp},
  });

  if (params.nlm_basis.size() > 0) {
    using gint::sturmian::atomic::Nlm;
    using gint::sturmian::atomic::NlmBasis;

    NlmBasis nlm_basis_conv;
    std::transform(params.nlm_basis.begin(), params.nlm_basis.end(),
                   std::back_inserter(nlm_basis_conv), [](std::array<int, 3> in) {
                     return Nlm{in[0], in[1], in[2]};
                   });

    intparams.update("nlm_basis", std::move(nlm_basis_conv));
  } else {
    const int n_max = params.n_max;
    const int l_max = params.l_max == Parameters::all ? n_max - 1 : params.l_max;
    const int m_max = params.m_max == Parameters::all ? l_max : params.m_max;

    intparams.update({
          {"n_max", n_max}, {"l_max", l_max}, {"m_max", m_max},
    });
  }
  return intparams;
}

krims::GenMap build_int_params(const Parameters& params, const MolecularSystem& system) {
  using gint::IntegralLookupKeys;
  assert_throw(!params.basis_type.empty(), ExcInvalidParameters("No basis type given."));

  const std::string gaussian_prefix = "gaussian/";
  const bool gaussian =
        params.basis_type.compare(0, gaussian_prefix.size(), gaussian_prefix) == 0;

  krims::GenMap intparams{
        {IntegralLookupKeys::basis_type, params.basis_type},
        {IntegralLookupKeys::structure, system.structure},
  };

  if (gaussian) {
    assert_throw(!params.basis_set.empty(),
                 ExcInvalidParameters("No gaussian basis set given."));

    intparams.update({
          {IntegralLookupKeys::orbital_type, gint::OrbitalType::REAL_MOLECULAR},
          {"basis_set", params.basis_set},
    });
  } else {
    intparams.update(build_int_params_sturmian(params, system));
  }

  return intparams;
}

krims::GenMap build_guess_params(const Parameters& params) {
  using linalgwrap::EigensystemSolverKeys;

  krims::GenMap guess_params{{ScfGuessKeys::method, params.guess_method}};
  guess_params.update(ScfGuessKeys::eigensolver_params,
                      {{EigensystemSolverKeys::method, params.guess_esolver}});
  return guess_params;
}

krims::GenMap build_scf_params(const Parameters& params) {
  using linalgwrap::EigensystemSolverKeys;

  ScfMsgType verbosity = ScfMsgType::Silent;
  if (params.print_iterations) verbosity |= ScfMsgType::IterationProcess;
  if (params.print_scf_summary) verbosity |= ScfMsgType::FinalSummary;

  krims::GenMap scfparams{
        // error
        {IopScfKeys::max_error_norm, params.error},
        {IopScfKeys::max_1e_energy_change, params.error * 100.},
        {IopScfKeys::max_tot_energy_change, params.error / 4.},
        //
        {IopScfKeys::max_iter, params.max_iter},
        {IopScfKeys::n_eigenpairs, params.n_eigenpairs},
        {IopScfKeys::verbosity, verbosity},
        {gscf::PulayDiisScfKeys::n_prev_steps, params.diis_size},
  };
  scfparams.update(IopScfKeys::eigensolver_params,
                   {{EigensystemSolverKeys::method, params.eigensolver}});

  return scfparams;
}

template <typename State>
HfResults build_restricted_hf_results(const krims::GenMap& scf_params,
                                      const gint::ERITensor_i<scalar_type>& eri,
                                      const State& state) {
  const auto& fbb = state.problem_matrix();
  const auto& soln = state.eigensolution();
  const size_t n_orbs_alpha = soln.n_ep();
  const size_t n_orbs = 2 * n_orbs_alpha;
  const size_t n_bas = fbb.n_rows();

  // Copy results back
  HfResults ret;
  ret.n_beta = fbb.indices_subspace(gscf::OrbitalSpace::OCC_BETA).length();
  ret.n_alpha = fbb.indices_subspace(gscf::OrbitalSpace::OCC_ALPHA).length();
  ret.n_bas = n_bas;

  ret.restricted = true;
  ret.n_orbs_alpha = n_orbs_alpha;
  ret.n_orbs_beta = n_orbs_alpha;

  ret.threshold = scf_params.at<double>(IopScfKeys::max_error_norm);

  // SCF results
  ret.n_iter = state.n_iter();
  ret.n_mtx_applies = state.n_mtx_applies();
  ret.final_error_norm = state.last_error_norm;
  ret.final_tot_energy_change = state.last_tot_energy_change;
  ret.final_1e_energy_change = state.last_1e_energy_change;

  // Energies:
  for (const auto& kv : fbb.energies()) {
    switch (kv.first.integral_type()) {
      case gint::IntegralType::kinetic:
        ret.energy_kinetic = kv.second;
        break;
      case gint::IntegralType::coulomb:
        ret.energy_coulomb = kv.second;
        break;
      case gint::IntegralType::exchange:
        ret.energy_exchange = kv.second;
        break;
      case gint::IntegralType::nuclear_attraction:
        ret.energy_nuclear_attraction = kv.second;
        break;
      default:
        continue;
    }
  }
  ret.energy_total = fbb.energy_total();
  ret.energy_nuclear_repulsion = fbb.energy_nuclear_repulsion();

  // Insert alpha and beta energies:
  std::copy(soln.evalues().begin(), soln.evalues().end(),
            std::back_inserter(ret.orbital_energies_f));
  std::copy(soln.evalues().begin(), soln.evalues().end(),
            std::back_inserter(ret.orbital_energies_f));
  assert_size(ret.orbital_energies_f.size(), n_orbs);

  // Compute the full fock matrix in MO space, i.e.  C^T * (F * C)
  // => Need a dot product here, so actually the dot of all vectors with another
  auto fock_ff = linalgwrap::dot(soln.evectors(), fbb * soln.evectors());
  assert_size(2 * fock_ff.n_rows(), n_orbs);
  assert_size(2 * fock_ff.n_cols(), n_orbs);
  ret.fock_ff.resize(n_orbs * n_orbs);

  for (size_t i = 0; i < fock_ff.n_rows(); ++i) {
    for (size_t j = 0; j < fock_ff.n_cols(); ++j) {
      const size_t ij = i * n_orbs + j;
      const size_t ij_beta = ij + n_orbs_alpha * n_orbs + n_orbs_alpha;
      // (i + ret.n_orbs_alpha) * n_orbs + (j + ret.n_orbs_alpha)
      //   = i * n_orbs + j + n_orbs_alpha * n_orbs + n_orbs_alpha
      //   = ij + n_orbs_alpha * n_orbs + n_orbs_alpha
      assert_internal(ij < ret.fock_ff.size());
      assert_internal(ij_beta < ret.fock_ff.size());
      ret.fock_ff[ij] = fock_ff(i, j);
      ret.fock_ff[ij_beta] = fock_ff(i, j);
    }
  }

  // Copy coefficients (also alpha and beta block is needed!)
  ret.coeff_fb.resize(n_orbs * n_bas);
  for (size_t f = 0; f < soln.evectors().n_vectors(); ++f) {
    for (size_t b = 0; b < n_bas; ++b) {
      const size_t fb = f * n_bas + b;
      const size_t fb_beta = (f + ret.n_orbs_alpha) * n_bas + b;
      assert_internal(fb < ret.coeff_fb.size());
      assert_internal(fb_beta < ret.coeff_fb.size());
      ret.coeff_fb[fb] = soln.evectors()[f][b];
      ret.coeff_fb[fb_beta] = soln.evectors()[f][b];
    }  // j
  }    // i

  // Form the alpha-alpha-alpha-alpha spin block by contraction
  // and copy it to all the other non-zero places
  // i.e. (aa|aa), (aa|bb), (bb|aa), (bb|bb)
  std::vector<double> eri_aaaa(n_orbs_alpha * n_orbs_alpha * n_orbs_alpha * n_orbs_alpha);
  eri.contract_with(soln.evectors(), soln.evectors(), soln.evectors(), soln.evectors(),
                    eri_aaaa);

  ret.repulsion_integrals_ffff.resize(n_orbs * n_orbs * n_orbs * n_orbs);
  for (size_t i = 0, ijkl = 0; i < n_orbs_alpha; ++i) {
    for (size_t j = 0; j < n_orbs_alpha; ++j) {
      for (size_t k = 0; k < n_orbs_alpha; ++k) {
        for (size_t l = 0; l < n_orbs_alpha; ++l, ++ijkl) {
          // Index of the element in the alpha-alpha-alpha-alpha
          // block to which we wish to copy our extracted values
          const size_t aaaa = ((i * n_orbs + j) * n_orbs + k) * n_orbs + l;

          // Index offset introduced by shifting the last shell
          // pair from the alpha-alpha to the beta-beta position
          const size_t bb_off = n_orbs * n_orbs_alpha + n_orbs_alpha;

          // The indices of the other elements in the full eri tensor
          // equal to the element at index aaaa by symmetry.
          const size_t aabb = aaaa + bb_off;
          const size_t bbaa = aaaa + n_orbs * n_orbs * bb_off;
          const size_t bbbb = aaaa + bb_off + n_orbs * n_orbs * bb_off;

          assert_internal(aaaa < n_orbs * n_orbs * n_orbs * n_orbs);
          assert_internal(aabb < n_orbs * n_orbs * n_orbs * n_orbs);
          assert_internal(bbaa < n_orbs * n_orbs * n_orbs * n_orbs);
          assert_internal(bbbb < n_orbs * n_orbs * n_orbs * n_orbs);
          assert_internal(ijkl < eri_aaaa.size());

          ret.repulsion_integrals_ffff[aaaa] = eri_aaaa[ijkl];
          ret.repulsion_integrals_ffff[aabb] = eri_aaaa[ijkl];
          ret.repulsion_integrals_ffff[bbaa] = eri_aaaa[ijkl];
          ret.repulsion_integrals_ffff[bbbb] = eri_aaaa[ijkl];
        }  // l
      }    // k
    }      // j
  }        // i

  return ret;
}

HfResults hartree_fock(const Parameters& params) {
  using gint::IntegralTypeKeys;
  typedef gint::IntegralLookup<matrix_type> integral_lookup_type;

  // Parse parameters
  const MolecularSystem system = build_molecular_system(params);
  const krims::GenMap int_params = build_int_params(params, system);
  const krims::GenMap guess_params = build_guess_params(params);
  const krims::GenMap scf_params = build_scf_params(params);

  // Lookup integral terms
  integral_lookup_type integrals(int_params);
  auto S_bb = integrals.lookup_integral(IntegralTypeKeys::overlap);
  auto T_bb = integrals.lookup_integral(IntegralTypeKeys::kinetic);
  auto V0_bb = integrals.lookup_integral(IntegralTypeKeys::nuclear_attraction);
  auto J_bb = integrals.lookup_integral(IntegralTypeKeys::coulomb);
  auto K_bb = integrals.lookup_integral(IntegralTypeKeys::exchange);

  // Checks about basis size:
  const size_t max_elec = std::max(system.n_alpha, system.n_beta);
  assert_throw(max_elec < S_bb.n_rows(), ExcTooSmallBasis(S_bb.n_rows(), max_elec));
  assert_throw(params.n_eigenpairs >= max_elec,
               krims::ExcTooLarge<size_t>(max_elec, params.n_eigenpairs));

  // Problem setup
  typedef RestrictedClosedIntegralOperator<matrix_type> fock_type;
  IntegralTermContainer<matrix_type> integral_container(
        {{std::move(T_bb), std::move(V0_bb)}}, std::move(J_bb), std::move(K_bb));
  fock_type fock_bb(integral_container, system);

  // Run solver
  auto guess = scf_guess(system, fock_bb, S_bb, guess_params);
  fock_bb.update(guess.evectors_ptr);

  auto result = run_scf(fock_bb, S_bb, guess, scf_params);

  return build_restricted_hf_results(scf_params, integrals.eri_tensor(), result);
}

}  // namespace iface
}  // namespace molsturm
