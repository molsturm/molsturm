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

#include "parse_parameters.hh"
#include "config.hh"
#include <gint/IntegralLookupKeys.hh>
#include <gint/OrbitalType.hh>
#include <gint/sturmian/atomic/NlmBasis.hh>
#include <gscf/PulayDiisScfKeys.hh>
#include <lazyten/EigensystemSolver.hh>
#include <molsturm/GuessAlgorithms.hh>
#include <molsturm/IopScfKeys.hh>
#include <molsturm/ScfMsgType.hh>
#include <molsturm/scf_guess.hh>

namespace molsturm {
namespace iface {

gint::Structure build_structure(const Parameters& params) {
  assert_throw(params.atom_numbers.size() > 0 || params.atoms.size() > 0,
               ExcInvalidParameters(
                     "At least one of atom_numbers, atoms needs to contain an entry."));
  assert_throw(params.atom_numbers.size() == 0 || params.atoms.size() == 0,
               ExcInvalidParameters(
                     "Exactly one of atom_numbers, atoms may contain atom definitions."));

  for (auto& num : params.atom_numbers) {
    assert_throw(num,
                 ExcInvalidParameters("Atomic numbers need to be larger than zero."));
  }

  size_t n_atoms = std::max(params.atom_numbers.size(), params.atoms.size());
  gint::Structure s;
  s.reserve(n_atoms);

  // The case of one atom is special ... deal with it first.
  if (n_atoms == 1) {
    if (params.atom_numbers.size() > 0) {
      const unsigned int atnum = static_cast<unsigned int>(params.atom_numbers[0]);
      s.push_back(gint::Atom{atnum, {{0., 0., 0.}}});
    } else {
      s.push_back(gint::Atom(params.atoms[0], {{0., 0., 0.}}));
    }
    return s;
  }

  // Check that we get the right number of coords
  assert_throw(params.coords.size() == 3 * n_atoms,
               ExcInvalidParameters("We expect the size of the coords vector to be "
                                    "exactly three times the number of atoms, i.e. " +
                                    std::to_string(3 * n_atoms) +
                                    " entries. But instead we got only " +
                                    std::to_string(params.coords.size()) + " entries."));

  for (size_t i = 0; i < n_atoms; ++i) {
    const std::array<double, 3> coord{
          {params.coords[3 * i], params.coords[3 * i + 1], params.coords[3 * i + 2]}};

    if (params.atom_numbers.size() > 0) {
      unsigned int atnum = static_cast<unsigned int>(params.atom_numbers[i]);
      s.emplace_back(atnum, std::move(coord));
    } else {
      s.emplace_back(params.atoms[i], std::move(coord));
    }
  }
  return s;
}

MolecularSystem build_molecular_system(const Parameters& params) {
  const gint::Structure s = build_structure(params);
  if (params.multiplicity > 0) {
    return MolecularSystem(std::move(s), params.charge, params.multiplicity);
  } else {
    return MolecularSystem(std::move(s), params.charge);
  }
}

bool parse_restricted(const Parameters& params, const MolecularSystem& system) {
  // Adjust value for restricted => Use automatically determined value
  // if the user did not override this.
  const bool restricted = params.internal_restricted_set_by_user
                                ? params.restricted
                                : system.n_alpha == system.n_beta;
  assert_throw(!restricted || system.n_alpha == system.n_beta,
               ExcInvalidParameters("Only systems with even electron count can be "
                                    "treated using restricted calculations. Use an "
                                    "unrestricted treatment for the other systems."));
  return restricted;
}

krims::GenMap build_int_params_sturmian(const Parameters& params,
                                        const MolecularSystem& system) {
  using gint::IntegralLookupKeys;
  assert_implemented(system.structure.n_atoms() == 1);

  assert_throw(params.k_exp > 0,
               ExcInvalidParameters("The parameter k_exp is mandatory for a sturmian "
                                    "basis with a value larger than zero."));
  krims::GenMap intparams{
        {IntegralLookupKeys::orbital_type, gint::OrbitalType::COMPLEX_ATOMIC},
        {"k_exponent", params.k_exp},
  };

  if (params.nlm_basis.size() > 0) {
    using gint::sturmian::atomic::Nlm;
    using gint::sturmian::atomic::NlmBasis;

    assert_throw(
          params.nlm_basis.size() % 3 == 0,
          ExcInvalidParameters(
                "The size of the nlm_basis array needs to be exactly divisible by 3"));

    NlmBasis nlm_basis_conv;
    nlm_basis_conv.reserve(params.nlm_basis.size() / 3);
    for (size_t i = 0; i < params.nlm_basis.size() / 3; ++i) {
      const auto& nlmbas = params.nlm_basis;
      nlm_basis_conv.push_back(Nlm{nlmbas[3 * i], nlmbas[3 * i + 1], nlmbas[3 * i + 2]});
    }
    intparams.update("nlm_basis", std::move(nlm_basis_conv));
  } else if (params.n_max > 0) {
    const int n_max = params.n_max;
    const int l_max = params.l_max == Parameters::all ? n_max - 1 : params.l_max;
    const int m_max = params.m_max == Parameters::all ? l_max : params.m_max;

    assert_throw(l_max < n_max,
                 ExcInvalidParameters("l_max needs to be smaller than n_max"));
    assert_throw(m_max <= l_max,
                 ExcInvalidParameters("m_max cannot be larger than l_max"));

    intparams.update({
          {"n_max", n_max}, {"l_max", l_max}, {"m_max", m_max},
    });
  } else {
    assert_throw(false,
                 ExcInvalidParameters("You need to provide either n_max or nlm_basis if "
                                      "you are using a sturmian basis set."));
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

krims::GenMap build_guess_params(const Parameters& params,
                                 const MolecularSystem& system) {
  using lazyten::EigensystemSolverKeys;
  const bool restricted = parse_restricted(params, system);

  krims::GenMap guess_params{{ScfGuessKeys::method, params.guess}};
  guess_params.update(ScfGuessKeys::eigensolver_params,
                      {{EigensystemSolverKeys::method, params.guess_esolver}});

  if (params.guess == "external") {
    // This is special, since we need to provide the guess externally
    auto& orben_f     = params.guess_external_orben_f;
    auto& orbcoeff_bf = params.guess_external_orbcoeff_bf;

    assert_throw(!orben_f.empty(),
                 ExcInvalidParameters("If guess == external the parameter "
                                      "guess_external_orben_f needs to be "
                                      "initialised."));
    assert_throw(!orbcoeff_bf.empty(),
                 ExcInvalidParameters("If guess == external the parameter "
                                      "guess_external_orbcoeff_bf needs to be "
                                      "initialised."));

    // Determine n_orbs and n_bas from the sizes provided:
    const size_t n_orbs = orben_f.size();
    const size_t n_bas  = orbcoeff_bf.size() / n_orbs;
    assert_throw(n_orbs * n_bas == orbcoeff_bf.size(),
                 ExcInvalidParameters("The sizes of guess_external_orben_f == " +
                                      std::to_string(orben_f.size()) +
                                      " and guess_external_orbcoeff_bf == " +
                                      std::to_string(orbcoeff_bf.size()) +
                                      " are not compatible, since the latter is not an "
                                      "integer multiple of the former."));

    // If we are unrestricted we need to build the block-diagonal structure inside
    // the eigensolution we pass onto the guess_external method and hence the
    // initial Fock matrix
    //
    // Hence the number of elements of each vector in the guess eigensolution
    // is twice the number of basis functions for unrestricted, and once the
    // number for restricted.
    assert_throw(
          restricted || n_orbs % 2 == 0,
          ExcInvalidParameters("For restricted calculations the number of orbitals in "
                               "guess_external_orben_f and guess_external_orbcoeff_bf "
                               "needs to be divisible by 2"));
    const size_t n_orbs_alpha = restricted ? n_orbs : n_orbs / 2;
    const size_t n_evec_elem  = restricted ? n_bas : 2 * n_bas;

    lazyten::Eigensolution<double, lazyten::SmallVector<double>> esolution;
    esolution.evalues() = orben_f;
    esolution.evectors() =
          lazyten::MultiVector<lazyten::SmallVector<double>>(n_evec_elem, n_orbs);

    for (size_t b = 0; b < n_bas; ++b) {
      for (size_t f = 0; f < n_orbs_alpha; ++f) {
        esolution.evectors()[f][b] = orbcoeff_bf[b * n_orbs + f];
      }

      // This loop will not be executed for external guesses for restricted hartree fock
      for (size_t f = n_orbs_alpha; f < n_orbs; ++f) {
        esolution.evectors()[f][b + n_bas] = orbcoeff_bf[b * n_orbs + f];
      }
    }  // b

    guess_params.update(GuessExternalKeys::eigensolution, std::move(esolution));
  }  // method == external

  return guess_params;
}

krims::GenMap build_scf_params(const Parameters& params, const MolecularSystem& system) {
  using lazyten::EigensystemSolverKeys;
  const bool restricted = parse_restricted(params, system);

  assert_throw(
        params.n_eigenpairs % 2 == 0,
        ExcInvalidParameters(
              "The n_eigenpairs parameter applies to the accumulated number of "
              "eigenpairs in the SCF calculations, i.e. the number of alpha plus "
              "the number of beta orbitals. This is the done even for restricted "
              "calculations. For now we further require this number to be even number."));
  const size_t n_eigenpairs = restricted ? params.n_eigenpairs / 2 : params.n_eigenpairs;

  const ScfMsgType verbosity =
        params.print_iterations ? ScfMsgType::IterationProcess : ScfMsgType::Silent;

  krims::GenMap scfparams{
        // error
        {IopScfKeys::max_error_norm, params.conv_tol},
        {IopScfKeys::max_1e_energy_change, params.conv_tol * 100.},
        {IopScfKeys::max_tot_energy_change, params.conv_tol / 4.},
        //
        {IopScfKeys::max_iter, params.max_iter},
        {IopScfKeys::n_eigenpairs, n_eigenpairs},
        {IopScfKeys::verbosity, verbosity},
        {gscf::PulayDiisScfKeys::n_prev_steps, params.diis_size},
  };
  scfparams.update(IopScfKeys::eigensolver_params,
                   {{EigensystemSolverKeys::method, params.eigensolver}});

  return scfparams;
}

}  // namespace iface
}  // namespace molsturm