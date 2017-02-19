#pragma once
#include <linalgwrap/TypeUtils.hh>
#include <linalgwrap/eigensystem.hh>
#include <linalgwrap/rescue.hh>
#include <string>

namespace molsturm {

// Idea:   Do Sturmian calculation with minimal nmax as a dense matrix
//         for a guess. Potentially even store the result *hard coded*
//         in the program (That's how ORCA does it).
//
//         max(n_alpha+n_beta)    nmax
//                          1     1
//                     <=   5     2
//                     <=  14     3
//                     <=  30     4
//                     <=  55     5
//                     <=  91     6

// TODO Code common interface
//      Ideally Integral operator -> guess
//      I.e. we would supply an (uninitialised) Fock matrix to compute the guess from
//
//      Check how ORCA and Q-Chem do this and mimic their behaviour

DefException1(ExcObtainingScfGuessFailed, std::string,
              << "Obtaining the guess for the SCF failed. Reason: " << arg1);

/** Obtain a core hamiltonian guess: Done by diagonalising the Hcore matrix */
template <typename IntegralOperator, typename OverlapMatrix>
linalgwrap::EigensolutionTypeFor<true, IntegralOperator> hcore_guess(
      const IntegralOperator& fock_bb, const OverlapMatrix& S_bb) {
  using namespace linalgwrap;

  // Setup Hcore (i.e. only take one electron terms from Fock)
  LazyMatrixSum<typename IntegralOperator::stored_matrix_type> hcore;
  for (const auto& term : fock_bb.terms_1e()) hcore += term;

  // Solve eigensystem for smallest real eigenvalues
  const size_t n_vectors = std::max(fock_bb.n_alpha(), fock_bb.n_beta());
  krims::GenMap params{{EigensystemSolverKeys::which, "SR"}};
  try {
    return eigensystem_hermitian(hcore, S_bb, n_vectors, params);
  } catch (const SolverException& e) {
    rescue::failed_eigenproblem(
          Eigenproblem<true, decltype(hcore), OverlapMatrix>(hcore, S_bb), params);
    assert_throw(false, ExcObtainingScfGuessFailed(
                              "Eigensolver for Hcore failed with message " + e.extra()));
    return EigensolutionTypeFor<true, IntegralOperator>{};
  }
}

/** Obtain a Löwdin guess: Done by diagonalising the overlap
 *
 * If obtaining the guess failed, the function throws
 * ExcObtainingScfGuessFailed.
 */
template <typename IntegralOperator, typename OverlapMatrix>
linalgwrap::EigensolutionTypeFor<true, IntegralOperator> loewdin_guess(
      const IntegralOperator& fock_bb, const OverlapMatrix& S_bb) {
  using namespace linalgwrap;

  // apply Löwdin normalisation to the basis functions
  //   - Diagonalise the overlap
  //   - Take 1/\sqrt{evals} at the diagonal
  //   - results in orthonormalised basis functions

  // Solve eigensystem for largest real eigenvalues
  krims::GenMap params{{EigensystemSolverKeys::which, "LR"}};
  const size_t n_vectors = std::max(fock_bb.n_alpha(), fock_bb.n_beta());

  try {
    auto sol = eigensystem_hermitian(S_bb, n_vectors, params);

    // Eigenvectors and eigenvalues.
    auto& evectors = sol.evectors();
    auto& evalues = sol.evalues();

    if (n_vectors != Constants<size_t>::all) {
      assert_dbg(evectors.n_vectors() == n_vectors, krims::ExcInternalError());
    }
    assert_dbg(evectors.n_elem() == S_bb.n_cols(), krims::ExcInternalError());

    for (size_t i = 0; i < evectors.n_vectors(); ++i) {
      evectors[i] *= 1. / sqrt(evalues[i]);
      evalues[i] = 1.;
    }
    return sol;
  } catch (const SolverException& e) {
    rescue::failed_eigenproblem(Eigenproblem<true, OverlapMatrix>(S_bb), params);
    assert_throw(false,
                 ExcObtainingScfGuessFailed(
                       "Eigensolver for overlap failed with message " + e.extra()));
  }
}

}  // namespace dummy_scf
