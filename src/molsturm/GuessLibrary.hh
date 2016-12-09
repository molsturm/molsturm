#pragma once
#include <linalgwrap/TypeUtils.hh>
#include <linalgwrap/eigensystem.hh>
#include <linalgwrap/rescue.hh>
#include <string>

namespace molsturm {

DefException1(ExcObtainingGuessFailed, std::string,
              << "Obtaining the guess for the SCF failed. Reason: " << arg1);

/** Obtain a Löwdin guess: Done by diagonalising the overlap
 *
 * Optionally one can tell the method that only a certain number of
 * guess vectors is required in order to reduce the effort.
 *
 * If obtaining the guess failed, the function throws
 * ExcObtainingGuessFailed.
 */
template <typename Matrix>
linalgwrap::MultiVector<typename linalgwrap::StoredTypeOf<Matrix>::type::vector_type>
loewdin_guess(const Matrix& overlap_bb,
              const size_t n_vectors = linalgwrap::Constants<size_t>::all) {
  using namespace linalgwrap;

  // apply Löwdin normalisation to the basis functions
  //   - Diagonalise the overlap
  //   - Take 1/\sqrt{evals} at the diagonal
  //   - results in orthonormalised basis functions

  // Solve eigensystem for smallest real eigenvalues
  krims::ParameterMap params{{EigensystemSolverKeys::which, "LR"}};

  try {
    auto sol = eigensystem_hermitian(overlap_bb, n_vectors, params);

    // Eigenvectors and eigenvalues.
    auto& evectors = sol.evectors();
    const auto& evalues = sol.evalues();

    if (n_vectors != Constants<typename Matrix::size_type>::all) {
      assert_dbg(evectors.n_vectors() == n_vectors, krims::ExcInternalError());
    }
    assert_dbg(evectors.n_elem() == overlap_bb.n_cols(), krims::ExcInternalError());

    for (size_t i = 0; i < evectors.n_vectors(); ++i) {
      evectors[i] *= 1. / sqrt(evalues[i]);
    }

    return std::move(sol.evectors());
  } catch (const SolverException& e) {
    rescue::failed_eigenproblem(Eigenproblem<true, Matrix>(overlap_bb), params);
    assert_throw(false,
                 ExcObtainingGuessFailed("Eigensolver for overlap failed with message " +
                                         e.extra()));
  }
}

}  // namespace dummy_scf
