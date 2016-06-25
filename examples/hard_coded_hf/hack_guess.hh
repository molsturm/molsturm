#pragma once
#include <linalgwrap/ArmadilloMatrix.hh>
namespace hard_coded_hf {
using namespace linalgwrap;

ArmadilloMatrix<double> hack_guess(const ArmadilloMatrix<double>& s_bb) {
  // apply Löwdin normalisation to the basis functions
  //   - Diagonalise the overlap
  //   - Take 1/\sqrt{evals} at the diagonal
  //   - results in orthonormalised basis functions

  const arma::mat& s_bb_data = s_bb.data();

  // Diagonalize s_bb
  arma::vec eval;
  arma::mat evec;
  arma::eig_sym(eval, evec, s_bb_data);

  // take 1/sqrt( . ) for each eigenvalue:
  std::transform(std::begin(eval), std::end(eval), std::begin(eval),
                 [](double elem) { return 1. / std::sqrt(elem); });

  // construct new Armadillo matrix of guess vectors:
  arma::mat guess = evec * arma::diagmat(eval);
  // If number of orbitals is different from number of basis functions,
  // guess should be rectangular.
  //      int b = s_bb.n_rows;
  //      guess = slice(evec * diagmat(eval), b, 5,0,0);

  // Return it properly enwrapped:
  // Note, that our matrices are row-major, but since armadillo
  // is column-major, we need to transpose it first
  return ArmadilloMatrix<double>(guess.t());
}
}
