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

#include "export_hf_results.hh"
#include <molsturm/FockOperator.hh>
#include <molsturm/IopScf.hh>
#include <molsturm/OverlapMatrix.hh>

namespace molsturm {
namespace iface {

namespace {

void export_vector(const bool restricted, const std::vector<scalar_type>& in,
                   std::vector<scalar_type>& out) {
  std::copy(in.begin(), in.end(), std::back_inserter(out));

  if (restricted) {
    // The above statement only inserts the alpha energies,
    // we replicate the statement to insert the same energies (as betas)
    // as well
    std::copy(in.begin(), in.end(), std::back_inserter(out));
  }
}

template <typename Matrix>
void export_ff_matrix(const bool restricted, const Matrix& in,
                      std::vector<scalar_type>& out) {
  assert_size(in.n_rows(), in.n_cols());
  const size_t n_orbs = restricted ? 2 * in.n_rows() : in.n_rows();

  // Shift from the alpha-alpha into the beta-beta block
  // if we are restricted (otherwise this has no meaning)
  const size_t n_orbs_alpha = in.n_cols();
  const size_t bb_shift = n_orbs_alpha * n_orbs + n_orbs_alpha;

  out.resize(n_orbs * n_orbs);
  for (size_t i = 0; i < in.n_rows(); ++i) {
    for (size_t j = 0; j < in.n_cols(); ++j) {
      const size_t ij = i * n_orbs + j;

      assert_internal(ij < out.size());
      out[ij] = in(i, j);

      if (restricted) {
        // Shift index and copy into beta-beta block as well.
        assert_internal(ij + bb_shift < out.size());
        out[ij + bb_shift] = in(i, j);
      }
    }  // j
  }    // i
}

template <typename Vector>
void export_coeff_fb(const bool restricted, const linalgwrap::MultiVector<Vector>& mv,
                     std::vector<scalar_type>& out) {
  const size_t n_orbs = restricted ? 2 * mv.n_vectors() : mv.n_vectors();
  const size_t n_bas = restricted ? mv.n_elem() : mv.n_elem() / 2;
  const size_t n_orbs_alpha = restricted ? mv.n_vectors() : mv.n_vectors() / 2;
  assert_internal(n_orbs_alpha * 2 == n_orbs);

  // Loop over alpha-alpha block
  out.resize(n_orbs * n_bas);
  for (size_t f = 0; f < n_orbs_alpha; ++f) {
    for (size_t b = 0; b < n_bas / 2; ++b) {
      const size_t fb = f * n_bas + b;
      const size_t fb_beta = (f + n_orbs_alpha) * n_bas + b;
      assert_internal(fb < out.size());
      assert_internal(fb_beta < out.size());

      // Copy the alpha-alpha block to the new location
      out[fb] = mv[f][b];

      // Also copy to the beta block. Either take another copy
      // of the alpha block or shift the f and b indices to get
      // into the beta block
      out[fb_beta] = restricted ? mv[f][b] : mv[f + n_orbs_alpha][b + n_bas];
    }  // b
  }    // f
}

template <typename Vector>
void export_eri_restricted(const gint::ERITensor_i<scalar_type>& eri,
                           const linalgwrap::MultiVector<Vector>& coeff_fb,
                           std::vector<scalar_type>& out) {
  const size_t n_orbs_alpha = coeff_fb.n_vectors();
  const size_t n_orbs = 2 * n_orbs_alpha;

  // Form the alpha-alpha-alpha-alpha spin block by contraction
  // and copy it to all the other non-zero places
  // i.e. (aa|aa), (aa|bb), (bb|aa), (bb|bb)
  std::vector<double> eri_aaaa(n_orbs_alpha * n_orbs_alpha * n_orbs_alpha * n_orbs_alpha);
  eri.contract_with(coeff_fb, coeff_fb, coeff_fb, coeff_fb, eri_aaaa);

  out.resize(n_orbs * n_orbs * n_orbs * n_orbs);
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
          // equal to the element at index aaaa by spin symmetry.
          const size_t aabb = aaaa + bb_off;
          const size_t bbaa = aaaa + n_orbs * n_orbs * bb_off;
          const size_t bbbb = aaaa + bb_off + n_orbs * n_orbs * bb_off;

          assert_internal(aaaa < out.size());
          assert_internal(aabb < out.size());
          assert_internal(bbaa < out.size());
          assert_internal(bbbb < out.size());
          assert_internal(ijkl < eri_aaaa.size());

          out[aaaa] = eri_aaaa[ijkl];
          out[aabb] = eri_aaaa[ijkl];
          out[bbaa] = eri_aaaa[ijkl];
          out[bbbb] = eri_aaaa[ijkl];
        }  // l
      }    // k
    }      // j
  }        // i
}

/** Split a block-diagonal coefficient multivector into its consituent blocks by copying
 */
template <typename Vector>
std::pair<linalgwrap::MultiVector<Vector>, linalgwrap::MultiVector<Vector>> coeff_blocks(
      const linalgwrap::MultiVector<Vector>& coeff_bf_full) {
  const size_t n_orbs_alpha = coeff_bf_full.n_vectors() / 2;  // == n_orbs_beta
  const size_t n_bas = coeff_bf_full.n_elem() / 2;
  assert_internal(n_orbs_alpha * 2 == coeff_bf_full.n_vectors());
  assert_internal(2 * n_bas == coeff_bf_full.n_elem());

  // Partition coeff_bf into alpha-alpha and beta-beta blocks
  linalgwrap::MultiVector<Vector> ca_bf(n_bas, n_orbs_alpha);
  for (size_t f = 0; f < n_orbs_alpha; ++f) {
    assert_internal(n_bas <= coeff_bf_full[f].size());
    assert_internal(n_bas <= ca_bf[f].size());
    std::copy(coeff_bf_full[f].begin(), coeff_bf_full[f].begin() + n_bas,
              ca_bf[f].begin());
  }

  linalgwrap::MultiVector<Vector> cb_bf(n_bas, n_orbs_alpha);
  for (size_t f = 0; f < n_orbs_alpha; ++f) {
    const size_t ff = f + n_orbs_alpha;
    for (size_t b = 0; b < n_bas; ++b) {
      assert_internal(std::distance(coeff_bf_full[ff].begin() + n_bas,
                                    coeff_bf_full[ff].end()) == cb_bf[f].size());
      std::copy(coeff_bf_full[ff].begin() + n_bas, coeff_bf_full[ff].end(),
                cb_bf[f].begin());
    }
  }

  return std::make_pair(std::move(ca_bf), std::move(cb_bf));
}

template <typename Vector>
void export_eri_unrestricted(const gint::ERITensor_i<scalar_type>& eri,
                             const linalgwrap::MultiVector<Vector>& coeff_fb,
                             std::vector<scalar_type>& out) {
  const size_t n_orbs_alpha = coeff_fb.n_vectors() / 2;  // == n_orbs_beta
  const size_t n_orbs = 2 * n_orbs_alpha;
  const size_t n_bas = coeff_fb.n_elem() / 2;

  assert_internal(n_orbs == coeff_fb.n_vectors());
  assert_internal(2 * n_bas == coeff_fb.n_elem());

  // Build the relevant blocks of the eri tensor
  std::vector<double> eri_aaaa(n_orbs_alpha * n_orbs_alpha * n_orbs_alpha * n_orbs_alpha);
  std::vector<double> eri_aabb(n_orbs_alpha * n_orbs_alpha * n_orbs_alpha * n_orbs_alpha);
  std::vector<double> eri_bbbb(n_orbs_alpha * n_orbs_alpha * n_orbs_alpha * n_orbs_alpha);
  std::vector<double> eri_bbaa(n_orbs_alpha * n_orbs_alpha * n_orbs_alpha * n_orbs_alpha);

  {
    // Partition into alpha-alpha and beta-beta blocks
    std::pair<MultiVector<Vector>, MultiVector<Vector>> cab_bf = coeff_blocks(coeff_bf);
    const auto& ca_bf = cab_bf.first;
    const auto& cb_bf = cab_bf.second;

    // Contract with the eri tensor
    eri.contract_with(ca_fb, ca_fb, ca_fb, ca_fb, eri_aaaa);
    eri.contract_with(ca_fb, ca_fb, cb_fb, cb_fb, eri_aabb);
    eri.contract_with(cb_fb, cb_fb, cb_fb, cb_fb, eri_bbbb);

    // TODO This contraction could be avoided if one uses the property that
    //      ( ij | kl ) = ( kl | ij ) to construct this thing from eri_aabb
    eri.contract_with(cb_fb, cb_fb, ca_fb, ca_fb, eri_bbaa);
  }

  out.resize(n_orbs * n_orbs * n_orbs * n_orbs);
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
          // which we fill here
          const size_t aabb = aaaa + bb_off;
          const size_t bbaa = aaaa + n_orbs * n_orbs * bb_off;
          const size_t bbbb = aaaa + bb_off + n_orbs * n_orbs * bb_off;

          assert_internal(aaaa < out.size());
          assert_internal(aabb < out.size());
          assert_internal(bbaa < out.size());
          assert_internal(bbbb < out.size());
          assert_internal(ijkl < eri_aaaa.size());

          out[aaaa] = eri_aaaa[ijkl];
          out[aabb] = eri_aabb[ijkl];
          // out[bbaa] = eri_aabb[ijkl];  // == eri_bbaa by symmetry
          out[bbaa] = eri_bbaa[ijkl];
          out[bbbb] = eri_bbbb[ijkl];
        }  // l
      }    // k
    }      // j
  }        // i
}

template <typename Vector>
void export_eri(const bool restricted, const gint::ERITensor_i<scalar_type>& eri,
                const linalgwrap::MultiVector<Vector>& coeff_bf,
                std::vector<scalar_type>& out) {
  if (restricted) {
    export_eri_restricted(eri, coeff_bf, out);
  } else {
    export_eri_unrestricted(eri, coeff_bf, out);
  }
}

template <typename Vector, typename OverlapMatrix>
typename Vector::scalar_type compute_spin_squared(
      const bool restricted, const OverlapMatrix& S_bb,
      const linalgwrap::MultiVector<Vector>& coeff_bf, const size_t n_alpha,
      const size_t n_beta) {
  using linalgwrap::MultiVector;
  using krims::range;

  // If restricted the <S^2> is always zeros
  if (restricted) return 0.;

  // Partition into alpha-alpha and beta-beta blocks
  // and filter out the occupied alpha and occupied beta coefficients
  std::pair<MultiVector<Vector>, MultiVector<Vector>> cab_bf = coeff_blocks(coeff_bf);
  const MultiVector<Vector> coa_bf = cab_bf.first.subview(range(n_alpha));
  const MultiVector<Vector> cob_bf = cab_bf.second.subview(range(n_beta));

  // Build the alpha-beta block of the overlap matrix in MO basis
  const auto& Sa_bb = S_bb.block_alpha();
  auto Sab_ff = linalgwrap::dot(cob_bf, Sa_bb * coa_bf);

  // Compute the exact value for <S^2> which we would expect for this
  // system if the determinant was an eigenfunction of S^2
  assert_internal(n_alpha >= n_beta);
  const double spin_total = (n_alpha - n_beta) / 2.;
  const double spin_squared_exact = spin_total * (1 + spin_total);

  // According to Szabo-Ostlund, p. 107 (2.271) this is the actual value for <S^2>
  return spin_squared_exact + n_beta - linalgwrap::norm_frobenius_squared(Sab_ff);
}

}  // namespace

template <typename State>
HfResults export_hf_results(const State& state, const gint::ERITensor_i<scalar_type>& eri,
                            const Parameters& params) {
  const auto& fbb = state.problem_matrix();
  const auto& soln = state.eigensolution();
  const bool restricted = fbb.restricted();

  const size_t n_bas = restricted ? fbb.n_rows() : fbb.n_rows() / 2;
  const size_t n_orbs_alpha = restricted ? soln.n_ep() : soln.n_ep() / 2;
  const size_t n_orbs_beta = restricted ? n_orbs_alpha : soln.n_ep() / 2;
#ifdef DEBUG
  const size_t n_orbs = n_orbs_alpha + n_orbs_beta;
#endif

  HfResults ret;

  // Size information
  ret.n_beta = fbb.indices_orbspace(gscf::OrbitalSpace::OCC_BETA).length();
  ret.n_alpha = fbb.indices_orbspace(gscf::OrbitalSpace::OCC_ALPHA).length();
  ret.n_bas = n_bas;
  ret.restricted = fbb.restricted();
  ret.n_orbs_alpha = n_orbs_alpha;
  ret.n_orbs_beta = n_orbs_beta;

  // SCF statistics
  ret.threshold = params.error;
  ret.n_iter = state.n_iter();
  ret.n_mtx_applies = state.n_mtx_applies();
  ret.final_error_norm = state.last_error_norm;
  ret.final_tot_energy_change = state.last_tot_energy_change;
  ret.final_1e_energy_change = state.last_1e_energy_change;

  // SCF analysis:
  ret.spin_squared = compute_spin_squared(restricted, state.overlap_matrix(),
                                          soln.evectors(), ret.n_alpha, ret.n_beta);

  // SCF energies:
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

  // Insert alpha and beta orbital energies and the coefficient matrices
  export_vector(restricted, soln.evalues(), ret.orben_f);
  export_coeff_fb(restricted, soln.evectors(), ret.orbcoeff_fb);
  assert_internal(ret.orben_f.size() == n_orbs);
  assert_internal(ret.orbcoeff_bf.size() == n_orbs * n_bas /* *2 */);

  if (params.export_overlap_matrix) {
    // Compute the overlap matrix in MO space, i.e. C^T (S * C)
    auto overlap_ff =
          linalgwrap::dot(soln.evectors(), state.overlap_matrix() * soln.evectors());
    export_ff_matrix(restricted, overlap_ff, ret.overlap_ff);
    assert_internal(ret.overlap_ff.size() == n_orbs * n_orbs);
  } else {
    // Empty vector objects won't proceed to the output python dictionary
    ret.overlap_ff.resize(0);
  }

  if (params.export_fock_matrix) {
    // Compute the full fock matrix in MO space, i.e.  C^T * (F * C)
    // => Need a dot product here, so actually the dot of all vectors with another
    auto fock_ff = linalgwrap::dot(soln.evectors(), fbb * soln.evectors());
    export_ff_matrix(restricted, fock_ff, ret.fock_ff);
    assert_internal(ret.fock_ff.size() == n_orbs * n_orbs);
  } else {
    // Empty vector objects won't proceed to the output python dictionary
    ret.fock_ff.resize(0);
  }

  if (params.export_hcore_matrix) {
    // Build the alpha-alpha block of the one electron terms in atomic basis function
    // space
    // (Note: This equals the beta-beta block for the one electron terms)
    linalgwrap::LazyMatrixSum<matrix_type> hcore_block;
    for (const auto& id_term : fbb.terms_1e()) hcore_block += id_term.second;

    // The full hcore matrix (alpha-alpha and beta-beta block along the diagonal)
    linalgwrap::BlockDiagonalMatrix<linalgwrap::LazyMatrixSum<matrix_type>, 2> hcore_bb{
          {{hcore_block, hcore_block}}};

    // Compute the full one electron matrix in MO space, i.e.  C^T * (Hcore_bb * C)
    // => Need a dot product here, so actually the dot of all vectors with another
    //    For restricted, where evectors() only runs over the alpha orbitals (betas are
    //    identical)
    //    we only need the alpha-alpha block. For unrestricted, where the alphas and betas
    //    might
    //    differ, we need both blocks.
    auto hcore_ff =
          restricted ? linalgwrap::dot(soln.evectors(), hcore_block * soln.evectors())
                     : linalgwrap::dot(soln.evectors(), hcore_bb * soln.evectors());
    export_ff_matrix(restricted, hcore_ff, ret.hcore_ff);
    assert_internal(ret.hcore_ff.size() == n_orbs * n_orbs);
  } else {
    ret.hcore_ff.resize(0);
  }

  if (params.export_repulsion_integrals) {
    export_eri(restricted, eri, soln.evectors(), ret.eri_ffff);
    assert_internal(ret.eri_ffff.size() == n_orbs * n_orbs * n_orbs * n_orbs);
  } else {
    ret.eri_ffff.resize(0);
  }

  return ret;
}

#define INSTANTIATE(STATE_TYPE)                                           \
  template HfResults export_hf_results(                                   \
        const IopScfState<FockOperator<matrix_type, STATE_TYPE>,          \
                          OverlapMatrix<matrix_type, STATE_TYPE>>& state, \
        const gint::ERITensor_i<scalar_type>& eri, const Parameters&)

INSTANTIATE(RestrictionType::Unrestricted);
INSTANTIATE(RestrictionType::RestrictedClosed);

#undef INSTANTIATE

}  // namespace iface
}  // namespace molsturm
