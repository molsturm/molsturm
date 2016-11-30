#pragma once
#include "detail/IopScfWrapper.hh"
#include <gscf/PlainScf.hh>
#include <gscf/PulayDiisScf.hh>

namespace molsturm {

// TODO Maybe a common base between this guy and linalgwrap's EigensystemSolver is
// sensible?

template <typename ProblemMatrix, typename OverlapMatrix>
struct IopScfState : public gscf::ScfStateBase<ProblemMatrix, OverlapMatrix> {
  typedef gscf::ScfStateBase<ProblemMatrix, OverlapMatrix> base_type;
  typedef typename base_type::probmat_type probmat_type;
  typedef typename base_type::overlap_type overlap_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::real_type real_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::vector_type vector_type;

  /** The orbital energies of the SCF */
  const std::vector<scalar_type>& orbital_energies() const {
    return base_type::eigensolution().evalues();
  }

  /** The orbital coefficients of the SCF */
  const linalgwrap::MultiVector<vector_type>& orbital_coeff() const {
    return base_type::eigensolution().evectors();
  }

  IopScfState(probmat_type probmat, const overlap_type& overlap_mat)
        : base_type{std::move(probmat), overlap_mat}, m_n_iter(0), m_n_mtx_applies(0) {}

  /** Transfer data from an arbitrary other state by copying
   * it into this one. This is needed, since we need to get the
   * data back here once the inner solvers are done with
   * (partially) solving the SCF.
   *
   * \note Do not use this to setup a guess state.
   * For this case the function obtain_guess_from exists.
   **/
  //@{
  template <typename DiagMat>
  void push_intermediate_results(
        gscf::ScfStateBase<probmat_type, overlap_type, DiagMat>&& other_state) {
    static_cast<linalgwrap::SolverStateBase&>(*this) = other_state;

    // Copy things from ScfStateBase
    // TODO This is all but good. The DiagMat template argument really has to
    //      disappear (or we need yet a layer ... a common base independent of
    //      DiagMat.
    base_type::last_error_norm = other_state.last_error_norm;
    base_type::problem_matrix_ptr = other_state.problem_matrix_ptr;
    base_type::push_new_eigensolution(other_state.previous_eigensolution(),
                                      other_state.eigenproblem_stats());
    base_type::push_new_eigensolution(other_state.eigensolution(),
                                      other_state.eigenproblem_stats());

    // Copy extra things from us:
    m_n_iter += other_state.n_iter();
    m_n_mtx_applies += other_state.n_mtx_applies();
  }

  void push_intermediate_results(base_type&& other_state) {
    static_cast<base_type&>(*this) = other_state;
    m_n_iter += other_state.n_iter();
    m_n_mtx_applies += other_state.n_mtx_applies();
  }
  //@}

  size_t n_iter() const override { return m_n_iter; }
  size_t n_mtx_applies() const override { return m_n_mtx_applies; }

 private:
  size_t m_n_iter;
  size_t m_n_mtx_applies;
};

/** Scf Solver which should be used for integral operators.
 *
 * TODO
 * Right now this is exactly the PulayDiisScf. Later we probably want
 * something more clever, which first does some DIIS, but later
 * switches it off (or for something else) when all of a sudden
 * the errors get to small (and the linear system to be solved in the
 * DIIS becomes ill-conditioned
 */
template <typename IntegralOperator, typename OverlapMatrix>
class IopScf : public gscf::ScfBase<IopScfState<IntegralOperator, OverlapMatrix>> {
 public:
  typedef IntegralOperator operator_type;
  typedef gscf::ScfBase<IopScfState<IntegralOperator, OverlapMatrix>> base_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::real_type real_type;
  typedef typename base_type::state_type state_type;
  typedef typename base_type::probmat_type probmat_type;
  typedef typename base_type::overlap_type overlap_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::vector_type vector_type;

  /** Construct scf from parameter map */
  IopScf(const krims::ParameterMap& map) { update_control_params(map); }

  /** \name Iteration control */
  ///@{
  //! Maximal total energy change between two cycles.
  real_type max_tot_energy_change = 1e-8;

  //! Maximal 1e energy change between two cycles.
  real_type max_1e_energy_change = 1e-5;

  //! How verbose should the solver be.
  ScfMsgType verbosity = ScfMsgType::Silent;

  // TODO introduce method parameter!

  /** Update control parameters from Parameter map */
  void update_control_params(const krims::ParameterMap& map) {
    base_type::update_control_params(map);

    max_tot_energy_change =
          map.at(IopScfKeys::max_tot_energy_change, max_tot_energy_change);
    max_1e_energy_change = map.at(IopScfKeys::max_1e_energy_change, max_1e_energy_change);
    verbosity = map.at(IopScfKeys::verbosity, verbosity);

    // Copy the map to the internal storage such that we
    // can pass it on to the actual eigensolvers.
    m_inner_params = map;
  }

  /** Get the current settings of all internal control parameters and
   *  update the ParameterMap accordingly.
   */
  void get_control_params(krims::ParameterMap& map) const {
    base_type::get_control_params(map);
    map.update(IopScfKeys::max_tot_energy_change, max_tot_energy_change);
    map.update(IopScfKeys::max_1e_energy_change, max_1e_energy_change);
    map.update(IopScfKeys::verbosity, verbosity);
  }
  ///@}

  /** Solve the state. Uses different SCF algorithms depending
   *  on certain properties of the iteration process */
  void solve_state(state_type& state) const override;

 private:
  /** Cache of the parameters which will be passed to the actual SCF.
   *
   * A mixture of the original parameters, which the user passed to us
   * and the current state which is reflected in the member variables
   * in this class and the subclasses. Should be updated with
   * get_control_params *before* the actual inner eigensolver invocation.
   */
  mutable krims::ParameterMap m_inner_params;
};

/** Try to find the self consistent field configuration for an integral operator
 *  and an overlap matrix
 */
template <typename IntegralOperator, typename OverlapMatrix>
typename IopScf<IntegralOperator, OverlapMatrix>::state_type run_scf(
      IntegralOperator iop, const OverlapMatrix& s, krims::ParameterMap map) {
  typedef IopScf<IntegralOperator, OverlapMatrix> scf;
  return scf{map}.solve(std::move(iop), s);
}

//
// -----------------------------------------------------------
//

template <typename IntegralOperator, typename OverlapMatrix>
void IopScf<IntegralOperator, OverlapMatrix>::solve_state(state_type& state) const {
  using namespace gscf;
  assert_dbg(!state.is_failed(), krims::ExcInvalidState("Cannot solve a failed state"));

  // For now always do DIIS:
  typedef detail::IopScfWrapper<PulayDiisScf<
        detail::IopScfStateWrapper<PulayDiisScfState<IntegralOperator, OverlapMatrix>>>>
        Solver;

  //
  //
  //

  typename Solver::state_type inner_state{state.problem_matrix(), state.overlap_matrix()};
  inner_state.obtain_guess_from(state);

  // Make sure that control parameters are up to date:
  get_control_params(m_inner_params);

  try {
    Solver{m_inner_params}.solve_state(inner_state);
    state.push_intermediate_results(std::move(inner_state));
  } catch (linalgwrap::SolverException& e) {
    // On exception still update the state reference
    state.push_intermediate_results(std::move(inner_state));
    throw;
  }
}

}  // namespace molsturm
