#pragma once
#include <linalgwrap/io.hh>

/** Wrap around an existing SCF and use the handler calls in order to
 *  extract status information, which is written to a file
 */
template <typename InnerScf>
class DebugScfWrapper : public InnerScf {
public:
  typedef InnerScf scf_type;
  typedef typename scf_type::probmat_type probmat_type;
  typedef typename scf_type::scalar_type scalar_type;
  typedef typename scf_type::state_type state_type;
  typedef typename scf_type::vector_type vector_type;

  DebugScfWrapper(linalgwrap::io::DataWriter_i<scalar_type>& writer,
                  const InnerScf& innerscf)
        : scf_type(innerscf), m_writer(writer) {}

protected:
  void before_iteration_step(state_type& s) const override {
    scf_type::before_iteration_step(s);

    if (s.n_iter_count() == 1) {
      write_fock_and_terms("guess", *s.problem_matrix_ptr());
    }
  }

  void on_update_eigenpairs(state_type& s) const override {
    scf_type::on_update_eigenpairs(s);
    assert_throw(m_writer, krims::ExcIO());
    m_writer.write("evals" + std::to_string(s.n_iter_count()),
                   linalgwrap::make_as_multivector<vector_type>(*s.eigenvalues_ptr()));
    m_writer.write("evecs" + std::to_string(s.n_iter_count()), *s.eigenvectors_ptr());
    assert_throw(m_writer, krims::ExcIO());
  }

  void on_update_problem_matrix(state_type& s) const override {
    scf_type::on_update_problem_matrix(s);
    std::string itstr = std::to_string(s.n_iter_count());
    write_fock_and_terms(itstr, *s.problem_matrix_ptr());
  }

private:
  /** Write both the inner terms of the fock matrix as well as the matrix itself to the
   *  DataWriter
   *
   * \param itstr    An identifier with is appended to the identifier of the terms or
   *                 "fock"
   **/
  void write_fock_and_terms(const std::string& itstr, const probmat_type& fock_bb) const {
    assert_throw(m_writer, krims::ExcIO());

    for (auto kv : fock_bb.terms_alpha()) {
      // Normalise the label: The id of the term may contain funny symbols
      std::string lala = m_writer.normalise_label(kv.first + "a" + itstr);
      m_writer.write(lala, kv.second);
    }
    for (auto kv : fock_bb.terms_beta()) {
      // Normalise the label: The id of the term may contain funny symbols
      std::string lalb = m_writer.normalise_label(kv.first + "b" + itstr);
      m_writer.write(lalb, kv.second);
    }
    m_writer.write("fock" + itstr, fock_bb);

    assert_throw(m_writer, krims::ExcIO());
  }

  linalgwrap::io::DataWriter_i<scalar_type>& m_writer;
};
