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

#pragma once
#include <lazyten/io.hh>

namespace molsturm {

/** Wrap around an existing SCF (Pulay, Iop, Plain)
 *  and use the handler calls in order to
 *  extract status information, which is written to a file
 */
template <typename InnerScf>
class ScfDebugWrapper : public InnerScf {
 public:
  typedef InnerScf scf_type;
  typedef typename scf_type::probmat_type probmat_type;
  typedef typename scf_type::scalar_type scalar_type;
  typedef typename scf_type::matrix_type matrix_type;
  typedef typename scf_type::state_type state_type;
  typedef typename scf_type::vector_type vector_type;

  ScfDebugWrapper(lazyten::io::DataWriter_i<scalar_type>& writer,
                  const InnerScf& innerscf)
        : scf_type(innerscf), m_writer(writer) {}

 protected:
  void before_iteration_step(state_type& s) const override {
    scf_type::before_iteration_step(s);

    if (s.n_iter() == 1) {
      write_fock_and_terms("guess", s.problem_matrix());
    }
  }

  matrix_type calculate_error(const state_type& s) const override {
    matrix_type error = scf_type::calculate_error(s);
    m_writer.write("error" + std::to_string(s.n_iter()), error);
    return error;
  }

  void on_update_eigenpairs(state_type& s) const override {
    scf_type::on_update_eigenpairs(s);
    assert_throw(m_writer, krims::ExcIO());
    m_writer.write(
          "evals" + std::to_string(s.n_iter()),
          lazyten::make_as_multivector<vector_type>(s.eigensolution().evalues()));

    auto tmp_copy =
          s.eigensolution().evectors();  // TODO right now we need this. Get rid later.
    m_writer.write("evecs" + std::to_string(s.n_iter()), tmp_copy);
    assert_throw(m_writer, krims::ExcIO());
  }

  void on_update_problem_matrix(state_type& s) const override {
    scf_type::on_update_problem_matrix(s);
    std::string itstr = std::to_string(s.n_iter());
    write_fock_and_terms(itstr, s.problem_matrix());
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
      std::string lala =
            m_writer.normalise_label(kv.first.integral_friendly_name() + "a" + itstr);
      m_writer.write(lala, kv.second);
    }
    for (auto kv : fock_bb.terms_beta()) {
      // Normalise the label: The id of the term may contain funny symbols
      std::string lalb =
            m_writer.normalise_label(kv.first.integral_friendly_name() + "b" + itstr);
      m_writer.write(lalb, kv.second);
    }
    m_writer.write("fock" + itstr, fock_bb);

    for (auto kv : fock_bb.energies()) {
      m_writer.write(m_writer.normalise_label(kv.first.integral_friendly_name() +
                                              "energy" + itstr),
                     kv.second);
    }

    assert_throw(m_writer, krims::ExcIO());
  }

  lazyten::io::DataWriter_i<scalar_type>& m_writer;
};

}  // namespace molsturm
