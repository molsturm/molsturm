#pragma once
#include <iostream>
#include <molsturm/MolecularSystem.hh>
#include <sstream>
#include <string>

namespace hf {
struct args_type {
  // The system we model
  molsturm::MolecularSystem system;

  // Basis
  std::string basis_type;

  // Sturmians
  bool sturmian = false;
  double k_exp = 1.0;
  size_t n_max = 3;
  size_t l_max = 2;
  size_t m_max = 2;

  // Gaussians
  bool gaussian = false;
  std::string basis_set = "<not avail>";

  // Convergence
  size_t max_iter = 25;
  double error = 5e-7;
  size_t diis_size = 4;
  size_t n_eigenpairs = 0;
  std::string eigensolver = "auto";
  std::string guess_esolver = "auto";
  std::string guess_method = "hcore";
};

/** Write the content of args_type to a stream */
std::ostream& operator<<(std::ostream& o, const args_type& args);

/** Quick and dirty function to parse a string to a different type.
 *  Return false if not possible */
template <typename T>
bool str_to_type(const std::string& in, T& out) {
  return static_cast<bool>(std::stringstream(in) >> out);
}

/** \brief Quick and dirty function to parse the commandline arguments
 *
 * \returns true if all is fine, else false
 */
bool parse_args(int argc, char** argv, args_type& parsed);

}  // namespace hf
