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

#include "parse_args.hh"
#include <algorithm>
#include <molsturm/read_xyz.hh>
#include <vector>

namespace hf {

std::ostream& operator<<(std::ostream& o, const args_type& args) {
  o << "basis_type:     " << args.basis_type << std::endl;

  if (args.sturmian) {
    o << "k_exp:          " << args.k_exp << '\n';

    if (args.nlm_basis.empty()) {
      o << "n_max:          " << args.n_max << '\n'
        << "l_max:          " << args.l_max << '\n'
        << "m_max:          " << args.m_max << '\n';
    } else {
      o << "nlm_basis:      " << args.nlm_basis << '\n';
    }
  }
  if (args.gaussian) {
    o << "basis_set:      " << args.basis_set << '\n';
  }

  o << "n_alpha:        " << args.system.n_alpha << '\n'
    << "n_beta:         " << args.system.n_beta << '\n'
    << "restricted      " << std::boolalpha << args.restricted << '\n'
    << "guess_method:   " << args.guess_method << '\n'
    << "guess_esolver:  " << args.guess_esolver << '\n'
    << "error:          " << args.error << '\n'
    << "max_iter        " << args.max_iter << '\n'
    << "diis_size:      " << args.diis_size << '\n'
    << "eigensolver     " << args.eigensolver << '\n'
    << "n_eigenpairs:   " << args.n_eigenpairs << '\n'
    << '\n'
    << "structure (distances in bohr):\n"
    << "------------------------\n"
    << args.system.structure << '\n'
    << "------------------------" << std::endl;
  return o;
}

/** Parse an nlm_basis string.
 *
 * Two formats are recognised:
 *    a) tuple format:
 *        (n,l,m)  (n,l,m), ...
 *    b) compressed format:
 *        - Words separated by spaces, first word gives lmax for n=1,
 *          second word for n=2, third for n=3, ..
 *        - Each word consists of an angular momentum letter (s,p,d, ...)
 *          and optionally a count, e.g. 3*s, 2*p, d
 *        - If a count is specified this many shells of increasing nmax,
 *          but of the same lmax are inserted.
 *        - E.g. 3*d 2*p s results in
 *            (n=1, lmax=2), (n=2, lmax=2), (n=3, lmax=2), (n=3, lmax=1),
 *            (n=4, lmax=1), (n=5, lmax=0)
 *          i.e. in the nlm tuples
 *            (1,0,0), (2,0,0), (2,1,-1),(2,1,0), (2,1,1),
 *            (3,0,0), (3,1,-1),(3,1,0), (3,1,1),
 *            (3,2,-2), (3,2,-1),(3,2,0), (3,2,1),(3,2,2),
 *            (4,0,0), (4,1,-1),(4,1,0), (4,1,1),
 *            (5,0,0)
 */
bool parse_nlm_basis(const std::string& str, NlmBasis& basis) {
  using gint::sturmian::atomic::Nlm;
  const std::vector<char> amletters = {'s', 'p', 'd', 'f', 'g', 'h', 'i',
                                       'j', 'k', 'l', 'm', 'n', 'o'};
  auto letter_to_am = [&amletters](char& c) {
    for (size_t l = 0; l < amletters.size(); ++l) {
      if (c == amletters[l]) return static_cast<int>(l);
    }
    return -1;
  };

  basis.clear();
  if (str.find('(') == std::string::npos) {
    // compressed format
    std::stringstream ss(str);
    std::string word;
    while (ss >> word) {
      if (word.size() == 1) {
        const int l = letter_to_am(word[0]);
        if (l == -1) {
          std::cerr << "parse_nlm_basis: Unknown am letter: " << word << std::endl;
          return false;
        }
        basis.add_shell(l);
      } else if (word.find('*') == std::string::npos) {
        std::cerr << "parse_nlm_basis: Unrecognised word: " << word << std::endl;
        return false;
      } else {
        const size_t starpos = word.find('*');
        if (starpos > word.size() - 1) {
          std::cerr << "parse_nlm_basis: Nothing follows a star: " << word << std::endl;
          return false;
        }
        const int l = letter_to_am(word[starpos + 1]);

        size_t times;
        std::stringstream st(word.substr(0, starpos));
        if (!(st >> times)) {
          std::cerr << "parse_nlm_basis: Invalid number before star: " << word
                    << std::endl;
          return false;
        }

        for (size_t i = 0; i < times; ++i) basis.add_shell(l);
      }
    }

    return true;
  }

  // Tuple format
  for (size_t pos = 0; pos < str.size(); ++pos) {
    if (str[pos] != '(') continue;
    const size_t start = pos;
    const size_t sep1 = str.find(',', start);
    const size_t sep2 = str.find(',', sep1);
    const size_t end = str.find(')', sep2);
    const size_t nextopen = str.find('(', sep1);

    // Errors
    if (sep1 == std::string::npos || sep2 == std::string::npos) {
      std::cerr << "parse_nlm_basis: Opening bracket does not start tuple." << std::endl;
      return false;
    }
    if (end == std::string::npos || end > nextopen) {
      std::cerr << "parse_nlm_basis: No closing bracket found." << std::endl;
      return false;
    }
    if (start == sep1 || sep1 == sep2 || sep2 == end) {
      std::cerr << "parse_nlm_basis: One tuple element is empty." << std::endl;
      return false;
    }

    // Parse nlm
    Nlm nlm;
    std::stringstream ssn(str.substr(start + 1, sep1 - start - 1));
    if (!(ssn >> nlm.n)) {
      std::cerr << "parse_nlm_basis: Invalid number for first tuple element."
                << std::endl;
      return false;
    }
    std::stringstream ssl(str.substr(sep1 + 1, sep2 - sep1 - 1));
    if (!(ssl >> nlm.l)) {
      std::cerr << "parse_nlm_basis: Invalid number for second tuple element."
                << std::endl;
      return false;
    }
    std::stringstream ssm(str.substr(sep2 + 1, end - sep2 - 1));
    if (!(ssm >> nlm.m)) {
      std::cerr << "parse_nlm_basis: Invalid number for third tuple element."
                << std::endl;
      return false;
    }
    basis.push_back(std::move(nlm));

    pos = end;  // Update string position
  }

  return true;
}

bool parse_args(int argc, char** argv, args_type& parsed) {
  // System
  bool had_Z_charge = false;
  bool had_xyz = false;
  bool had_alpha = false;
  bool had_beta = false;
  bool had_charge = false;
  bool had_multiplicity = false;
  bool had_atomic_xyz = false;
  bool had_restricted = false;

  // Basis
  bool had_basis_type = false;
  bool had_basis_set = false;
  bool had_k_exp = false;
  bool had_n_max = false;
  bool had_l_max = false;
  bool had_m_max = false;
  bool had_nlm_basis = false;

  // Convergence
  bool had_error = false;
  bool had_max_iter = false;
  bool had_diis_size = false;
  bool had_n_eigenpairs = false;

  // Temporary
  gint::Structure structure;
  size_t n_alpha = 0;
  size_t n_beta = 0;
  size_t multiplicity = 0;
  size_t charge = 0;

  // Parsing
  for (int i = 1; i < argc; ++i) {
    std::string flag(argv[i]);

    // Check that an argument exists:
    if (++i >= argc) {
      std::cerr << "Flag " << flag << " needs an argument." << std::endl;
      return false;
    }

    // Store the argument:
    std::string argument(argv[i]);

    // Interpret argument
    if (flag == std::string("--Z_charge")) {
      had_Z_charge = true;
      if (had_xyz) {
        std::cerr << "Can only have one of --Z_charge or --xyz" << std::endl;
        return false;
      }
      double Z;
      if (!str_to_type<double>(argument, Z)) {
        std::cerr << "Invalid float provided to --Z_charge: " << argument << std::endl;
        return false;
      }
      structure = gint::Structure{{static_cast<unsigned int>(Z), {{0, 0, 0}}}};
    } else if (flag == std::string("--xyz")) {
      had_xyz = true;
      if (had_Z_charge) {
        std::cerr << "Can only have one of --Z_charge or --xyz" << std::endl;
        return false;
      }
      std::ifstream xyz(argument);
      if (!xyz) {
        std::cerr << "Could not open xyz file: " << argument << std::endl;
        return false;
      }
      if (had_atomic_xyz) {
        structure = molsturm::read_xyz(xyz, 1.);
      } else {
        structure = molsturm::read_xyz(xyz);
      }
    } else if (flag == std::string("--atomic_units_xyz")) {
      had_atomic_xyz = true;
      if (had_xyz) {
        std::cerr << "--atomic_units_xyz needs to be specified before --xyz" << std::endl;
        return false;
      }
    } else if (flag == std::string("--restricted")) {
      had_restricted = true;
      if (argument == "true") {
        parsed.restricted = true;
      } else if (argument == "false") {
        parsed.restricted = false;
      } else {
        std::cerr << "--restricted only allows true or false as arguments" << std::endl;
        return false;
      }
    } else if (flag == std::string("--alpha")) {
      had_alpha = true;
      if (had_charge || had_multiplicity) {
        std::cerr << "--charge/--multiplicity and --alpha/--beta are exclusive"
                  << std::endl;
      }
      if (!str_to_type<size_t>(argument, n_alpha)) {
        std::cerr << "Invalid int provided to --alpha: " << argument << std::endl;
        return false;
      }
    } else if (flag == std::string("--beta")) {
      had_beta = true;
      if (had_charge || had_multiplicity) {
        std::cerr << "--charge/--multiplicity and --alpha/--beta are exclusive"
                  << std::endl;
      }
      if (!str_to_type<size_t>(argument, n_beta)) {
        std::cerr << "Invalid int provided to --beta: " << argument << std::endl;
        return false;
      }
    } else if (flag == std::string("--charge")) {
      had_charge = true;
      if (had_alpha || had_beta) {
        std::cerr << "--charge/--multiplicity and --alpha/--beta are exclusive"
                  << std::endl;
      }
      if (!str_to_type<size_t>(argument, charge)) {
        std::cerr << "Invalid int provided to --charge: " << argument << std::endl;
        return false;
      }
    } else if (flag == std::string("--multiplicity")) {
      had_multiplicity = true;
      if (had_alpha || had_beta) {
        std::cerr << "--charge/--multiplicity and --alpha/--beta are exclusive"
                  << std::endl;
      }
      if (!str_to_type<size_t>(argument, multiplicity)) {
        std::cerr << "Invalid int provided to --multiplicity: " << argument << std::endl;
        return false;
      }
    } else if (flag == std::string("--kexp")) {
      had_k_exp = true;
      if (!str_to_type<double>(argument, parsed.k_exp)) {
        std::cerr << "Invalid double provided to --kexp: " << argument << std::endl;
        return false;
      }
    } else if (flag == std::string("--n_max")) {
      if (had_nlm_basis) {
        std::cerr << "--n_max, --l_max, --m_max and --nlm_basis are exclusive."
                  << std::endl;
        return false;
      }
      had_n_max = true;
      if (!str_to_type<int>(argument, parsed.n_max)) {
        std::cerr << "Invalid int provided to --n_max: " << argument << std::endl;
        return false;
      }
    } else if (flag == std::string("--l_max")) {
      if (had_nlm_basis) {
        std::cerr << "--n_max, --l_max, --m_max and --nlm_basis are exclusive."
                  << std::endl;
        return false;
      }
      had_l_max = true;
      if (!str_to_type<int>(argument, parsed.l_max)) {
        std::cerr << "Invalid int provided to --l_max: " << argument << std::endl;
        return false;
      }
    } else if (flag == std::string("--m_max")) {
      if (had_nlm_basis) {
        std::cerr << "--n_max, --l_max, --m_max and --nlm_basis are exclusive."
                  << std::endl;
        return false;
      }
      had_m_max = true;
      if (!str_to_type<int>(argument, parsed.m_max)) {
        std::cerr << "Invalid int provided to --m_max: " << argument << std::endl;
        return false;
      }
    } else if (flag == std::string("--nlm_basis")) {
      if (had_n_max || had_l_max || had_m_max) {
        std::cerr << "--n_max, --l_max, --m_max and --nlm_basis are exclusive."
                  << std::endl;
        return false;
      }
      had_nlm_basis = true;

      const bool res = parse_nlm_basis(argument, parsed.nlm_basis);
      if (!res) {
        std::cerr << "Error when parsing argument string given to --nlm_basis."
                  << std::endl;
        return false;
      }
    } else if (flag == std::string("--basis_set")) {
      had_basis_set = true;
      parsed.basis_set = argument;
    } else if (flag == std::string("--basis_type")) {
      had_basis_type = true;

      const std::vector<std::string> valid_atomic_basis{
            "cs_static14", "cs_dummy", "cs_naive", "cs_reference", "cs_reference_pc"};
      const std::vector<std::string> valid_gaussian_basis{"libint"};
      const auto res_atomic = std::find(std::begin(valid_atomic_basis),
                                        std::end(valid_atomic_basis), argument);
      const auto res_gaussian = std::find(std::begin(valid_gaussian_basis),
                                          std::end(valid_gaussian_basis), argument);

      if (res_atomic != std::end(valid_atomic_basis)) {
        parsed.sturmian = true;
        parsed.gaussian = false;
        parsed.basis_type = "sturmian/atomic/" + argument;
      } else if (res_gaussian != std::end(valid_gaussian_basis)) {
        parsed.gaussian = true;
        parsed.sturmian = false;
        parsed.basis_type = "gaussian/" + argument;
      } else {
        std::cerr << "Invalid argument provided to --basis_type:   " << argument
                  << "valid are:  ";
        for (auto& s : valid_atomic_basis) std::cerr << s << "  ";
        for (auto& s : valid_gaussian_basis) std::cerr << s << "  ";
        std::cerr << std::endl;
        return false;
      }
    } else if (flag == std::string("--error")) {
      had_error = true;
      if (!str_to_type<double>(argument, parsed.error)) {
        std::cerr << "Invalid double provided to --error: " << argument << std::endl;
        return false;
      }
    } else if (flag == std::string("--max_iter")) {
      had_max_iter = true;
      if (!str_to_type<size_t>(argument, parsed.max_iter)) {
        std::cerr << "Invalid double provided to --max_iter: " << argument << std::endl;
        return false;
      }
    } else if (flag == std::string("--diis_size")) {
      had_diis_size = true;
      if (!str_to_type<size_t>(argument, parsed.diis_size)) {
        std::cerr << "Invalid size_t provided to --diis_size: " << argument << std::endl;
        return false;
      }
    } else if (flag == std::string("--n_eigenpairs")) {
      had_n_eigenpairs = true;
      if (!str_to_type<size_t>(argument, parsed.n_eigenpairs)) {
        std::cerr << "Invalid size_t provided to --n_eigenpairs: " << argument
                  << std::endl;
        return false;
      }
      if (parsed.n_eigenpairs % 2 != 0) {
        std::cerr << "n_eigenpairs needs an even number." << std::endl;
        return false;
      }
    } else if (flag == std::string("--eigensolver")) {
      parsed.eigensolver = argument;
    } else if (flag == std::string("--guess_method")) {
      parsed.guess_method = argument;
    } else if (flag == std::string("--guess_esolver")) {
      parsed.guess_esolver = argument;
    } else {
      std::cerr << "Unknown flag: " << flag << std::endl;
      std::cerr << "Valid are: --basis_type, --n_max, --l_max, --n_max, --kexp, "
                   "--Z_charge, --alpha, --beta, --error, --max_iter, --diis_size, "
                   "--n_eigenpairs, --basis_set, --guess_method, --xyz, --charge, "
                   "--multiplicity, --atomic_units_xyz, --eigensolver, --guess_esolver, "
                   "--nlm_basis, --restricted"
                << std::endl;
      return false;
    }
  }

  //
  // Error handling and system setup
  //
  bool error_encountered = false;
  if (!had_xyz && !had_Z_charge) {
    error_encountered = true;
    std::cerr << "Need flag --xyz <file> to supply an xyz file or --Z_charge <double>"
              << std::endl;
  }
  if (had_alpha && !had_beta) {
    n_beta = n_alpha;
  }
  if (had_beta && !had_alpha) {
    error_encountered = true;
    std::cerr << "Need flag --alpha <int> to supply number of alpha electrons."
              << std::endl;
  }

  if (had_charge && !had_multiplicity) {
    parsed.system = molsturm::MolecularSystem(std::move(structure), charge);
  } else if (had_charge && had_multiplicity) {
    parsed.system = molsturm::MolecularSystem(std::move(structure), charge, multiplicity);
  } else if (had_multiplicity) {
    parsed.system = molsturm::MolecularSystem(std::move(structure), 0, multiplicity);
  } else if (had_charge) {
    parsed.system = molsturm::MolecularSystem(std::move(structure), charge);
  } else if (had_alpha) {
    parsed.system = molsturm::MolecularSystem(std::move(structure), {{n_alpha, n_beta}});
  } else {
    parsed.system = molsturm::MolecularSystem(std::move(structure));
  }

  if (!had_basis_type) {
    error_encountered = true;
    std::cerr << "Need flag --basis_type <string> to supply basis type to use."
              << std::endl;
  }
  if (parsed.gaussian) {
    if (!had_basis_set) {
      error_encountered = true;
      std::cerr << "Need flag --basis_set <string> to supply basis set to use."
                << std::endl;
    }
  } else if (parsed.sturmian) {
    if (!had_k_exp) {
      error_encountered = true;
      std::cerr << "Need flag --kexp <double> to supply k exponent." << std::endl;
    }
    if (!had_n_max && !had_nlm_basis) {
      error_encountered = true;
      std::cerr
            << "Need flag --n_max <int> or --nlm_basis <basisstring> to supply a basis."
            << std::endl;
    }
    if (parsed.system.structure.size() != 1) {
      error_encountered = true;
      std::cerr << "Need a molecule with exactly one atom for sturmian calculations."
                << std::endl;
    }
  }

  //
  // Set defaults
  //
  if (!had_error) {
    parsed.error = 5e-7;
  }

  if (!had_max_iter) {
    parsed.max_iter = 25;
  }
  if (!had_diis_size) {
    parsed.diis_size = 4;
  }

  if (!had_n_eigenpairs) {
    parsed.n_eigenpairs = 2 * std::max(parsed.system.n_alpha, parsed.system.n_beta);
  }

  if (!had_l_max && !had_nlm_basis) {
    parsed.l_max = parsed.n_max - 1;
  }
  if (!had_m_max && !had_nlm_basis) {
    parsed.m_max = parsed.l_max;
  }

  if (!had_restricted) {
    if (parsed.system.n_alpha == parsed.system.n_beta) {
      parsed.restricted = true;
    } else {
      parsed.restricted = false;
    }
  }

  return !error_encountered;
}
}  // namespace hf
