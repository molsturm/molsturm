#include "hack_guess.hh"
#include <gint/IntegralLookup.hh>
#include <gint/version.hh>
#include <gscf/version.hh>
#include <iostream>
#include <linalgwrap/io.hh>
#include <linalgwrap/version.hh>
#include <molsturm/IopPlainScfSolver.hh>
#include <molsturm/RestrictedClosedIntegralOperator.hh>
#include <molsturm/version.hh>

namespace hard_coded_hf {
using namespace molsturm;
using namespace linalgwrap;

/** Run an (atomic) SCF based on the integral data Sturmian14.
 *
 * \param Z       Number of nuclei
 * \param k_exp   Sturmian exponent
 * \param n_alpha Number of alpha electrons
 * \param n_beta  Number of beta electrons
 */
void run_rhf_sturmian_debug(double Z, double k_exp, size_t n_alpha,
                            size_t n_beta) {
  //
  // Types and settings
  //
  // Types of scalar and matrix
  typedef double scalar_type;
  typedef SmallMatrix<scalar_type> stored_matrix_type;

  // Types of the integrals we use:
  const gint::OrbitalType otype = gint::REAL_ATOMIC;
  const std::string basis_type("Sturmians");

  // The lookup class type to get the actual integrals
  typedef gint::IntegralLookup<stored_matrix_type, otype> int_lookup_type;

  // The type of the integral terms:
  typedef typename int_lookup_type::integral_matrix_type int_term_type;

  //
  // Integral terms
  //
  // Generate integral lookup object
  int_lookup_type integrals(basis_type);

  // Get the integral as actual objects.
  int_term_type S_bb = integrals("overlap");
  int_term_type T_bb = integrals("kinetic");
  int_term_type V0_bb = integrals("nuclear_attraction");
  int_term_type J_bb = integrals("coulomb");
  int_term_type K_bb = integrals("exchange");

  // Combine 1e terms:
  std::vector<int_term_type> terms_1e{std::move(T_bb), std::move(V0_bb)};

  //
  // Debug output
  //
  std::ofstream mathematicafile("/tmp/debug_molsturm_rhf_sturmian.m");
  auto debugout = linalgwrap::io::make_formatted_stream_writer<
        linalgwrap::io::Mathematica, scalar_type>(mathematicafile, 1e-5);

  //
  // Problem setup
  //
  // TODO hack a guess:
  stored_matrix_type guess_bf =
        hack_guess(static_cast<stored_matrix_type>(S_bb));
  debugout.write("sbb", S_bb);
  debugout.write("guess", guess_bf);

  // The term container for the fock operator matrix
  IntegralTermContainer<stored_matrix_type> integral_container(
        std::move(terms_1e), std::move(J_bb), std::move(K_bb));

  RestrictedClosedIntegralOperator<stored_matrix_type> m_fock(
        integral_container, guess_bf, n_alpha, n_beta);

  IopPlainScfSolver<decltype(m_fock)> solver(debugout);
  solver.solve(m_fock, static_cast<stored_matrix_type>(S_bb));
}

struct args_type {
  double Z = 4.0;
  double k_exp = 1.0;
  size_t n_alpha = 2;
  size_t n_beta = n_alpha;
};

/** Quick and dirty function to parse a string to a different type.
 *  Return false if not possible */
template <typename T>
bool str_to_type(const std::string& in, T& out) {
  std::stringstream ss(in);
  if (!(ss >> out)) {
    return false;
  }
  return true;
}

/** \brief Quick and dirty function to parse the commandline arguments
 *
 * \returns true if all is fine, else false
 */
bool parse_args(int argc, char** argv, args_type& parsed) {
  bool had_Z = false;
  bool had_alpha = false;
  bool had_beta = false;
  bool had_k_exp = false;

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
    if (flag == std::string("--Z")) {
      had_Z = true;
      if (!str_to_type<double>(argument, parsed.Z)) {
        std::cerr << "Invalid double provided to --Z: " << argument
                  << std::endl;
        return false;
      }
    } else if (flag == std::string("--alpha")) {
      had_alpha = true;
      if (!str_to_type<size_t>(argument, parsed.n_alpha)) {
        std::cerr << "Invalid int provided to --alpha: " << argument
                  << std::endl;
        return false;
      }
    } else if (flag == std::string("--beta")) {
      had_beta = true;
      if (!str_to_type<size_t>(argument, parsed.n_beta)) {
        std::cerr << "Invalid int provided to --beta: " << argument
                  << std::endl;
        return false;
      }
    } else if (flag == std::string("--kexp")) {
      had_k_exp = true;
      if (!str_to_type<double>(argument, parsed.k_exp)) {
        std::cerr << "Invalid double provided to --kexp: " << argument
                  << std::endl;
        return false;
      }
    } else {
      std::cerr << "Unknown flag: " << flag << std::endl;
      std::cerr << "Valid are: --Z, --alpha, --beta, --kexp" << std::endl;
      return false;
    }
  }

  // Error handling
  if (!had_Z) {
    std::cerr << "Need flag --Z <double> to supply nuclear charge."
              << std::endl;
  }
  if (!had_alpha) {
    std::cerr << "Need flag --alpha <int> to supply number of alpha electrons."
              << std::endl;
  }
  if (!had_beta) {
    std::cerr << "Need flag --beta <int> to supply number of beta electrons."
              << std::endl;
  }
  if (!had_k_exp) {
    std::cerr << "Need flag --kexp <double> to supply k exponent." << std::endl;
  }

  if (had_Z && had_alpha && had_beta && had_k_exp) {
    return true;
  }
  return false;
}

int main(int argc, char** argv) {
  std::cout << "molsturm version:   " << molsturm::version::version_string()
            << std::endl
            << "gscf version:       " << gscf::version::version_string()
            << std::endl
            << "gint version:       " << gint::version::version_string()
            << std::endl
            << "linalgwrap version: " << linalgwrap::version::version_string()
            << std::endl;

  args_type args;
  if (!parse_args(argc, argv, args)) {
    return 1;
  }

  std::cout << "The following configuration was read:" << std::endl
            << "Z:        " << args.Z << std::endl
            << "k_exp:    " << args.k_exp << std::endl
            << "n_alpha:  " << args.n_alpha << std::endl
            << "n_beta:   " << args.n_beta << std::endl
            << std::endl;

  run_rhf_sturmian_debug(args.Z, args.k_exp, args.n_alpha, args.n_beta);
  return 0;
}

}  // namespace hard_coded_hf

int main(int argc, char** argv) { return hard_coded_hf::main(argc, argv); }
