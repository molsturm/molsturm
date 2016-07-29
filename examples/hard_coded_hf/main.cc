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
void run_rhf_sturmian_debug(double k_exp, size_t n_max, size_t l_max, double Z,
                            size_t n_alpha, size_t n_beta) {
  //
  // Types and settings
  //
  // Types of scalar and matrix
  typedef double scalar_type;
  typedef SmallMatrix<scalar_type> stored_matrix_type;

  // Parameters for the basis
  ParameterMap intparams;
  intparams.update_copy("k_exponent", k_exp);
  intparams.update_copy("Z_charge", Z);
  intparams.update_copy("n_max", static_cast<int>(n_max));
  intparams.update_copy("l_max", static_cast<int>(l_max));

  // Types of the integrals we use:
  const gint::OrbitalType otype = gint::COMPLEX_ATOMIC;
  const std::string basis_type("atomic/cs_dummy");

  // The lookup class type to get the actual integrals
  typedef gint::IntegralLookup<stored_matrix_type, otype> int_lookup_type;

  // The type of the integral terms:
  typedef typename int_lookup_type::integral_matrix_type int_term_type;

  //
  // Integral terms
  //
  // Generate integral lookup object
  int_lookup_type integrals(basis_type, intparams);

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
        linalgwrap::io::Mathematica, scalar_type>(mathematicafile, 1e-12);

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

  RestrictedClosedIntegralOperator<stored_matrix_type> fock_bb(
        integral_container, guess_bf, n_alpha, n_beta);

  debugout.write("guessfock", fock_bb);

  IopPlainScfSolver<decltype(fock_bb)> solver(debugout);
  solver.scf_control().max_iter = 10;
  solver.solve(fock_bb, static_cast<stored_matrix_type>(S_bb));
}

struct args_type {
  double Z = 4.0;
  double k_exp = 1.0;
  size_t n_alpha = 2;
  size_t n_beta = n_alpha;
  size_t n_max = 3;
  size_t l_max = n_max - 1;
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
  bool had_n_max = false;
  bool had_l_max = false;

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
    } else if (flag == std::string("--nmax")) {
      had_n_max = true;
      if (!str_to_type<size_t>(argument, parsed.n_max)) {
        std::cerr << "Invalid int provided to --nmax: " << argument
                  << std::endl;
        return false;
      }
    } else if (flag == std::string("--lmax")) {
      had_l_max = true;
      if (!str_to_type<size_t>(argument, parsed.n_max)) {
        std::cerr << "Invalid int provided to --lmax: " << argument
                  << std::endl;
        return false;
      }
    } else {
      std::cerr << "Unknown flag: " << flag << std::endl;
      std::cerr << "Valid are: --Z, --alpha, --beta, --kexp, --lmax, --nmax"
                << std::endl;
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
  if (!had_n_max) {
    std::cerr << "Need flag --nmax <int> to supply maximal principle quantum "
                 "number."
              << std::endl;
  }

  if (!had_l_max) {
    parsed.l_max = parsed.n_max - 1;
  }

  if (had_Z && had_alpha && had_beta && had_k_exp && had_n_max) {
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
            << "k_exp:    " << args.k_exp << std::endl
            << "n_max:    " << args.n_max << std::endl
            << "l_max:    " << args.l_max << std::endl
            << "Z:        " << args.Z << std::endl
            << "n_alpha:  " << args.n_alpha << std::endl
            << "n_beta:   " << args.n_beta << std::endl
            << std::endl;

  run_rhf_sturmian_debug(args.k_exp, args.n_max, args.l_max, args.Z,
                         args.n_alpha, args.n_beta);
  return 0;
}

}  // namespace hard_coded_hf

int main(int argc, char** argv) { return hard_coded_hf::main(argc, argv); }
