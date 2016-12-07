#include <gint/IntegralLookup.hh>
#include <gint/version.hh>
#include <gscf/version.hh>
#include <iostream>
#include <linalgwrap/io.hh>
#include <linalgwrap/version.hh>
#include <molsturm/GuessLibrary.hh>
#include <molsturm/IopScf.hh>
#include <molsturm/IopScfKeys.hh>
#include <molsturm/RestrictedClosedIntegralOperator.hh>
#include <molsturm/ScfDebugWrapper.hh>
#include <molsturm/version.hh>

DefException1(
      ExcTooSmallBasis, size_t,
      << "A basis of size " << arg1
      << " is too small to incorporate this many electrons. Choose a larger basis.");

namespace hf {
using namespace molsturm;
using namespace linalgwrap;
using namespace krims;

/** Print the resulting eigenstates */
template <typename ProblemMatrix, typename OverlapMatrix>
void print_res(const IopScfState<ProblemMatrix, OverlapMatrix>& res) {
  io::OstreamState oldstate(std::cout);

  const std::string indt = "    ";
  std::cout << "Obtained Hartree-Fock orbitals: " << std::endl;
  std::cout << indt << "a b" << std::endl;

  size_t n_alpha = res.problem_matrix().n_alpha();
  size_t n_beta = res.problem_matrix().n_beta();
  const auto& orben = res.orbital_energies();
  for (const auto& val : orben) {
    std::cout << indt;
    if (n_alpha > 0) {
      --n_alpha;
      std::cout << "* ";
    } else {
      std::cout << "  ";
    }
    if (n_beta > 0) {
      --n_beta;
      std::cout << "* ";
    } else {
      std::cout << "  ";
    }
    std::cout << std::right << std::setw(17) << std::setprecision(12) << val << std::endl;
  }
}

/** Run an (atomic) SCF based on the integral data Sturmian14.
 *
 * \param Z       Number of nuclei
 * \param k_exp   Sturmian exponent
 * \param n_alpha Number of alpha electrons
 * \param n_beta  Number of beta electrons
 */
void run_rhf_sturmian(double k_exp, size_t n_max, size_t l_max, double Z, size_t n_alpha,
                      size_t n_beta, bool debug, double error) {
  //
  // Types and settings
  //
  // Types of scalar and matrix
  typedef double scalar_type;
  typedef SmallMatrix<scalar_type> stored_matrix_type;
  typedef stored_matrix_type::vector_type vector_type;

  // Types of the integrals we use:
  const gint::OrbitalType otype = gint::COMPLEX_ATOMIC;

  // The lookup class type to get the actual integrals
  typedef gint::IntegralLookup<otype> int_lookup_type;

  // The type of the integral terms:
  typedef typename int_lookup_type::integral_type integral_type;

  //
  // Integral terms
  //
  // Generate integral lookup object
  krims::ParameterMap intparams{{"basis_type", "atomic/static14"},
                                // {"basis_type", "atomic/cs_dummy"},
                                {"k_exponent", k_exp},
                                {"Z_charge", Z},
                                {"n_max", static_cast<int>(n_max)},
                                {"l_max", static_cast<int>(l_max)}};
  int_lookup_type integrals{intparams};

  // Get the integral as actual objects.
  integral_type S_bb = integrals("overlap");
  integral_type T_bb = integrals("kinetic");
  integral_type V0_bb = integrals("nuclear_attraction");
  integral_type J_bb = integrals("coulomb");
  integral_type K_bb = integrals("exchange");

  std::cout << "Basis size:  " << S_bb.n_rows() << std::endl << std::endl;

  // Combine 1e terms:
  std::vector<integral_type> terms_1e{std::move(T_bb), std::move(V0_bb)};

  //
  // Problem setup
  //
  // TODO Make this configurable
  //   It is important to note, that the value should not be too small
  //   (not enough virtuals) and not too large
  //   (numerical problems with eigensolvers)
  size_t n_eigenpairs = std::max(n_alpha, n_beta) * 2 + 2;
  n_eigenpairs = std::min(S_bb.n_rows() / 2 - 2, n_eigenpairs);
  assert_throw(std::max(n_alpha, n_beta) <= n_eigenpairs,
               ExcTooSmallBasis(S_bb.n_rows()));

  auto guess_bf_ptr = std::make_shared<linalgwrap::MultiVector<vector_type>>(
        loewdin_guess(S_bb, n_eigenpairs));

  // The term container for the fock operator matrix
  IntegralTermContainer<stored_matrix_type> integral_container(
        std::move(terms_1e), std::move(J_bb), std::move(K_bb));

  RestrictedClosedIntegralOperator<stored_matrix_type> fock_bb(
        integral_container, guess_bf_ptr, n_alpha, n_beta);

  // TODO make threshold configurable
  krims::ParameterMap params{
        {IopScfKeys::max_error_norm, error},
        {IopScfKeys::max_iter, 20ul},
        {IopScfKeys::n_eigenpairs, n_eigenpairs},
        {IopScfKeys::verbosity, ScfMsgType::FinalSummary | ScfMsgType::IterationProcess},
        {IopScfKeys::n_prev_steps, size_t(4)},
  };

  if (debug) {
    std::ofstream mathematicafile("/tmp/debug_molsturm_rhf_sturmian.m");
    auto debugout =
          linalgwrap::io::make_formatted_stream_writer<linalgwrap::io::Mathematica,
                                                       scalar_type>(mathematicafile,
                                                                    1e-12);
    debugout.write("guess", *guess_bf_ptr);
    debugout.write("sbb", S_bb);
    IopScf<decltype(fock_bb), decltype(S_bb)> solver(params);
    ScfDebugWrapper<decltype(solver)> solwrap(debugout, solver);
    auto res = solwrap.solve(fock_bb, S_bb);
    print_res(res);
  } else {
    auto res = run_scf(fock_bb, S_bb, params);
    print_res(res);
  }
}

struct args_type {
  double Z = 4.0;
  double k_exp = 1.0;
  size_t n_alpha = 2;
  size_t n_beta = n_alpha;
  size_t n_max = 3;
  size_t l_max = n_max - 1;
  double error = 5e-7;
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
  bool had_error = false;

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
        std::cerr << "Invalid double provided to --Z: " << argument << std::endl;
        return false;
      }
    } else if (flag == std::string("--alpha")) {
      had_alpha = true;
      if (!str_to_type<size_t>(argument, parsed.n_alpha)) {
        std::cerr << "Invalid int provided to --alpha: " << argument << std::endl;
        return false;
      }
    } else if (flag == std::string("--beta")) {
      had_beta = true;
      if (!str_to_type<size_t>(argument, parsed.n_beta)) {
        std::cerr << "Invalid int provided to --beta: " << argument << std::endl;
        return false;
      }
    } else if (flag == std::string("--kexp")) {
      had_k_exp = true;
      if (!str_to_type<double>(argument, parsed.k_exp)) {
        std::cerr << "Invalid double provided to --kexp: " << argument << std::endl;
        return false;
      }
    } else if (flag == std::string("--nmax")) {
      had_n_max = true;
      if (!str_to_type<size_t>(argument, parsed.n_max)) {
        std::cerr << "Invalid int provided to --nmax: " << argument << std::endl;
        return false;
      }
    } else if (flag == std::string("--lmax")) {
      had_l_max = true;
      if (!str_to_type<size_t>(argument, parsed.l_max)) {
        std::cerr << "Invalid int provided to --lmax: " << argument << std::endl;
        return false;
      }
    } else if (flag == std::string("--error")) {
      had_error = true;
      if (!str_to_type<double>(argument, parsed.error)) {
        std::cerr << "Invalid double provided to --error: " << argument << std::endl;
        return false;
      }
    } else {
      std::cerr << "Unknown flag: " << flag << std::endl;
      std::cerr << "Valid are: --Z, --alpha, --beta, --kexp, --lmax, --nmax, --error"
                << std::endl;
      return false;
    }
  }

  // Error handling
  if (!had_Z) {
    std::cerr << "Need flag --Z <double> to supply nuclear charge." << std::endl;
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

  if (!had_error) {
    parsed.error = 5e-7;
  }

  if (had_Z && had_alpha && had_beta && had_k_exp && had_n_max) {
    return true;
  }
  return false;
}

int main(int argc, char** argv) {
  std::cout << "molsturm version:   " << molsturm::version::version_string() << std::endl
            << "gscf version:       " << gscf::version::version_string() << std::endl
            << "gint version:       " << gint::version::version_string() << std::endl
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
            << "error:   " << args.error << std::endl
            << std::endl;

  const bool debug_mode = false;
  run_rhf_sturmian(args.k_exp, args.n_max, args.l_max, args.Z, args.n_alpha, args.n_beta,
                   debug_mode, args.error);
  return 0;
}

}  // namespace hf

int main(int argc, char** argv) { return hf::main(argc, argv); }
