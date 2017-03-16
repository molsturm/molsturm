#pragma once
#include <gscf/FocklikeMatrix_i.hh>
#include <linalgwrap/LazyMatrix_i.hh>

namespace molsturm {

template <typename StoredMatrix>
class IntegralOperatorBase
      : public linalgwrap::LazyMatrix_i<StoredMatrix>,
        public gscf::FocklikeMatrix_i<typename StoredMatrix::scalar_type> {
 public:
  typedef linalgwrap::LazyMatrix_i<StoredMatrix> base_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::stored_matrix_type stored_matrix_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::lazy_matrix_expression_ptr_type
        lazy_matrix_expression_ptr_type;

  // TODO implement general integral operator stuff in here
};

//@{
/** \brief struct representing a type (std::true_type, std::false_type) which
 *  indicates whether T is an IntegralOperator
 *
 * The definition is done using SFINAE, such that even for types not having a
 * typedef scalar_type this expression is valid.
 *  */
template <typename T, typename = void>
struct IsIntegralOperator : public std::false_type {};

template <typename T>
struct IsIntegralOperator<T, krims::VoidType<typename T::stored_matrix_type>>
      : public std::is_base_of<IntegralOperatorBase<typename T::stored_matrix_type>, T> {
};
//@}

}  // namespace molsturm
