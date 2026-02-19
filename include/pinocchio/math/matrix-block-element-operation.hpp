//
// Copyright (c) 2026 INRIA
//

#ifndef __pinocchio_math_matrix_block_element_operation_hpp__
#define __pinocchio_math_matrix_block_element_operation_hpp__

namespace pinocchio
{

  template<typename Derived>
  struct MatrixBlockElementPlain;

  template<typename Derived>
  struct MatrixBlockElementOperation : MatrixBlockElementBase<Derived>
  {

    typedef MatrixBlockElementBase<Derived> Base;
    using Base::derived;

    template<typename OtherDerived>
    void evalTo(MatrixBlockElementPlain<OtherDerived> & res) const
    {
      derived().evalTo(res.derived());
    }
  };

} // namespace pinocchio

#endif // #ifndef __pinocchio_math_matrix_block_element_operation_hpp__
