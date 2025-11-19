//
// Copyright (c) 2023-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraint_data_base_hpp__
#define __pinocchio_algorithm_constraint_data_base_hpp__

#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/common/data-entity.hpp"

namespace pinocchio
{

  template<typename Derived>
  struct ConstraintDataBase
  : NumericalBase<Derived>
  , DataEntity<Derived>
  {
    typedef typename traits<Derived>::Scalar Scalar;
    typedef typename traits<Derived>::ConstraintModel ConstraintModel;

    typedef
      typename traits<ConstraintModel>::ActiveComplianceVectorTypeRef ActiveComplianceVectorTypeRef;
    typedef typename traits<ConstraintModel>::ActiveComplianceVectorTypeConstRef
      ActiveComplianceVectorTypeConstRef;

    typedef DataEntity<Derived> Base;

    /// \brief Cast to Derived
    Derived & derived()
    {
      return static_cast<Derived &>(*this);
    }

    /// \brief Const cast to Derived
    const Derived & derived() const
    {
      return static_cast<const Derived &>(*this);
    }

    /// \brief Returns the name of the underlying constraint if this is a variant.
    std::string shortname() const
    {
      return derived().shortnameImpl();
    }

    /// \brief Returns the name of the class.
    static std::string classname()
    {
      return Derived::classnameImpl();
    }

    /// \brief Print the name of the class
    void disp(std::ostream & os) const
    {
      using namespace std;
      os << shortname() << endl;
    }

    /// \brief Print the name of the class
    friend std::ostream &
    operator<<(std::ostream & os, const ConstraintDataBase<Derived> & constraint)
    {
      constraint.disp(os);
      return os;
    }

    /// \brief Comparison operator
    template<typename OtherDerived>
    bool operator==(const ConstraintDataBase<OtherDerived> &) const
    {
      return true;
    }

    /// \brief Comparison operator
    template<typename OtherDerived>
    bool operator!=(const ConstraintDataBase<OtherDerived> & other) const
    {
      return !(*this == other);
    }

  protected:
    /// \brief Default constructor
    ConstraintDataBase()
    {
    }
  };
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraint_data_base_hpp__
