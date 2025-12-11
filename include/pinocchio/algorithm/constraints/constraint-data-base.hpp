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

    typedef DataEntity<Derived> Base;

    // -------------------------------
    // METHODS SPECIFIC TO CLASS
    // -------------------------------

    // CRTP related ------------------

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

    /// \brief Cast to base.
    ConstraintDataBase & base()
    {
      return *this;
    }

    /// \brief Const cast to base.
    const ConstraintDataBase & base() const
    {
      return *this;
    }

    // Constructors ------------------

  protected:
    /// \brief Default constructor
    ConstraintDataBase()
    {
    }

    // Operators ---------------------

  public:
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

    /// \brief Returns the name of the class.
    static std::string classname()
    {
      return Derived::classnameImpl();
    }

    /// \brief Returns the name of the underlying constraint if this is a variant.
    std::string shortname() const
    {
      return derived().shortnameImpl();
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
  };
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraint_data_base_hpp__
