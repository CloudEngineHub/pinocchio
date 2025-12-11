//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_constraint_model_common_parameters_hpp__
#define __pinocchio_algorithm_constraints_constraint_model_common_parameters_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"

namespace pinocchio
{

  template<typename _Scalar>
  struct BaumgarteCorrectorParametersTpl;

  template<typename Derived>
  struct ConstraintModelCommonParameters
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    template<typename OtherDerived>
    friend struct ConstraintModelCommonParameters;

    typedef ConstraintModelCommonParameters<Derived> Self;
    typedef typename traits<Derived>::Scalar Scalar;
    typedef typename traits<Derived>::ComplianceVectorType ComplianceVectorType;
    typedef typename traits<Derived>::ComplianceVectorTypeRef ComplianceVectorTypeRef;
    typedef typename traits<Derived>::ComplianceVectorTypeConstRef ComplianceVectorTypeConstRef;

    typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;

    // -------------------------------
    // METHODS SPECIFIC TO CLASS
    // -------------------------------

    // Constructors ------------------

  protected:
    /// \brief Default constructor - protected so that the class cannot be instanciated on its own.
    ConstraintModelCommonParameters()
    {
    }

    // Operators ---------------------

  public:
    /// \brief Cast to NewScalar
    template<typename NewScalar, typename OtherDerived>
    void cast(ConstraintModelCommonParameters<OtherDerived> & other) const
    {
      other.m_compliance = m_compliance.template cast<NewScalar>();
      other.m_baumgarte_parameters = m_baumgarte_parameters.template cast<NewScalar>();
    }

    /// \brief Comparison operator
    bool operator==(const Self & other) const
    {
      return m_compliance == other.m_compliance
             && m_baumgarte_parameters == other.m_baumgarte_parameters;
    }

    /// \brief Comparison operator
    bool operator!=(const Self & other) const
    {
      return !(*this == other);
    }

    /// \brief Returns the compliance internally stored in the constraint model.
    ComplianceVectorTypeConstRef compliance_impl() const
    {
      return m_compliance;
    }

    /// \brief Returns the compliance internally stored in the constraint model.
    ComplianceVectorTypeRef compliance_impl()
    {
      return m_compliance;
    }

    /// \brief Returns the Baumgarte parameters internally stored in the constraint model
    const BaumgarteCorrectorParameters & baumgarte_corrector_parameters_impl() const
    {
      return m_baumgarte_parameters;
    }

    /// \brief Returns the Baumgarte parameters internally stored in the constraint model
    BaumgarteCorrectorParameters & baumgarte_corrector_parameters_impl()
    {
      return m_baumgarte_parameters;
    }

  protected:
    // ------------------------------
    // MEMBERS
    // ------------------------------

    ComplianceVectorType m_compliance;
    BaumgarteCorrectorParameters m_baumgarte_parameters;
  };

} // namespace pinocchio

#endif // __pinocchio_algorithm_constraints_constraint_model_common_parameters_hpp__
