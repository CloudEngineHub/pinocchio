//
// Copyright (c) 2023-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_constraint_model_visitor_hpp__
#define __pinocchio_algorithm_constraints_constraint_model_visitor_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-data-base.hpp"
#include "pinocchio/multibody/visitor/fusion.hpp"

namespace pinocchio
{

  namespace visitors
  {

    namespace bf = boost::fusion;
    using fusion::NoArg;

    namespace internal
    {
      template<typename T>
      struct NoRun
      {
        static T run()
        {
          assert(false && "Should never happened.");
          // Hacky way to not have to return something real. The system should throw before.
          const typename std::remove_reference<T>::type * null_ptr = NULL;
          return *null_ptr;
        }
      };

      // Specialization for reference types
      template<typename T>
      struct NoRun<T &>
      {
        static T & run()
        {
          assert(false && "Should never happen.");
          // Hacky way to not have to return something real. The system should throw before.
          T * null_ptr = nullptr;
          return *null_ptr;
        }
      };

      template<>
      struct NoRun<void>
      {
        static void run()
        {
          return;
        }
      };
    } // namespace internal

    ///
    /// \brief Base structure for \b Unary visitation of a ConstraintModel.
    ///        This structure provides runners to call the right visitor according to the number of
    ///        arguments.
    ///

    template<typename ConstraintModelVisitorDerived, typename ReturnType = void>
    struct ConstraintUnaryVisitorBase
    {
      template<
        typename Scalar,
        int Options,
        template<typename, int> class ConstraintCollectionTpl,
        typename ArgsTmp>
      static ReturnType run(
        ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
        ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata,
        ArgsTmp args)
      {
        typedef ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintModel;
        typedef ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintData;

        ModelAndDataVisitor<ConstraintModel, ConstraintData, ArgsTmp> visitor(cdata, args);

        return boost::apply_visitor(visitor, cmodel);
      }

      template<
        typename Scalar,
        int Options,
        template<typename, int> class ConstraintCollectionTpl,
        typename ArgsTmp>
      static ReturnType run(
        const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
        ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata,
        ArgsTmp args)
      {
        typedef ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintModel;
        typedef ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintData;

        ModelAndDataVisitor<ConstraintModel, ConstraintData, ArgsTmp> visitor(cdata, args);

        return boost::apply_visitor(visitor, cmodel);
      }

      template<
        typename Scalar,
        int Options,
        template<typename, int> class ConstraintCollectionTpl,
        typename ArgsTmp>
      static ReturnType run(
        const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
        const ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata,
        ArgsTmp args)
      {
        typedef ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintModel;
        typedef ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintData;

        ModelAndDataVisitor<ConstraintModel, const ConstraintData, ArgsTmp> visitor(cdata, args);

        return boost::apply_visitor(visitor, cmodel);
      }

      template<typename Scalar, int Options, template<typename, int> class ConstraintCollectionTpl>
      static ReturnType run(
        const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
        const ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata)
      {
        typedef ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintModel;
        typedef ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintData;

        ModelAndDataVisitor<ConstraintModel, const ConstraintData, NoArg> visitor(cdata);

        return boost::apply_visitor(visitor, cmodel);
      }

      template<
        typename Scalar,
        int Options,
        template<typename, int> class ConstraintCollectionTpl,
        typename ArgsTmp>
      static ReturnType
      run(const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel, ArgsTmp args)
      {
        ModelVisitor<Scalar, Options, ConstraintCollectionTpl, ArgsTmp> visitor(args);
        return boost::apply_visitor(visitor, cmodel);
      }

      template<typename Scalar, int Options, template<typename, int> class ConstraintCollectionTpl>
      static ReturnType
      run(const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
      {
        ModelVisitor<Scalar, Options, ConstraintCollectionTpl, NoArg> visitor;
        return boost::apply_visitor(visitor, cmodel);
      }

      template<
        typename Scalar,
        int Options,
        template<typename, int> class ConstraintCollectionTpl,
        typename ArgsTmp>
      static ReturnType
      run(ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel, ArgsTmp args)
      {
        ModelVisitor<Scalar, Options, ConstraintCollectionTpl, ArgsTmp> visitor(args);
        return boost::apply_visitor(visitor, cmodel);
      }

      template<typename Scalar, int Options, template<typename, int> class ConstraintCollectionTpl>
      static ReturnType run(ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
      {
        ModelVisitor<Scalar, Options, ConstraintCollectionTpl, NoArg> visitor;
        return boost::apply_visitor(visitor, cmodel);
      }

    private:
      template<
        typename Scalar,
        int Options,
        template<typename, int> class ConstraintCollectionTpl,
        typename ArgsTmp>
      struct ModelVisitor : public boost::static_visitor<ReturnType>
      {
        typedef ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintModel;
        typedef ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintData;

        ModelVisitor(ArgsTmp args)
        : args(args)
        {
        }

        template<typename ConstraintModelDerived>
        ReturnType operator()(const ConstraintModelBase<ConstraintModelDerived> & cmodel) const
        {
          return bf::invoke(
            &ConstraintModelVisitorDerived::template algo<ConstraintModelDerived>,
            bf::append(boost::ref(cmodel.derived()), args));
        }

        template<typename ConstraintModelDerived>
        ReturnType operator()(ConstraintModelBase<ConstraintModelDerived> & cmodel) const
        {
          return bf::invoke(
            &ConstraintModelVisitorDerived::template algo<ConstraintModelDerived>,
            bf::append(boost::ref(cmodel.derived()), args));
        }

        template<typename ConstraintDataDerived>
        ReturnType operator()(const ConstraintDataBase<ConstraintDataDerived> & cdata) const
        {
          return bf::invoke(
            &ConstraintModelVisitorDerived::template algo<ConstraintDataDerived>,
            bf::append(boost::ref(cdata.derived()), args));
        }

        ReturnType operator()(const BlankConstraintModel &) const
        {
          PINOCCHIO_THROW_PRETTY(
            std::invalid_argument, "The constraint model is of type BlankConstraintModel.");
          return internal::NoRun<ReturnType>::run();
        }

        ArgsTmp args;
      };

      template<typename Scalar, int Options, template<typename, int> class ConstraintCollectionTpl>
      struct ModelVisitor<Scalar, Options, ConstraintCollectionTpl, NoArg>
      : public boost::static_visitor<ReturnType>
      {
        typedef ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintModel;
        typedef ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintData;

        ModelVisitor()
        {
        }

        template<typename ConstraintModelDerived>
        ReturnType operator()(const ConstraintModelBase<ConstraintModelDerived> & cmodel) const
        {
          return ConstraintModelVisitorDerived::template algo<ConstraintModelDerived>(
            cmodel.derived());
        }

        template<typename ConstraintModelDerived>
        ReturnType operator()(ConstraintModelBase<ConstraintModelDerived> & cmodel) const
        {
          return ConstraintModelVisitorDerived::template algo<ConstraintModelDerived>(
            cmodel.derived());
        }

        template<typename ConstraintDataDerived>
        ReturnType operator()(const ConstraintDataBase<ConstraintDataDerived> & cdata) const
        {
          return ConstraintModelVisitorDerived::template algo<ConstraintDataDerived>(
            cdata.derived());
        }

        ReturnType operator()(const BlankConstraintModel &) const
        {
          PINOCCHIO_THROW_PRETTY(
            std::invalid_argument, "The constraint model is of type BlankConstraintModel.");
          return internal::NoRun<ReturnType>::run();
        }

      }; // struct ModelVisitor

      template<typename ConstraintModel, typename ConstraintData, typename ArgsTmp>
      struct ModelAndDataVisitor : public boost::static_visitor<ReturnType>
      {

        ModelAndDataVisitor(ConstraintData & cdata, ArgsTmp args)
        : cdata(cdata)
        , args(args)
        {
        }

        template<typename ConstraintModelDerived>
        ReturnType operator()(ConstraintModelBase<ConstraintModelDerived> & cmodel) const
        {
          typedef typename ConstraintModelBase<ConstraintModelDerived>::ConstraintData
            ConstraintDataDerived;
          using ConstraintDataGet = typename std::conditional<
            std::is_const<ConstraintData>::value, const ConstraintDataDerived,
            ConstraintDataDerived>::type;

          return bf::invoke(
            &ConstraintModelVisitorDerived::template algo<ConstraintModelDerived>,
            bf::append(
              boost::ref(cmodel.derived()), boost::ref(boost::get<ConstraintDataGet>(cdata)),
              args));
        }

        template<typename ConstraintModelDerived>
        ReturnType operator()(const ConstraintModelBase<ConstraintModelDerived> & cmodel) const
        {
          typedef typename ConstraintModelBase<ConstraintModelDerived>::ConstraintData
            ConstraintDataDerived;
          using ConstraintDataGet = typename std::conditional<
            std::is_const<ConstraintData>::value, const ConstraintDataDerived,
            ConstraintDataDerived>::type;

          return bf::invoke(
            &ConstraintModelVisitorDerived::template algo<ConstraintModelDerived>,
            bf::append(
              boost::ref(cmodel.derived()), boost::ref(boost::get<ConstraintDataGet>(cdata)),
              args));
        }

        ReturnType operator()(const BlankConstraintModel &) const
        {
          PINOCCHIO_THROW_PRETTY(
            std::invalid_argument, "The constraint model is of type BlankConstraintModel.");
          return internal::NoRun<ReturnType>::run();
        }

        ConstraintData & cdata;
        ArgsTmp args;
      }; // struct ModelAndDataVisitor

      template<typename ConstraintModel, typename ConstraintData>
      struct ModelAndDataVisitor<ConstraintModel, ConstraintData, NoArg>
      : public boost::static_visitor<ReturnType>
      {

        ModelAndDataVisitor(ConstraintData & cdata)
        : cdata(cdata)
        {
        }

        template<typename ConstraintModelDerived>
        ReturnType operator()(ConstraintModelBase<ConstraintModelDerived> & cmodel) const
        {
          typedef typename ConstraintModelBase<ConstraintModelDerived>::ConstraintData
            ConstraintDataDerived;
          using ConstraintDataGet = typename std::conditional<
            std::is_const<ConstraintData>::value, const ConstraintDataDerived,
            ConstraintDataDerived>::type;

          return bf::invoke(
            &ConstraintModelVisitorDerived::template algo<ConstraintModelDerived>,
            bf::make_vector(
              boost::ref(cmodel.derived()), boost::ref(boost::get<ConstraintDataGet>(cdata))));
        }

        template<typename ConstraintModelDerived>
        ReturnType operator()(const ConstraintModelBase<ConstraintModelDerived> & cmodel) const
        {
          typedef typename ConstraintModelBase<ConstraintModelDerived>::ConstraintData
            ConstraintDataDerived;
          using ConstraintDataGet = typename std::conditional<
            std::is_const<ConstraintData>::value, const ConstraintDataDerived,
            ConstraintDataDerived>::type;

          return bf::invoke(
            &ConstraintModelVisitorDerived::template algo<ConstraintModelDerived>,
            bf::make_vector(
              boost::ref(cmodel.derived()), boost::ref(boost::get<ConstraintDataGet>(cdata))));
        }

        ReturnType operator()(const BlankConstraintModel &) const
        {
          PINOCCHIO_THROW_PRETTY(
            std::invalid_argument, "The constraint model is of type BlankConstraintModel.");
          return internal::NoRun<ReturnType>::run();
        }

        ConstraintData & cdata;
      }; // struct ModelAndDataVisitor
    };

    // ----------------------------------------------------------------------
    // Implementation of the visitors
    // ----------------------------------------------------------------------

    /**
     * @brief      ConstraintDataComparisonOperatorVisitor fusion visitor
     */
    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl,
      typename ConstraintDataDerived>
    struct ConstraintDataComparisonOperatorVisitor
    : visitors::ConstraintUnaryVisitorBase<
        ConstraintDataComparisonOperatorVisitor<
          Scalar,
          Options,
          ConstraintCollectionTpl,
          ConstraintDataDerived>,
        bool>
    {
      typedef boost::fusion::vector<const ConstraintDataDerived &> ArgsType;

      template<typename ConstraintData>
      static bool algo(
        const ConstraintDataBase<ConstraintData> & cdata_lhs,
        const ConstraintDataDerived & cdata_rhs)
      {
        return cdata_lhs.derived() == cdata_rhs;
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl,
      typename ConstraintDataDerived>
    bool isEqual(
      const ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata_generic,
      const ConstraintDataBase<ConstraintDataDerived> & cdata)
    {
      typedef ConstraintDataComparisonOperatorVisitor<
        Scalar, Options, ConstraintCollectionTpl, ConstraintDataDerived>
        Algo;
      return Algo::run(cdata_generic, typename Algo::ArgsType(boost::ref(cdata.derived())));
    }

    /**
     * @brief      ConstraintModelShortnameVisitor visitor
     */
    struct ConstraintModelShortnameVisitor : boost::static_visitor<std::string>
    {
      template<typename ConstraintModelDerived>
      std::string operator()(const ConstraintModelBase<ConstraintModelDerived> & cmodel) const
      {
        return cmodel.shortname();
      }
      std::string operator()(const BlankConstraintModel &) const
      {
        PINOCCHIO_THROW_PRETTY(
          std::invalid_argument, "The constraint model is of type BlankConstraintModel.");
        return internal::NoRun<std::string>::run();
      }

      template<
        typename Scalar,
        int Options,
        template<typename S, int O> class ConstraintCollectionTpl>
      static std::string
      run(const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
      {
        return boost::apply_visitor(ConstraintModelShortnameVisitor(), cmodel);
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    inline std::string
    shortname(const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
    {
      return ConstraintModelShortnameVisitor::run(cmodel);
    }

    /**
     * @brief      ConstraintDataShortnameVisitor visitor
     */
    struct ConstraintDataShortnameVisitor : boost::static_visitor<std::string>
    {
      template<typename ConstraintDataDerived>
      std::string operator()(const ConstraintDataBase<ConstraintDataDerived> & cdata) const
      {
        return cdata.shortname();
      }
      std::string operator()(const BlankConstraintData &) const
      {
        PINOCCHIO_THROW_PRETTY(
          std::invalid_argument, "The constraint data is of type BlankConstraintData.");
        return internal::NoRun<std::string>::run();
      }

      template<
        typename Scalar,
        int Options,
        template<typename S, int O> class ConstraintCollectionTpl>
      static std::string
      run(const ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata)
      {
        return boost::apply_visitor(ConstraintDataShortnameVisitor(), cdata);
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    inline std::string
    shortname(const ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata)
    {
      return ConstraintDataShortnameVisitor::run(cdata);
    }

    /**
     * @brief      ConstraintModelCreateDataVisitor fusion visitor
     */
    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    struct ConstraintModelCreateDataVisitor
    : visitors::ConstraintUnaryVisitorBase<
        ConstraintModelCreateDataVisitor<Scalar, Options, ConstraintCollectionTpl>,
        typename ConstraintCollectionTpl<Scalar, Options>::ConstraintDataVariant>
    {
      typedef NoArg ArgsType;
      typedef ConstraintCollectionTpl<Scalar, Options> ConstraintCollection;
      typedef typename ConstraintCollection::ConstraintModelVariant ConstraintModelVariant;
      typedef typename ConstraintCollection::ConstraintDataVariant ConstraintDataVariant;

      template<typename ConstraintModel>
      static ConstraintDataVariant
      algo(const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel)
      {
        return cmodel.createData();
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl>
    createData(const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
    {
      typedef ConstraintModelCreateDataVisitor<Scalar, Options, ConstraintCollectionTpl> Algo;
      return Algo::run(cmodel);
    }

    /**
     * @brief      ConstraintModelMaxResidualSizeVisitor visitor
     */
    template<typename Scalar, int Options>
    struct ConstraintModelMaxResidualSizeVisitor
    : visitors::
        ConstraintUnaryVisitorBase<ConstraintModelMaxResidualSizeVisitor<Scalar, Options>, int>
    {
      typedef NoArg ArgsType;

      template<typename ConstraintModel>
      static int algo(const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel)
      {
        return cmodel.maxResidualSize();
      }
    };

    template<typename Scalar, int Options, template<typename, int> class ConstraintCollectionTpl>
    int maxResidualSize(const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
    {
      typedef ConstraintModelMaxResidualSizeVisitor<Scalar, Options> Algo;
      return Algo::run(cmodel);
    }

    /**
     * @brief      ConstraintModelResidualSizeVisitor visitor
     */
    template<typename Scalar, int Options>
    struct ConstraintModelResidualSizeVisitor
    : visitors::ConstraintUnaryVisitorBase<ConstraintModelResidualSizeVisitor<Scalar, Options>, int>
    {
      typedef NoArg ArgsType;

      template<typename ConstraintModel>
      static int algo(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const typename ConstraintModel::ConstraintData & cdata)
      {
        return cmodel.residualSize(cdata);
      }
    };

    template<typename Scalar, int Options, template<typename, int> class ConstraintCollectionTpl>
    int residualSize(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata)
    {
      typedef ConstraintModelResidualSizeVisitor<Scalar, Options> Algo;
      return Algo::run(cmodel, cdata);
    }

    /**
     * @brief      ConstraintModelRetriveComplianceVisitor visitor
     */
    template<typename VectorLike>
    struct ConstraintModelRetriveComplianceVisitor
    : ConstraintUnaryVisitorBase<ConstraintModelRetriveComplianceVisitor<VectorLike>>
    {
      typedef boost::fusion::vector<VectorLike &> ArgsType;

      template<typename ConstraintModel>
      static void algo(
        const ::pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const typename ConstraintModel::ConstraintData & cdata,
        const Eigen::MatrixBase<VectorLike> & res)
      {
        return cmodel.retrieveCompliance(cdata.derived(), res.const_cast_derived());
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl,
      typename VectorLike>
    void retrieveCompliance(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata,
      const Eigen::MatrixBase<VectorLike> & res)
    {
      typedef ConstraintModelRetriveComplianceVisitor<VectorLike> Algo;
      typename Algo::ArgsType args(res.const_cast_derived());
      return Algo::run(cmodel, cdata, args);
    }

    /**
     * @brief      ConstraintModelgetRowSparsityPatternVisitor visitor
     */
    template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
    struct ConstraintModelgetRowSparsityPatternVisitor
    : visitors::ConstraintUnaryVisitorBase<
        ConstraintModelgetRowSparsityPatternVisitor<Scalar, Options, JointCollectionTpl>,
        const Eigen::Matrix<bool, Eigen::Dynamic, 1, Options> &>
    {
      typedef const Eigen::Matrix<bool, Eigen::Dynamic, 1, Options> & ReturnType;
      typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      typedef boost::fusion::vector<const Model &, const Data &, const Eigen::Index> ArgsType;

      template<typename ConstraintModel>
      static ReturnType algo(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const typename ConstraintModel::ConstraintData & cdata,
        const Model & model,
        const Data & data,
        const Eigen::Index row_id)
      {
        return cmodel.getRowSparsityPattern(model, data, cdata.derived(), row_id);
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class JointCollectionTpl,
      template<typename S, int O> class ConstraintCollectionTpl>
    const Eigen::Matrix<bool, Eigen::Dynamic, 1, Options> & getRowSparsityPattern(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata,
      const Eigen::Index row_id)
    {
      typedef ConstraintModelgetRowSparsityPatternVisitor<Scalar, Options, JointCollectionTpl> Algo;
      return Algo::run(cmodel, cdata, typename Algo::ArgsType(model, data, row_id));
    }

    /**
     * @brief      ConstraintModelgetRowIndexesVisitor visitor
     */
    template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
    struct ConstraintModelgetRowIndexesVisitor
    : visitors::ConstraintUnaryVisitorBase<
        ConstraintModelgetRowIndexesVisitor<Scalar, Options, JointCollectionTpl>,
        const std::vector<Eigen::Index> &>
    {
      typedef const std::vector<Eigen::Index> & ReturnType;
      typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      typedef boost::fusion::vector<const Model &, const Data &, const Eigen::Index> ArgsType;

      template<typename ConstraintModel>
      static ReturnType algo(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const typename ConstraintModel::ConstraintData & cdata,
        const Model & model,
        const Data & data,
        const Eigen::Index row_id)
      {
        return cmodel.getRowIndexes(model, data, cdata.derived(), row_id);
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class JointCollectionTpl,
      template<typename S, int O> class ConstraintCollectionTpl>
    const std::vector<Eigen::Index> & getRowIndexes(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata,
      const Eigen::Index row_id)
    {
      typedef ConstraintModelgetRowIndexesVisitor<Scalar, Options, JointCollectionTpl> Algo;
      return Algo::run(cmodel, cdata, typename Algo::ArgsType(model, data, row_id));
    }

    /**
     * @brief      ConstraintModelCalcVisitor fusion visitor
     */
    template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
    struct ConstraintModelCalcVisitor
    : visitors::ConstraintUnaryVisitorBase<
        ConstraintModelCalcVisitor<Scalar, Options, JointCollectionTpl>>
    {
      typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      typedef boost::fusion::vector<const Model &, const Data &> ArgsType;

      template<typename ConstraintModel>
      static void algo(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        typename ConstraintModel::ConstraintData & cdata,
        const Model & model,
        const Data & data)
      {
        cmodel.calc(model, data, cdata.derived());
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class JointCollectionTpl,
      template<typename S, int O> class ConstraintCollectionTpl>
    void calc(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata)
    {
      typedef ConstraintModelCalcVisitor<Scalar, Options, JointCollectionTpl> Algo;
      Algo::run(cmodel, cdata, typename Algo::ArgsType(model, data));
    }

    /**
     * @brief      ConstraintModelJacobianVisitor fusion visitor
     */
    template<
      typename Scalar,
      int Options,
      template<typename, int> class JointCollectionTpl,
      typename JacobianMatrix>
    struct ConstraintModelJacobianVisitor
    : visitors::ConstraintUnaryVisitorBase<
        ConstraintModelJacobianVisitor<Scalar, Options, JointCollectionTpl, JacobianMatrix>>
    {
      typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      typedef boost::fusion::vector<const Model &, const Data &, JacobianMatrix &> ArgsType;

      template<typename ConstraintModel>
      static void algo(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const typename ConstraintModel::ConstraintData & cdata,
        const Model & model,
        const Data & data,
        const Eigen::MatrixBase<JacobianMatrix> & jacobian_matrix)
      {
        cmodel.jacobian(model, data, cdata.derived(), jacobian_matrix.const_cast_derived());
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class JointCollectionTpl,
      template<typename S, int O> class ConstraintCollectionTpl,
      typename JacobianMatrix>
    void jacobian(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata,
      const Eigen::MatrixBase<JacobianMatrix> & jacobian_matrix)
    {
      typedef ConstraintModelJacobianVisitor<Scalar, Options, JointCollectionTpl, JacobianMatrix>
        Algo;
      Algo::run(
        cmodel, cdata, typename Algo::ArgsType(model, data, jacobian_matrix.const_cast_derived()));
    }

    /**
     * @brief      ConstraintModelJacobianMatrixProductVisitor visitor
     */
    template<
      typename Scalar,
      int Options,
      template<typename, int> class JointCollectionTpl,
      typename InputMatrix,
      typename OutputMatrix,
      AssignmentOperatorType op>
    struct ConstraintModelJacobianMatrixProductVisitor
    : visitors::ConstraintUnaryVisitorBase<ConstraintModelJacobianMatrixProductVisitor<
        Scalar,
        Options,
        JointCollectionTpl,
        InputMatrix,
        OutputMatrix,
        op>>
    {
      typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      typedef boost::fusion::vector<
        const Model &,
        const Data &,
        const InputMatrix &,
        OutputMatrix &,
        AssignmentOperatorTag<op>>
        ArgsType;

      template<typename ConstraintModel>
      static void algo(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const typename ConstraintModel::ConstraintData & cdata,
        const Model & model,
        const Data & data,
        const Eigen::MatrixBase<InputMatrix> & input_matrix,
        const Eigen::MatrixBase<OutputMatrix> & result_matrix,
        AssignmentOperatorTag<op> aot)
      {
        cmodel.jacobianMatrixProduct(
          model, data, cdata.derived(), input_matrix.derived(), result_matrix.const_cast_derived(),
          aot);
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class JointCollectionTpl,
      template<typename S, int O> class ConstraintCollectionTpl,
      typename InputMatrix,
      typename OutputMatrix,
      AssignmentOperatorType op = SETTO>
    void jacobianMatrixProduct(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata,
      const Eigen::MatrixBase<InputMatrix> & input_matrix,
      const Eigen::MatrixBase<OutputMatrix> & result_matrix,
      AssignmentOperatorTag<op> aot = SetTo())
    {
      typedef ConstraintModelJacobianMatrixProductVisitor<
        Scalar, Options, JointCollectionTpl, InputMatrix, OutputMatrix, op>
        Algo;

      typename Algo::ArgsType args(
        model, data, input_matrix.derived(), result_matrix.const_cast_derived(), aot);
      Algo::run(cmodel, cdata, args);
    }

    /**
     * @brief      ConstraintModelJacobianTransposeMatrixProductVisitor visitor
     */
    template<
      typename Scalar,
      int Options,
      template<typename, int> class JointCollectionTpl,
      typename InputMatrix,
      typename OutputMatrix,
      AssignmentOperatorType op>
    struct ConstraintModelJacobianTransposeMatrixProductVisitor
    : visitors::ConstraintUnaryVisitorBase<ConstraintModelJacobianTransposeMatrixProductVisitor<
        Scalar,
        Options,
        JointCollectionTpl,
        InputMatrix,
        OutputMatrix,
        op>>
    {
      typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      typedef boost::fusion::vector<
        const Model &,
        const Data &,
        const InputMatrix &,
        OutputMatrix &,
        AssignmentOperatorTag<op>>
        ArgsType;

      template<typename ConstraintModel>
      static void algo(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const typename ConstraintModel::ConstraintData & cdata,
        const Model & model,
        const Data & data,
        const Eigen::MatrixBase<InputMatrix> & input_matrix,
        const Eigen::MatrixBase<OutputMatrix> & result_matrix,
        AssignmentOperatorTag<op> aot)
      {
        cmodel.jacobianTransposeMatrixProduct(
          model, data, cdata.derived(), input_matrix.derived(), result_matrix.const_cast_derived(),
          aot);
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class JointCollectionTpl,
      template<typename S, int O> class ConstraintCollectionTpl,
      typename InputMatrix,
      typename OutputMatrix,
      AssignmentOperatorType op = SETTO>
    void jacobianTransposeMatrixProduct(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata,
      const Eigen::MatrixBase<InputMatrix> & input_matrix,
      const Eigen::MatrixBase<OutputMatrix> & result_matrix,
      AssignmentOperatorTag<op> aot = SetTo())
    {
      typedef ConstraintModelJacobianTransposeMatrixProductVisitor<
        Scalar, Options, JointCollectionTpl, InputMatrix, OutputMatrix, op>
        Algo;

      typename Algo::ArgsType args(
        model, data, input_matrix.derived(), result_matrix.const_cast_derived(), aot);
      Algo::run(cmodel, cdata, args);
    }

    /**
     * @brief      ConstraintModelMapConstraintForceToJointSpace visitor
     */
    template<
      typename Scalar,
      int Options,
      template<typename, int> class JointCollectionTpl,
      typename ConstraintForceLike,
      class ForceAllocator,
      typename JointTorquesLike,
      ReferenceFrame rf>
    struct ConstraintModelMapConstraintForceToJointSpaceVisitor
    : visitors::ConstraintUnaryVisitorBase<ConstraintModelMapConstraintForceToJointSpaceVisitor<
        Scalar,
        Options,
        JointCollectionTpl,
        ConstraintForceLike,
        ForceAllocator,
        JointTorquesLike,
        rf>>
    {
      typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      typedef std::vector<ForceTpl<Scalar, Options>, ForceAllocator> ForceVector;
      typedef boost::fusion::vector<
        const Model &,
        const Data &,
        const ConstraintForceLike &,
        ForceVector &,
        JointTorquesLike &,
        const ReferenceFrameTag<rf>>
        ArgsType;

      template<typename ConstraintModel>
      static void algo(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const typename ConstraintModel::ConstraintData & cdata,
        const Model & model,
        const Data & data,
        const Eigen::MatrixBase<ConstraintForceLike> & constraint_forces,
        std::vector<ForceTpl<Scalar, Options>, ForceAllocator> & joint_forces,
        const Eigen::MatrixBase<JointTorquesLike> & joint_torques,
        const ReferenceFrameTag<rf> reference_frame)
      {
        cmodel.mapConstraintForceToJointSpace(
          model, data, cdata, constraint_forces, joint_forces, joint_torques.const_cast_derived(),
          reference_frame);
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename, int> class JointCollectionTpl,
      template<typename S, int O> class ConstraintCollectionTpl,
      typename ConstraintForceLike,
      class ForceAllocator,
      typename JointTorquesLike,
      ReferenceFrame rf>
    void mapConstraintForceToJointSpace(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<ConstraintForceLike> & constraint_forces,
      std::vector<ForceTpl<Scalar, Options>, ForceAllocator> & joint_forces,
      const Eigen::MatrixBase<JointTorquesLike> & joint_torques,
      const ReferenceFrameTag<rf> reference_frame)
    {
      typedef ConstraintModelMapConstraintForceToJointSpaceVisitor<
        Scalar, Options, JointCollectionTpl, ConstraintForceLike, ForceAllocator, JointTorquesLike,
        rf>
        Algo;

      typename Algo::ArgsType args(
        model, data, constraint_forces.derived(), joint_forces, joint_torques.const_cast_derived(),
        reference_frame);
      Algo::run(cmodel, cdata, args);
    }

    /**
     * @brief      ConstraintModelMapJointSpaceToConstraintMotionVisitor visitor
     */
    template<
      typename Scalar,
      int Options,
      template<typename, int> class JointCollectionTpl,
      class MotionAllocator,
      typename GeneralizedVelocityLike,
      typename ConstraintMotionLike,
      ReferenceFrame rf>
    struct ConstraintModelMapJointSpaceToConstraintMotionVisitor
    : visitors::ConstraintUnaryVisitorBase<ConstraintModelMapJointSpaceToConstraintMotionVisitor<
        Scalar,
        Options,
        JointCollectionTpl,
        MotionAllocator,
        GeneralizedVelocityLike,
        ConstraintMotionLike,
        rf>>
    {
      typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      typedef std::vector<MotionTpl<Scalar, Options>, MotionAllocator> MotionVector;
      typedef boost::fusion::vector<
        const Model &,
        const Data &,
        const MotionVector &,
        const GeneralizedVelocityLike &,
        ConstraintMotionLike &,
        const ReferenceFrameTag<rf>>
        ArgsType;

      template<typename ConstraintModel>
      static void algo(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const typename ConstraintModel::ConstraintData & cdata,
        const Model & model,
        const Data & data,
        const std::vector<MotionTpl<Scalar, Options>, MotionAllocator> & joint_motions,
        const Eigen::MatrixBase<GeneralizedVelocityLike> & generalized_velocity,
        const Eigen::MatrixBase<ConstraintMotionLike> & constraint_motions,
        const ReferenceFrameTag<rf> reference_frame)
      {
        cmodel.mapJointSpaceToConstraintMotion(
          model, data, cdata, joint_motions, generalized_velocity.derived(),
          constraint_motions.const_cast_derived(), reference_frame);
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename, int> class JointCollectionTpl,
      template<typename S, int O> class ConstraintCollectionTpl,
      class MotionAllocator,
      typename GeneralizedVelocityLike,
      typename ConstraintMotionLike,
      ReferenceFrame rf>
    void mapJointSpaceToConstraintMotion(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const std::vector<MotionTpl<Scalar, Options>, MotionAllocator> & joint_motions,
      const Eigen::MatrixBase<GeneralizedVelocityLike> & generalized_velocity,
      const Eigen::MatrixBase<ConstraintMotionLike> & constraint_motions,
      const ReferenceFrameTag<rf> reference_frame)
    {
      typedef ConstraintModelMapJointSpaceToConstraintMotionVisitor<
        Scalar, Options, JointCollectionTpl, MotionAllocator, GeneralizedVelocityLike,
        ConstraintMotionLike, rf>
        Algo;

      typename Algo::ArgsType args(
        model, data, joint_motions, generalized_velocity.derived(),
        constraint_motions.const_cast_derived(), reference_frame);
      Algo::run(cmodel, cdata, args);
    }

    /**
     * @brief      ConstraintModelAppendCouplingConstraintInertiasVisitor visitor
     */
    template<
      typename Scalar,
      int Options,
      template<typename, int> class JointCollectionTpl,
      typename VectorNLike,
      ReferenceFrame rf>
    struct ConstraintModelAppendCouplingConstraintInertiasVisitor
    : visitors::ConstraintUnaryVisitorBase<ConstraintModelAppendCouplingConstraintInertiasVisitor<
        Scalar,
        Options,
        JointCollectionTpl,
        VectorNLike,
        rf>>
    {
      typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      typedef boost::fusion::
        vector<const Model &, Data &, const VectorNLike &, ReferenceFrameTag<rf>>
          ArgsType;

      template<typename ConstraintModel>
      static void algo(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const typename ConstraintModel::ConstraintData & cdata,
        const Model & model,
        Data & data,
        const Eigen::MatrixBase<VectorNLike> & diagonal_constraint_inertia,
        const ReferenceFrameTag<rf> reference_frame)
      {
        cmodel.appendCouplingConstraintInertias(
          model, data, cdata.derived(), diagonal_constraint_inertia.derived(), reference_frame);
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class JointCollectionTpl,
      template<typename S, int O> class ConstraintCollectionTpl,
      typename VectorNLike,
      ReferenceFrame rf>
    void appendCouplingConstraintInertias(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> & cdata,
      const Eigen::MatrixBase<VectorNLike> & diagonal_constraint_inertia,
      const ReferenceFrameTag<rf> reference_frame)
    {
      typedef ConstraintModelAppendCouplingConstraintInertiasVisitor<
        Scalar, Options, JointCollectionTpl, VectorNLike, rf>
        Algo;

      typename Algo::ArgsType args(
        model, data, diagonal_constraint_inertia.derived(), reference_frame);
      Algo::run(cmodel, cdata, args);
    }

    /**
     * @brief      ConstraintModelComplianceMemberVisitor visitor
     */
    /// \brief ComplianceMemberGetter - default behavior for false for
    /// HasComplianceMember
    template<bool HasComplianceMember, typename ReturnType>
    struct ComplianceMemberGetter
    {
      template<typename ConstraintModelDerived>
      static ReturnType run(const ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        std::stringstream ss;
        ss << cmodel.shortname() << " does not have the compliance member.\n";
        PINOCCHIO_THROW(std::invalid_argument, ss.str());
        return internal::NoRun<ReturnType>::run();
      }
      template<typename ConstraintModelDerived>
      static ReturnType run(ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        std::stringstream ss;
        ss << cmodel.shortname() << " does not have the compliance member.\n";
        PINOCCHIO_THROW(std::invalid_argument, ss.str());
        return internal::NoRun<ReturnType>::run();
      }
    };

    /// \brief ComplianceMemberGetter - partial specialization for true for
    /// HasComplianceMember
    template<typename ReturnType>
    struct ComplianceMemberGetter<true, ReturnType>
    {
      template<typename ConstraintModelDerived>
      static ReturnType run(const ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        return cmodel.compliance();
      }
      template<typename ConstraintModelDerived>
      static ReturnType run(ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        return cmodel.compliance();
      }
    };

    template<typename ReturnType>
    struct ConstraintModelComplianceMemberVisitor
    : ConstraintUnaryVisitorBase<ConstraintModelComplianceMemberVisitor<ReturnType>, ReturnType>
    {
      typedef NoArg ArgsType;

      template<typename ConstraintModelDerived>
      static ReturnType algo(const ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        static constexpr bool has_compliance_member =
          traits<ConstraintModelDerived>::has_compliance_member;
        return ComplianceMemberGetter<has_compliance_member, ReturnType>::run(cmodel);
      }

      template<typename ConstraintModelDerived>
      static ReturnType algo(ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        static constexpr bool has_compliance_member =
          traits<ConstraintModelDerived>::has_compliance_member;
        return ComplianceMemberGetter<has_compliance_member, ReturnType>::run(cmodel);
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    typename traits<
      ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl>>::ComplianceVectorTypeConstRef
    getCompliance(const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
    {
      typedef typename traits<ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl>>::
        ComplianceVectorTypeConstRef ReturnType;
      typedef ConstraintModelComplianceMemberVisitor<ReturnType> Algo;
      return Algo::run(cmodel);
    }

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    typename traits<
      ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl>>::ComplianceVectorTypeRef
    getCompliance(ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
    {
      typedef typename traits<ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl>>::
        ComplianceVectorTypeRef ReturnType;
      typedef ConstraintModelComplianceMemberVisitor<ReturnType> Algo;
      return Algo::run(cmodel);
    }

    /**
     * @brief      ConstraintModelBaumgarteCorrectorParametersVisitor visitor
     */
    /// \brief BaumgarteCorrectorParametersGetter - default behavior for false for
    /// HasBaumgarteCorrector
    template<bool HasBaumgarteCorrector, typename BaumgarteReturnType>
    struct BaumgarteCorrectorParametersGetter
    {
      template<typename ConstraintModelDerived>
      static BaumgarteReturnType run(const ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        std::stringstream ss;
        ss << cmodel.shortname() << " does not have baumgarte corrector parameters.\n";
        PINOCCHIO_THROW(std::invalid_argument, ss.str());
        return internal::NoRun<BaumgarteReturnType>::run();
      }
      template<typename ConstraintModelDerived>
      static BaumgarteReturnType run(ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        std::stringstream ss;
        ss << cmodel.shortname() << " does not have baumgarte corrector parameters.\n";
        PINOCCHIO_THROW(std::invalid_argument, ss.str());
        return internal::NoRun<BaumgarteReturnType>::run();
      }
    };

    /// \brief BaumgarteCorrectorParametersGetter - partial specialization for true for
    /// HasBaumgarteCorrector
    template<typename BaumgarteReturnType>
    struct BaumgarteCorrectorParametersGetter<true, BaumgarteReturnType>
    {
      template<typename ConstraintModelDerived>
      static BaumgarteReturnType run(const ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        return cmodel.baumgarte_corrector_parameters();
      }
      template<typename ConstraintModelDerived>
      static BaumgarteReturnType run(ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        return cmodel.baumgarte_corrector_parameters();
      }
    };

    template<typename BaumgarteReturnType>
    struct ConstraintModelBaumgarteCorrectorParametersVisitor
    : ConstraintUnaryVisitorBase<
        ConstraintModelBaumgarteCorrectorParametersVisitor<BaumgarteReturnType>,
        BaumgarteReturnType>
    {
      typedef NoArg ArgsType;

      template<typename ConstraintModelDerived>
      static BaumgarteReturnType algo(const ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        static constexpr bool has_baumgarte_corrector =
          traits<ConstraintModelDerived>::has_baumgarte_corrector;
        return BaumgarteCorrectorParametersGetter<
          has_baumgarte_corrector, BaumgarteReturnType>::run(cmodel);
      }

      template<typename ConstraintModelDerived>
      static BaumgarteReturnType algo(ConstraintModelBase<ConstraintModelDerived> & cmodel)
      {
        static constexpr bool has_baumgarte_corrector =
          traits<ConstraintModelDerived>::has_baumgarte_corrector;
        return BaumgarteCorrectorParametersGetter<
          has_baumgarte_corrector, BaumgarteReturnType>::run(cmodel);
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    const BaumgarteCorrectorParametersTpl<Scalar> & getBaumgarteCorrectorParameters(
      const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
    {
      typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;
      typedef const BaumgarteCorrectorParameters & ReturnType;
      typedef ConstraintModelBaumgarteCorrectorParametersVisitor<ReturnType> Algo;
      return Algo::run(cmodel);
    }

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    BaumgarteCorrectorParametersTpl<Scalar> & getBaumgarteCorrectorParameters(
      ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel)
    {
      typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;
      typedef BaumgarteCorrectorParameters & ReturnType;
      typedef ConstraintModelBaumgarteCorrectorParametersVisitor<ReturnType> Algo;
      return Algo::run(cmodel);
    }

    /**
     * @brief      ConstraintModelSetComplianceVisitor visitor
     */
    template<typename InputVector>
    struct ConstraintModelSetComplianceVisitor
    : visitors::ConstraintUnaryVisitorBase<ConstraintModelSetComplianceVisitor<InputVector>>
    {
      typedef boost::fusion::vector<const InputVector &> ArgsType;

      template<typename ConstraintModel>
      static void algo(
        pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const Eigen::MatrixBase<InputVector> & input_vector)
      {
        cmodel.setCompliance(input_vector.derived());
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl,
      typename InputVector>
    void setCompliance(
      ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const Eigen::MatrixBase<InputVector> & input_vector)
    {
      typedef ConstraintModelSetComplianceVisitor<InputVector> Algo;

      typename Algo::ArgsType args(input_vector.derived());
      Algo::run(cmodel, args);
    }

    /**
     * @brief      ConstraintModelSetBaumgarteCorrectorParameters visitor
     */
    template<typename _BaumgarteCorrectorParameters>
    struct ConstraintModelSetBaumgarteCorrectorParameters
    : visitors::ConstraintUnaryVisitorBase<
        ConstraintModelSetBaumgarteCorrectorParameters<_BaumgarteCorrectorParameters>>
    {
      typedef boost::fusion::vector<const _BaumgarteCorrectorParameters &> ArgsType;

      template<typename ConstraintModel>
      static void algo(
        pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        const _BaumgarteCorrectorParameters & baumgarte_corrector_parameters_in)
      {
        cmodel.setBaumgarteCorrectorParameters(baumgarte_corrector_parameters_in);
      }
    };

    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    void setBaumgarteCorrectorParameters(
      ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const BaumgarteCorrectorParametersTpl<Scalar> & baumgarte_corrector_parameters_in)
    {
      typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;
      typedef ConstraintModelSetBaumgarteCorrectorParameters<BaumgarteCorrectorParameters> Algo;

      typename Algo::ArgsType args(baumgarte_corrector_parameters_in);
      Algo::run(cmodel, args);
    }

  } // namespace visitors

} // namespace pinocchio

#endif // ifdef __pinocchio_algorithm_constraints_constraint_model_visitor_hpp__
