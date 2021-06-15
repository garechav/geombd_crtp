#ifndef GEOMBD_FUNCTORS_HXX
#define GEOMBD_FUNCTORS_HXX

#define EIGEN_NO_DEBUG
#define EIGEN_MPL2_ONLY
//#define EIGEN_UNROLLING_LIMIT 30

#include "Eigen/Core"
#include "geombd/XprHelper.h"

namespace helper
{
  template<typename T> struct argument_type;
  template<typename T, typename U> struct argument_type<T(U)> { typedef U type; };
}

#define GEOMBD_EIGEN_PLAIN_TYPE(T) Eigen::internal::plain_matrix_type< typename helper::argument_type<void(T)>::type >::type


//! Overriden/Overloading operators for static expressions
//! see https://eigen.tuxfamily.org/dox/TopicNewExpressionType.html
//! and http://library.isr.ist.utl.pt/docs/roswiki/eigen(2f)Cookbook.html
//!------------------------------------------------------------------------------!//

//! Skew Symmetric Matrix
//!------------------------------------------------------------------------------!//
template <class ArgType> class SkewSym;

namespace Eigen {
  namespace internal {
    template <class ArgType>
    struct traits<SkewSym<ArgType> >
    {
      typedef Eigen::Dense StorageKind;
      typedef Eigen::MatrixXpr XprKind;
      typedef typename ArgType::StorageIndex StorageIndex;
      typedef typename ArgType::Scalar Scalar;
      enum {
        Flags = Eigen::ColMajor,
        RowsAtCompileTime = 3,//ArgType::RowsAtCompileTime,
        ColsAtCompileTime = 3,//ArgType::RowsAtCompileTime,
        MaxRowsAtCompileTime = 3,//ArgType::MaxRowsAtCompileTime,
        MaxColsAtCompileTime = 3,//ArgType::MaxRowsAtCompileTime
      };
    };
  }
}

template <class ArgType>
class SkewSym : public Eigen::MatrixBase<SkewSym<ArgType> >
{
public:
  SkewSym(const ArgType& arg)
    : m_arg(arg)
  {
    EIGEN_STATIC_ASSERT(ArgType::ColsAtCompileTime == 1,
                        YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_A_MATRIX);
  }

  typedef typename Eigen::internal::ref_selector<SkewSym>::type Nested;

  typedef Eigen::Index Index;
  Index rows() const { return 3; }//m_arg.rows(); }
  Index cols() const { return 3; }//m_arg.rows(); }

  typedef typename Eigen::internal::ref_selector<ArgType>::type ArgTypeNested;
  ArgTypeNested m_arg;
};

namespace Eigen {
  namespace internal {
    template<typename ArgType>
    struct evaluator<SkewSym<ArgType> >
        : evaluator_base<SkewSym<ArgType> >
    {
      typedef SkewSym<ArgType> XprType;
      typedef typename nested_eval<ArgType, XprType::ColsAtCompileTime>::type ArgTypeNested;
      typedef typename remove_all<ArgTypeNested>::type ArgTypeNestedCleaned;
      typedef typename XprType::CoeffReturnType CoeffReturnType;

      enum {
        CoeffReadCost = evaluator<ArgTypeNestedCleaned>::CoeffReadCost,
        Flags = Eigen::ColMajor
      };

      evaluator(const XprType& xpr)
        : m_argImpl(xpr.m_arg)
      { }

      CoeffReturnType coeff(Index row, Index col) const
      {
        CoeffReturnType out = 0;

        if(row == 0){
            if(col == 1) out = -m_argImpl.coeff(2);
            if(col == 2) out =  m_argImpl.coeff(1);
          }else if(row == 1){
            if(col == 0) out =  m_argImpl.coeff(2);
            if(col == 2) out = -m_argImpl.coeff(0);
          }else if(row == 2){
            if(col == 0) out = -m_argImpl.coeff(1);
            if(col == 1) out =  m_argImpl.coeff(0);
          }

        return out;
      }

      evaluator<ArgTypeNestedCleaned> m_argImpl;
      const Index m_rows = 3;
    };
  }
}

template <class ArgType>
SkewSym<ArgType> Skew(const Eigen::MatrixBase<ArgType>& arg)
{
  return SkewSym<ArgType>(arg.derived());
}


//! My Static Cross Product
//!------------------------------------------------------------------------------!//
template <class ArgType> class CrossPro;

namespace Eigen {
  namespace internal {
    template <class ArgType>
    struct traits<CrossPro<ArgType> >
    {
      typedef Eigen::Dense StorageKind;
      typedef Eigen::MatrixXpr XprKind;
      typedef typename ArgType::StorageIndex StorageIndex;
      typedef typename ArgType::Scalar Scalar;
      enum {
        Flags = Eigen::ColMajor,
        RowsAtCompileTime = 3,
        ColsAtCompileTime = 1,
        MaxRowsAtCompileTime = 3,
        MaxColsAtCompileTime = 1,
      };
    };
  }
}

template <class ArgType>
class CrossPro : public Eigen::MatrixBase<CrossPro<ArgType> >
{
public:
  CrossPro(const ArgType& arg1, const ArgType& arg2)
    : m_arg1(arg1), m_arg2(arg2)
  {
    EIGEN_STATIC_ASSERT(ArgType::ColsAtCompileTime == 1,
                        YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_A_MATRIX);
  }

  typedef typename Eigen::internal::ref_selector<CrossPro>::type Nested;

  typedef Eigen::Index Index;
  Index rows() const { return 3; }
  Index cols() const { return 1; }

  typedef typename Eigen::internal::ref_selector<ArgType>::type ArgTypeNested;
  ArgTypeNested m_arg1, m_arg2;
};

namespace Eigen {
  namespace internal {
    template<typename ArgType>
    struct evaluator<CrossPro<ArgType> >
        : evaluator_base<CrossPro<ArgType> >
    {
      typedef CrossPro<ArgType> XprType;
      typedef typename nested_eval<ArgType, XprType::ColsAtCompileTime>::type ArgTypeNested;
      typedef typename remove_all<ArgTypeNested>::type ArgTypeNestedCleaned;
      typedef typename XprType::CoeffReturnType CoeffReturnType;

      enum {
        CoeffReadCost = evaluator<ArgTypeNestedCleaned>::CoeffReadCost,
        Flags = Eigen::ColMajor
      };

      evaluator(const XprType& xpr)
        : m_argImpl1(xpr.m_arg1), m_argImpl2(xpr.m_arg2)
      { }

      CoeffReturnType coeff(Index row, Index col) const
      {
        CoeffReturnType out;

        if(row == 0) out = m_argImpl1.coeff(1)*m_argImpl2.coeff(2) - m_argImpl1.coeff(2)*m_argImpl2.coeff(1);
        if(row == 0) out = m_argImpl1.coeff(2)*m_argImpl2.coeff(0) - m_argImpl1.coeff(0)*m_argImpl2.coeff(2);
        if(row == 0) out = m_argImpl1.coeff(0)*m_argImpl2.coeff(1) - m_argImpl1.coeff(1)*m_argImpl2.coeff(0);

        return out;
      }

      evaluator<ArgTypeNestedCleaned> m_argImpl1, m_argImpl2;
      const Index m_rows = 3;
    };
  }
}

template <class ArgType>
CrossPro<ArgType> iCross(const Eigen::MatrixBase<ArgType>& arg1, const Eigen::MatrixBase<ArgType>& arg2)
{
  return CrossPro<ArgType>(arg1.derived(), arg2.derived());
}


//! overload operators from
//! https://eigen.tuxfamily.org/dox/UserManual_CustomizingEigen.html


//! Overloading Spatial Matrix Operations
//!------------------------------------------------------------------------------!//
//EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
typedef Eigen::Matrix<double, 6, 6, 0, 6, 6> SpatialMatrixBase;

class MyMatrixType : public SpatialMatrixBase
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  MyMatrixType(void):SpatialMatrixBase() {}

  // This constructor allows you to construct MyVectorType from Eigen expressions
  template<typename OtherDerived>
  MyMatrixType(const Eigen::MatrixBase<OtherDerived>& other)
    : SpatialMatrixBase(other)
  { }

  // This method allows you to assign Eigen expressions to MyVectorType
  template<typename OtherDerived>
  MyMatrixType& operator=(const Eigen::MatrixBase <OtherDerived>& other)
  {
    this->SpatialMatrixBase::operator=(other);
    return *this;
  }
};


//https://stackoverflow.com/questions/65341405/eigen-derived-object-cannot-make-operation-with-scalar


namespace geo{

  //! Performing R*A*R.transpose() statically for Rx joint type
  //!------------------------------------------------------------------------------!//
  // Forward declaration
  template<typename ScalarType, typename Matrix3TypeIn, typename Matrix3TypeOut> struct Mat3ProjRxAlgo;

  template<typename ScalarType, typename Matrix3TypeIn, typename Matrix3TypeOut>
  void Mat3ProjRx(const ScalarType & sq,
                  const ScalarType & cq,
                  const Eigen::MatrixBase<Matrix3TypeIn> & MatIn,
                  Eigen::MatrixBase<Matrix3TypeOut> & MatOut)
  { Mat3ProjRxAlgo<ScalarType,Matrix3TypeIn,Matrix3TypeOut>::run(sq, cq, MatIn, MatOut); }

  template<typename ScalarType, typename Matrix3TypeIn, typename Matrix3TypeOut>
  struct Mat3ProjRxAlgo
  {
  public :
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inline static void run(const ScalarType & sq,
                           const ScalarType & cq,
                           const Eigen::MatrixBase<Matrix3TypeIn> & MatIn,
                           Eigen::MatrixBase<Matrix3TypeOut> & MatOut)
    {
      typename Matrix3TypeIn::Scalar a, b, c, d, e, f;

      a = MatIn.coeff(0,0);  b = MatIn.coeff(0,1);  c = MatIn.coeff(0,2);
      d = MatIn.coeff(1,1);  e = MatIn.coeff(1,2);  f = MatIn.coeff(2,2);

      ScalarType s2, c2, sc, c2s2, esc_;

      s2 = sq*sq;  c2 = cq*cq;  sc = sq*cq;  c2s2 = c2 - s2;

      esc_ = 2*e*sc;
      MatOut.coeffRef(0,0) = a;
      MatOut.coeffRef(1,0) = b*cq - c*sq;
      MatOut.coeffRef(2,0) = c*cq + b*sq;

      MatOut.coeffRef(0,1) = MatOut.coeff(1,0);
      MatOut.coeffRef(1,1) = d*c2 - esc_ + f*s2;
      MatOut.coeffRef(2,1) = e*c2s2 + (d-f)*sc;

      MatOut.coeffRef(0,2) = MatOut.coeff(2,0);
      MatOut.coeffRef(1,2) = MatOut.coeff(2,1);
      MatOut.coeffRef(2,2) = f*c2 + esc_ + d*s2;
    }
  };


  //! Performing Ad_dual*M*Ad statically for Rx joint type
  //!------------------------------------------------------------------------------!//
  // Forward declaration
  template<typename Vector3Type, typename Matrix3Type, typename Matrix6TypeIn, typename Matrix6TypeOut> struct Mat6ProjRxAlgo;

  template<typename Vector3Type, typename Matrix3Type, typename Matrix6TypeIn, typename Matrix6TypeOut>
  void Mat6ProjRx(bool P_z,
                  Eigen::MatrixBase<Vector3Type> & P_,
                  Eigen::MatrixBase<Matrix3Type> & R_,
                  const Eigen::MatrixBase<Matrix6TypeIn> & MatIn,
                  Eigen::MatrixBase<Matrix6TypeOut> & MatOut)
  { Mat6ProjRxAlgo<Vector3Type, Matrix3Type, Matrix6TypeIn, Matrix6TypeOut>::run(P_z, P_, R_, MatIn, MatOut); }

  template<typename Vector3Type, typename Matrix3Type, typename Matrix6TypeIn, typename Matrix6TypeOut>
  struct Mat6ProjRxAlgo
  {
  public :
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inline static void run(bool P_z,
                           Eigen::MatrixBase<Vector3Type> & P_,
                           Eigen::MatrixBase<Matrix3Type> & R_,
                           const Eigen::MatrixBase<Matrix6TypeIn> & MatIn,
                           Eigen::MatrixBase<Matrix6TypeOut> & MatOut)
    {
      typedef const Eigen::Block<Matrix6TypeIn,3,3> constBlock3;
      typedef Eigen::Block<Matrix6TypeOut,3,3> Block3;

      typename Matrix6TypeIn::PlainObject & MatIn_ = const_cast<Matrix6TypeIn &>(MatIn.derived());
      constBlock3 & Ai = MatIn_.template block<3,3>(0,0);
      constBlock3 & Bi = MatIn_.template block<3,3>(0,3);
      constBlock3 & Di = MatIn_.template block<3,3>(3,3);

      Block3 Ao = MatOut.template block<3,3>(0,0);
      Block3 Bo = MatOut.template block<3,3>(0,3);
      Block3 Co = MatOut.template block<3,3>(3,0);
      Block3 Do = MatOut.template block<3,3>(3,3);

      typedef typename Matrix3Type::Scalar MyScalar;
      MyScalar sq, cq;
      sq = R_.coeff(2,1);  cq = R_.coeff(1,1);

      Mat3ProjRx(sq, cq, Ai, Ao);

      Do.noalias() = R_*Bi; // tmp variable
      Co.noalias() = Do*R_.transpose();

      Mat3ProjRx(sq, cq, Di, Do);

      Bo = Co;
      if(P_z){
          //! Linear
          typename Matrix3Type::PlainObject SkP = Skew(P_);
          Bo.noalias() -= Ao*SkP;

          typename Matrix3Type::PlainObject Dtmp1;
          Dtmp1.noalias() = SkP*Bo;
          Do.noalias() += Dtmp1;
          Do.noalias() -= Co.transpose()*SkP;
        }

      Co = Bo.transpose();
    }
  };


  //! Performing R*A*R.transpose() statically for Ry joint type
  //!------------------------------------------------------------------------------!//
  // Forward declaration
  template<typename ScalarType, typename Matrix3TypeIn, typename Matrix3TypeOut> struct Mat3ProjRyAlgo;

  template<typename ScalarType, typename Matrix3TypeIn, typename Matrix3TypeOut>
  void Mat3ProjRy(const ScalarType & sq,
                  const ScalarType & cq,
                  const Eigen::MatrixBase<Matrix3TypeIn> & MatIn,
                  Eigen::MatrixBase<Matrix3TypeOut> & MatOut)
  { Mat3ProjRyAlgo<ScalarType,Matrix3TypeIn,Matrix3TypeOut>::run(sq, cq, MatIn, MatOut); }

  template<typename ScalarType, typename Matrix3TypeIn, typename Matrix3TypeOut>
  struct Mat3ProjRyAlgo
  {
  public :
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inline static void run(const ScalarType & sq,
                           const ScalarType & cq,
                           const Eigen::MatrixBase<Matrix3TypeIn> & MatIn,
                           Eigen::MatrixBase<Matrix3TypeOut> & MatOut)
    {
      typename Matrix3TypeIn::Scalar a, b, c, d, e, f;

      a = MatIn.coeff(0,0);  b = MatIn.coeff(0,1);  c = MatIn.coeff(0,2);
      d = MatIn.coeff(1,1);  e = MatIn.coeff(1,2);  f = MatIn.coeff(2,2);

      ScalarType s2, c2, sc, c2s2, csc;

      s2 = sq*sq;  c2 = cq*cq;  sc = sq*cq;  c2s2 = c2 - s2;

      csc = 2*c*sc;
      MatOut.coeffRef(0,0) = a*c2 + csc + f*s2;
      MatOut.coeffRef(1,0) = b*cq + e*sq;
      MatOut.coeffRef(2,0) = (f-a)*sc + c*c2s2;

      MatOut.coeffRef(0,1) = MatOut.coeff(1,0);
      MatOut.coeffRef(1,1) = d;
      MatOut.coeffRef(2,1) = cq*e - b*sq;

      MatOut.coeffRef(0,2) = MatOut.coeff(2,0);
      MatOut.coeffRef(1,2) = MatOut.coeff(2,1);
      MatOut.coeffRef(2,2) = f*c2 - csc + a*s2;
    }
  };


  //! Performing Ad_dual*M*Ad statically for Ry joint type
  //!------------------------------------------------------------------------------!//
  // Forward declaration
  template<typename Vector3Type, typename Matrix3Type, typename Matrix6TypeIn, typename Matrix6TypeOut> struct Mat6ProjRyAlgo;

  template<typename Vector3Type, typename Matrix3Type, typename Matrix6TypeIn, typename Matrix6TypeOut>
  void Mat6ProjRy(bool P_z,
                  Eigen::MatrixBase<Vector3Type> & P_,
                  Eigen::MatrixBase<Matrix3Type> & R_,
                  const Eigen::MatrixBase<Matrix6TypeIn> & MatIn,
                  Eigen::MatrixBase<Matrix6TypeOut> & MatOut)
  { Mat6ProjRyAlgo<Vector3Type, Matrix3Type, Matrix6TypeIn, Matrix6TypeOut>::run(P_z, P_, R_, MatIn, MatOut); }

  template<typename Vector3Type, typename Matrix3Type, typename Matrix6TypeIn, typename Matrix6TypeOut>
  struct Mat6ProjRyAlgo
  {
  public :
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inline static void run(bool P_z,
                           Eigen::MatrixBase<Vector3Type> & P_,
                           Eigen::MatrixBase<Matrix3Type> & R_,
                           const Eigen::MatrixBase<Matrix6TypeIn> & MatIn,
                           Eigen::MatrixBase<Matrix6TypeOut> & MatOut)
    {
      typedef const Eigen::Block<Matrix6TypeIn,3,3> constBlock3;
      typedef Eigen::Block<Matrix6TypeOut,3,3> Block3;

      typename Matrix6TypeIn::PlainObject & MatIn_ = const_cast<Matrix6TypeIn &>(MatIn.derived());
      constBlock3 & Ai = MatIn_.template block<3,3>(0,0);
      constBlock3 & Bi = MatIn_.template block<3,3>(0,3);
      constBlock3 & Di = MatIn_.template block<3,3>(3,3);

      Block3 Ao = MatOut.template block<3,3>(0,0);
      Block3 Bo = MatOut.template block<3,3>(0,3);
      Block3 Co = MatOut.template block<3,3>(3,0);
      Block3 Do = MatOut.template block<3,3>(3,3);

      typedef typename Matrix3Type::Scalar MyScalar;
      MyScalar sq, cq;
      sq = R_.coeff(0,2);  cq = R_.coeff(0,0);

      Mat3ProjRy(sq, cq, Ai, Ao);

      Do.noalias() = R_*Bi; // tmp variable
      Co.noalias() = Do*R_.transpose();

      Mat3ProjRy(sq, cq, Di, Do);

      Bo = Co;
      if(P_z){
          //! Linear
          typename Matrix3Type::PlainObject SkP = Skew(P_);
          Bo.noalias() -= Ao*SkP;

          typename Matrix3Type::PlainObject Dtmp1;
          Dtmp1.noalias() = SkP*Bo;
          Do.noalias() += Dtmp1;
          Do.noalias() -= Co.transpose()*SkP;
        }

      Co = Bo.transpose();
    }
  };


  //! Performing R*A*R.transpose() statically for Rz joint type
  //!------------------------------------------------------------------------------!//
  // Forward declaration
  template<typename ScalarType, typename Matrix3TypeIn, typename Matrix3TypeOut> struct Mat3ProjRzAlgo;

  template<typename ScalarType, typename Matrix3TypeIn, typename Matrix3TypeOut>
  void Mat3ProjRz(const ScalarType & sq,
                  const ScalarType & cq,
                  const Eigen::MatrixBase<Matrix3TypeIn> & MatIn,
                  Eigen::MatrixBase<Matrix3TypeOut> & MatOut)
  { Mat3ProjRzAlgo<ScalarType,Matrix3TypeIn,Matrix3TypeOut>::run(sq, cq, MatIn, MatOut); }

  template<typename ScalarType, typename Matrix3TypeIn, typename Matrix3TypeOut>
  struct Mat3ProjRzAlgo
  {
  public :
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inline static void run(const ScalarType & sq,
                           const ScalarType & cq,
                           const Eigen::MatrixBase<Matrix3TypeIn> & MatIn,
                           Eigen::MatrixBase<Matrix3TypeOut> & MatOut)
    {
      typename Matrix3TypeIn::Scalar a, b, c, d, e, f;

      a = MatIn.coeff(0,0);  b = MatIn.coeff(0,1);  c = MatIn.coeff(0,2);
      d = MatIn.coeff(1,1);  e = MatIn.coeff(1,2);  f = MatIn.coeff(2,2);

      ScalarType s2, c2, sc, c2s2, bsc;

      s2 = sq*sq;  c2 = cq*cq;  sc = sq*cq;  c2s2 = c2 - s2;

      bsc = 2*b*sc;
      MatOut.coeffRef(0,0) = a*c2 - bsc + d*s2;
      MatOut.coeffRef(1,0) = (a-d)*sc + b*c2s2;
      MatOut.coeffRef(2,0) = c*cq - e*sq;

      MatOut.coeffRef(0,1) = MatOut.coeff(1,0);
      MatOut.coeffRef(1,1) = a*s2 + bsc + d*c2;
      MatOut.coeffRef(2,1) = c*sq + e*cq;

      MatOut.coeffRef(0,2) = MatOut.coeff(2,0);
      MatOut.coeffRef(1,2) = MatOut.coeff(2,1);
      MatOut.coeffRef(2,2) = f;
    }
  };


  //! Performing Ad_dual*M*Ad statically for Rz joint type
  //!------------------------------------------------------------------------------!//
  // Forward declaration
  template<typename Vector3Type, typename Matrix3Type, typename Matrix6TypeIn, typename Matrix6TypeOut> struct Mat6ProjRzAlgo;

  template<typename Vector3Type, typename Matrix3Type, typename Matrix6TypeIn, typename Matrix6TypeOut>
  void Mat6ProjRz(bool P_z,
                  Eigen::MatrixBase<Vector3Type> & P_,
                  Eigen::MatrixBase<Matrix3Type> & R_,
                  const Eigen::MatrixBase<Matrix6TypeIn> & MatIn,
                  Eigen::MatrixBase<Matrix6TypeOut> & MatOut)
  { Mat6ProjRzAlgo<Vector3Type, Matrix3Type, Matrix6TypeIn, Matrix6TypeOut>::run(P_z, P_, R_, MatIn, MatOut); }

  template<typename Vector3Type, typename Matrix3Type, typename Matrix6TypeIn, typename Matrix6TypeOut>
  struct Mat6ProjRzAlgo
  {
  public :
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inline static void run(bool P_z,
                           Eigen::MatrixBase<Vector3Type> & P_,
                           Eigen::MatrixBase<Matrix3Type> & R_,
                           const Eigen::MatrixBase<Matrix6TypeIn> & MatIn,
                           Eigen::MatrixBase<Matrix6TypeOut> & MatOut)
    {
      typedef const Eigen::Block<Matrix6TypeIn,3,3> constBlock3;
      typedef Eigen::Block<Matrix6TypeOut,3,3> Block3;

      typename Matrix6TypeIn::PlainObject & MatIn_ = const_cast<Matrix6TypeIn &>(MatIn.derived());
      constBlock3 & Ai = MatIn_.template block<3,3>(0,0);
      constBlock3 & Bi = MatIn_.template block<3,3>(0,3);
      constBlock3 & Di = MatIn_.template block<3,3>(3,3);

      Block3 Ao = MatOut.template block<3,3>(0,0);
      Block3 Bo = MatOut.template block<3,3>(0,3);
      Block3 Co = MatOut.template block<3,3>(3,0);
      Block3 Do = MatOut.template block<3,3>(3,3);

      typedef typename Matrix3Type::Scalar MyScalar;
      MyScalar sq, cq;
      sq = R_.coeff(1,0);  cq = R_.coeff(0,0);

      Mat3ProjRz(sq, cq, Ai, Ao);

      Do.noalias() = R_*Bi; // tmp variable
      Co.noalias() = Do*R_.transpose();

      Mat3ProjRz(sq, cq, Di, Do);

      Bo = Co;
      if(P_z){
          //! Linear
          typename Matrix3Type::PlainObject SkP = Skew(P_);
          Bo.noalias() -= Ao*SkP;

          typename Matrix3Type::PlainObject Dtmp1;
          Dtmp1.noalias() = SkP*Bo;
          Do.noalias() += Dtmp1;
          Do.noalias() -= Co.transpose()*SkP;
        }

      Co = Bo.transpose();
    }
  };

}


//! First partial differentiation structs
//!------------------------------------------------------------------------------!//
//!------------------------------------------------------------------------------!//
//!------------------------------------------------------------------------------!//
namespace geo{

  //! Performing the permutation operator
  //!------------------------------------------------------------------------------!//
  // Forward declaration
  template<typename Vector6Type, typename Matrix6TypeIn, typename Matrix6TypeOut> struct PermutatorAlgo;

  template<typename Vector6Type, typename Matrix6TypeIn, typename Matrix6TypeOut>
  void permutator(const Vector6Type & Twist,
                  const Eigen::MatrixBase<Matrix6TypeIn> & Inertia,
                  Eigen::MatrixBase<Matrix6TypeOut> & D_p_aux)
  { PermutatorAlgo<Vector6Type, Matrix6TypeIn, Matrix6TypeOut>::run(Twist, Inertia, D_p_aux); }

  template<typename Vector6Type, typename Matrix6TypeIn, typename Matrix6TypeOut>
  struct PermutatorAlgo
  {
  public :
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inline static void run(const Vector6Type & T_,
                           const Eigen::MatrixBase<Matrix6TypeIn> & Inertia,
                           Eigen::MatrixBase<Matrix6TypeOut> & D_p_)
    {
      typedef Eigen::Block<Vector6Type,3,1> Segment3;
      typedef const Eigen::Block<Vector6Type,3,1> constSegment3;

      Vector6Type & T = const_cast<Vector6Type &>(T_.derived());
      constSegment3 & T_v = T.template segment<3>(0);
      constSegment3 & T_w = T.template segment<3>(3);

      typedef Eigen::Block<Matrix6TypeOut,3,3> Block3;
      Block3 UL = D_p_.template block<3,3>(0,0);  Block3 UR = D_p_.template block<3,3>(0,3);
      Block3 LL = D_p_.template block<3,3>(3,0);  Block3 LR = D_p_.template block<3,3>(3,3);

      //! adDual
      UL = Skew(T_w);  UR.setZero();
      LL = Skew(T_v);  LR = UL;
      typename Matrix6TypeOut::PlainObject D_aux;
      D_aux.noalias() = D_p_*Inertia;  // tmp

      //! adBar
      typename Vector6Type::PlainObject momentum;
      momentum.noalias() = Inertia*T_;
      Segment3 Mu_f = momentum.template segment<3>(0);
      Segment3 Mu_m = momentum.template segment<3>(3);

      UL.setZero();  UR = -Skew(Mu_f);
      LL = UR;       LR = -Skew(Mu_m);

      D_p_ += D_aux;
    }
  };


  //! Adjoint dual operator
  //!------------------------------------------------------------------------------!//
  // Forward declaration
  template<typename Vector3Type, typename Matrix3Type, typename Matrix6Type> struct AdDualAlgo;

  template<typename Vector3Type, typename Matrix3Type, typename Matrix6Type>
  void AdDual(const Vector3Type & P_,
              const Eigen::MatrixBase<Matrix3Type> & R_,
              Eigen::MatrixBase<Matrix6Type> & Adjoint)
  { AdDualAlgo<Vector3Type, Matrix3Type, Matrix6Type>::run(P_, R_, Adjoint); }

  template<typename Vector3Type, typename Matrix3Type, typename Matrix6Type>
  struct AdDualAlgo
  {
  public :
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inline static void run(const Vector3Type & P_,
                           const Eigen::MatrixBase<Matrix3Type> & R_,
                           Eigen::MatrixBase<Matrix6Type> & Adjoint)
    {
      typedef Eigen::Block<Matrix6Type,3,3> Block3;
      Block3 UL = Adjoint.template block<3,3>(0,0);  Block3 UR = Adjoint.template block<3,3>(0,3);
      Block3 LL = Adjoint.template block<3,3>(3,0);  Block3 LR = Adjoint.template block<3,3>(3,3);

      UL = R_;                    UR = Matrix3Type::Zero();
      LL.noalias() = Skew(P_)*R_;  LR = R_;
    }
  };


  //! Performing Ad_dual*D_M*Ad statically for Rx joint type
  //!------------------------------------------------------------------------------!//
  // Forward declaration
  template<typename Vector3Type, typename Matrix3Type, typename Matrix6TypeIn, typename Matrix6TypeOut> struct D_Mat6ProjRxAlgo;

  template<typename Vector3Type, typename Matrix3Type, typename Matrix6TypeIn, typename Matrix6TypeOut>
  void D_Mat6ProjRx(bool P_z,
                    int nS,
                    Eigen::MatrixBase<Vector3Type> & P_,
                    Eigen::MatrixBase<Matrix3Type> & R_,
                    const Eigen::MatrixBase<Matrix6TypeIn> & MatIn,
                    Eigen::MatrixBase<Matrix6TypeOut> & MatOut)
  { D_Mat6ProjRxAlgo<Vector3Type, Matrix3Type, Matrix6TypeIn, Matrix6TypeOut>::run(P_z, nS, P_, R_, MatIn, MatOut); }

  template<typename Vector3Type, typename Matrix3Type, typename Matrix6TypeIn, typename Matrix6TypeOut>
  struct D_Mat6ProjRxAlgo
  {
  public :
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inline static void run(bool P_z,
                           int nS,
                           Eigen::MatrixBase<Vector3Type> & P_,
                           Eigen::MatrixBase<Matrix3Type> & R_,
                           const Eigen::MatrixBase<Matrix6TypeIn> & MatIn,
                           Eigen::MatrixBase<Matrix6TypeOut> & MatOut)
    {
      typedef const Eigen::Block<Matrix6TypeIn,3,3> constBlock3;
      typedef Eigen::Block<Matrix6TypeOut,3,3> Block3;
      //! Linear
      Matrix3Type SkP, Dtmp1;

      typedef typename Matrix3Type::Scalar MyScalar;
      MyScalar sq, cq;
      sq = R_.coeffRef(2,1);  cq = R_.coeffRef(1,1);

      //! Linear
      if(P_z)  SkP = Skew(P_);

      Matrix6TypeIn & MatIn_ = const_cast<Matrix6TypeIn &>(MatIn.derived());

      for(int iter = 0; iter < nS*6; iter += 6){
          constBlock3 & Ai = MatIn_.template block<3,3>(iter,0);
          constBlock3 & Bi = MatIn_.template block<3,3>(iter,3);
          constBlock3 & Di = MatIn_.template block<3,3>(iter+3,3);

          Block3 Ao = MatOut.template block<3,3>(iter,0);
          Block3 Bo = MatOut.template block<3,3>(iter,3);
          Block3 Co = MatOut.template block<3,3>(iter+3,0);
          Block3 Do = MatOut.template block<3,3>(iter+3,3);

          Mat3ProjRx(sq, cq, Ai, Ao);

          Do.noalias() = R_*Bi; // tmp variable
          Co.noalias() = Do*R_.transpose();

          Mat3ProjRx(sq, cq, Di, Do);

          Bo = Co;
          if(P_z){
              //! Linear
              Bo.noalias() -= Ao*SkP;

              Dtmp1.noalias() = SkP*Bo;
              Do.noalias() += Dtmp1;
              Do.noalias() -= Co.transpose()*SkP;
            }
          Co = Bo.transpose();
        }
    }
  };


  //! Performing Ad_dual*D_M*Ad statically for Ry joint type
  //!------------------------------------------------------------------------------!//
  // Forward declaration
  template<typename Vector3Type, typename Matrix3Type, typename Matrix6TypeIn, typename Matrix6TypeOut> struct D_Mat6ProjRyAlgo;

  template<typename Vector3Type, typename Matrix3Type, typename Matrix6TypeIn, typename Matrix6TypeOut>
  void D_Mat6ProjRy(bool P_z,
                    int nS,
                    Eigen::MatrixBase<Vector3Type> & P_,
                    Eigen::MatrixBase<Matrix3Type> & R_,
                    const Eigen::MatrixBase<Matrix6TypeIn> & MatIn,
                    Eigen::MatrixBase<Matrix6TypeOut> & MatOut)
  { D_Mat6ProjRyAlgo<Vector3Type, Matrix3Type, Matrix6TypeIn, Matrix6TypeOut>::run(P_z, nS, P_, R_, MatIn, MatOut); }

  template<typename Vector3Type, typename Matrix3Type, typename Matrix6TypeIn, typename Matrix6TypeOut>
  struct D_Mat6ProjRyAlgo
  {
  public :
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inline static void run(bool P_z,
                           int nS,
                           Eigen::MatrixBase<Vector3Type> & P_,
                           Eigen::MatrixBase<Matrix3Type> & R_,
                           const Eigen::MatrixBase<Matrix6TypeIn> & MatIn,
                           Eigen::MatrixBase<Matrix6TypeOut> & MatOut)
    {
      typedef const Eigen::Block<Matrix6TypeIn,3,3> constBlock3;
      typedef Eigen::Block<Matrix6TypeOut,3,3> Block3;
      //! Linear
      Matrix3Type SkP, Dtmp1;

      typedef typename Matrix3Type::Scalar MyScalar;
      MyScalar sq, cq;
      sq = R_.coeffRef(0,2);  cq = R_.coeffRef(0,0);

      //! Linear
      if(P_z)  SkP = Skew(P_);

      Matrix6TypeIn & MatIn_ = const_cast<Matrix6TypeIn &>(MatIn.derived());

      for(int iter = 0; iter < nS*6; iter += 6){
          constBlock3 & Ai = MatIn_.template block<3,3>(iter,0);
          constBlock3 & Bi = MatIn_.template block<3,3>(iter,3);
          constBlock3 & Di = MatIn_.template block<3,3>(iter+3,3);

          Block3 Ao = MatOut.template block<3,3>(iter,0);
          Block3 Bo = MatOut.template block<3,3>(iter,3);
          Block3 Co = MatOut.template block<3,3>(iter+3,0);
          Block3 Do = MatOut.template block<3,3>(iter+3,3);

          Mat3ProjRy(sq, cq, Ai, Ao);

          Do.noalias() = R_*Bi; // tmp variable
          Co.noalias() = Do*R_.transpose();

          Mat3ProjRy(sq, cq, Di, Do);

          Bo = Co;
          if(P_z){
              //! Linear
              Bo.noalias() -= Ao*SkP;

              Dtmp1.noalias() = SkP*Bo;
              Do.noalias() += Dtmp1;
              Do.noalias() -= Co.transpose()*SkP;
            }
          Co = Bo.transpose();
        }
    }
  };


  //! Performing Ad_dual*D_M*Ad statically for Rz joint type
  //!------------------------------------------------------------------------------!//
  // Forward declaration
  template<typename Vector3Type, typename Matrix3Type, typename Matrix6TypeIn, typename Matrix6TypeOut> struct D_Mat6ProjRzAlgo;

  template<typename Vector3Type, typename Matrix3Type, typename Matrix6TypeIn, typename Matrix6TypeOut>
  void D_Mat6ProjRz(bool P_z,
                    int nS,
                    Eigen::MatrixBase<Vector3Type> & P_,
                    Eigen::MatrixBase<Matrix3Type> & R_,
                    const Eigen::MatrixBase<Matrix6TypeIn> & MatIn,
                    Eigen::MatrixBase<Matrix6TypeOut> & MatOut)
  { D_Mat6ProjRzAlgo<Vector3Type, Matrix3Type, Matrix6TypeIn, Matrix6TypeOut>::run(P_z, nS, P_, R_, MatIn, MatOut); }

  template<typename Vector3Type, typename Matrix3Type, typename Matrix6TypeIn, typename Matrix6TypeOut>
  struct D_Mat6ProjRzAlgo
  {
  public :
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inline static void run(bool P_z,
                           int nS,
                           Eigen::MatrixBase<Vector3Type> & P_,
                           Eigen::MatrixBase<Matrix3Type> & R_,
                           const Eigen::MatrixBase<Matrix6TypeIn> & MatIn,
                           Eigen::MatrixBase<Matrix6TypeOut> & MatOut)
    {
      typedef const Eigen::Block<Matrix6TypeIn,3,3> constBlock3;
      typedef Eigen::Block<Matrix6TypeOut,3,3> Block3;
      //! Linear
      Matrix3Type SkP, Dtmp1;

      typedef typename Matrix3Type::Scalar MyScalar;
      MyScalar sq, cq;
      sq = R_.coeffRef(1,0);  cq = R_.coeffRef(0,0);

      //! Linear
      if(P_z)  SkP = Skew(P_);

      Matrix6TypeIn & MatIn_ = const_cast<Matrix6TypeIn &>(MatIn.derived());

      for(int iter = 0; iter < nS*6; iter += 6){
          constBlock3 & Ai = MatIn_.template block<3,3>(iter,0);
          constBlock3 & Bi = MatIn_.template block<3,3>(iter,3);
          constBlock3 & Di = MatIn_.template block<3,3>(iter+3,3);

          Block3 Ao = MatOut.template block<3,3>(iter,0);
          Block3 Bo = MatOut.template block<3,3>(iter,3);
          Block3 Co = MatOut.template block<3,3>(iter+3,0);
          Block3 Do = MatOut.template block<3,3>(iter+3,3);

          Mat3ProjRz(sq, cq, Ai, Ao);

          Do.noalias() = R_*Bi; // tmp variable
          Co.noalias() = Do*R_.transpose();

          Mat3ProjRz(sq, cq, Di, Do);

          Bo = Co;
          if(P_z){
              //! Linear
              Bo.noalias() -= Ao*SkP;

              Dtmp1.noalias() = SkP*Bo;
              Do.noalias() += Dtmp1;
              Do.noalias() -= Co.transpose()*SkP;
            }
          Co = Bo.transpose();
        }
    }
  };


  //! Performing Ad_dual*D_P statically
  //!------------------------------------------------------------------------------!//
  // Forward declaration
  template<typename Vector3Type, typename Matrix3Type, typename D_Vector6TypeIn, typename D_Vector6TypeOut> struct D_Vec6ProjAlgo;

  template<typename Vector3Type, typename Matrix3Type, typename D_Vector6TypeIn, typename D_Vector6TypeOut>
  void D_Vec6Proj(bool P_z,
                  Eigen::MatrixBase<Vector3Type> & P_,
                  Eigen::MatrixBase<Matrix3Type> & R_,
                  const Eigen::MatrixBase<D_Vector6TypeIn> & VecIn,
                  Eigen::MatrixBase<D_Vector6TypeOut> & VecOut)
  { D_Vec6ProjAlgo<Vector3Type, Matrix3Type, D_Vector6TypeIn, D_Vector6TypeOut>::run(P_z, P_, R_, VecIn, VecOut); }

  template<typename Vector3Type, typename Matrix3Type, typename D_Vector6TypeIn, typename D_Vector6TypeOut>
  struct D_Vec6ProjAlgo
  {
  public :
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inline static void run(bool P_z,
                           Eigen::MatrixBase<Vector3Type> & P_,
                           Eigen::MatrixBase<Matrix3Type> & R_,
                           const Eigen::MatrixBase<D_Vector6TypeIn> & VecIn,
                           Eigen::MatrixBase<D_Vector6TypeOut> & VecOut)
    {
      Matrix3Type skP;
      VecOut.template topRows<3>() = R_*VecIn.template topRows<3>();
      VecOut.template bottomRows<3>() = R_*VecIn.template bottomRows<3>();

      if(P_z) {
          skP = Skew(P_);
          VecOut.template bottomRows<3>() += skP*VecIn.template topRows<3>();
        }
    }
  };

}


//! -ad(Sx) effect
//!------------------------------------------------------------------------------!//
template <class ArgType> class ad_Sx;

namespace Eigen {
  namespace internal {
    template <class ArgType>
    struct traits<ad_Sx<ArgType> >
    {
      typedef Eigen::Dense StorageKind;
      typedef Eigen::MatrixXpr XprKind;
      typedef typename ArgType::StorageIndex StorageIndex;
      typedef typename ArgType::Scalar Scalar;
      enum {
        Flags = Eigen::ColMajor,
        RowsAtCompileTime = 6,//ArgType::RowsAtCompileTime,
        ColsAtCompileTime = ArgType::ColsAtCompileTime,
        MaxRowsAtCompileTime = 6,//ArgType::MaxRowsAtCompileTime,
        MaxColsAtCompileTime = ArgType::MaxColsAtCompileTime
      };
    };
  }
}

template <class ArgType>
class ad_Sx : public Eigen::MatrixBase<ad_Sx<ArgType> >
{
public:
  ad_Sx(const ArgType& arg)
    : m_arg(arg)
  {
    EIGEN_STATIC_ASSERT(ArgType::RowsAtCompileTime == 6,
                        YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_A_MATRIX);
  }

  typedef typename Eigen::internal::ref_selector<ad_Sx>::type Nested;

  typedef Eigen::Index Index;
  Index rows() const { return m_arg.rows(); }
  Index cols() const { return m_arg.cols(); }

  typedef typename Eigen::internal::ref_selector<ArgType>::type ArgTypeNested;
  ArgTypeNested m_arg;
};

namespace Eigen {
  namespace internal {
    template<typename ArgType>
    struct evaluator<ad_Sx<ArgType> >
        : evaluator_base<ad_Sx<ArgType> >
    {
      typedef ad_Sx<ArgType> XprType;
      typedef typename nested_eval<ArgType, XprType::ColsAtCompileTime>::type ArgTypeNested;
      typedef typename remove_all<ArgTypeNested>::type ArgTypeNestedCleaned;
      typedef typename XprType::CoeffReturnType CoeffReturnType;

      enum {
        CoeffReadCost = evaluator<ArgTypeNestedCleaned>::CoeffReadCost,
        Flags = Eigen::ColMajor
      };

      evaluator(const XprType& xpr)
        : m_argImpl(xpr.m_arg)
      { }

      CoeffReturnType coeff(Index row, Index col) const
      {
        CoeffReturnType out = 0;

        if(row==1)       out =  m_argImpl.coeff(2,col);
        else if(row==2)  out = -m_argImpl.coeff(1,col);
        else if(row==4)  out =  m_argImpl.coeff(5,col);
        else if(row==5)  out = -m_argImpl.coeff(4,col);

        return out;
      }

      evaluator<ArgTypeNestedCleaned> m_argImpl;
      const Index m_rows = 6;
    };
  }
}

template <class ArgType>
ad_Sx<ArgType> adSx(const Eigen::MatrixBase<ArgType>& arg)
{
  return ad_Sx<ArgType>(arg.derived());
}


//! -ad(Sy) effect
//!------------------------------------------------------------------------------!//
template <class ArgType> class ad_Sy;

namespace Eigen {
  namespace internal {
    template <class ArgType>
    struct traits<ad_Sy<ArgType> >
    {
      typedef Eigen::Dense StorageKind;
      typedef Eigen::MatrixXpr XprKind;
      typedef typename ArgType::StorageIndex StorageIndex;
      typedef typename ArgType::Scalar Scalar;
      enum {
        Flags = Eigen::ColMajor,
        RowsAtCompileTime = 6,//ArgType::RowsAtCompileTime,
        ColsAtCompileTime = ArgType::ColsAtCompileTime,
        MaxRowsAtCompileTime = 6,//ArgType::MaxRowsAtCompileTime,
        MaxColsAtCompileTime = ArgType::MaxColsAtCompileTime
      };
    };
  }
}

template <class ArgType>
class ad_Sy : public Eigen::MatrixBase<ad_Sy<ArgType> >
{
public:
  ad_Sy(const ArgType& arg)
    : m_arg(arg)
  {
    EIGEN_STATIC_ASSERT(ArgType::RowsAtCompileTime == 6,
                        YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_A_MATRIX);
  }

  typedef typename Eigen::internal::ref_selector<ad_Sy>::type Nested;

  typedef Eigen::Index Index;
  Index rows() const { return m_arg.rows(); }
  Index cols() const { return m_arg.cols(); }

  typedef typename Eigen::internal::ref_selector<ArgType>::type ArgTypeNested;
  ArgTypeNested m_arg;
};

namespace Eigen {
  namespace internal {
    template<typename ArgType>
    struct evaluator<ad_Sy<ArgType> >
        : evaluator_base<ad_Sy<ArgType> >
    {
      typedef ad_Sy<ArgType> XprType;
      typedef typename nested_eval<ArgType, XprType::ColsAtCompileTime>::type ArgTypeNested;
      typedef typename remove_all<ArgTypeNested>::type ArgTypeNestedCleaned;
      typedef typename XprType::CoeffReturnType CoeffReturnType;

      enum {
        CoeffReadCost = evaluator<ArgTypeNestedCleaned>::CoeffReadCost,
        Flags = Eigen::ColMajor
      };

      evaluator(const XprType& xpr)
        : m_argImpl(xpr.m_arg)
      { }

      CoeffReturnType coeff(Index row, Index col) const
      {
        CoeffReturnType out = 0;

        if(row==0)       out = -m_argImpl.coeff(2,col);
        else if(row==2)  out =  m_argImpl.coeff(0,col);
        else if(row==3)  out = -m_argImpl.coeff(5,col);
        else if(row==5)  out =  m_argImpl.coeff(3,col);

        return out;
      }

      evaluator<ArgTypeNestedCleaned> m_argImpl;
      const Index m_rows = 6;
    };
  }
}

template <class ArgType>
ad_Sy<ArgType> adSy(const Eigen::MatrixBase<ArgType>& arg)
{
  return ad_Sy<ArgType>(arg.derived());
}


//! -ad(Sz) effect
//!------------------------------------------------------------------------------!//
template <class ArgType> class ad_Sz;

namespace Eigen {
  namespace internal {
    template <class ArgType>
    struct traits<ad_Sz<ArgType> >
    {
      typedef Eigen::Dense StorageKind;
      typedef Eigen::MatrixXpr XprKind;
      typedef typename ArgType::StorageIndex StorageIndex;
      typedef typename ArgType::Scalar Scalar;
      enum {
        Flags = Eigen::ColMajor,
        RowsAtCompileTime = 6,//ArgType::RowsAtCompileTime,
        ColsAtCompileTime = ArgType::ColsAtCompileTime,
        MaxRowsAtCompileTime = 6,//ArgType::MaxRowsAtCompileTime,
        MaxColsAtCompileTime = ArgType::MaxColsAtCompileTime
      };
    };
  }
}

template <class ArgType>
class ad_Sz : public Eigen::MatrixBase<ad_Sz<ArgType> >
{
public:
  ad_Sz(const ArgType& arg)
    : m_arg(arg)
  {
    EIGEN_STATIC_ASSERT(ArgType::RowsAtCompileTime == 6,
                        YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_A_MATRIX);
  }

  typedef typename Eigen::internal::ref_selector<ad_Sz>::type Nested;

  typedef Eigen::Index Index;
  Index rows() const { return m_arg.rows(); }
  Index cols() const { return m_arg.cols(); }

  typedef typename Eigen::internal::ref_selector<ArgType>::type ArgTypeNested;
  ArgTypeNested m_arg;
};

namespace Eigen {
  namespace internal {
    template<typename ArgType>
    struct evaluator<ad_Sz<ArgType> >
        : evaluator_base<ad_Sz<ArgType> >
    {
      typedef ad_Sz<ArgType> XprType;
      typedef typename nested_eval<ArgType, XprType::ColsAtCompileTime>::type ArgTypeNested;
      typedef typename remove_all<ArgTypeNested>::type ArgTypeNestedCleaned;
      typedef typename XprType::CoeffReturnType CoeffReturnType;

      enum {
        CoeffReadCost = evaluator<ArgTypeNestedCleaned>::CoeffReadCost,
        Flags = Eigen::ColMajor
      };

      evaluator(const XprType& xpr)
        : m_argImpl(xpr.m_arg)
      { }

      CoeffReturnType coeff(Index row, Index col) const
      {
        CoeffReturnType out = 0;

        if(row==0)  out = m_argImpl.coeff(1,col);
        else if(row==1)  out = -m_argImpl.coeff(0,col);
        else if(row==3)  out =  m_argImpl.coeff(4,col);
        else if(row==4)  out = -m_argImpl.coeff(3,col);

        return out;
      }

      evaluator<ArgTypeNestedCleaned> m_argImpl;
      const Index m_rows = 6;
    };
  }
}

template <class ArgType>
ad_Sz<ArgType> adSz(const Eigen::MatrixBase<ArgType>& arg)
{
  return ad_Sz<ArgType>(arg.derived());
}


#endif // GEOMBD_FUNCTORS_HXX
