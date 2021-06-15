#ifndef GEOMBD_DIFFERENTIATION_JOINT_TYPE_REVOLUTE_Z_HXX
#define GEOMBD_DIFFERENTIATION_JOINT_TYPE_REVOLUTE_Z_HXX

#define EIGEN_NO_DEBUG
#define EIGEN_MPL2_ONLY
#define EIGEN_UNROLLING_LIMIT 30

#include "Eigen/Core"
#include "Eigen/Geometry"

//#include <iostream>
//#include <memory>

namespace geo{


  //! Derived -> Differentiation of Joint Type Revolute On Z Axis
  //!------------------------------------------------------------------------------!//
  struct D_JointTypeRz : public D_CRTPInterface<D_JointTypeRz>
  {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! FK Forward Declaration.
    template<typename ScalarType, typename Matrix3Type>
    inline static void
    runD_FK(const ScalarType & qi,
            typename Eigen::MatrixBase<Matrix3Type> & R);


    //! TwCbPb at Root Forward Declaration.
    template<typename ScalarType, typename Vector6Type, typename Matrix6Type, typename D_Vector6Type>
    inline static void
    runD_TCP_root(const ScalarType & vi,
                  Eigen::MatrixBase<Vector6Type> & S_i,
                  Eigen::MatrixBase<Vector6Type> & p_,
                  const Eigen::MatrixBase<Matrix6Type> & M_,
                  Eigen::MatrixBase<D_Vector6Type> & D_dq_p_);


    //! TwCbPb Forward Declaration.
    template<typename ScalarType, typename Matrix3Type, typename Matrix6Type,
             typename Vector3Type, typename Vector6Type, typename D_Vector6Type>
    inline static void
    runD_TwCbPb(bool zeroFlag,
                const ScalarType & vi,
                const Eigen::MatrixBase<Matrix3Type> & R_,
                const Eigen::MatrixBase<Vector3Type> & P_,
                const Eigen::MatrixBase<Vector6Type> & S_l,
                const Eigen::MatrixBase<Matrix6Type> & M_,
                Eigen::MatrixBase<Vector6Type> & S_i,
                Eigen::MatrixBase<Vector6Type> & c_,
                Eigen::MatrixBase<Vector6Type> & p_,
                Eigen::MatrixBase<D_Vector6Type> & D_q_V_,
                Eigen::MatrixBase<D_Vector6Type> & D_dq_V_,
                const Eigen::MatrixBase<D_Vector6Type> & D_q_Vj_,
                const Eigen::MatrixBase<D_Vector6Type> & D_dq_Vj_,
                Eigen::MatrixBase<D_Vector6Type> & D_q_c_,
                Eigen::MatrixBase<D_Vector6Type> & D_dq_c_,
                Eigen::MatrixBase<D_Vector6Type> & D_q_p_,
                Eigen::MatrixBase<D_Vector6Type> & D_dq_p_);


    //! Forward declaration of inertial expressions at leaf.
    template<typename ScalarType, typename Vector3Type, typename Matrix3Type, typename Vector6Type,
             typename Matrix6Type, typename RowVectorXType, typename D_Vector6Type, typename D_Matrix6Type>
    inline static void
    runD_InertiaLeaf(ScalarType & u,
                     ScalarType & iD,
                     const ScalarType tau,
                     Eigen::MatrixBase<Vector6Type> & U_,
                     Eigen::MatrixBase<Vector6Type> & c_,
                     Eigen::MatrixBase<Vector6Type> & P_a_,
                     Eigen::MatrixBase<Matrix6Type> & M_a_,
                     Eigen::MatrixBase<Vector6Type> & P_A_,
                     Eigen::MatrixBase<Matrix6Type> & M_A_,
                     bool P_z,
                     Eigen::MatrixBase<Vector3Type> & P_,
                     Eigen::MatrixBase<Matrix3Type> & R_,
                     Eigen::MatrixBase<Vector6Type> & P_Aj_,
                     Eigen::MatrixBase<Matrix6Type> & M_Aj_,
                     Eigen::MatrixBase<D_Matrix6Type> & D_M_Aj_,
                     Eigen::MatrixBase<RowVectorXType> & D_q_u_,
                     Eigen::MatrixBase<RowVectorXType> & D_dq_u_,
                     Eigen::MatrixBase<D_Vector6Type> & D_q_p_,
                     Eigen::MatrixBase<D_Vector6Type> & D_dq_p_,
                     Eigen::MatrixBase<D_Vector6Type> & D_q_Pa_,
                     Eigen::MatrixBase<D_Vector6Type> & D_dq_Pa_,
                     Eigen::MatrixBase<D_Vector6Type> & D_q_PA_,
                     Eigen::MatrixBase<D_Vector6Type> & D_dq_PA_,
                     Eigen::MatrixBase<D_Vector6Type> & D_q_PAj_,
                     Eigen::MatrixBase<D_Vector6Type> & D_dq_PAj_,
                     Eigen::MatrixBase<D_Vector6Type> & D_q_c_,
                     Eigen::MatrixBase<D_Vector6Type> & D_dq_c_);


    //! Forward declaration of inertial expressions.
    template<typename ScalarType, typename Vector6iType, typename Vector3Type, typename Matrix3Type, typename Vector6Type,
             typename Matrix6Type, typename VectorXType, typename RowVectorXType, typename D_Vector6Type, typename D_Matrix6Type>
    inline static void
    runD_Inertial(bool rootFlag,
                  Eigen::MatrixBase<Vector6iType> & iVec,
                  ScalarType & u,
                  ScalarType & iD,
                  const ScalarType tau,
                  Eigen::MatrixBase<Vector6Type> & U_,
                  Eigen::MatrixBase<Vector6Type> & c_,
                  Eigen::MatrixBase<Vector6Type> & P_a_,
                  Eigen::MatrixBase<Matrix6Type> & M_a_,
                  Eigen::MatrixBase<Vector6Type> & P_A_,
                  Eigen::MatrixBase<Matrix6Type> & M_A_,
                  bool P_z,
                  Eigen::MatrixBase<Vector3Type> & P_,
                  Eigen::MatrixBase<Matrix3Type> & R_,
                  Eigen::MatrixBase<Vector6Type> & P_Aj_,
                  Eigen::MatrixBase<Matrix6Type> & M_Aj_,
                  Eigen::MatrixBase<D_Vector6Type> & D_U_h_,
                  Eigen::MatrixBase<VectorXType> & D_U_v_,
                  Eigen::MatrixBase<RowVectorXType> & D_invD_,
                  Eigen::MatrixBase<D_Matrix6Type> & D_M_A_,
                  Eigen::MatrixBase<D_Matrix6Type> & D_M_Aj_,
                  Eigen::MatrixBase<RowVectorXType> & D_q_u_,
                  Eigen::MatrixBase<RowVectorXType> & D_dq_u_,
                  Eigen::MatrixBase<D_Vector6Type> & D_q_Pa_,
                  Eigen::MatrixBase<D_Vector6Type> & D_dq_Pa_,
                  Eigen::MatrixBase<D_Vector6Type> & D_q_PA_,
                  Eigen::MatrixBase<D_Vector6Type> & D_dq_PA_,
                  Eigen::MatrixBase<D_Vector6Type> & D_q_PAj_,
                  Eigen::MatrixBase<D_Vector6Type> & D_dq_PAj_,
                  Eigen::MatrixBase<D_Vector6Type> & D_q_c_,
                  Eigen::MatrixBase<D_Vector6Type> & D_dq_c_);


    //! Forward Declaration for Acceleration at Root.
    template<typename ScalarType, typename Vector3Type, typename Matrix3Type, typename Vector6Type,
             typename D_Vector6Type, typename RowVectorXType, typename MatrixXType>
    inline static void
    runD_AccelRoot(ScalarType u,
                   ScalarType iD,
                   ScalarType* ddq,
                   const Eigen::MatrixBase<Vector3Type> & P_r,
                   const Eigen::MatrixBase<Matrix3Type> & R_r,
                   const Eigen::MatrixBase<Vector6Type> & U_r,
                   Eigen::MatrixBase<Vector6Type> & Acc_i_r,
                   Eigen::MatrixBase<D_Vector6Type> & D_U_h_,
                   Eigen::MatrixBase<RowVectorXType> & D_invD_,
                   Eigen::MatrixBase<RowVectorXType> & D_q_u_,
                   Eigen::MatrixBase<RowVectorXType> & D_dq_u_,
                   Eigen::MatrixBase<D_Vector6Type> & D_q_A_,
                   Eigen::MatrixBase<D_Vector6Type> & D_dq_A_,
                   Eigen::MatrixBase<MatrixXType> & D_ddq_);


    //! Forward Declaration for Acceleration.
    template<typename ScalarType, typename IndexType, typename Vector3Type, typename Matrix3Type,
             typename Vector6Type, typename D_Vector6Type, typename RowVectorXType, typename MatrixXType>
    inline static void
    runD_Accel(IndexType ID,
               bool zeroFlag,
               ScalarType u,
               ScalarType iD,
               ScalarType* ddq,
               const Eigen::MatrixBase<Vector3Type> & P_,
               const Eigen::MatrixBase<Matrix3Type> & R_,
               const Eigen::MatrixBase<Vector6Type> & c_,
               const Eigen::MatrixBase<Vector6Type> & U_,
               Eigen::MatrixBase<Vector6Type> & A_,
               const Eigen::MatrixBase<Vector6Type> & Aj_,
               bool isLeaf,
               std::vector< IndexType >* Pre_,
               std::vector< IndexType >* Suc_,
               std::vector< IndexType >* PreSuc_,
               Eigen::MatrixBase<D_Vector6Type> & D_U_h_,
               Eigen::MatrixBase<RowVectorXType> & D_invD_,
               Eigen::MatrixBase<MatrixXType> & D_ddq_,
               Eigen::MatrixBase<RowVectorXType> & D_q_u_,
               Eigen::MatrixBase<RowVectorXType> & D_dq_u_,
               Eigen::MatrixBase<D_Vector6Type> & D_q_c_,
               Eigen::MatrixBase<D_Vector6Type> & D_dq_c_,
               Eigen::MatrixBase<D_Vector6Type> & D_q_A_,
               Eigen::MatrixBase<D_Vector6Type> & D_dq_A_,
               Eigen::MatrixBase<D_Vector6Type> & D_q_Aj_,
               Eigen::MatrixBase<D_Vector6Type> & D_dq_Aj_);


  };


  //! Member-function declaration
  //!------------------------------------------------------------------------------!//
  //!------------------------------------------------------------------------------!//
  //!------------------------------------------------------------------------------!//


  //! Forward Kinematics Declaration
  //! EIGEN_ASM_COMMENT("MyBegin");
  //!------------------------------------------------------------------------------!//
  template<typename ScalarType, typename Matrix3Type>
  inline void
  D_JointTypeRz::runD_FK(const ScalarType & qi,
                         typename Eigen::MatrixBase<Matrix3Type> & R) {
    static ScalarType sqi, cqi;
    SINCOS<ScalarType>(qi, &sqi, &cqi);

    R.coeffRef(0,0) = cqi;  R.coeffRef(0,1) = -sqi;
    R.coeffRef(1,0) = sqi;  R.coeffRef(1,1) =  cqi;
  }


  //! TwCbPb at Root Declaration.
  template<typename ScalarType, typename Vector6Type, typename Matrix6Type, typename D_Vector6Type>
  inline void
  D_JointTypeRz::runD_TCP_root(const ScalarType & vi,
                               Eigen::MatrixBase<Vector6Type> & S_i,
                               Eigen::MatrixBase<Vector6Type> & p_,
                               const Eigen::MatrixBase<Matrix6Type> & M_,
                               Eigen::MatrixBase<D_Vector6Type> & D_dq_p_) {

    S_i.template tail<1>()(0) = vi;

    Matrix6Type & M = const_cast<Matrix6Type &>(M_.derived());

    ScalarType vi2 = vi * vi;

    p_.coeffRef(0) = -M.coeff(1,5)*vi2;
    p_.coeffRef(1) =  M.coeff(0,5)*vi2;
    p_.coeffRef(3) = -M.coeff(4,5)*vi2;
    p_.coeffRef(4) =  M.coeff(3,5)*vi2;

    //! P bias differentiation
    //!------------------------------------------------------------------------------!//
    typedef Eigen::Block<Matrix6Type,3,3> Block3;
    Matrix6Type D_p_;

    Block3 UR_ = D_p_.template block<3,3>(0,3);
    Block3 LL_ = D_p_.template block<3,3>(3,0);  Block3 LR_ = D_p_.template block<3,3>(3,3);

    //! adDual
//    D_p_.setZero();
    D_p_.template row(0) = -M.template row(1);  // adDual(Sz) effect
    D_p_.template row(1) =  M.template row(0);
    D_p_.template row(3) = -M.template row(4);
    D_p_.template row(4) =  M.template row(3);

    Vector6Type mu;
    mu.noalias() = M.template rightCols<1>()*vi;

    //! adBar
    UR_.noalias() -= Skew(mu.segment(0,3));
    LL_.noalias() -= Skew(mu.segment(0,3));  LR_.noalias() -= Skew(mu.segment(3,3));

    D_dq_p_.noalias() = D_p_.template rightCols<1>();
  }


  //! Apply the Ad operator on se(3) element to get the Twist and bias
  //! TODO: Twist is only used here, thus consider to store it only when multiple children
  //!------------------------------------------------------------------------------!//
  template<typename ScalarType, typename Matrix3Type, typename Matrix6Type,
           typename Vector3Type, typename Vector6Type, typename D_Vector6Type>
  inline void
  D_JointTypeRz::runD_TwCbPb(bool zeroFlag,
                             const ScalarType & vi,
                             const Eigen::MatrixBase<Matrix3Type> & R_,
                             const Eigen::MatrixBase<Vector3Type> & P_,
                             const Eigen::MatrixBase<Vector6Type> & S_l,
                             const Eigen::MatrixBase<Matrix6Type> & M_,
                             Eigen::MatrixBase<Vector6Type> & S_i,
                             Eigen::MatrixBase<Vector6Type> & c_,
                             Eigen::MatrixBase<Vector6Type> & p_,
                             Eigen::MatrixBase<D_Vector6Type> & D_q_V_,
                             Eigen::MatrixBase<D_Vector6Type> & D_dq_V_,
                             const Eigen::MatrixBase<D_Vector6Type> & D_q_Vj_,
                             const Eigen::MatrixBase<D_Vector6Type> & D_dq_Vj_,
                             Eigen::MatrixBase<D_Vector6Type> & D_q_c_,
                             Eigen::MatrixBase<D_Vector6Type> & D_dq_c_,
                             Eigen::MatrixBase<D_Vector6Type> & D_q_p_,
                             Eigen::MatrixBase<D_Vector6Type> & D_dq_p_) {

    EIGEN_STATIC_ASSERT(Vector6Type::ColsAtCompileTime == 1,
                        YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_A_MATRIX);

    //! Twist
    //!------------------------------------------------------------------------------!//
    typedef const Eigen::Block<Vector6Type,3,1> constSegment3;
    typedef Eigen::Block<Vector6Type,3,1> Segment3;

    Vector6Type & S_ = const_cast<Vector6Type &>(S_l.derived());
    constSegment3 & Sup   = S_.template segment<3>(0);
    constSegment3 & Sdown = S_.template segment<3>(3);

    Segment3 Up_   = S_i.template segment<3>(0);
    Segment3 Down_ = S_i.template segment<3>(3);

    short int icol = D_q_V_.cols();
    if(zeroFlag) {
        typename Matrix3Type::PlainObject SkP = Skew(P_);
        Down_.noalias() = Sup - SkP*Sdown;  //tmp
        Up_.noalias()   = R_.transpose()*Down_;
        //!-------------------------------------------------------
        D_q_V_.block(3,0,3,icol-1).noalias() = D_q_Vj_.template topRows<3>() - SkP*D_q_Vj_.template bottomRows<3>();  //tmp
        D_q_V_.block(0,0,3,icol-1).noalias() = R_.transpose()*D_q_V_.block(3,0,3,icol-1);

        D_dq_V_.block(3,0,3,icol-1).noalias() = D_dq_Vj_.template topRows<3>() - SkP*D_dq_Vj_.template bottomRows<3>();  //tmp
        D_dq_V_.block(0,0,3,icol-1).noalias() = R_.transpose()*D_dq_V_.block(3,0,3,icol-1);
      } else {
        Up_.noalias() = R_.transpose()*Sup;
        //!-------------------------------------------------------
        D_q_V_.block(0,0,3,icol-1).noalias()  = R_.transpose()*D_q_Vj_.template topRows<3>();
        D_dq_V_.block(0,0,3,icol-1).noalias() = R_.transpose()*D_dq_Vj_.template topRows<3>();
      }

    Down_.noalias() = R_.transpose()*Sdown;
    Down_.coeffRef(2) += vi;
    //!-------------------------------------------------------
    D_q_V_.block(3,0,3,icol-1).noalias()  = R_.transpose()*D_q_Vj_.template bottomRows<3>();
    D_dq_V_.block(3,0,3,icol-1).noalias() = R_.transpose()*D_dq_Vj_.template bottomRows<3>();

    //! C bias
    //!------------------------------------------------------------------------------!//
    c_ = adSz(S_i);

    //! Twist differentiation -> taking advantage of c bias computation.
    //!------------------------------------------------------------------------------!//
    D_q_V_.template rightCols<1>() = c_;
    //!------------------------------------------------------------------------------!//

    c_ *= vi;

    //! C bias differentiation -> taking advantage of twist computation.
    //!------------------------------------------------------------------------------!//
    D_Vector6Type D_c_aux;
    D_c_aux.noalias() = D_q_V_*vi;
    D_q_c_ = adSz(D_c_aux);

    D_c_aux.noalias() = D_dq_V_*vi;
    D_c_aux.template rightCols<1>() += S_i;
    D_dq_c_ = adSz(D_c_aux);
    //!------------------------------------------------------------------------------!//

    //! P bias
    //!------------------------------------------------------------------------------!//
    typedef const Eigen::Block<Matrix6Type,3,3> constBlock3;
    typedef Eigen::Block<Matrix6Type,3,3> Block3;

    Matrix6Type & M = const_cast<Matrix6Type &>(M_.derived());
    constBlock3 & UL = M.template block<3,3>(0,0);
    constBlock3 & UR = M.template block<3,3>(0,3);
    constBlock3 & LL = M.template block<3,3>(3,0);
    constBlock3 & LR = M.template block<3,3>(3,3);

    typename Vector3Type::PlainObject muUp, muDo;
    //    typedef typename GEOMBD_EIGEN_PLAIN_TYPE(Vector3Type) PlainVector3Type;
    //    PlainVector3Type muUp, muDo;

    Segment3 p_Up = p_.template segment<3>(0);
    Segment3 p_Do = p_.template segment<3>(3);

    //! M is symmetric and its shape is known
    muUp.noalias() = UL*Up_;  muUp.noalias() += UR*Down_;
    muDo.noalias() = LL*Up_;  muDo.noalias() += LR*Down_;

    //! TODO: Cross is very expensive, consider enable staticness with traits expressions
    p_Up.noalias() = Down_.cross(muUp);
    p_Do.noalias() = Up_.cross(muUp);
    p_Do.noalias() += Down_.cross(muDo);

    //! P bias differentiation
    //!------------------------------------------------------------------------------!//
    Matrix6Type D_p_, D_aux;

    Block3 UL_ = D_p_.template block<3,3>(0,0);  Block3 UR_ = D_p_.template block<3,3>(0,3);
    Block3 LL_ = D_p_.template block<3,3>(3,0);  Block3 LR_ = D_p_.template block<3,3>(3,3);

    //! adDual
    UL_.noalias() = Skew(Down_);  UR_.setZero();
    LL_.noalias() = Skew(Up_);    LR_.noalias() = UL_;

    D_aux.noalias() = D_p_*M.derived();  // tmp

    D_p_.noalias() = D_aux;

    //! adBar
    UR_.noalias() -= Skew(muUp);
    LL_.noalias() -= Skew(muUp);  LR_.noalias() -= Skew(muDo);

    D_q_p_.noalias()  = D_p_*D_q_V_.derived();
    D_dq_p_.noalias() = D_p_*D_dq_V_.derived();
  }


  //! Implementation for inertial expressions at leaf.
  //! ToDo: Currently a leaf's parent isn't endowed to be branched.
  template<typename ScalarType, typename Vector3Type, typename Matrix3Type, typename Vector6Type,
           typename Matrix6Type, typename RowVectorXType, typename D_Vector6Type, typename D_Matrix6Type>
  inline void
  D_JointTypeRz::runD_InertiaLeaf(ScalarType & u,
                                  ScalarType & iD,
                                  const ScalarType tau,
                                  Eigen::MatrixBase<Vector6Type> & U_,
                                  Eigen::MatrixBase<Vector6Type> & c_,
                                  Eigen::MatrixBase<Vector6Type> & P_a_,
                                  Eigen::MatrixBase<Matrix6Type> & M_a_,
                                  Eigen::MatrixBase<Vector6Type> & P_A_,
                                  Eigen::MatrixBase<Matrix6Type> & M_A_,
                                  bool P_z,
                                  Eigen::MatrixBase<Vector3Type> & P_,
                                  Eigen::MatrixBase<Matrix3Type> & R_,
                                  Eigen::MatrixBase<Vector6Type> & P_Aj_,
                                  Eigen::MatrixBase<Matrix6Type> & M_Aj_,
                                  Eigen::MatrixBase<D_Matrix6Type> & D_M_Aj_,
                                  Eigen::MatrixBase<RowVectorXType> & D_q_u_,
                                  Eigen::MatrixBase<RowVectorXType> & D_dq_u_,
                                  Eigen::MatrixBase<D_Vector6Type> & D_q_p_,
                                  Eigen::MatrixBase<D_Vector6Type> & D_dq_p_,
                                  Eigen::MatrixBase<D_Vector6Type> & D_q_Pa_,
                                  Eigen::MatrixBase<D_Vector6Type> & D_dq_Pa_,
                                  Eigen::MatrixBase<D_Vector6Type> & D_q_PA_,
                                  Eigen::MatrixBase<D_Vector6Type> & D_dq_PA_,
                                  Eigen::MatrixBase<D_Vector6Type> & D_q_PAj_,
                                  Eigen::MatrixBase<D_Vector6Type> & D_dq_PAj_,
                                  Eigen::MatrixBase<D_Vector6Type> & D_q_c_,
                                  Eigen::MatrixBase<D_Vector6Type> & D_dq_c_) {
    //! Solve U, inverse of D and u.
    //!------------------------------------------------------------------------------!//
    U_ = M_A_.template rightCols<1>();  iD = 1 / U_.coeff(5);  u = tau - P_A_.coeff(5);


    //! Solve partial derivative of u as D_q_u and D_dq_u.
    //!------------------------------------------------------------------------------!//
    D_q_u_ = -D_q_PA_.template bottomRows<1>();  D_dq_u_ = -D_dq_PA_.template bottomRows<1>();


    //! Prepare inertial expressions.
    //!------------------------------------------------------------------------------!//
    //! Ma
    Vector6Type UD_ = U_*iD;
    M_a_.noalias() = M_A_ - UD_*U_.transpose();

    //! Pa
    P_a_.noalias() = P_A_ + M_a_*c_ + u*UD_;


    //! Inertial back projection.
    //!------------------------------------------------------------------------------!//
    //    typename GEOMBD_EIGEN_PLAIN_TYPE(Matrix6Type) Mtmp;
    Matrix6Type Mtmp, AdjointDual;

    Mat6ProjRz(P_z, P_, R_, M_a_, Mtmp);

    M_Aj_ += Mtmp;

    //! Ad*Xi*
    typedef Eigen::Block<Vector6Type,3,1> Segment3;

    Segment3 Pa_up   = P_a_.template segment<3>(0);
    Segment3 Pa_down = P_a_.template segment<3>(3);

    typename Vector6Type::PlainObject P_A_i;

    Segment3 Up_   = P_A_i.template segment<3>(0);
    Segment3 Down_ = P_A_i.template segment<3>(3);

    Up_.noalias()   = R_*Pa_up;
    Down_.noalias() = R_*Pa_down;

    if(P_z) Down_.noalias() += P_.cross(Up_);  // If P != 0

    P_Aj_.noalias() += P_A_i;


    //! Inertial back-projection differentiation.
    //!------------------------------------------------------------------------------!//
    //! M_a differentiation -> Single-motion-revolute screws enable skew-symmetric properties.
    //!------------------------------------------------------------------------------!//
    Matrix6Type M_a_S = Matrix6Type::Zero();  Matrix6Type M_a_S_;
    M_a_S.template row(0) = -M_a_.template row(1);  // mimicking effect ad_dual(-Sz)*Ma
    M_a_S.template row(1) =  M_a_.template row(0);
    M_a_S.template row(3) = -M_a_.template row(4);
    M_a_S.template row(4) =  M_a_.template row(3);

    M_a_S_.noalias() = M_a_S + M_a_S.transpose();   // mimicking effect ad_dual(-Sz)*Ma - Ma*ad(Sz)

    //! Transform M_a_S_ since it is now symmetric.
    Mat6ProjRz(P_z, P_, R_, M_a_S_, Mtmp);

    D_M_Aj_ = Mtmp;

    //! P_a differentiation -> D_q_Pa and D_dq_Pa.
    //!------------------------------------------------------------------------------!//
    D_q_Pa_ = D_q_p_;
    D_q_Pa_.noalias() += M_a_*D_q_c_;
    D_q_Pa_.noalias() += UD_*D_q_u_;

    D_dq_Pa_ = D_dq_p_;
    D_dq_Pa_.noalias() += M_a_*D_dq_c_;
    D_dq_Pa_.noalias() += UD_*D_dq_u_;


    //! D_q_Pa back projection.
    //!------------------------------------------------------------------------------!//
    P_A_i << -P_a_.coeff(1), P_a_.coeff(0), 0, -P_a_.coeff(4), P_a_.coeff(3), 0;
    AdDual(P_, R_, AdjointDual);
    D_q_PAj_.template rightCols<1>() = AdjointDual*P_A_i;
    D_q_PAj_.noalias() += AdjointDual*D_q_Pa_;


    //! D_dq_Pa back projection.
    //!------------------------------------------------------------------------------!//
    D_dq_PAj_.noalias() += AdjointDual*D_dq_Pa_;
  }


  //! Implementation for inertial expressions.
  template<typename ScalarType, typename Vector6iType, typename Vector3Type, typename Matrix3Type, typename Vector6Type,
           typename Matrix6Type, typename VectorXType, typename RowVectorXType, typename D_Vector6Type, typename D_Matrix6Type>
  inline void
  D_JointTypeRz::runD_Inertial(bool rootFlag,
                               Eigen::MatrixBase<Vector6iType> & iVec,  // ID, j, nP, nS, nPj, #children
                               ScalarType & u,
                               ScalarType & iD,
                               const ScalarType tau,
                               Eigen::MatrixBase<Vector6Type> & U_,
                               Eigen::MatrixBase<Vector6Type> & c_,
                               Eigen::MatrixBase<Vector6Type> & P_a_,
                               Eigen::MatrixBase<Matrix6Type> & M_a_,
                               Eigen::MatrixBase<Vector6Type> & P_A_,
                               Eigen::MatrixBase<Matrix6Type> & M_A_,
                               bool P_z,
                               Eigen::MatrixBase<Vector3Type> & P_,
                               Eigen::MatrixBase<Matrix3Type> & R_,
                               Eigen::MatrixBase<Vector6Type> & P_Aj_,
                               Eigen::MatrixBase<Matrix6Type> & M_Aj_,
                               Eigen::MatrixBase<D_Vector6Type> & D_U_h_,
                               Eigen::MatrixBase<VectorXType> & D_U_v_,
                               Eigen::MatrixBase<RowVectorXType> & D_invD_,
                               Eigen::MatrixBase<D_Matrix6Type> & D_M_A_,
                               Eigen::MatrixBase<D_Matrix6Type> & D_M_Aj_,
                               Eigen::MatrixBase<RowVectorXType> & D_q_u_,
                               Eigen::MatrixBase<RowVectorXType> & D_dq_u_,
                               Eigen::MatrixBase<D_Vector6Type> & D_q_Pa_,
                               Eigen::MatrixBase<D_Vector6Type> & D_dq_Pa_,
                               Eigen::MatrixBase<D_Vector6Type> & D_q_PA_,
                               Eigen::MatrixBase<D_Vector6Type> & D_dq_PA_,
                               Eigen::MatrixBase<D_Vector6Type> & D_q_PAj_,
                               Eigen::MatrixBase<D_Vector6Type> & D_dq_PAj_,
                               Eigen::MatrixBase<D_Vector6Type> & D_q_c_,
                               Eigen::MatrixBase<D_Vector6Type> & D_dq_c_) {
    //! Solve U and its partial derivative.
    //!------------------------------------------------------------------------------!
    U_ = M_A_.template rightCols<1>();
    //!-------------------------------------------------------
    D_U_v_ = D_M_A_.template rightCols<1>();
    short int nS = iVec.coeffRef(3) - 1;
    D_U_h_ = Eigen::Map<D_Vector6Type>(D_U_v_.derived().data(), 6, nS);


    //! Solve the inverse of D and its partial derivative.
    //!------------------------------------------------------------------------------!//
    iD = 1 / U_.coeffRef(5);
    //!-------------------------------------------------------
    D_invD_.noalias() = -pow(iD,2)*D_U_h_.template bottomRows<1>();


    //! Solve u and its partial derivative.
    //!------------------------------------------------------------------------------!//
    u = tau - P_A_.coeffRef(5);
    //!-------------------------------------------------------
    D_q_u_ = -D_q_PA_.template bottomRows<1>();  D_dq_u_ = -D_dq_PA_.template bottomRows<1>();


    //! Differentiating inertial back-projection to parent. Only if parent is not the root.
    //!------------------------------------------------------------------------------!//
    if (rootFlag) {

        //! Prepare inertial expression M_a.
        //!------------------------------------------------------------------------------!//
        Vector6Type UD_ = U_*iD;
        M_a_ = M_A_;
        M_a_.noalias() -= UD_*U_.transpose();

        //! Prepare inertial expression P_a.
        //!------------------------------------------------------------------------------!//
        P_a_ = P_A_;
        P_a_.noalias() += M_a_*c_;
        P_a_.noalias() += u*UD_;

        //! Prepare inertial expression differentiation D_M_a.
        //!------------------------------------------------------------------------------!//
        Matrix6Type Mtmp, AdjointDual;  VectorXType kronAux;  D_Matrix6Type D_MtmpI, D_MtmpII;

        short int ID, j, nP, nPj, nChild;
        ID = iVec.coeff(0);  j = iVec.coeff(1);  nP = iVec.coeff(2);  nPj = iVec.coeff(4);  nChild = iVec.coeff(5);

        //!------------------------------------------------------------------------------!//
        D_MtmpI.noalias() = D_U_v_*UD_.transpose();
        D_MtmpII = D_M_A_ - D_MtmpI;
        Mtmp.noalias() = U_*U_.transpose();

        AdDual(P_, R_, AdjointDual);
        int ite = 0;
        for(int iter = 0; iter < nS*6; iter += 6){
            D_MtmpII.middleRows(iter, 6).noalias() -= D_MtmpI.middleRows(iter, 6).transpose() + D_invD_.coeffRef(ite)*Mtmp;
            ite++;
          }

        //! then D_MtmpII = D_M_a and D_MtmpI is AdjointDual*D_M_a*Adjoint
        D_Mat6ProjRz(P_z, nS, P_, R_, D_MtmpII, D_MtmpI);

        //! Prepare inertial expression differentiation D_q_Pa.
        //!------------------------------------------------------------------------------!//
        kronAux = D_MtmpII*c_;
        D_Vector6Type D_Vec, D_VecII;
        D_Vec = Eigen::Map<D_Vector6Type>(kronAux.data(),6,nS);

        D_Vec.noalias() += (u*iD)*D_U_h_ + (u*U_)*D_invD_;

        D_q_Pa_.leftCols(nP).noalias() = M_a_*D_q_c_;
        D_q_Pa_.rightCols(nS) = D_Vec;

        D_q_Pa_.noalias() += D_q_PA_ + UD_*D_q_u_;

        //! Prepare inertial expression differentiation D_dq_Pa.
        //!------------------------------------------------------------------------------!//
        D_dq_Pa_.noalias() = D_dq_PA_ + UD_*D_dq_u_;
        D_dq_Pa_.leftCols(nP).noalias() += M_a_*D_dq_c_;

        //!------------------------------------------------------------------------------!//
        //!                         INERTIAL BACK PROJECTION                             !//
        //!------------------------------------------------------------------------------!//
        //! Back projection of M_a.
        //!------------------------------------------------------------------------------!//
        Mat6ProjRz(P_z, P_, R_, M_a_, Mtmp);
        M_Aj_.noalias() += Mtmp;


        //! Back projection of P_a.
        //!------------------------------------------------------------------------------!//
        typedef Eigen::Block<Vector6Type,3,1> Segment3;

        Segment3 Pa_up   = P_a_.template segment<3>(0);
        Segment3 Pa_down = P_a_.template segment<3>(3);

        Vector6Type P_A_i;

        Segment3 Up_   = P_A_i.template segment<3>(0);
        Segment3 Down_ = P_A_i.template segment<3>(3);

        Up_.noalias()   = R_*Pa_up;
        Down_.noalias() = R_*Pa_down;

        if(P_z) Down_.noalias() += P_.cross(Up_);  // If P != 0

        P_Aj_.noalias() += P_A_i;


        //! Back projection of the partial derivative D_M_a.
        //! M_a differentiation -> Single-motion-revolute screws enable skew-symmetric properties.
        //!------------------------------------------------------------------------------!//
        Matrix6Type M_a_S;  Matrix6Type M_a_S_;
        M_a_S.setZero();
        M_a_S.template row(0) = -M_a_.template row(1);  // mimicking effect ad_dual(-Sz)*Ma
        M_a_S.template row(1) =  M_a_.template row(0);
        M_a_S.template row(3) = -M_a_.template row(4);
        M_a_S.template row(4) =  M_a_.template row(3);

        M_a_S_.noalias() = M_a_S + M_a_S.transpose();   // mimicking effect ad_dual(-Sz)*Ma - Ma*ad(Sz)

        //! Transform M_a_S_ since it is now symmetric.
        Mat6ProjRz(P_z, P_, R_, M_a_S_, Mtmp);

        //!------------------------------------------------------------------------------!//
        //!                                  BRANCHING                                   !//
        //!------------------------------------------------------------------------------!//
        P_A_i << -P_a_.coeff(1), P_a_.coeff(0), 0, -P_a_.coeff(4), P_a_.coeff(3), 0;
        if (nChild == 1) {  //! Serial Body
            //!-------------------------------------------------------------------------  D_MA
            D_M_Aj_.template topRows<6>() = Mtmp;
            D_M_Aj_.middleRows(6, 6*nS) = D_MtmpI;
            //!-------------------------------------------------------------------------  D_q_PA
            D_q_PAj_.noalias() += AdjointDual*D_q_Pa_;

            D_q_PAj_.col(nP-1) += AdjointDual*P_A_i;
            //!-------------------------------------------------------------------------  D_dq_PA
            D_dq_PAj_.noalias() += AdjointDual*D_dq_Pa_;
          } else {  //! Branched Body
            //!-------------------------------------------------------------------------  D_MA
            D_M_Aj_.template middleRows<6>((ID-j)*6) = Mtmp;
            D_M_Aj_.middleRows((ID-j+1)*6, 6*nS) = D_MtmpI;
            //!-------------------------------------------------------------------------  D_q_PA
            D_Vec.noalias() = AdjointDual*D_q_Pa_;
            D_Vec.leftCols(nPj) += D_q_PAj_.leftCols(nPj);

            D_VecII = D_q_PAj_.rightCols(D_q_PAj_.cols()-nPj);

            D_q_PAj_.derived().conservativeResize(Eigen::NoChange, D_Vec.cols()+D_VecII.cols());
            D_q_PAj_ << D_Vec, D_VecII;

            D_q_PAj_.col(nP-1) += AdjointDual*P_A_i;
            //!-------------------------------------------------------------------------  D_dq_PA
            D_Vec.noalias() = AdjointDual*D_dq_Pa_;
            D_Vec.leftCols(nPj) += D_dq_PAj_.leftCols(nPj);

            D_VecII = D_dq_PAj_.rightCols(D_dq_PAj_.cols()-nPj);

            D_dq_PAj_.derived().conservativeResize(Eigen::NoChange, D_Vec.cols()+D_VecII.cols());
            D_dq_PAj_ << D_Vec, D_VecII;
          }

      }

  }


  //! Implementation for acceleration expressions at root.
  template<typename ScalarType, typename Vector3Type, typename Matrix3Type, typename Vector6Type,
           typename D_Vector6Type, typename RowVectorXType, typename MatrixXType>
  inline void
  D_JointTypeRz::runD_AccelRoot(ScalarType u,
                                ScalarType iD,
                                ScalarType* ddq,
                                const Eigen::MatrixBase<Vector3Type> & P_r,
                                const Eigen::MatrixBase<Matrix3Type> & R_r,
                                const Eigen::MatrixBase<Vector6Type> & U_r,
                                Eigen::MatrixBase<Vector6Type> & Acc_i_r,
                                Eigen::MatrixBase<D_Vector6Type> & D_U_h_,
                                Eigen::MatrixBase<RowVectorXType> & D_invD_,
                                Eigen::MatrixBase<RowVectorXType> & D_q_u_,
                                Eigen::MatrixBase<RowVectorXType> & D_dq_u_,
                                Eigen::MatrixBase<D_Vector6Type> & D_q_A_,
                                Eigen::MatrixBase<D_Vector6Type> & D_dq_A_,
                                Eigen::MatrixBase<MatrixXType> & D_ddq_) {

    EIGEN_STATIC_ASSERT(Vector6Type::ColsAtCompileTime == 1,
                        YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_A_MATRIX);

    //! Acceleration bias and its derivative.
    //!------------------------------------------------------------------------------!//
    //! g = [ 0, 0, 9.81, 0, 0, 0 ]^T
    //! Acc_a = Ad(G[Sz])*g = [ 0 0 9.81 0 0 0 ]^T
    //! D_Acc_a = -ad(Sz)*Ad(G)*g = 0


    //! Joint acceleration and its derivatives D_q_ddq and D_dq_ddq.
    //!------------------------------------------------------------------------------!//
    ScalarType ddq_;
    ddq_ = u - 9.81*U_r.template segment<1>(2)(0);
    (*ddq) = iD*ddq_;
    //!-------------------------------------------------------
    typename RowVectorXType::PlainObject D_q_ddq, D_dq_ddq;
    D_q_ddq.noalias() = iD*D_q_u_;  D_dq_ddq.noalias() = iD*D_dq_u_;
    D_q_ddq.segment(1,D_U_h_.cols()) += ddq_*D_invD_ - (iD*9.81)*D_U_h_.row(2);


    //! Update spatial acceleration and its derivative.
    //!------------------------------------------------------------------------------!//
    Acc_i_r.coeffRef(2) = 9.81;
    Acc_i_r.coeffRef(5) = iD*ddq_;
    //!-------------------------------------------------------
    D_q_A_.template bottomRows<1>() = D_q_ddq;  D_dq_A_.template bottomRows<1>() = D_dq_ddq;


    //! Fill up D_ddq.
    //!------------------------------------------------------------------------------!//
    D_ddq_.row(0) << D_q_ddq, D_dq_ddq;
  }


  //! Implementation for acceleration expressions.
  template<typename ScalarType, typename IndexType, typename Vector3Type, typename Matrix3Type,
           typename Vector6Type, typename D_Vector6Type, typename RowVectorXType, typename MatrixXType>
  inline void
  D_JointTypeRz::runD_Accel(IndexType ID,
                            bool zeroFlag,
                            ScalarType u,
                            ScalarType iD,
                            ScalarType* ddq,
                            const Eigen::MatrixBase<Vector3Type> & P_,
                            const Eigen::MatrixBase<Matrix3Type> & R_,
                            const Eigen::MatrixBase<Vector6Type> & c_,
                            const Eigen::MatrixBase<Vector6Type> & U_,
                            Eigen::MatrixBase<Vector6Type> & A_,
                            const Eigen::MatrixBase<Vector6Type> & Aj_,
                            bool isLeaf,
                            std::vector< IndexType >* Pre_,
                            std::vector< IndexType >* Suc_,
                            std::vector< IndexType >* PreSuc_,
                            Eigen::MatrixBase<D_Vector6Type> & D_U_h_,
                            Eigen::MatrixBase<RowVectorXType> & D_invD_,
                            Eigen::MatrixBase<MatrixXType> & D_ddq_,
                            Eigen::MatrixBase<RowVectorXType> & D_q_u_,
                            Eigen::MatrixBase<RowVectorXType> & D_dq_u_,
                            Eigen::MatrixBase<D_Vector6Type> & D_q_c_,
                            Eigen::MatrixBase<D_Vector6Type> & D_dq_c_,
                            Eigen::MatrixBase<D_Vector6Type> & D_q_A_,
                            Eigen::MatrixBase<D_Vector6Type> & D_dq_A_,
                            Eigen::MatrixBase<D_Vector6Type> & D_q_Aj_,
                            Eigen::MatrixBase<D_Vector6Type> & D_dq_Aj_) {

    EIGEN_STATIC_ASSERT(Vector6Type::ColsAtCompileTime == 1,
                        YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_A_MATRIX);

    //! Acceleration bias and its derivative.
    //!------------------------------------------------------------------------------!//
    typedef const Eigen::Block<Vector6Type,3,1> constSegment3;
    typedef Eigen::Block<Vector6Type,3,1> Segment3;

    Vector6Type & Aj__ = const_cast<Vector6Type &>(Aj_.derived());
    constSegment3 & AjUp = Aj__.template segment<3>(0);
    constSegment3 & AjDown = Aj__.template segment<3>(3);

    typename Vector6Type::PlainObject Aa_, AdAj;

    Segment3 AaUp = Aa_.template segment<3>(0);
    Segment3 AaDown = Aa_.template segment<3>(3);

    if(zeroFlag) {
        //! delete cross product
        AaDown.noalias() = AjUp - P_.cross(AjDown);  //tmp
        AaUp.noalias()   = R_.transpose()*AaDown;
      } else {
        AaUp.noalias()   = R_.transpose()*AjUp;
      }

    AaDown.noalias() = R_.transpose()*AjDown;

    //!-------------------------------------------------------
    //! mimicking effect -ad(Sz)*Ad*Aj
    AdAj << Aa_.coeff(1), -Aa_.coeff(0), 0, Aa_.coeff(4), -Aa_.coeff(3), 0;
    //!-------------------------------------------------------

    Aa_.noalias() += c_;
    //!------------------------------------------------------------------------------!//
    Eigen::Matrix<ScalarType, 6, 6> AdjointDual;
    AdDual(P_, R_, AdjointDual);
    typename D_Vector6Type::PlainObject D_q_Aa_, D_dq_Aa_;

    D_q_Aa_.noalias() = AdjointDual.transpose()*D_q_Aj_;
    D_q_c_.template rightCols<1>() += AdAj;

    D_q_Aa_(Eigen::all, *Pre_) += D_q_c_;

    D_dq_Aa_.noalias() = AdjointDual.transpose()*D_dq_Aj_;
    D_dq_Aa_(Eigen::all, *Pre_) += D_dq_c_;


    //! Joint acceleration and its derivatives D_q_ddq and D_dq_ddq.
    //!------------------------------------------------------------------------------!//
    ScalarType ddq__, ddq_;
    ddq__ = u - U_.transpose()*Aa_;  ddq_ = iD*ddq__;
    (*ddq) = ddq_;
    //!-------------------------------------------------------
    typename RowVectorXType::PlainObject D_q_ddq, D_dq_ddq;
    D_q_ddq = -iD*U_.transpose()*D_q_Aa_;
    D_q_ddq(Eigen::all, *PreSuc_) += iD*D_q_u_;

    D_dq_ddq = -iD*U_.transpose()*D_dq_Aa_;
    D_dq_ddq(Eigen::all, *PreSuc_) += iD*D_dq_u_;


    //! Update spatial acceleration and its derivatives D_q_A and D_dq_A.
    //!------------------------------------------------------------------------------!//
    A_ = Aa_;
    A_.coeffRef(5) += ddq_;


    //! Special treatment when body is not a leaf.
    //!------------------------------------------------------------------------------!//
    if (isLeaf) {
        D_q_ddq(Eigen::all, *Suc_) += ddq__*D_invD_ - iD*Aa_.transpose()*D_U_h_;
        //!------------------------------------------------------------------------------!//
        D_q_A_ = D_q_Aa_;
        D_q_A_.template bottomRows<1>() += D_q_ddq;  //! TODO: D_q_A_ can be tmp used as D_q_Aa_; saving an assignation

        D_dq_A_ = D_dq_Aa_;
        D_dq_A_.template bottomRows<1>() += D_dq_ddq;
      }


    //! Fill up D_ddq.
    //!------------------------------------------------------------------------------!//
    D_ddq_.row(ID) << D_q_ddq, D_dq_ddq;
  }

} // end namespace geo

#endif // GEOMBD_DIFFERENTIATION_JOINT_TYPE_REVOLUTE_Z_HXX
