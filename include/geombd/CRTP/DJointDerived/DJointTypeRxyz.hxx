/**
 *    \file include/geombd/CRTP/DJointDerived/DJointTypeRxyz.hxx
 *    \author Alvaro Paz, Gustavo Arechavaleta
 *    \version 1.0
 *    \date 2021
 *
 *    Derived class for Rxyz joint type
 *    Copyright (c) 2021 Cinvestav
 *    This library is distributed under the MIT License.
 */

#ifndef GEOMBD_DIFFERENTIATION_JOINT_TYPE_REVOLUTE_XYZ_HXX
#define GEOMBD_DIFFERENTIATION_JOINT_TYPE_REVOLUTE_XYZ_HXX

#define EIGEN_NO_DEBUG
#define EIGEN_MPL2_ONLY
#define EIGEN_UNROLLING_LIMIT 30

#include "Eigen/Core"
#include "Eigen/Geometry"

//#include <iostream>
//#include <memory>

namespace geo{


  //! Derived -> Differentiation of Joint Type Revolute On XYZ Axis
  //!------------------------------------------------------------------------------!//
  struct D_JointTypeRxyz : public D_CRTPInterface<D_JointTypeRxyz>
  {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! FK Forward Declaration.
    template<typename ScalarType, typename Vector3Type, typename Matrix3Type>
    inline static void
    runD_FK(const ScalarType & qi,
            const Eigen::MatrixBase<Vector3Type> & S,
            typename Eigen::MatrixBase<Matrix3Type> & R);


    //! TwCbPb at Root Forward Declaration.
    template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename Matrix6Type, typename D_Vector6Type>
    inline static void
    runD_TCP_root(const ScalarType & vi,
                  const Eigen::MatrixBase<Vector3Type> & S,
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
                const Eigen::MatrixBase<Vector3Type> & S,
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
                     const Eigen::MatrixBase<Vector3Type> & S,
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
                  const Eigen::MatrixBase<Vector3Type> & S,
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
                   const Eigen::MatrixBase<Vector3Type> & S,
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
               const Eigen::MatrixBase<Vector3Type> & S,
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
  template<typename ScalarType, typename Vector3Type, typename Matrix3Type>
  inline void
  D_JointTypeRxyz::runD_FK(const ScalarType & qi,
                           const Eigen::MatrixBase<Vector3Type> & S,
                           typename Eigen::MatrixBase<Matrix3Type> & R) {
    static ScalarType sqi, cqi;
    ::geo::SINCOS<ScalarType>(qi, &sqi, &cqi);

    ScalarType k0, k1, k2, k3, k4, k5, k6;
    k0 = cqi-1;
    k1 = k0*S[2]*S[2];    k4 = S[1]*S[2]*k0;
    k2 = k0*S[1]*S[1];    k5 = S[0]*S[2]*k0;
    k3 = k0*S[0]*S[0];    k6 = S[0]*S[1]*k0;

    R.coeffRef(0,0) = k2+k1+1;       R.coeffRef(0,1) = -S[2]*sqi-k6;  R.coeffRef(0,2) = S[1]*sqi-k5;
    R.coeffRef(1,0) = S[2]*sqi-k6;   R.coeffRef(1,1) = k3+k1+1;       R.coeffRef(1,2) = -S[0]*sqi-k4;
    R.coeffRef(2,0) = -S[1]*sqi-k5;  R.coeffRef(2,1) = S[0]*sqi-k4;   R.coeffRef(2,2) = k3+k2+1;

  }


  //! TwCbPb at Root Declaration.
  template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename Matrix6Type, typename D_Vector6Type>
  inline void
  D_JointTypeRxyz::runD_TCP_root(const ScalarType & vi,
                                 const Eigen::MatrixBase<Vector3Type> & S,
                                 Eigen::MatrixBase<Vector6Type> & S_i,
                                 Eigen::MatrixBase<Vector6Type> & p_,
                                 const Eigen::MatrixBase<Matrix6Type> & M_,
                                 Eigen::MatrixBase<D_Vector6Type> & D_dq_p_) {

    typename Vector3Type::PlainObject Twist_w;
    Twist_w = S*vi;
    S_i.template tail<3>() = Twist_w;

    //! P bias
    //!------------------------------------------------------------------------------!//
    typedef Eigen::Block<Vector6Type,3,1> Segment3;

    typedef const Eigen::Block<Matrix6Type,3,3> constBlock3;
    typedef Eigen::Block<Matrix6Type,3,3> Block3;

    Matrix6Type & M = const_cast<Matrix6Type &>(M_.derived());
    constBlock3 & UR = M.template block<3,3>(0,3);
    constBlock3 & LR = M.template block<3,3>(3,3);

    typename Vector3Type::PlainObject muUp, muDo;

    Segment3 p_Up = p_.template segment<3>(0);
    Segment3 p_Do = p_.template segment<3>(3);

    //! M is symmetric and its shape is known
    muUp.noalias() = UR*Twist_w;
    muDo.noalias() = LR*Twist_w;

    //! TODO: Cross is very expensive, consider enable staticness with traits expressions
    p_Up.noalias() = Skew(Twist_w)*muUp;
    p_Do.noalias() = Skew(Twist_w)*muDo;

    //! P bias differentiation
    //!------------------------------------------------------------------------------!//
    typedef Eigen::Block<Matrix6Type,3,3> Block3;
    Matrix6Type D_p_;

    Block3 UR_ = D_p_.template block<3,3>(0,3);
    Block3 LL_ = D_p_.template block<3,3>(3,0);  Block3 LR_ = D_p_.template block<3,3>(3,3);

    //! adDual
//    D_p_.setZero();
    D_p_.template topRows<3>()    = -Skew(S) * M.template topRows<3>();  //! why S instead Twist_w?
    D_p_.template bottomRows<3>() = -Skew(S) * M.template bottomRows<3>();


    Vector6Type mu;
    mu.noalias() = M.template rightCols<3>()*Twist_w;

    //! adBar
    UR_.noalias() -= Skew(mu.segment(0,3));
    LL_.noalias() -= Skew(mu.segment(0,3));  LR_.noalias() -= Skew(mu.segment(3,3));

    D_dq_p_.noalias() = D_p_.template rightCols<3>()*S;
  }


  //! Apply the Ad operator on se(3) element to get the Twist and bias
  //! TODO: Twist is only used here, thus consider to store it only when multiple children
  //!------------------------------------------------------------------------------!//
  template<typename ScalarType, typename Matrix3Type, typename Matrix6Type,
           typename Vector3Type, typename Vector6Type, typename D_Vector6Type>
  inline void
  D_JointTypeRxyz::runD_TwCbPb(bool zeroFlag,
                               const ScalarType & vi,
                               const Eigen::MatrixBase<Vector3Type> & S,
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

    typename Vector3Type::PlainObject Twist_w;
    Twist_w = S*vi;

    Down_ += Twist_w;
    //!-------------------------------------------------------
    D_q_V_.block(3,0,3,icol-1).noalias()  = R_.transpose()*D_q_Vj_.template bottomRows<3>();
    D_dq_V_.block(3,0,3,icol-1).noalias() = R_.transpose()*D_dq_Vj_.template bottomRows<3>();

    //! C bias
    //!------------------------------------------------------------------------------!//
    Segment3 Cup   = c_.template segment<3>(0);
    Segment3 Cdown = c_.template segment<3>(3);

    Cup = Up_.cross(S);
    Cdown = Down_.cross(S);

    //! Twist differentiation -> taking advantage of c bias computation.
    //!------------------------------------------------------------------------------!//
    D_q_V_.template rightCols<1>() = c_;
    //!------------------------------------------------------------------------------!//

    c_ *= vi;

    //! C bias differentiation -> taking advantage of twist computation.
    //!------------------------------------------------------------------------------!//
    D_Vector6Type D_c_aux;
    D_c_aux.noalias() = D_q_V_*vi;
    typename Matrix3Type::PlainObject skew_S = Skew(S);
    D_q_c_.template topRows<3>()    = -skew_S * D_c_aux.template topRows<3>();  // -ad(Sxyz)* effect
    D_q_c_.template bottomRows<3>() = -skew_S * D_c_aux.template bottomRows<3>();

    D_c_aux.noalias() = D_dq_V_*vi;
    D_c_aux.template rightCols<1>() += S_i;
    D_dq_c_.template topRows<3>()    = -skew_S * D_c_aux.template topRows<3>();  // -ad(Sxyz)* effect
    D_dq_c_.template bottomRows<3>() = -skew_S * D_c_aux.template bottomRows<3>();
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
  D_JointTypeRxyz::runD_InertiaLeaf(ScalarType & u,
                                    ScalarType & iD,
                                    const ScalarType tau,
                                    const Eigen::MatrixBase<Vector3Type> & S,
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
    U_.noalias() = M_A_.template rightCols<3>()*S;
    ScalarType D;  D = S.dot( U_.template segment<3>(3) );  iD = 1 / D;
    u = tau - S.dot( P_A_.template segment<3>(3) );


    //! Solve partial derivative of u as D_q_u and D_dq_u.
    //!------------------------------------------------------------------------------!//
    D_q_u_  = - S.transpose() * D_q_PA_.template bottomRows<3>();
    D_dq_u_ = - S.transpose() * D_dq_PA_.template bottomRows<3>();


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

    Mat6ProjRxyz(P_z, P_, R_, M_a_, Mtmp);

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
    //! M_a differentiation -> Screws due to pure-revolute-motions enable skew-symmetric properties.
    //!------------------------------------------------------------------------------!//
    Matrix6Type M_a_S = Matrix6Type::Zero();  Matrix6Type M_a_S_;
    typename Matrix3Type::PlainObject skew_S = Skew(S);

    M_a_S.template topRows<3>()    = skew_S * M_a_.template topRows<3>();  // mimicking ad_dual(-Sxyz)*Ma effect
    M_a_S.template bottomRows<3>() = skew_S * M_a_.template bottomRows<3>();

    M_a_S_.noalias() = M_a_S + M_a_S.transpose();   // mimicking ad_dual(-Sxyz)*Ma - Ma*ad(Sxyz) effect

    //! Transform M_a_S_ since it is now symmetric.
    Mat6ProjRxyz(P_z, P_, R_, M_a_S_, Mtmp);

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
    P_A_i.template segment<3>(0) = S.cross(P_a_.template segment<3>(0));  // mimicking ad_dual(-Sxyz)*Pa effect
    P_A_i.template segment<3>(3) = S.cross(P_a_.template segment<3>(3));


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
  D_JointTypeRxyz::runD_Inertial(bool rootFlag,
                                 Eigen::MatrixBase<Vector6iType> & iVec,  // ID, j, nP, nS, nPj, #children
                                 ScalarType & u,
                                 ScalarType & iD,
                                 const ScalarType tau,
                                 const Eigen::MatrixBase<Vector3Type> & S,
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
    U_.noalias() = M_A_.template rightCols<3>()*S;
    //!-------------------------------------------------------
    D_U_v_ = D_M_A_.template rightCols<3>()*S;
    short int nS = iVec.coeffRef(3) - 1;
    D_U_h_ = Eigen::Map<D_Vector6Type>(D_U_v_.derived().data(), 6, nS);


    //! Solve the inverse of D and its partial derivative.
    //!------------------------------------------------------------------------------!//
    ScalarType D;
    D = S.dot( U_.template segment<3>(3) );
    iD = 1 / D;
    //!-------------------------------------------------------
    D_invD_.noalias() = -pow(iD,2)*S.transpose()*D_U_h_.template bottomRows<3>();


    //! Solve u and its partial derivative.
    //!------------------------------------------------------------------------------!//
    u = tau - S.dot( P_A_.template segment<3>(3) );
    //!-------------------------------------------------------
    D_q_u_ = -S.transpose()*D_q_PA_.template bottomRows<3>();
    D_dq_u_ = -S.transpose()*D_dq_PA_.template bottomRows<3>();


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
        D_Mat6ProjRxyz(P_z, nS, P_, R_, D_MtmpII, D_MtmpI);

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
        Mat6ProjRxyz(P_z, P_, R_, M_a_, Mtmp);
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
        typename Matrix3Type::PlainObject skew_S = Skew(S);

        M_a_S.template topRows<3>()    = skew_S * M_a_.template topRows<3>();  // mimicking ad_dual(-Sxyz)*Ma effect
        M_a_S.template bottomRows<3>() = skew_S * M_a_.template bottomRows<3>();

        M_a_S_.noalias() = M_a_S + M_a_S.transpose();   // mimicking ad_dual(-Sxyz)*Ma - Ma*ad(Sxyz) effect

        //! Transform M_a_S_ since it is now symmetric.
        Mat6ProjRxyz(P_z, P_, R_, M_a_S_, Mtmp);

        //!------------------------------------------------------------------------------!//
        //!                                  BRANCHING                                   !//
        //!------------------------------------------------------------------------------!//
        P_A_i.template segment<3>(0) = S.cross(P_a_.template segment<3>(0));  // mimicking ad_dual(-Sxyz)*Pa effect
        P_A_i.template segment<3>(3) = S.cross(P_a_.template segment<3>(3));

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
  D_JointTypeRxyz::runD_AccelRoot(ScalarType u,
                                  ScalarType iD,
                                  ScalarType* ddq,
                                  const Eigen::MatrixBase<Vector3Type> & S,
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
    //! Acc_a = Ad(G[Sxyz])*g = [ 9.81*R_r.row(2), 0, 0, 0 ]^T
    //! D_Acc_a = -ad(Sz)*Ad(G[Sxyz])*g = 0

    Eigen::Matrix<ScalarType, 1, 3> R_row;
    R_row = 9.81*R_r.row(2);
    Acc_i_r.setZero();  D_q_A_.setZero();
    Acc_i_r.template segment<3>(0) = R_row.transpose();

    D_q_A_.template block<3,1>(0,0) = -S.cross(R_row.transpose());


    //! Joint acceleration and its derivatives D_q_ddq and D_dq_ddq.
    //!------------------------------------------------------------------------------!//
    ScalarType ddq_;
    ddq_ = u - U_r.transpose()*Acc_i_r;
    (*ddq) = ddq_*iD;


    //!-------------------------------------------------------
    RowVectorXType D_q_ddq, D_dq_ddq;
    D_q_ddq.noalias() = iD*D_q_u_;  D_dq_ddq.noalias() = iD*D_dq_u_;
    ScalarType U_D_A = U_r.transpose()*D_q_A_.template leftCols<1>();
    D_q_ddq.coeffRef(0) -= iD*U_D_A;
    D_q_ddq.segment(1,D_U_h_.cols()).noalias() += ddq_*D_invD_ - iD*R_row*D_U_h_.template topRows<3>();


    //! Update spatial acceleration and its derivative.
    //!------------------------------------------------------------------------------!//
    Acc_i_r.template segment<3>(3) = S*iD*ddq_;
    //!-------------------------------------------------------
    D_q_A_.template bottomRows<3>() = S*D_q_ddq;  D_dq_A_.template bottomRows<3>() = S*D_dq_ddq;


    //! Fill up D_ddq.
    //!------------------------------------------------------------------------------!//
    D_ddq_.row(0) << D_q_ddq, D_dq_ddq;
  }


  //! Implementation for acceleration expressions.
  template<typename ScalarType, typename IndexType, typename Vector3Type, typename Matrix3Type,
           typename Vector6Type, typename D_Vector6Type, typename RowVectorXType, typename MatrixXType>
  inline void
  D_JointTypeRxyz::runD_Accel(IndexType ID,
                              bool zeroFlag,
                              ScalarType u,
                              ScalarType iD,
                              ScalarType* ddq,
                              const Eigen::MatrixBase<Vector3Type> & S,
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
    //! mimicking effect -ad(Sxyz)*Ad*Aj
    AdAj.template segment<3>(0) = - S.cross( Aa_.template segment<3>(0));
    AdAj.template segment<3>(3) = - S.cross( Aa_.template segment<3>(3));

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
    A_.template segment<3>(3) += S*ddq_;


    //! Special treatment when body is not a leaf.
    //!------------------------------------------------------------------------------!//
    if (isLeaf) {
        D_q_ddq(Eigen::all, *Suc_) += ddq__*D_invD_ - iD*Aa_.transpose()*D_U_h_;
        //!------------------------------------------------------------------------------!//
        D_q_A_ = D_q_Aa_;
        D_q_A_.template bottomRows<3>() += S*D_q_ddq;  //! TODO: D_q_A_ can be tmp used as D_q_Aa_; saving an assignation

        D_dq_A_ = D_dq_Aa_;
        D_dq_A_.template bottomRows<3>() += S*D_dq_ddq;
      }


    //! Fill up D_ddq.
    //!------------------------------------------------------------------------------!//
    D_ddq_.row(ID) << D_q_ddq, D_dq_ddq;
  }

} // end namespace geo

#endif // GEOMBD_DIFFERENTIATION_JOINT_TYPE_REVOLUTE_XYZ_HXX
