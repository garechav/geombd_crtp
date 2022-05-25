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
#include "geombd/CRTP/traits/GeoOperator"
#include "geombd/CRTP/utils/functors.hxx"

//#include <iostream>
//#include <memory>

namespace geo{


  //! Derived -> Differentiation of Joint Type Revolute On XYZ Axis
  //!------------------------------------------------------------------------------!//
  struct D_JointTypeRxyz : public D_CRTPInterface<D_JointTypeRxyz>
  {
  public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW    

    //! Forward Kinematics Declaration
    //! EIGEN_ASM_COMMENT("MyBegin");
    //!------------------------------------------------------------------------------!//
    template<typename ScalarType, typename Vector3Type, typename Matrix3Type>
    inline static void
    runD_FK(const ScalarType & qi,
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


    //! TCP01 Declaration
    //!------------------------------------------------------------------------------!//
    template<typename ScalarType, typename Vector3Type, typename Vector6Type>
    inline static void
    runD_TCP01(const ScalarType & vi,
               typename Eigen::MatrixBase<Vector3Type> & Sw_,
               typename Eigen::MatrixBase<Vector6Type> & S_,
               typename Eigen::MatrixBase<Vector6Type> & C_) {
      //!------------------------------------------------------------------------------!//
      typedef Eigen::Block<Vector6Type,3,1> Segment3;

      Segment3 S_up = S_.template segment<3>(0);
      Segment3 S_dw = S_.template segment<3>(3);

      Segment3 C_up = C_.template segment<3>(0);
      Segment3 C_dw = C_.template segment<3>(3);
      //!------------------------------------------------------------------------------!//
      S_dw.noalias() += Sw_*vi;
      //!------------------------------------------------------------------------------!//
      C_up = S_up.cross(Sw_);      // ad(Twist)*Sxyz effect
      C_dw = S_dw.cross(Sw_);
      //!------------------------------------------------------------------------------!//
    }


    //! TCP02 Declaration
    //!------------------------------------------------------------------------------!//
    template<typename Vector3Type, typename D_Vector6Type>
    inline static void
    runD_TCP02(typename Eigen::MatrixBase<Vector3Type> & Sw_,
               typename Eigen::MatrixBase<D_Vector6Type> & D_q_c_,
               typename Eigen::MatrixBase<D_Vector6Type> & D_dq_c_,
               typename Eigen::MatrixBase<D_Vector6Type> & D_q_c_aux_,
               typename Eigen::MatrixBase<D_Vector6Type> & D_dq_c_aux_) {

      typedef typename Vector3Type::Scalar ScalarType;
      typedef Eigen::Matrix<ScalarType , 3 , 3> Matrix3Type;
      //!------------------------------------------------------------------------------!//
      Matrix3Type skew_S = Skew(Sw_);
      D_q_c_.template topRows<3>()    = -skew_S * D_q_c_aux_.template topRows<3>();     // -ad(Sxyz)* effect
      D_q_c_.template bottomRows<3>() = -skew_S * D_q_c_aux_.template bottomRows<3>();
      //!------------------------------------------------------------------------------!//
      D_dq_c_.template topRows<3>()    = -skew_S * D_dq_c_aux_.template topRows<3>();   // -ad(Sxyz)* effect
      D_dq_c_.template bottomRows<3>() = -skew_S * D_dq_c_aux_.template bottomRows<3>();

    }


    //! TCProot Declaration.
    //!------------------------------------------------------------------------------!//
    template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename Matrix6Type, typename D_Vector6Type>
    inline static void
    runD_TCProot(const ScalarType & vi,
                 typename Eigen::MatrixBase<Vector3Type> & S_,
                 typename Eigen::MatrixBase<Vector6Type> & T_,
                 typename Eigen::MatrixBase<Vector6Type> & p_,
                 typename Eigen::MatrixBase<Matrix6Type> & M_,
                 typename Eigen::MatrixBase<D_Vector6Type> & D_dq_p_) {

      EIGEN_STATIC_ASSERT(Vector6Type::ColsAtCompileTime == 1,
                          YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_A_MATRIX);

      typename Vector3Type::PlainObject Twist_w;
      Twist_w = S_*vi;
      T_.template tail<3>() = Twist_w;

      //! P bias
      //!------------------------------------------------------------------------------!//
      typedef Eigen::Block<Vector6Type,3,1> Segment3;

      typedef const Eigen::Block<Matrix6Type,3,3> constBlock3;
      typedef Eigen::Block<Matrix6Type,3,3> Block3;

      constBlock3 & UR = M_.template block<3,3>(0,3);
      constBlock3 & LR = M_.template block<3,3>(3,3);

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
      Matrix6Type D_p_;

      Block3 UR_ = D_p_.template block<3,3>(0,3);
      Block3 LR_ = D_p_.template block<3,3>(3,3);

      //! adDual
      //    D_p_.setZero();
      UR_.noalias() = Skew(Twist_w) * UR;  //! why S instead Twist_w? and negative?
      LR_.noalias() = Skew(Twist_w) * LR;

      //! adBar
      UR_.noalias() -= Skew(muUp);
      LR_.noalias() -= Skew(muDo);

      D_dq_p_.noalias() = D_p_.template rightCols<3>()*S_;

    }


    //! Inertial Terms 01 Declaration.
    //!------------------------------------------------------------------------------!//
    template<typename Vector3Type, typename Vector6Type, typename Matrix6Type, typename VectorXType, typename D_Matrix6Type>
    inline static void
    runInertia01(typename Eigen::MatrixBase<Vector3Type> & Sw_,
                 typename Eigen::MatrixBase<Vector6Type> & U_,
                 typename Eigen::MatrixBase<Matrix6Type> & M_A_,
                 typename Eigen::MatrixBase<VectorXType> & D_U_v_,
                 typename Eigen::MatrixBase<D_Matrix6Type> & D_M_A_i_) {

      //! Solve U and its partial derivative.
      //!------------------------------------------------------------------------------!//
      U_.noalias() = M_A_.template rightCols<3>()*Sw_;
      //!------------------------------------------------------------------------------!//
      D_U_v_ = D_M_A_i_.template rightCols<3>()*Sw_;

    }


    //! Inertial Terms 02 Declaration.
    //!------------------------------------------------------------------------------!//
    template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename RowVectorXType, typename D_Vector6Type>
    inline static void
    runInertia02(ScalarType & invD_,
                 ScalarType & u_,
                 typename Eigen::MatrixBase<Vector3Type> & Sw_,
                 typename Eigen::MatrixBase<Vector6Type> & U_,
                 typename Eigen::MatrixBase<Vector6Type> & P_A_,
                 typename Eigen::MatrixBase<RowVectorXType> & D_invD_,
                 typename Eigen::MatrixBase<RowVectorXType> & D_q_u_,
                 typename Eigen::MatrixBase<RowVectorXType> & D_dq_u_,
                 typename Eigen::MatrixBase<D_Vector6Type> & D_U_h_,
                 typename Eigen::MatrixBase<D_Vector6Type> & D_q_PA_,
                 typename Eigen::MatrixBase<D_Vector6Type> & D_dq_PA_) {

      //! Solve the inverse of D and its partial derivative.
      //!------------------------------------------------------------------------------!//
      ScalarType D;
      D = Sw_.dot( U_.template segment<3>(3) );
      invD_ = 1 / D;
      //!------------------------------------------------------------------------------!//
      D_invD_.noalias() = -pow(invD_,2)*Sw_.transpose()*D_U_h_.template bottomRows<3>();

      //! Solve u and its partial derivative.
      //!------------------------------------------------------------------------------!//
      u_ -= Sw_.dot( P_A_.template segment<3>(3) );
      //!------------------------------------------------------------------------------!//
      D_q_u_ = -Sw_.transpose()*D_q_PA_.template bottomRows<3>();
      D_dq_u_ = -Sw_.transpose()*D_dq_PA_.template bottomRows<3>();

    }


    //! Inertial Terms 03 Declaration.
    //!------------------------------------------------------------------------------!//
    template<typename IndexType, typename Vector3Type, typename Matrix3Type, typename Vector6Type, typename Matrix6Type, typename D_Matrix6Type>
    inline static void
    runInertia03(bool P_z_,
                 IndexType nS_,
                 typename Eigen::MatrixBase<Vector3Type> & Sw_,
                 typename Eigen::MatrixBase<Vector3Type> & P_,
                 typename Eigen::MatrixBase<Matrix3Type> & R_,
                 typename Eigen::MatrixBase<Vector6Type> & P_a_,
                 typename Eigen::MatrixBase<Vector6Type> & P_A_i_,
                 typename Eigen::MatrixBase<Matrix6Type> & M_a_,
                 typename Eigen::MatrixBase<Matrix6Type> & Mtmp_,
                 typename Eigen::MatrixBase<D_Matrix6Type> & D_M_A_i_) {

      //! Back projection of M_a.
      //!------------------------------------------------------------------------------!//
      Mat6ProjRxyz(P_z_, P_.derived(), R_.derived(), M_a_.derived(), Mtmp_.derived());

      //! Back projection of D_M_a.
      //!------------------------------------------------------------------------------!//
      D_Mat6ProjRxyz(P_z_, nS_, P_.derived(), R_.derived(), D_M_A_i_.derived());

      //! Back projection of the partial derivative D_M_a.
      //! M_a differentiation -> Single-motion-revolute screws enable skew-symmetric properties.
      //!------------------------------------------------------------------------------!//
      Matrix6Type M_a_S, M_a_S_;
      typename Matrix3Type::PlainObject skew_S = Skew(Sw_);

      M_a_S.template topRows<3>()    = skew_S * M_a_.template topRows<3>();  // mimicking ad_dual(-Sxyz)*Ma effect
      M_a_S.template bottomRows<3>() = skew_S * M_a_.template bottomRows<3>();

      M_a_S_.noalias() = M_a_S + M_a_S.transpose();   // mimicking ad_dual(-Sxyz)*Ma - Ma*ad(Sxyz) effect

      //! Transform M_a_S_ since it is now symmetric.
      //! Back projection of the partial derivative D_M_a.
      //!------------------------------------------------------------------------------!//
      Mat6ProjRxyz(P_z_, P_.derived(), R_.derived(), M_a_S_.derived(), M_a_S.derived());

      D_M_A_i_.template topRows<6>() = M_a_S;

      P_A_i_.template segment<3>(0) = Sw_.cross(P_a_.template segment<3>(0));  // mimicking ad_dual(-Sxyz)*Pa effect
      P_A_i_.template segment<3>(3) = Sw_.cross(P_a_.template segment<3>(3));

    }


    //! Leaf Terms 01 Declaration.
    //!------------------------------------------------------------------------------!//
    template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename Matrix6Type, typename RowVectorXType, typename D_Vector6Type>
    inline static void
    runLeaf01(ScalarType & invD_,
              ScalarType & u_,
              typename Eigen::MatrixBase<Vector3Type> & Sw_,
              typename Eigen::MatrixBase<Vector6Type> & U_,
              typename Eigen::MatrixBase<Vector6Type> & P_A_,
              typename Eigen::MatrixBase<Matrix6Type> & M_A_,
              typename Eigen::MatrixBase<RowVectorXType> & D_q_u_,
              typename Eigen::MatrixBase<RowVectorXType> & D_dq_u_,
              typename Eigen::MatrixBase<D_Vector6Type> & D_q_PA_,
              typename Eigen::MatrixBase<D_Vector6Type> & D_dq_PA_) {

      //! Solve U, inverse of D and u.
      //!------------------------------------------------------------------------------!//
      U_.noalias() = M_A_.template rightCols<3>()*Sw_;
      ScalarType D;  D = Sw_.dot( U_.template segment<3>(3) );  invD_ = 1 / D;
      u_ -= Sw_.dot( P_A_.template segment<3>(3) );

      //! Solve partial derivative of u as D_q_u and D_dq_u.
      //!------------------------------------------------------------------------------!//
      D_q_u_  = - Sw_.transpose() * D_q_PA_.template bottomRows<3>();
      D_dq_u_ = - Sw_.transpose() * D_dq_PA_.template bottomRows<3>();

    }


    //! Leaf Terms 02 Declaration.
    //!------------------------------------------------------------------------------!//
    template<typename Vector3Type, typename Matrix3Type, typename Vector6Type, typename Matrix6Type, typename D_Matrix6Type>
    inline static void
    runLeaf02(bool P_z_,
              typename Eigen::MatrixBase<Vector3Type> & Sw_,
              typename Eigen::MatrixBase<Vector3Type> & P_,
              typename Eigen::MatrixBase<Matrix3Type> & R_,
              typename Eigen::MatrixBase<Vector6Type> & P_A_i_,
              typename Eigen::MatrixBase<Vector6Type> & P_a_,
              typename Eigen::MatrixBase<Matrix6Type> & M_a_,
              typename Eigen::MatrixBase<Matrix6Type> & M_A_j_,
              typename Eigen::MatrixBase<D_Matrix6Type> & D_M_A_j_) {

      //! Inertial back projection.
      //!------------------------------------------------------------------------------!//
      //      typename GEOMBD_EIGEN_PLAIN_TYPE(Matrix6Type) Mtmp_;
      Matrix6Type Mtmp_;

      //! Back projection of M_a.
      //!------------------------------------------------------------------------------!//
      Mat6ProjRxyz(P_z_, P_.derived(), R_.derived(), M_a_.derived(), Mtmp_.derived());
      M_A_j_ += Mtmp_;

      //! Inertial back-projection differentiation.
      //!------------------------------------------------------------------------------!//
      //! M_a differentiation -> Single-motion-revolute screws enable skew-symmetric properties.
      //!------------------------------------------------------------------------------!//
      Matrix6Type M_a_S = Matrix6Type::Zero();  Matrix6Type M_a_S_;
      typename Matrix3Type::PlainObject skew_S = Skew(Sw_);

      M_a_S.template topRows<3>()    = skew_S * M_a_.template topRows<3>();  // mimicking ad_dual(-Sxyz)*Ma effect
      M_a_S.template bottomRows<3>() = skew_S * M_a_.template bottomRows<3>();

      M_a_S_.noalias() = M_a_S + M_a_S.transpose();   // mimicking ad_dual(-Sxyz)*Ma - Ma*ad(Sxyz) effect

      //! Transform M_a_S_ since it is now symmetric.
      //!------------------------------------------------------------------------------!//
      Mat6ProjRxyz(P_z_, P_.derived(), R_.derived(), M_a_S_.derived(), M_a_S.derived()); /// this will return as MtmpTop

      D_M_A_j_ = M_a_S;

      P_A_i_.template segment<3>(0) = Sw_.cross(P_a_.template segment<3>(0));  // mimicking ad_dual(-Sxyz)*Pa effect
      P_A_i_.template segment<3>(3) = Sw_.cross(P_a_.template segment<3>(3));

    }


    //! Spatial Acceleration 01 Declaration.
    //!------------------------------------------------------------------------------!//
    template<typename Vector3Type, typename Vector6Type>
    inline static void
    runAccel01(typename Eigen::MatrixBase<Vector3Type> & Sw_,
               typename Eigen::MatrixBase<Vector6Type> & AdAj_,
               typename Eigen::MatrixBase<Vector6Type> & Aa_) {

      //! Mimicking effect -ad(Sxyz)*Ad*Aj
      //!------------------------------------------------------------------------------!//
      AdAj_.template segment<3>(0) = - Sw_.cross( Aa_.template segment<3>(0));
      AdAj_.template segment<3>(3) = - Sw_.cross( Aa_.template segment<3>(3));

    }


    //! Spatial Acceleration 02 Declaration.
    //!------------------------------------------------------------------------------!//
    template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename RowVectorXType, typename D_Vector6Type>
    inline static void
    runAccel02(ScalarType* ddq_,
               typename Eigen::MatrixBase<Vector3Type> & Sw_,
               typename Eigen::MatrixBase<Vector6Type> & A_,
               typename Eigen::MatrixBase<D_Vector6Type> & D_q_A_,
               typename Eigen::MatrixBase<D_Vector6Type> & D_dq_A_,
               typename Eigen::MatrixBase<RowVectorXType> & D_q_ddq_,
               typename Eigen::MatrixBase<RowVectorXType> & D_dq_ddq_) {

      //! Update spatial acceleration and its derivatives D_q_A_ and D_dq_A_.
      //!------------------------------------------------------------------------------!//
      A_.template segment<3>(3) += Sw_*(*ddq_);
      //!------------------------------------------------------------------------------!//
      D_q_A_.template bottomRows<3>() += Sw_*D_q_ddq_;

      D_dq_A_.template bottomRows<3>() += Sw_*D_dq_ddq_;

    }


    //! Implementation for acceleration expressions at root.
    template<typename ScalarType, typename Vector3Type, typename Matrix3Type, typename Vector6Type,
             typename D_Vector6Type, typename RowVectorXType, typename MatrixXType>
    inline static void
    runAccelRoot(ScalarType u,
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
      //! D_Acc_a = -ad(Sxyz)*Ad(G[Sxyz])*g = 0

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

  };

} // end namespace geo

#endif // GEOMBD_DIFFERENTIATION_JOINT_TYPE_REVOLUTE_XYZ_HXX
