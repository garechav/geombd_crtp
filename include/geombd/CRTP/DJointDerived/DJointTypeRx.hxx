/**
 *    \file include/geombd/CRTP/DJointDerived/DJointTypeRx.hxx
 *    \author Alvaro Paz, Gustavo Arechavaleta
 *    \version 1.0
 *    \date 2021
 *
 *    Derived class for Rx joint type
 *    Copyright (c) 2021 Cinvestav
 *    This library is distributed under the MIT License.
 */

#ifndef GEOMBD_DIFFERENTIATION_JOINT_TYPE_REVOLUTE_X_HXX
#define GEOMBD_DIFFERENTIATION_JOINT_TYPE_REVOLUTE_X_HXX

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


  //! Derived -> Differentiation of Joint Type Revolute On X Axis
  //!------------------------------------------------------------------------------!//
  struct D_JointTypeRx : public D_CRTPInterface<D_JointTypeRx>
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
      SINCOS<ScalarType>(qi, &sqi, &cqi);

      R.coeffRef(1,1) = cqi;  R.coeffRef(1,2) = -sqi;
      R.coeffRef(2,1) = sqi;  R.coeffRef(2,2) =  cqi;
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
      S_dw.coeffRef(0) += vi;
      //!------------------------------------------------------------------------------!//      
      C_up.coeffRef(1) =  S_up.coeff(2);  // -ad(Sx) effect
      C_up.coeffRef(2) = -S_up.coeff(1);

      C_dw.coeffRef(1) =  S_dw.coeff(2);
      C_dw.coeffRef(2) = -S_dw.coeff(1);
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
      //!------------------------------------------------------------------------------!//
      D_q_c_.noalias() = apply_adSx(D_q_c_aux_);
      //!------------------------------------------------------------------------------!//
      D_dq_c_.noalias() = apply_adSx(D_dq_c_aux_);
      //!------------------------------------------------------------------------------!//
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

      T_.coeffRef(3) = vi;

      ScalarType vi2 = vi*vi;

      p_.coeffRef(1) = -M_.coeff(2,3);  // adDual(Sx) effect
      p_.coeffRef(2) =  M_.coeff(1,3);
      p_.coeffRef(4) = -M_.coeff(5,3);
      p_.coeffRef(5) =  M_.coeff(4,3);

      //! P bias differentiation
      //!------------------------------------------------------------------------------!//
      D_dq_p_.noalias() = p_;  // adDual(Sx)*M*S

      Vector6Type mu;
      mu.noalias() = M_.template col(3)*vi;

      D_dq_p_.coeffRef(1) -= mu.coeff(2);  // -adBar(mu)*Sx effect
      D_dq_p_.coeffRef(2) += mu.coeff(1);
      D_dq_p_.coeffRef(4) -= mu.coeff(5);
      D_dq_p_.coeffRef(5) += mu.coeff(4);

      p_ *= vi2;

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
      U_ = M_A_.template col(3);
      //!------------------------------------------------------------------------------!//
      D_U_v_ = D_M_A_i_.template col(3);

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
      invD_ = 1 / U_.coeffRef(3);
      //!------------------------------------------------------------------------------!//
      D_invD_.noalias() = -pow(invD_,2)*D_U_h_.template row(3);

      //! Solve u and its partial derivative.
      //!------------------------------------------------------------------------------!//
      u_ -= P_A_.coeffRef(3);
      //!------------------------------------------------------------------------------!//
      D_q_u_  = -D_q_PA_.template row(3);
      D_dq_u_ = -D_dq_PA_.template row(3);

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
      Mat6ProjRx(P_z_, P_.derived(), R_.derived(), M_a_.derived(), Mtmp_.derived());

      //! Back projection of D_M_a.
      //!------------------------------------------------------------------------------!//
      D_Mat6ProjRx(P_z_, nS_, P_.derived(), R_.derived(), D_M_A_i_.derived());

      //! Creation of D_M_a_i.
      //! M_a differentiation -> Single-motion-revolute screws enable skew-symmetric properties.
      //!------------------------------------------------------------------------------!//
      Matrix6Type M_a_S, M_a_S_;
      M_a_S.setZero();
      M_a_S.template row(1) = -M_a_.template row(2);  // mimicking effect ad_dual(-Sx)*Ma
      M_a_S.template row(2) =  M_a_.template row(1);
      M_a_S.template row(4) = -M_a_.template row(5);
      M_a_S.template row(5) =  M_a_.template row(4);

      M_a_S_.noalias() = M_a_S + M_a_S.transpose();   // mimicking effect ad_dual(-Sx)*Ma - Ma*ad(Sx)

      //! Transform M_a_S_ since it is now symmetric.
      //! Back projection of the partial derivative D_M_a.
      //!------------------------------------------------------------------------------!//
      Mat6ProjRx(P_z_, P_.derived(), R_.derived(), M_a_S_.derived(), M_a_S.derived());

      D_M_A_i_.template topRows<6>() = M_a_S;

      P_A_i_ << 0, -P_a_.coeff(2), P_a_.coeff(1), 0, -P_a_.coeff(5), P_a_.coeff(4);

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
      U_ = M_A_.template col(3);  invD_ = 1 / U_.coeff(3);  u_ -= P_A_.coeff(3);

      //! Solve partial derivative of u as D_q_u and D_dq_u.
      //!------------------------------------------------------------------------------!//
      D_q_u_ = -D_q_PA_.template row(3);  D_dq_u_ = -D_dq_PA_.template row(3);

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
      Mat6ProjRx(P_z_, P_.derived(), R_.derived(), M_a_.derived(), Mtmp_.derived());
      M_A_j_ += Mtmp_;

      //! Inertial back-projection differentiation.
      //!------------------------------------------------------------------------------!//
      //! M_a differentiation -> Single-motion-revolute screws enable skew-symmetric properties.
      //!------------------------------------------------------------------------------!//
      Matrix6Type M_a_S, M_a_S_;
      M_a_S.setZero();
      M_a_S.template row(1) = -M_a_.template row(2);  // mimicking effect ad_dual(-Sx)*Ma
      M_a_S.template row(2) =  M_a_.template row(1);
      M_a_S.template row(4) = -M_a_.template row(5);
      M_a_S.template row(5) =  M_a_.template row(4);

      M_a_S_.noalias() = M_a_S + M_a_S.transpose();   // mimicking effect ad_dual(-Sx)*Ma - Ma*ad(Sx)

      //! Transform M_a_S_ since it is now symmetric.
      //!------------------------------------------------------------------------------!//
      Mat6ProjRx(P_z_, P_.derived(), R_.derived(), M_a_S_.derived(), M_a_S.derived()); /// this will return as MtmpTop

      D_M_A_j_ = M_a_S;

      P_A_i_ << 0, -P_a_.coeff(2), P_a_.coeff(1), 0, -P_a_.coeff(5), P_a_.coeff(4);

    }


    //! Spatial Acceleration 01 Declaration.
    //!------------------------------------------------------------------------------!//
    template<typename Vector3Type, typename Vector6Type>
    inline static void
    runAccel01(typename Eigen::MatrixBase<Vector3Type> & Sw_,
               typename Eigen::MatrixBase<Vector6Type> & AdAj_,
               typename Eigen::MatrixBase<Vector6Type> & Aa_) {

      //! Mimicking effect -ad(Sx)*Ad*Aj
      //!------------------------------------------------------------------------------!//
      AdAj_ << 0, Aa_.coeff(2), -Aa_.coeff(1), 0, Aa_.coeff(5), -Aa_.coeff(4);

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
      A_.coeffRef(3) += (*ddq_);
      //!------------------------------------------------------------------------------!//
      D_q_A_.template row(3) += D_q_ddq_;

      D_dq_A_.template row(3) += D_dq_ddq_;

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

      //! Acceleration bias and its derivative.
      //!------------------------------------------------------------------------------!//
      //! g = [ 0, 0, 9.81, 0, 0, 0 ]^T
      //! Acc_a = Ad(G[Sx])*g = [ 0, 9.81*sq, 9.81*cq, 0, 0, 0 ]^T
      //! D_Acc_a = -ad(Sx)*Ad(G)*g = [ 0, 9.81*cq, -9.81*sq, 0, 0, 0 ]^T

      ScalarType ddq_, sq, cq;
      sq = R_r.coeff(2,1);  cq = R_r.coeff(1,1);

      Acc_i_r.setZero();  D_q_A_.setZero();
      Acc_i_r.coeffRef(1) = 9.81*sq;  Acc_i_r.coeffRef(2) = 9.81*cq;
      D_q_A_.coeffRef(1) = 9.81*cq;   D_q_A_.coeffRef(2) = -9.81*sq;


      //! Joint acceleration and its derivatives D_q_ddq and D_dq_ddq.
      //!------------------------------------------------------------------------------!//
      Vector6Type & U_ = const_cast<Vector6Type &>(U_r.derived());
      ddq_ = u - U_.transpose()*Acc_i_r;
      (*ddq) = iD*ddq_;


      //!-------------------------------------------------------
      RowVectorXType D_q_ddq, D_dq_ddq;
      D_q_ddq.noalias() = iD*D_q_u_;  D_dq_ddq.noalias() = iD*D_dq_u_;
      ScalarType U_D_A = U_.transpose()*D_q_A_.template leftCols<1>();
      D_q_ddq.coeffRef(0) -= iD*U_D_A;
      D_q_ddq.segment(1,D_U_h_.cols()).noalias() += ddq_*D_invD_ - iD*Acc_i_r.transpose()*D_U_h_;


      //! Update spatial acceleration and its derivative.
      //!------------------------------------------------------------------------------!//
      Acc_i_r.coeffRef(3) = iD*ddq_;
      //!-------------------------------------------------------
      D_q_A_.template row(3) = D_q_ddq;  D_dq_A_.template row(3) = D_dq_ddq;


      //! Fill up D_ddq.
      //!------------------------------------------------------------------------------!//
      D_ddq_.row(0) << D_q_ddq, D_dq_ddq;

    }

  };

} // end namespace geo

#endif // GEOMBD_DIFFERENTIATION_JOINT_TYPE_REVOLUTE_X_HXX
