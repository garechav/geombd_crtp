/**
 *    \file include/geombd/CRTP/DJointDerived/DJointTypeRz.hxx
 *    \author Alvaro Paz, Gustavo Arechavaleta
 *    \version 1.0
 *    \date 2021
 *
 *    Derived class for Rz joint type
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
      SINCOS<ScalarType>(qi, &sqi, &cqi);

      R.coeffRef(0,0) = cqi;  R.coeffRef(0,1) = -sqi;
      R.coeffRef(1,0) = sqi;  R.coeffRef(1,1) =  cqi;
    }


    //! TCP01 Declaration
    //!------------------------------------------------------------------------------!//
    template<typename ScalarType, typename Vector6Type>
    inline static void
    runD_TCP01(const ScalarType & vi,
                              typename Eigen::MatrixBase<Vector6Type> & S_,
                              typename Eigen::MatrixBase<Vector6Type> & C_) {
      //!------------------------------------------------------------------------------!//
      typedef Eigen::Block<Vector6Type,3,1> Segment3;

      Segment3 S_up = S_.template segment<3>(0);
      Segment3 S_dw = S_.template segment<3>(3);

      Segment3 C_up = C_.template segment<3>(0);
      Segment3 C_dw = C_.template segment<3>(3);
      //!------------------------------------------------------------------------------!//
      S_dw.coeffRef(2) += vi;
      //!------------------------------------------------------------------------------!//
      C_up.coeffRef(0) =  S_up.coeff(1);  // -ad(Sz) effect
      C_up.coeffRef(1) = -S_up.coeff(0);

      C_dw.coeffRef(0) =  S_dw.coeff(1);
      C_dw.coeffRef(1) = -S_dw.coeff(0);
      //!------------------------------------------------------------------------------!//
    }


    //! TCP02 Declaration
    //!------------------------------------------------------------------------------!//
    template<typename D_Vector6Type>
    inline static void
    runD_TCP02(typename Eigen::MatrixBase<D_Vector6Type> & D_q_c_,
               typename Eigen::MatrixBase<D_Vector6Type> & D_dq_c_,
               typename Eigen::MatrixBase<D_Vector6Type> & D_q_c_aux_,
               typename Eigen::MatrixBase<D_Vector6Type> & D_dq_c_aux_) {
      //!------------------------------------------------------------------------------!//
      D_q_c_.noalias() = apply_adSz(D_q_c_aux_);
      //!------------------------------------------------------------------------------!//
      D_dq_c_.noalias() = apply_adSz(D_dq_c_aux_);
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

      T_.template tail<1>()(0) = vi;

      ScalarType vi2 = vi*vi;

      p_.coeffRef(0) = -M_.coeff(1,5);  // adDual(Sz) effect
      p_.coeffRef(1) =  M_.coeff(0,5);
      p_.coeffRef(3) = -M_.coeff(4,5);
      p_.coeffRef(4) =  M_.coeff(3,5);

      //! P bias differentiation
      //!------------------------------------------------------------------------------!//
      D_dq_p_.noalias() = p_;  // adDual(Sz)*M*S

      Vector6Type mu;
      mu.noalias() = M_.template rightCols<1>()*vi;

      D_dq_p_.coeffRef(3) -= mu.coeff(4);  // adBar(mu) effect
      D_dq_p_.coeffRef(4) += mu.coeff(3);

      p_ *= vi2;

    }


    //! Inertial Terms 01 Declaration.
    //!------------------------------------------------------------------------------!//
    template<typename Vector6Type, typename Matrix6Type, typename VectorXType, typename D_Matrix6Type>
    inline static void
    runInertia01(typename Eigen::MatrixBase<Vector6Type> & U_,
                 typename Eigen::MatrixBase<Matrix6Type> & M_A_,
                 typename Eigen::MatrixBase<VectorXType> & D_U_v_,
                 typename Eigen::MatrixBase<D_Matrix6Type> & D_M_A_i_) {

      //! Solve U and its partial derivative.
      //!------------------------------------------------------------------------------!//
      U_ = M_A_.template rightCols<1>();
      //!------------------------------------------------------------------------------!//
      D_U_v_ = D_M_A_i_.template rightCols<1>();

    }


    //! Inertial Terms 02 Declaration.
    //!------------------------------------------------------------------------------!//
    template<typename ScalarType, typename Vector6Type, typename RowVectorXType, typename D_Vector6Type>
    inline static void
    runInertia02(ScalarType & invD_,
                 ScalarType & u_,
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
      invD_ = 1 / U_.coeffRef(5);
      //!------------------------------------------------------------------------------!//
      D_invD_.noalias() = -pow(invD_,2)*D_U_h_.template bottomRows<1>();

      //! Solve u and its partial derivative.
      //!------------------------------------------------------------------------------!//
      u_ -= P_A_.coeffRef(5);
      //!------------------------------------------------------------------------------!//
      D_q_u_  = -D_q_PA_.template bottomRows<1>();
      D_dq_u_ = -D_dq_PA_.template bottomRows<1>();

    }


    //! Inertial Terms 03 Declaration.
    //!------------------------------------------------------------------------------!//
    template<typename IndexType, typename Vector3Type, typename Matrix3Type, typename Vector6Type, typename Matrix6Type, typename D_Matrix6Type>
    inline static void
    runInertia03(bool P_z_,
                 IndexType nS_,
                 typename Eigen::MatrixBase<Vector3Type> & P_,
                 typename Eigen::MatrixBase<Matrix3Type> & R_,
                 typename Eigen::MatrixBase<Vector6Type> & P_a_,
                 typename Eigen::MatrixBase<Vector6Type> & P_A_i_,
                 typename Eigen::MatrixBase<Matrix6Type> & M_a_,
                 typename Eigen::MatrixBase<Matrix6Type> & Mtmp_,
                 typename Eigen::MatrixBase<D_Matrix6Type> & D_M_A_i_) {

      //! Back projection of M_a.
      //!------------------------------------------------------------------------------!//
      Mat6ProjRz(P_z_, P_.derived(), R_.derived(), M_a_.derived(), Mtmp_.derived());

      //! Back projection of D_M_a.
      //!------------------------------------------------------------------------------!//
      D_Mat6ProjRz(P_z_, nS_, P_.derived(), R_.derived(), D_M_A_i_.derived());

      //! Creation of D_M_a_i.
      //! M_a differentiation -> Single-motion-revolute screws enable skew-symmetric properties.
      //!------------------------------------------------------------------------------!//
      Matrix6Type M_a_S, M_a_S_;
      M_a_S.setZero();
      M_a_S.template row(0) = -M_a_.template row(1);  // mimicking effect ad_dual(-Sz)*Ma
      M_a_S.template row(1) =  M_a_.template row(0);
      M_a_S.template row(3) = -M_a_.template row(4);
      M_a_S.template row(4) =  M_a_.template row(3);

      M_a_S_.noalias() = M_a_S + M_a_S.transpose();   // mimicking effect ad_dual(-Sz)*Ma - Ma*ad(Sz)

      //! Transform M_a_S_ since it is now symmetric.
      //! Back projection of the partial derivative D_M_a.
      //!------------------------------------------------------------------------------!//
      Mat6ProjRz(P_z_, P_.derived(), R_.derived(), M_a_S_.derived(), M_a_S.derived()); /// this will return as MtmpTop

      D_M_A_i_.template topRows<6>() = M_a_S;

      P_A_i_ << -P_a_.coeff(1), P_a_.coeff(0), 0, -P_a_.coeff(4), P_a_.coeff(3), 0;

    }


  };

} // end namespace geo

#endif // GEOMBD_DIFFERENTIATION_JOINT_TYPE_REVOLUTE_XYZ_HXX
