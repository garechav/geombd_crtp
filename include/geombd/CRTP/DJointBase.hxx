/**
 *    \file include/geombd/CRTP/DJointBase.hxx
 *    \author Alvaro Paz, Gustavo Arechavaleta
 *    \version 1.0
 *    \date 2021
 *
 *    Class to implement the CRTP base (interface)
 *    Copyright (c) 2021 Cinvestav
 *    This library is distributed under the MIT License.
 */

#ifdef EIGEN_VECTORIZE

#endif

#ifndef GEOMBD_JOINT_BASE_DIFFERENTIATION_CRTP_HXX
#define GEOMBD_JOINT_BASE_DIFFERENTIATION_CRTP_HXX

//! Include Eigen Library
//!--------------------------------------------------------------------------------!//
#define EIGEN_NO_DEBUG
#define EIGEN_MPL2_ONLY
#define EIGEN_UNROLLING_LIMIT 30

#include "Eigen/Core"
//#include "Eigen/../unsupported/Eigen/KroneckerProduct"
//!--------------------------------------------------------------------------------!//

namespace geo {

  //! Joint Type Base
  //!------------------------------------------------------------------------------!//
  template<typename Derived>
  struct D_CRTPInterface{
  public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW


    //! Recurring Pattern for the Forward Kinematics.
    //!------------------------------------------------------------------------------!//
    template<typename ScalarType, typename Vector3Type, typename Matrix3Type>
    EIGEN_ALWAYS_INLINE void
    D_FwdKin(const ScalarType & qi,
             const Eigen::MatrixBase<Vector3Type> & S_,
             typename Eigen::MatrixBase<Matrix3Type> & R) {
      static_cast<Derived*>(this)->runD_FK(qi, S_, R);
      // static_cast<Derived&>(*this).runD_FK(qi, S_, R);
    }


    //! Recurring Pattern for the TCP 01.
    //!------------------------------------------------------------------------------!//
    template<typename ScalarType, typename Vector3Type, typename Vector6Type>
    EIGEN_ALWAYS_INLINE void
    D_TCP01(const ScalarType & vi,
            typename Eigen::MatrixBase<Vector3Type> & Sw_,
            typename Eigen::MatrixBase<Vector6Type> & S_,
            typename Eigen::MatrixBase<Vector6Type> & C_) {
      static_cast<Derived*>(this)->runD_TCP01(vi, Sw_, S_, C_);
      // static_cast<Derived&>(*this).runD_TCP01(vi, S_, C_);
    }


    //! Recurring Pattern for the TCP 02.
    //!------------------------------------------------------------------------------!//
    template<typename Vector3Type, typename D_Vector6Type>
    EIGEN_ALWAYS_INLINE void
    D_TCP02(typename Eigen::MatrixBase<Vector3Type> & Sw_,
            typename Eigen::MatrixBase<D_Vector6Type> & D_q_c_,
            typename Eigen::MatrixBase<D_Vector6Type> & D_dq_c_,
            typename Eigen::MatrixBase<D_Vector6Type> & D_q_c_aux_,
            typename Eigen::MatrixBase<D_Vector6Type> & D_dq_c_aux_) {
      static_cast<Derived*>(this)->runD_TCP02(Sw_.derived(), D_q_c_.derived(), D_dq_c_.derived(), D_q_c_aux_.derived(), D_dq_c_aux_.derived());
    }


    //! Recurring Pattern for the TCP root.
    template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename Matrix6Type, typename D_Vector6Type>
    EIGEN_ALWAYS_INLINE void
    D_TCProot(const ScalarType & vi,
              typename Eigen::MatrixBase<Vector3Type> & S_,
              typename Eigen::MatrixBase<Vector6Type> & T_,
              typename Eigen::MatrixBase<Vector6Type> & p_,
              typename Eigen::MatrixBase<Matrix6Type> & M_,
              typename Eigen::MatrixBase<D_Vector6Type> & D_dq_p_) {
      static_cast<Derived*>(this)->runD_TCProot(vi, S_.derived(), T_.derived(), p_.derived(), M_.derived(), D_dq_p_.derived());
    }


    //! Recurring Pattern for the Inertial Terms 01.
    template<typename Vector3Type, typename Vector6Type, typename Matrix6Type, typename VectorXType, typename D_Matrix6Type>
    EIGEN_ALWAYS_INLINE void
    Inertia01(typename Eigen::MatrixBase<Vector3Type> & Sw_,
              typename Eigen::MatrixBase<Vector6Type> & U_,
              typename Eigen::MatrixBase<Matrix6Type> & M_A_,
              typename Eigen::MatrixBase<VectorXType> & D_U_v_,
              typename Eigen::MatrixBase<D_Matrix6Type> & D_M_A_i_) {
      static_cast<Derived*>(this)->runInertia01(Sw_.derived(), U_.derived(), M_A_.derived(), D_U_v_.derived(), D_M_A_i_.derived());
    }


    //! Recurring Pattern for the Inertial Terms 02.
    template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename RowVectorXType, typename D_Vector6Type>
    EIGEN_ALWAYS_INLINE void
    Inertia02(ScalarType & invD_,
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
      static_cast<Derived*>(this)->runInertia02(invD_, u_, Sw_.derived(), U_.derived(), P_A_.derived(), D_invD_.derived(), D_q_u_.derived(),
                                                D_dq_u_.derived(), D_U_h_.derived(), D_q_PA_.derived(), D_dq_PA_.derived());
    }


    //! Recurring Pattern for the Inertial Terms 03.
    template<typename IndexType, typename Vector3Type, typename Matrix3Type, typename Vector6Type, typename Matrix6Type, typename D_Matrix6Type>
    EIGEN_ALWAYS_INLINE void
    Inertia03(bool P_z_,
              IndexType nS_,
              typename Eigen::MatrixBase<Vector3Type> & Sw_,
              typename Eigen::MatrixBase<Vector3Type> & P_,
              typename Eigen::MatrixBase<Matrix3Type> & R_,
              typename Eigen::MatrixBase<Vector6Type> & P_a_,
              typename Eigen::MatrixBase<Vector6Type> & P_A_i_,
              typename Eigen::MatrixBase<Matrix6Type> & M_a_,
              typename Eigen::MatrixBase<Matrix6Type> & Mtmp_,
              typename Eigen::MatrixBase<D_Matrix6Type> & D_M_A_i_) {
      static_cast<Derived*>(this)->runInertia03(P_z_, nS_, Sw_.derived(), P_.derived(), R_.derived(), P_a_.derived(), P_A_i_.derived(),
                                                M_a_.derived(), Mtmp_.derived(), D_M_A_i_.derived());
    }


    //! Recurring Pattern for the Leaf body 01.
    template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename Matrix6Type, typename RowVectorXType, typename D_Vector6Type>
    EIGEN_ALWAYS_INLINE void
    Leaf01(ScalarType & invD_,
           ScalarType & u_,
           typename Eigen::MatrixBase<Vector3Type> & Sw_,
           typename Eigen::MatrixBase<Vector6Type> & U_,
           typename Eigen::MatrixBase<Vector6Type> & P_A_,
           typename Eigen::MatrixBase<Matrix6Type> & M_A_,
           typename Eigen::MatrixBase<RowVectorXType> & D_q_u_,
           typename Eigen::MatrixBase<RowVectorXType> & D_dq_u_,
           typename Eigen::MatrixBase<D_Vector6Type> & D_q_PA_,
           typename Eigen::MatrixBase<D_Vector6Type> & D_dq_PA_) {
      static_cast<Derived*>(this)->runLeaf01(invD_, u_, Sw_.derived(), U_.derived(), P_A_.derived(), M_A_.derived(),
                                             D_q_u_.derived(), D_dq_u_.derived(), D_q_PA_.derived(), D_dq_PA_.derived());
    }


    //! Recurring Pattern for the Leaf body 02.
    template<typename Vector3Type, typename Matrix3Type, typename Vector6Type, typename Matrix6Type, typename D_Matrix6Type>
    EIGEN_ALWAYS_INLINE void
    Leaf02(bool P_z_,
           typename Eigen::MatrixBase<Vector3Type> & Sw_,
           typename Eigen::MatrixBase<Vector3Type> & P_,
           typename Eigen::MatrixBase<Matrix3Type> & R_,
           typename Eigen::MatrixBase<Vector6Type> & P_A_i_,
           typename Eigen::MatrixBase<Vector6Type> & P_a_,
           typename Eigen::MatrixBase<Matrix6Type> & M_a_,
           typename Eigen::MatrixBase<Matrix6Type> & M_A_j_,
           typename Eigen::MatrixBase<D_Matrix6Type> & D_M_A_j_) {
      static_cast<Derived*>(this)->runLeaf02(P_z_, Sw_.derived(), P_.derived(), R_.derived(), P_A_i_.derived(), P_a_.derived(),
                                             M_a_.derived(), M_A_j_.derived(), D_M_A_j_.derived());
    }


    //! Recurring Pattern for the spatial acceleration 01.
    template<typename Vector3Type, typename Vector6Type>
    EIGEN_ALWAYS_INLINE void
    Accel01(typename Eigen::MatrixBase<Vector3Type> & Sw_,
            typename Eigen::MatrixBase<Vector6Type> & AdAj_,
            typename Eigen::MatrixBase<Vector6Type> & Aa_) {
      static_cast<Derived*>(this)->runAccel01(Sw_.derived(), AdAj_.derived(), Aa_.derived());
    }


    //! Recurring Pattern for the spatial acceleration 02.
    template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename RowVectorXType, typename D_Vector6Type>
    EIGEN_ALWAYS_INLINE void
    Accel02(ScalarType* ddq_,
            Eigen::MatrixBase<Vector3Type> & Sw_,
            Eigen::MatrixBase<Vector6Type> & A_,
            Eigen::MatrixBase<D_Vector6Type> & D_q_A_,
            Eigen::MatrixBase<D_Vector6Type> & D_dq_A_,
            Eigen::MatrixBase<RowVectorXType> & D_q_ddq_,
            Eigen::MatrixBase<RowVectorXType> & D_dq_ddq_) {
      static_cast<Derived*>(this)->runAccel02(ddq_, Sw_.derived(), A_.derived(), D_q_A_.derived(), D_dq_A_.derived(), D_q_ddq_.derived(), D_dq_ddq_.derived());
    }


    //! Recurring Pattern for the spatial acceleration 02 at root.
    template<typename ScalarType, typename Vector3Type, typename Matrix3Type, typename Vector6Type,
             typename D_Vector6Type, typename RowVectorXType, typename MatrixXType>
    EIGEN_ALWAYS_INLINE void
    AccelRoot(ScalarType u,
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
      static_cast<Derived*>(this)->runAccelRoot(u, iD, ddq, S, P_r, R_r, U_r, Acc_i_r, D_U_h_, D_invD_,
                                                D_q_u_, D_dq_u_, D_q_A_, D_dq_A_, D_ddq_);
    }


  };


}

#include "DJointDerived/DJointDerived.hpp"
#include "DJointDerived/DVisitors.hxx"

#endif // GEOMBD_JOINT_BASE_DIFFERENTIATION_CRTP_HXX
