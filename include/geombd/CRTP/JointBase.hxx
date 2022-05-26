/**
 *    \file include/geombd/CRTP/JointBase.hxx
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

#ifndef GEOMBD_JOINT_BASE_CRTP_HXX
#define GEOMBD_JOINT_BASE_CRTP_HXX

//! Include Eigen Library
//!--------------------------------------------------------------------------------!//
#define EIGEN_NO_DEBUG
#include "Eigen/Core"
//!--------------------------------------------------------------------------------!//

namespace geoCRTP{

  //! Joint Type Base
  //!------------------------------------------------------------------------------!//
  template<typename Derived>
  struct CRTPInterface{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Recurring Pattern for the Forward Kinematics.
    template<typename ScalarType, typename Vector3Type, typename Matrix3Type>
    EIGEN_ALWAYS_INLINE void
    FwdKin(const ScalarType & qi,
           typename Eigen::MatrixBase<Vector3Type> & S,
           typename Eigen::MatrixBase<Matrix3Type> & R) {
      static_cast<Derived*>(this)->runFK(qi, S, R);
      // static_cast<Derived&>(*this).runFK(qi, S, R);
    }


    //! Recurring Pattern for Twist, C bias and P bias at root.
    template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename Matrix6Type>
    EIGEN_ALWAYS_INLINE void
    TCP_root(const ScalarType & vi,
             Eigen::MatrixBase<Vector3Type> & S_,
             Eigen::MatrixBase<Vector6Type> & S_i,
             Eigen::MatrixBase<Vector6Type> & p_,
             const Eigen::MatrixBase<Matrix6Type> & M_) {
      static_cast<Derived*>(this)->runTCP_root(vi, S_, S_i, p_.derived(), M_.derived());
    }


    //! Recurring Pattern for Twist, C bias and P bias.
    template<typename ScalarType, typename Matrix3Type, typename Matrix6Type,
             typename Vector3Type, typename Vector6Type>
    EIGEN_ALWAYS_INLINE void
    TwCbPb(bool zeroFlag,
           const ScalarType & vi,
           const Eigen::MatrixBase<Vector3Type> & S_,
           const Eigen::MatrixBase<Matrix3Type> & R_,
           const Eigen::MatrixBase<Vector3Type> & P_,
           const Eigen::MatrixBase<Vector6Type> & S_l,
           const Eigen::MatrixBase<Matrix6Type> & M_,
           Eigen::MatrixBase<Vector6Type> & S_i,
           Eigen::MatrixBase<Vector6Type> & c_,
           Eigen::MatrixBase<Vector6Type> & p_) {
      static_cast<Derived*>(this)->runTwCbPb(zeroFlag, vi, S_.derived(), R_.derived(), P_.derived(), S_l.derived(), M_.derived(),
                                             S_i.derived(), c_.derived(), p_.derived());
    }


    //! Recurring Pattern for U, u & invD expressions.
    template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename Matrix6Type>
    EIGEN_ALWAYS_INLINE void
    UuiD(ScalarType & u,
         ScalarType & iD,
         const ScalarType tau,
         const Eigen::MatrixBase<Vector3Type> & S_,
         Eigen::MatrixBase<Vector6Type> & U_r,
         const Eigen::MatrixBase<Vector6Type> & P_A_r,
         const Eigen::MatrixBase<Matrix6Type> & M_A_r) {
      static_cast<Derived*>(this)->runUuiD(u, iD, tau, S_, U_r, P_A_r, M_A_r);
    }


    //! Recurring pattern for preparing inertial expressions.
    template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename Matrix6Type>
    EIGEN_ALWAYS_INLINE void
    PreIner(ScalarType u,
            ScalarType iD,
            const Eigen::MatrixBase<Vector3Type> & S_,
            const Eigen::MatrixBase<Vector6Type> & U_r,
            const Eigen::MatrixBase<Vector6Type> & c_r,
            Eigen::MatrixBase<Vector6Type> & P_a_r,
            Eigen::MatrixBase<Matrix6Type> & M_a_r,
            const Eigen::MatrixBase<Vector6Type> & P_A_r,
            const Eigen::MatrixBase<Matrix6Type> & M_A_r) {
      static_cast<Derived*>(this)->runPreIner(u, iD, S_, U_r, c_r, P_a_r, M_a_r, P_A_r, M_A_r);
    }


    //! Recurring pattern for back-projection of inertial expressions.
    template<typename Vector3Type, typename Matrix3Type, typename Vector6Type, typename Matrix6Type>
    EIGEN_ALWAYS_INLINE void
    InerProj(bool P_z,
             const Eigen::MatrixBase<Vector3Type> & S_,
             const Eigen::MatrixBase<Vector3Type> & P_r,
             const Eigen::MatrixBase<Matrix3Type> & R_r,
             const Eigen::MatrixBase<Vector6Type> & P_a_r,
             Eigen::MatrixBase<Vector6Type> & P_A_r,
             const Eigen::MatrixBase<Matrix6Type> & M_a_r,
             Eigen::MatrixBase<Matrix6Type> & M_A_r) {
      static_cast<Derived*>(this)->runInerProj(P_z, S_, P_r, R_r, P_a_r, P_A_r, M_a_r, M_A_r);
    }


    //! Recurring pattern for acceleration expression.
    template<typename ScalarType, typename Vector3Type, typename Matrix3Type, typename Vector6Type>
    EIGEN_ALWAYS_INLINE void
    Accel(bool zeroFlag,
          ScalarType u,
          ScalarType iD,
          ScalarType* ddq,
          const Eigen::MatrixBase<Vector3Type> & S_,
          const Eigen::MatrixBase<Vector3Type> & P_r,
          const Eigen::MatrixBase<Matrix3Type> & R_r,
          const Eigen::MatrixBase<Vector6Type> & c_r,
          const Eigen::MatrixBase<Vector6Type> & U_r,
          Eigen::MatrixBase<Vector6Type> & Acc_i_r,
          const Eigen::MatrixBase<Vector6Type> & Acc_j_r) {
      static_cast<Derived*>(this)->runAccel(zeroFlag, u, iD, ddq, S_, P_r, R_r, c_r, U_r, Acc_i_r, Acc_j_r);
    }


    //! Recurring pattern for acceleration expression at root.
    template<typename ScalarType, typename Vector3Type, typename Matrix3Type, typename Vector6Type>
    EIGEN_ALWAYS_INLINE void
    AccelRoot(ScalarType u,
              ScalarType iD,
              ScalarType* ddq,
              const Eigen::MatrixBase<Vector3Type> & S_,
              const Eigen::MatrixBase<Vector3Type> & P_r,
              const Eigen::MatrixBase<Matrix3Type> & R_r,
              const Eigen::MatrixBase<Vector6Type> & U_r,
              Eigen::MatrixBase<Vector6Type> & Acc_i_r) {
      static_cast<Derived*>(this)->runAccelRoot(u, iD, ddq, S_, P_r, R_r, U_r, Acc_i_r);
    }


  };


}

#include "JointDerived/JointDerived.hpp"
#include "JointDerived/Visitors.hxx"

#endif // GEOMBD_JOINT_BASE_CRTP_HXX
