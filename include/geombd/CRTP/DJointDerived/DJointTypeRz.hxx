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

#ifndef GEOMBD_DIFFERENTIATION_JOINT_TYPE_REVOLUTE_Z_HXX
#define GEOMBD_DIFFERENTIATION_JOINT_TYPE_REVOLUTE_Z_HXX

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


  //! Derived -> Differentiation of Joint Type Revolute On Z Axis
  //!------------------------------------------------------------------------------!//
  struct D_JointTypeRz : public D_CRTPInterface<D_JointTypeRz>
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


  };

} // end namespace geo

#endif // GEOMBD_DIFFERENTIATION_JOINT_TYPE_REVOLUTE_Z_HXX
