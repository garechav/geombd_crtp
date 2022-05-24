/**
 *    \file include/geombd/CRTP/DJointDerived/DVisitors.hxx
 *    \author Alvaro Paz, Gustavo Arechavaleta
 *    \version 1.0
 *    \date 2021
 *
 *    Classes for visitors
 *    Copyright (c) 2021 Cinvestav
 *    This library is distributed under the MIT License.
 */

#ifndef GEOMBD_DIFFERENTIATION_VISITORS_HXX
#define GEOMBD_DIFFERENTIATION_VISITORS_HXX

#define EIGEN_NO_DEBUG
#define EIGEN_MPL2_ONLY
#define EIGEN_UNROLLING_LIMIT 30

#include "Eigen/Core"

//#include <boost/foreach.hpp>
//#include <boost/shared_ptr.hpp>
#include <boost/variant.hpp>

namespace geo{

  //! Forward Kinematics Visitor
  //!------------------------------------------------------------------------------!//
  template<typename ScalarType, typename Vector3Type, typename Matrix3Type>
  class D_FwdKin_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    D_FwdKin_visitor( ) : R(nullptr), S(nullptr) {}

    ScalarType qi;
    Vector3Type* S;
    Matrix3Type* R;

    template<typename Derived>
    int operator()(D_CRTPInterface<Derived> & BaseType) const {
      BaseType.D_FwdKin(qi, (*S).derived(), (*R).derived());
      return 0;
    }

  };


  //! TCP 01 Visitor
  //!------------------------------------------------------------------------------!//
  template<typename ScalarType, typename Vector6Type>
  class D_TCP01_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    D_TCP01_visitor( ) : S_(nullptr), C_(nullptr) {}

    ScalarType vi;
    Vector6Type* S_;
    Vector6Type* C_;

    template<typename Derived>
    int operator()(D_CRTPInterface<Derived> & BaseType) const {
      BaseType.D_TCP01(vi, (*S_).derived(), (*C_).derived());
      return 0;
    }

  };


  //! TCP 02 Visitor
  //!------------------------------------------------------------------------------!//
  template<typename D_Vector6Type>
  class D_TCP02_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    D_TCP02_visitor( ) : D_q_c_(nullptr), D_dq_c_(nullptr), D_q_c_aux_(nullptr), D_dq_c_aux_(nullptr) {}

    D_Vector6Type* D_q_c_;    D_Vector6Type* D_q_c_aux_;
    D_Vector6Type* D_dq_c_;   D_Vector6Type* D_dq_c_aux_;

    template<typename Derived>
    int operator()(D_CRTPInterface<Derived> & BaseType) const {
      BaseType.D_TCP02((*D_q_c_).derived(), (*D_dq_c_).derived(), (*D_q_c_aux_).derived(), (*D_dq_c_aux_).derived());
      return 0;
    }

  };


  //! TCP root Visitor
  //!------------------------------------------------------------------------------!//
  template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename Matrix6Type, typename D_Vector6Type>
  class D_TCProot_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor
    D_TCProot_visitor( ) : S_(nullptr), T_(nullptr), p_(nullptr), M_(nullptr), D_dq_p_(nullptr) {}

    //! Members
    ScalarType vi;  Vector3Type* S_;  Vector6Type* T_;  Vector6Type* p_;  Matrix6Type* M_;  D_Vector6Type* D_dq_p_;

    //! Method -> member function
    template<typename Derived>
    int operator()(D_CRTPInterface<Derived> & BaseType) const {
      BaseType.D_TCProot(vi, (*S_).derived(), (*T_).derived(), (*p_).derived(), (*M_).derived(), (*D_dq_p_).derived());
      return 0;
    }

  };


  //! Inertia 01 Visitor
  //!------------------------------------------------------------------------------!//
  template<typename Vector6Type, typename Matrix6Type, typename VectorXType, typename D_Matrix6Type>
  class Inertia01_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor
    Inertia01_visitor( ) : U_(nullptr), M_A_(nullptr), D_U_v_(nullptr), D_M_A_i_(nullptr) {}

    //! Members
    Vector6Type* U_;  Matrix6Type* M_A_;  VectorXType* D_U_v_;  D_Matrix6Type* D_M_A_i_;

    //! Method -> member function
    template<typename Derived>
    int operator()(D_CRTPInterface<Derived> & BaseType) const {
      BaseType.Inertia01((*U_).derived(), (*M_A_).derived(), (*D_U_v_).derived(), (*D_M_A_i_).derived());
      return 0;
    }

  };


  //! Inertia 02 Visitor
  //!------------------------------------------------------------------------------!//
  template<typename ScalarType, typename Vector6Type, typename RowVectorXType, typename D_Vector6Type>
  class Inertia02_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor
    Inertia02_visitor( ) : U_(nullptr), P_A_(nullptr), D_invD_(nullptr), D_q_u_(nullptr), D_dq_u_(nullptr),
      D_U_h_(nullptr), D_q_PA_(nullptr), D_dq_PA_(nullptr) {}

    //! Members
    ScalarType* invD_;  ScalarType* u_;  Vector6Type* U_;  Vector6Type* P_A_;
    RowVectorXType* D_invD_;  RowVectorXType* D_q_u_;  RowVectorXType* D_dq_u_;
    D_Vector6Type* D_U_h_;  D_Vector6Type* D_q_PA_;  D_Vector6Type* D_dq_PA_;

    //! Method -> member function
    template<typename Derived>
    int operator()(D_CRTPInterface<Derived> & BaseType) const {
      ScalarType & _u_ = (*u_);  ScalarType & _invD_ = (*invD_);
      BaseType.Inertia02(_invD_, _u_, (*U_).derived(), (*P_A_).derived(), (*D_invD_).derived(), (*D_q_u_).derived(),
                         (*D_dq_u_).derived(), (*D_U_h_).derived(), (*D_q_PA_).derived(), (*D_dq_PA_).derived());
      return 0;
    }

  };


  //! Inertia 03 Visitor
  //!------------------------------------------------------------------------------!//
  template<typename IndexType, typename Vector3Type, typename Matrix3Type, typename Vector6Type, typename Matrix6Type, typename D_Matrix6Type>
  class Inertia03_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor
    Inertia03_visitor( ) : P_(nullptr), R_(nullptr), P_a_(nullptr), P_A_i_(nullptr), M_a_(nullptr),
      Mtmp_(nullptr), D_M_A_i_(nullptr) {}

    //! Members
    bool P_z_;  IndexType nS_;  Vector3Type* P_;  Matrix3Type* R_;
    Vector6Type* P_a_;  Vector6Type* P_A_i_;  Matrix6Type* M_a_;  Matrix6Type* Mtmp_;  D_Matrix6Type* D_M_A_i_;

    //! Method -> member function
    template<typename Derived>
    int operator()(D_CRTPInterface<Derived> & BaseType) const {
      BaseType.Inertia03(P_z_, nS_, (*P_).derived(), (*R_).derived(), (*P_a_).derived(), (*P_A_i_).derived(),
                         (*M_a_).derived(), (*Mtmp_).derived(), (*D_M_A_i_).derived());
      return 0;
    }

  };


  //! Leaf 01 Visitor
  //!------------------------------------------------------------------------------!//
  template<typename ScalarType, typename Vector6Type, typename Matrix6Type, typename RowVectorXType, typename D_Vector6Type>
  class Leaf01_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor
    Leaf01_visitor( ) : U_(nullptr), P_A_(nullptr), M_A_(nullptr), D_q_u_(nullptr), D_dq_u_(nullptr), D_q_PA_(nullptr), D_dq_PA_(nullptr) {}

    //! Members
    ScalarType* invD_;  ScalarType* u_;  Vector6Type* U_;  Vector6Type* P_A_;  Matrix6Type* M_A_;
    RowVectorXType* D_q_u_;  RowVectorXType* D_dq_u_;  D_Vector6Type* D_q_PA_;  D_Vector6Type* D_dq_PA_;

    //! Method -> member function
    template<typename Derived>
    int operator()(D_CRTPInterface<Derived> & BaseType) const {
      ScalarType & _u_ = (*u_);  ScalarType & _invD_ = (*invD_);
      BaseType.Leaf01(_invD_, _u_, (*U_).derived(), (*P_A_).derived(), (*M_A_).derived(),
                      (*D_q_u_).derived(), (*D_dq_u_).derived(), (*D_q_PA_).derived(), (*D_dq_PA_).derived());
      return 0;
    }

  };


  //! Leaf 02 Visitor
  //!------------------------------------------------------------------------------!//
  template<typename Vector3Type, typename Matrix3Type, typename Vector6Type, typename Matrix6Type, typename D_Matrix6Type>
  class Leaf02_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor
    Leaf02_visitor( ) : P_(nullptr), R_(nullptr), P_A_i_(nullptr), P_a_(nullptr), M_a_(nullptr),
      M_A_j_(nullptr), D_M_A_j_(nullptr) {}

    //! Members
    bool P_z_;  Vector3Type* P_;  Matrix3Type* R_;  Vector6Type* P_A_i_;  Vector6Type* P_a_;
    Matrix6Type* M_a_;  Matrix6Type* M_A_j_;  D_Matrix6Type* D_M_A_j_;

    //! Method -> member function
    template<typename Derived>
    int operator()(D_CRTPInterface<Derived> & BaseType) const {
      BaseType.Leaf02(P_z_, (*P_).derived(), (*R_).derived(), (*P_A_i_).derived(), (*P_a_).derived(),
                      (*M_a_).derived(), (*M_A_j_).derived(), (*D_M_A_j_).derived());
      return 0;
    }

  };


  //! Spatial Acceleration 01 Visitor
  //!------------------------------------------------------------------------------!//
  template<typename Vector6Type>
  class Accel01_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor
    Accel01_visitor( ) : AdAj_(nullptr), Aa_(nullptr) {}

    //! Members
    Vector6Type* AdAj_;  Vector6Type* Aa_;

    //! Method -> member function
    template<typename Derived>
    int operator()(D_CRTPInterface<Derived> & BaseType) const {
      BaseType.Accel01((*AdAj_).derived(), (*Aa_).derived());
      return 0;
    }

  };


  //! Visitor for computing the acceleration.
  //!------------------------------------------------------------------------------!//
  template<typename ScalarType, typename Vector6Type, typename RowVectorXType, typename D_Vector6Type>
  class Accel02_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor
    Accel02_visitor( ) : ddq_(nullptr), A_(nullptr), D_q_A_(nullptr), D_dq_A_(nullptr), D_q_ddq_(nullptr), D_dq_ddq_(nullptr) {}

    //! Members
    ScalarType* ddq_;  Vector6Type* A_;  D_Vector6Type* D_q_A_;   D_Vector6Type* D_dq_A_;  RowVectorXType* D_q_ddq_;  RowVectorXType* D_dq_ddq_;

    //! Method -> member function
    template<typename Derived>
    int operator()(D_CRTPInterface<Derived> & BaseType) const {
      BaseType.Accel02(ddq_, (*A_).derived(), (*D_q_A_).derived(), (*D_dq_A_).derived(), (*D_q_ddq_).derived(), (*D_dq_ddq_).derived());
      return 0;
    }
  };


  //! Visitor for computing the acceleration at root.
  //!------------------------------------------------------------------------------!//
  template<typename ScalarType, typename Vector3Type, typename Matrix3Type, typename Vector6Type,
           typename D_Vector6Type, typename RowVectorXType, typename MatrixXType>
  class AccelRoot_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor
    AccelRoot_visitor( ) : ddq(nullptr), S(nullptr), P_(nullptr), R_(nullptr), U_(nullptr), Acc_i_(nullptr), D_U_h_(nullptr),
      D_invD_(nullptr), D_q_u_(nullptr), D_dq_u_(nullptr), D_q_A_(nullptr), D_dq_A_(nullptr), D_ddq_(nullptr){}

    //! Members
    ScalarType u, iD;  ScalarType* ddq;
    Vector3Type* S;  Vector3Type* P_;  Matrix3Type* R_;  Vector6Type* U_;  Vector6Type* Acc_i_;
    //!-------------------------------------------------------
    D_Vector6Type* D_U_h_;  RowVectorXType* D_invD_;  RowVectorXType* D_q_u_;  RowVectorXType* D_dq_u_;
    D_Vector6Type* D_q_A_;  D_Vector6Type* D_dq_A_;   MatrixXType* D_ddq_;

    //! Method -> member function
    template<typename Derived>
    int operator()(D_CRTPInterface<Derived> & BaseType) const {
      BaseType.AccelRoot(u, iD, ddq, (*S).derived(), (*P_).derived(), (*R_).derived(), (*U_).derived(), (*Acc_i_).derived(),
                         (*D_U_h_).derived(), (*D_invD_).derived(), (*D_q_u_).derived(), (*D_dq_u_).derived(),
                         (*D_q_A_).derived(), (*D_dq_A_).derived(), (*D_ddq_).derived());
      return 0;
    }

  };

}

#endif // GEOMBD_DIFFERENTIATION_VISITORS_HXX
