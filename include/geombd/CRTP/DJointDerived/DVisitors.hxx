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























//!------------------------------------------------------------------------------!//
//!------------------------------------------------------------------------------!//
//!------------------------------------------------------------------------------!//







//  //! Visitor for Twist and bias elements at root
//  //!------------------------------------------------------------------------------!//
//  template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename Matrix6Type, typename D_Vector6Type>
//  class D_TCP_root_visitor : public boost::static_visitor<int> {
//  public:
//    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

//    //! Constructor
//    D_TCP_root_visitor( ) : S(nullptr), S_i(nullptr), p_(nullptr), M_(nullptr), D_dq_p_(nullptr) {}

//    //! Members
//    ScalarType vi;  Vector3Type* S;  Vector6Type* S_i;  Vector6Type* p_;  Matrix6Type* M_;  D_Vector6Type* D_dq_p_;

//    //! Method -> member function
//    template<typename Derived>
//    int operator()(D_CRTPInterface<Derived> & BaseType) const {
//      BaseType.D_TCP_root(vi, (*S).derived(), (*S_i).derived(), (*p_).derived(), (*M_).derived(), (*D_dq_p_).derived());
//      return 0;
//    }

//  };


//  //! Visitor for Twist and bias elements.
//  //!------------------------------------------------------------------------------!//
//  template<typename ScalarType, typename Matrix3Type, typename Matrix6Type,
//           typename Vector3Type, typename Vector6Type, typename D_Vector6Type>
//  class D_TwCbPb_visitor : public boost::static_visitor<int> {
//  public:
//    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

//    //! Constructor
//    D_TwCbPb_visitor( ) : S(nullptr), R_(nullptr), P_(nullptr), S_j(nullptr), M_(nullptr), S_i(nullptr), c_(nullptr), p_(nullptr),
//      D_q_V_(nullptr), D_dq_V_(nullptr), D_q_Vj_(nullptr), D_dq_Vj_(nullptr), D_q_c_(nullptr), D_dq_c_(nullptr),
//      D_q_p_(nullptr), D_dq_p_(nullptr){}

//    //! Members
//    bool zeroFlag;
//    ScalarType vi;    Matrix3Type* R_;   Vector3Type* P_;  Vector3Type* S;  Vector6Type* S_j;
//    Matrix6Type* M_;  Vector6Type* S_i;  Vector6Type* c_;  Vector6Type* p_;
//    //!-------------------------------------------------------
//    D_Vector6Type* D_q_V_;    D_Vector6Type* D_dq_V_;  D_Vector6Type* D_q_Vj_;  D_Vector6Type* D_q_p_;
//    D_Vector6Type* D_dq_Vj_;  D_Vector6Type* D_q_c_;   D_Vector6Type* D_dq_c_;  D_Vector6Type* D_dq_p_;

//    //! Method -> member function
//    template<typename Derived>
//    int operator()(D_CRTPInterface<Derived> & BaseType) const {
//      BaseType.D_TwCbPb(zeroFlag, vi, (*S).derived(), (*R_).derived(),  (*P_).derived(), (*S_j).derived(), (*M_).derived(), (*S_i).derived(),
//                        (*c_).derived(), (*p_).derived(), (*D_q_V_).derived(), (*D_dq_V_).derived(), (*D_q_Vj_).derived(),
//                        (*D_dq_Vj_).derived(), (*D_q_c_).derived(), (*D_dq_c_).derived(), (*D_q_p_).derived(), (*D_dq_p_).derived());
//      return 0;
//    }

//  };


//  //! Visitor for inertial expressions at leaf.
//  //!------------------------------------------------------------------------------!//
//  template<typename ScalarType, typename Vector3Type, typename Matrix3Type, typename Vector6Type,
//           typename Matrix6Type, typename RowVectorXType, typename D_Vector6Type, typename D_Matrix6Type>
//  class D_InertiaLeaf_visitor : public boost::static_visitor<int> {
//  public:
//    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

//    //! Constructor
//    D_InertiaLeaf_visitor( ) : u(nullptr), iD(nullptr), S(nullptr), U_(nullptr), c_(nullptr), P_a_(nullptr), M_a_(nullptr),
//      P_A_(nullptr), M_A_(nullptr), P_(nullptr), R_(nullptr), P_Aj_(nullptr), M_Aj_(nullptr), D_M_Aj_(nullptr),
//      D_q_u_(nullptr), D_dq_u_(nullptr), D_q_p_(nullptr), D_dq_p_(nullptr), D_q_Pa_(nullptr), D_dq_Pa_(nullptr),
//      D_q_PA_(nullptr), D_dq_PA_(nullptr), D_q_PAj_(nullptr), D_dq_PAj_(nullptr), D_q_c_(nullptr), D_dq_c_(nullptr) {}

//    //! Members
//    ScalarType* u;  ScalarType* iD;  ScalarType tau;  Vector3Type* S;  Vector6Type* U_;  Vector6Type* c_;
//    Vector6Type* P_a_;  Matrix6Type* M_a_;  Vector6Type* P_A_;  Matrix6Type* M_A_;
//    bool P_z;  Vector3Type* P_;  Matrix3Type* R_;  Vector6Type* P_Aj_;  Matrix6Type* M_Aj_;
//    //!-------------------------------------------------------
//    D_Matrix6Type* D_M_Aj_;  RowVectorXType* D_q_u_;   RowVectorXType* D_dq_u_;
//    D_Vector6Type* D_q_p_;   D_Vector6Type* D_dq_p_;   D_Vector6Type* D_q_Pa_;   D_Vector6Type* D_dq_Pa_;   D_Vector6Type* D_q_c_;
//    D_Vector6Type* D_q_PA_;  D_Vector6Type* D_dq_PA_;  D_Vector6Type* D_q_PAj_;  D_Vector6Type* D_dq_PAj_;  D_Vector6Type* D_dq_c_;

//    //! Method -> member function
//    template<typename Derived>
//    int operator()(D_CRTPInterface<Derived> & BaseType) const {
//      ScalarType & u_ = (*u);  ScalarType & iD_ = (*iD);
//      BaseType.D_InertiaLeaf(u_, iD_, tau, (*S).derived(), (*U_).derived(), (*c_).derived(), (*P_a_).derived(), (*M_a_).derived(),
//                             (*P_A_).derived(), (*M_A_).derived(), P_z, (*P_).derived(), (*R_).derived(), (*P_Aj_).derived(),
//                             (*M_Aj_).derived(), (*D_M_Aj_).derived(), (*D_q_u_).derived(), (*D_dq_u_).derived(),
//                             (*D_q_p_).derived(), (*D_dq_p_).derived(), (*D_q_Pa_).derived(), (*D_dq_Pa_).derived(),
//                             (*D_q_PA_).derived(), (*D_dq_PA_).derived(), (*D_q_PAj_).derived(), (*D_dq_PAj_).derived(),
//                             (*D_q_c_).derived(), (*D_dq_c_).derived());
//      return 0;
//    }

//  };


//  //! Visitor for inertial expressions.
//  //!------------------------------------------------------------------------------!//
//  template<typename ScalarType, typename Vector6iType, typename Vector3Type, typename Matrix3Type, typename Vector6Type,
//           typename Matrix6Type, typename VectorXType, typename RowVectorXType, typename D_Vector6Type, typename D_Matrix6Type>
//  class D_Inertial_visitor : public boost::static_visitor<int> {
//  public:
//    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

//    //! Constructor
//    D_Inertial_visitor( ) : u(nullptr), iD(nullptr), S(nullptr), U_(nullptr), c_(nullptr), P_a_(nullptr), M_a_(nullptr),
//      P_A_(nullptr), M_A_(nullptr), P_(nullptr), R_(nullptr), P_Aj_(nullptr), M_Aj_(nullptr), D_U_h_(nullptr),
//      D_U_v_(nullptr), D_invD_(nullptr), D_M_A_(nullptr), D_M_Aj_(nullptr){}

//    //! Members
//    Vector6iType* indexVec;  bool rootFlag;
//    ScalarType* u;  ScalarType* iD;  ScalarType tau;  Vector3Type* S;  Vector6Type* U_;  Vector6Type* c_;
//    Vector6Type* P_a_;  Matrix6Type* M_a_;  Vector6Type* P_A_;  Matrix6Type* M_A_;
//    bool P_z;  Vector3Type* P_;  Matrix3Type* R_;  Vector6Type* P_Aj_;  Matrix6Type* M_Aj_;
//    //!-------------------------------------------------------
//    D_Vector6Type* D_U_h_;  VectorXType* D_U_v_;  RowVectorXType* D_invD_;
//    D_Matrix6Type* D_M_A_;  D_Matrix6Type* D_M_Aj_;
//    //!-------------------------------------------------------
//    RowVectorXType* D_q_u_;   RowVectorXType* D_dq_u_;
//    D_Vector6Type* D_q_Pa_;   D_Vector6Type* D_dq_Pa_;  D_Vector6Type* D_q_c_;   D_Vector6Type* D_q_PA_;
//    D_Vector6Type* D_dq_PA_;  D_Vector6Type* D_q_PAj_;  D_Vector6Type* D_dq_PAj_;  D_Vector6Type* D_dq_c_;


//    //! Method -> member function
//    template<typename Derived>
//    int operator()(D_CRTPInterface<Derived> & BaseType) const {
//      ScalarType & u_ = (*u);  ScalarType & iD_ = (*iD);
//      BaseType.D_Inertial(rootFlag, (*indexVec).derived(), u_, iD_, tau, (*S).derived(), (*U_).derived(), (*c_).derived(), (*P_a_).derived(),
//                          (*M_a_).derived(), (*P_A_).derived(), (*M_A_).derived(), P_z, (*P_).derived(), (*R_).derived(),
//                          (*P_Aj_).derived(), (*M_Aj_).derived(), (*D_U_h_).derived(), (*D_U_v_).derived(),
//                          (*D_invD_).derived(), (*D_M_A_).derived(), (*D_M_Aj_).derived(), (*D_q_u_).derived(), (*D_dq_u_).derived(),
//                          (*D_q_Pa_).derived(), (*D_dq_Pa_).derived(), (*D_q_PA_).derived(), (*D_dq_PA_).derived(), (*D_q_PAj_).derived(),
//                          (*D_dq_PAj_).derived(), (*D_q_c_).derived(), (*D_dq_c_).derived());
//      return 0;
//    }

//  };


//  //! Visitor for computing the acceleration at root.
//  //!------------------------------------------------------------------------------!//
//  template<typename ScalarType, typename Vector3Type, typename Matrix3Type, typename Vector6Type,
//           typename D_Vector6Type, typename RowVectorXType, typename MatrixXType>
//  class D_Accel_root_visitor : public boost::static_visitor<int> {
//  public:
//    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

//    //! Constructor
//    D_Accel_root_visitor( ) : ddq(nullptr), S(nullptr), P_(nullptr), R_(nullptr), U_(nullptr), Acc_i_(nullptr), D_U_h_(nullptr),
//      D_invD_(nullptr), D_q_u_(nullptr), D_dq_u_(nullptr), D_q_A_(nullptr), D_dq_A_(nullptr), D_ddq_(nullptr){}

//    //! Members
//    ScalarType u, iD;  ScalarType* ddq;
//    Vector3Type* S;  Vector3Type* P_;  Matrix3Type* R_;  Vector6Type* U_;  Vector6Type* Acc_i_;
//    //!-------------------------------------------------------
//    D_Vector6Type* D_U_h_;  RowVectorXType* D_invD_;  RowVectorXType* D_q_u_;  RowVectorXType* D_dq_u_;
//    D_Vector6Type* D_q_A_;  D_Vector6Type* D_dq_A_;  MatrixXType* D_ddq_;

//    //! Method -> member function
//    template<typename Derived>
//    int operator()(D_CRTPInterface<Derived> & BaseType) const {
//      BaseType.D_AccelRoot(u, iD, ddq, (*S).derived(), (*P_).derived(), (*R_).derived(), (*U_).derived(), (*Acc_i_).derived(),
//                           (*D_U_h_).derived(), (*D_invD_).derived(), (*D_q_u_).derived(), (*D_dq_u_).derived(),
//                           (*D_q_A_).derived(), (*D_dq_A_).derived(), (*D_ddq_).derived());
//      return 0;
//    }

//  };


//  //! Visitor for computing the acceleration.
//  //!------------------------------------------------------------------------------!//
//  template<typename ScalarType, typename IndexType, typename Vector3Type, typename Matrix3Type,
//           typename Vector6Type, typename D_Vector6Type, typename RowVectorXType, typename MatrixXType, typename VectorXiType>
//  class D_Accel_visitor : public boost::static_visitor<int> {
//  public:
//    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

//    //! Constructor
//    D_Accel_visitor( ) : ddq(nullptr), S(nullptr), P_(nullptr), R_(nullptr), c_(nullptr), U_(nullptr), Acc_i_(nullptr), Acc_j_(nullptr),
//      D_U_h_(nullptr), D_invD_(nullptr), D_ddq_(nullptr), D_q_u_(nullptr), D_dq_u_(nullptr), D_q_c_(nullptr), D_dq_c_(nullptr),
//      D_q_A_(nullptr), D_dq_A_(nullptr), D_q_Aj_(nullptr), D_dq_Aj_(nullptr), Pre_(nullptr), Suc_(nullptr), PreSuc_(nullptr) {}

//    //! Members
//    bool zeroFlag;  ScalarType u, iD;  ScalarType* ddq;
//    Vector3Type* S;  Vector3Type* P_;  Matrix3Type* R_;  Vector6Type* c_;
//    Vector6Type* U_;  Vector6Type* Acc_i_;  Vector6Type* Acc_j_;
//    //!-------------------------------------------------------
//    bool isLeaf;  IndexType ID;
//    VectorXiType* Pre_;      VectorXiType* Suc_;       VectorXiType* PreSuc_;
//    D_Vector6Type* D_U_h_;   RowVectorXType* D_invD_;  MatrixXType* D_ddq_;
//    RowVectorXType* D_q_u_;  RowVectorXType* D_dq_u_;  D_Vector6Type* D_q_c_;  D_Vector6Type* D_dq_c_;
//    D_Vector6Type* D_q_A_;   D_Vector6Type* D_dq_A_;   D_Vector6Type* D_q_Aj_;  D_Vector6Type* D_dq_Aj_;

//    //! Method -> member function
//    template<typename Derived>
//    int operator()(D_CRTPInterface<Derived> & BaseType) const {
//      BaseType.D_Accel(ID, zeroFlag, u, iD, ddq, (*S).derived(), (*P_).derived(), (*R_).derived(), (*c_).derived(), (*U_).derived(),
//                       (*Acc_i_).derived(), (*Acc_j_).derived(), isLeaf, (*Pre_).derived(), (*Suc_).derived(), (*PreSuc_).derived(), (*D_U_h_).derived(), (*D_invD_).derived(),
//                       (*D_ddq_).derived(), (*D_q_u_).derived(), (*D_dq_u_).derived(), (*D_q_c_).derived(), (*D_dq_c_).derived(),
//                       (*D_q_A_).derived(), (*D_dq_A_).derived(), (*D_q_Aj_).derived(), (*D_dq_Aj_).derived());
//      return 0;
//    }

//  };

}

#endif // GEOMBD_DIFFERENTIATION_VISITORS_HXX
