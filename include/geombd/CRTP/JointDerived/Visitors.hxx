/**
 *    \file include/geombd/CRTP/JointDerived/Visitors.hxx
 *    \author Alvaro Paz, Gustavo Arechavaleta
 *    \version 1.0
 *    \date 2021
 *
 *    Classes for visitors
 *    Copyright (c) 2021 Cinvestav
 *    This library is distributed under the MIT License.
 */

#ifndef GEOMBD_VISITORS_HXX
#define GEOMBD_VISITORS_HXX

#define EIGEN_NO_DEBUG
#include "Eigen/Core"

//#include <boost/foreach.hpp>
//#include <boost/shared_ptr.hpp>
#include <boost/variant.hpp>


namespace geoCRTP{

  //! Forward Kinematics Visitor
  //!------------------------------------------------------------------------------!//
  template<typename ScalarType, typename Vector3Type, typename Matrix3Type>
  class FwdKin_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    FwdKin_visitor( ) : R(nullptr), S(nullptr) {}

    ScalarType qi;
    Vector3Type* S;
    Matrix3Type* R;

    template<typename Derived>
    int operator()(CRTPInterface<Derived> & BaseType) const {
        BaseType.FwdKin(qi, (*S).derived(), (*R).derived());
      return 0;
    }

  };


  //! Visitor for Twist and bias elements at root
  //!------------------------------------------------------------------------------!//
  template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename Matrix6Type>
  class TCP_root_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor
    TCP_root_visitor( ) : S(nullptr), S_i(nullptr), p_(nullptr), M_(nullptr) {}

    //! Members
    ScalarType vi;  Vector3Type* S;  Vector6Type* S_i;  Vector6Type* p_;  Matrix6Type* M_;

    //! Method -> member function
    template<typename Derived>
    int operator()(CRTPInterface<Derived> & BaseType) const {
        BaseType.TCP_root(vi, (*S).derived(), (*S_i).derived(), (*p_).derived(), (*M_).derived());
      return 0;
    }

  };


  //! Visitor for Twist and bias elements.
  //!------------------------------------------------------------------------------!//
  template<typename ScalarType, typename Matrix3Type, typename Matrix6Type,
           typename Vector3Type, typename Vector6Type>
  class TwCbPb_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor
    TwCbPb_visitor( ) : S(nullptr), R_(nullptr), P_(nullptr), S_j(nullptr),
                        M_(nullptr), S_i(nullptr), c_(nullptr), p_(nullptr) {}

    //! Members
    bool zeroFlag;
    ScalarType vi;    Vector3Type* S;    Matrix3Type* R_;   Vector3Type* P_;  Vector6Type* S_j;
    Matrix6Type* M_;  Vector6Type* S_i;  Vector6Type* c_;   Vector6Type* p_;

    //! Method -> member function
    template<typename Derived>
    int operator()(CRTPInterface<Derived> & BaseType) const {
      BaseType.TwCbPb(zeroFlag, vi, (*S).derived(), (*R_).derived(),  (*P_).derived(), (*S_j).derived(), (*M_).derived(),
                      (*S_i).derived(), (*c_).derived(), (*p_).derived());
      return 0;
    }

  };


  //! Visitor for U, u & invD expressions
  //!------------------------------------------------------------------------------!//
  template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename Matrix6Type>
  class UuiD_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor
    UuiD_visitor( ) : u(nullptr), iD(nullptr), S(nullptr), U_(nullptr), P_A_(nullptr), M_A_(nullptr) {}

    //! Members
    ScalarType* u;    ScalarType* iD;     ScalarType tau;     Vector3Type* S;
    Vector6Type* U_;  Vector6Type* P_A_;  Matrix6Type* M_A_;

    //! Method -> member function
    template<typename Derived>
    int operator()(CRTPInterface<Derived> & BaseType) const {
      ScalarType & u_ = (*u);  ScalarType & iD_ = (*iD);
      BaseType.UuiD(u_, iD_, tau, (*S).derived(), (*U_).derived(), (*P_A_).derived(), (*M_A_).derived());
      return 0;
    }

  };


  //! Visitor for preparing inertial expressions
  //!------------------------------------------------------------------------------!//
  template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename Matrix6Type>
  class PreIner_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor
    PreIner_visitor( ) : S(nullptr), U_(nullptr), c_(nullptr), P_a_(nullptr),
                         M_a_(nullptr), P_A_(nullptr), M_A_(nullptr) {}

    //! Members
    ScalarType u, iD;   Vector3Type* S;
    Vector6Type* U_;    Vector6Type* c_;    Vector6Type* P_a_;
    Matrix6Type* M_a_;  Vector6Type* P_A_;  Matrix6Type* M_A_;

    //! Method -> member function
    template<typename Derived>
    int operator()(CRTPInterface<Derived> & BaseType) const {
      BaseType.PreIner(u, iD, (*S).derived(), (*U_).derived(), (*c_).derived(), (*P_a_).derived(),
                       (*M_a_).derived(), (*P_A_).derived(), (*M_A_).derived());
      return 0;
    }

  };


  //! Visitor for inertial back projection
  //!------------------------------------------------------------------------------!//
  template<typename Vector3Type, typename Matrix3Type, typename Vector6Type, typename Matrix6Type>
  class InerProj_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor
    InerProj_visitor( ) : S(nullptr), P_(nullptr), R_(nullptr), P_a_(nullptr), P_A_(nullptr), M_a_(nullptr), M_A_(nullptr) {}

    //! Members
    bool P_z;           Vector3Type* S;
    Vector3Type* P_;    Matrix3Type* R_;
    Vector6Type* P_a_;  Vector6Type* P_A_;
    Matrix6Type* M_a_;  Matrix6Type* M_A_;

    //! Method -> member function
    template<typename Derived>
    int operator()(CRTPInterface<Derived> & BaseType) const {
      BaseType.InerProj(P_z, (*S).derived(), (*P_).derived(), (*R_).derived(), (*P_a_).derived(), (*P_A_).derived(),
                        (*M_a_).derived(), (*M_A_).derived());
      return 0;
    }

  };


  //! Visitor for computing the acceleration
  //!------------------------------------------------------------------------------!//
  template<typename ScalarType, typename Vector3Type, typename Matrix3Type, typename Vector6Type>
  class Accel_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor
    Accel_visitor( ) : ddq(nullptr), S_(nullptr), P_(nullptr), R_(nullptr), c_(nullptr),
                       U_(nullptr), Acc_i_(nullptr), Acc_j_(nullptr) {}

    //! Members
    bool zeroFlag;
    ScalarType u, iD;
    ScalarType* ddq;

    Vector3Type* S_;  Vector3Type* P_;  Matrix3Type* R_;  Vector6Type* c_;
    Vector6Type* U_;  Vector6Type* Acc_i_;  Vector6Type* Acc_j_;

    //! Method -> member function
    template<typename Derived>
    int operator()(CRTPInterface<Derived> & BaseType) const {
      BaseType.Accel(zeroFlag, u, iD, ddq, (*S_).derived(), (*P_).derived(), (*R_).derived(), (*c_).derived(),
                     (*U_).derived(), (*Acc_i_).derived(), (*Acc_j_).derived());
      return 0;
    }

  };


  //! Visitor for computing the acceleration at root
  //!------------------------------------------------------------------------------!//
  template<typename ScalarType, typename Vector3Type, typename Matrix3Type, typename Vector6Type>
  class Accel_root_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor
    Accel_root_visitor( ) : ddq(nullptr), S_(nullptr), P_(nullptr), R_(nullptr), U_(nullptr), Acc_i_(nullptr) {}

    //! Members
    ScalarType u, iD;
    ScalarType* ddq;

    Vector3Type* S_;  Vector3Type* P_;  Matrix3Type* R_;  Vector6Type* U_;  Vector6Type* Acc_i_;

    //! Method -> member function
    template<typename Derived>
    int operator()(CRTPInterface<Derived> & BaseType) const {
      BaseType.AccelRoot(u, iD, ddq, (*S_).derived(), (*P_).derived(), (*R_).derived(), (*U_).derived(), (*Acc_i_).derived());
      return 0;
    }

  };


  //!------------------------------------------------------------------------------!//
  //!----------------------Inverse Inertia Matrix Visitors-------------------------!//
  //!------------------------------------------------------------------------------!//


  //! Visitor for U & invD expressions
  //!------------------------------------------------------------------------------!//
  template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename Matrix6Type>
  class UinvD_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor
    UinvD_visitor( ) : iD(nullptr), S(nullptr), U_(nullptr), M_A_(nullptr) {}

    //! Members
    ScalarType* iD;  Vector3Type* S;  Vector6Type* U_;  Matrix6Type* M_A_;

    //! Method -> member function
    template<typename Derived>
    int operator()(CRTPInterface<Derived> & BaseType) const {
      ScalarType & iD_ = (*iD);
      BaseType.UinvD(iD_, (*S).derived(), (*U_).derived(), (*M_A_).derived());
      return 0;
    }

  };


  //! Visitor for inverse inertia expressions
  //!------------------------------------------------------------------------------!//
  template<typename IntType, typename ScalarType, typename Vector3Type, typename D_Vector6Type, typename RowVectorXrType>
  class invInertia_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor
    invInertia_visitor( ) : iD(nullptr), S(nullptr), Fcrb_(nullptr), iHrow_(nullptr) {}

    //! Members
    IntType n_ID;   ScalarType* iD;   Vector3Type* S;  D_Vector6Type* Fcrb_;  RowVectorXrType* iHrow_;

    //! Method -> member function
    template<typename Derived>
    int operator()(CRTPInterface<Derived> & BaseType) const {
      ScalarType & iD_ = (*iD);
      BaseType.invInertia(n_ID, iD_, (*S).derived(), (*Fcrb_).derived(), (*iHrow_).derived());
      return 0;
    }

  };


  //! Visitor for H selector
  //!------------------------------------------------------------------------------!//
  template<typename IntType, typename Vector3Type, typename D_Vector6Type, typename MatrixXrType>
  class HSelector_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor
    HSelector_visitor( ) : S_(nullptr), Pc_(nullptr), iH_total_(nullptr) {}

    //! Members
    IntType n, ID;   Vector3Type* S_;  D_Vector6Type* Pc_;  MatrixXrType* iH_total_;

    //! Method -> member function
    template<typename Derived>
    int operator()(CRTPInterface<Derived> & BaseType) const {
      BaseType.invInertia(n, ID, (*S_).derived(), (*Pc_).derived(), (*iH_total_).derived());
      return 0;
    }

  };

}

#endif // GEOMBD_VISITORS_HXX
