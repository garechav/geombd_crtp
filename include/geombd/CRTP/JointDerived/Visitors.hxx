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
  template<typename ScalarType, typename Matrix3Type>
  class FwdKin_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    FwdKin_visitor( ) : R(nullptr) {}

    ScalarType qi;
    Matrix3Type* R;

    template<typename Derived>
    int operator()(CRTPInterface<Derived> & BaseType) const {
        BaseType.FwdKin(qi, (*R).derived());
      return 0;
    }

  };


  //! Visitor for Twist and bias elements at root
  //!------------------------------------------------------------------------------!//
  template<typename ScalarType, typename Vector6Type, typename Matrix6Type>
  class TCP_root_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor
    TCP_root_visitor( ) : S_i(nullptr), p_(nullptr), M_(nullptr) {}

    //! Members
    ScalarType vi;  Vector6Type* S_i;  Vector6Type* p_;  Matrix6Type* M_;

    //! Method -> member function
    template<typename Derived>
    int operator()(CRTPInterface<Derived> & BaseType) const {
        BaseType.TCP_root(vi, (*S_i).derived(), (*p_).derived(), (*M_).derived());
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
    TwCbPb_visitor( ) : R_(nullptr), P_(nullptr), S_j(nullptr), M_(nullptr),
                        S_i(nullptr), c_(nullptr), p_(nullptr) {}

    //! Members
    bool zeroFlag;
    ScalarType vi;    Matrix3Type* R_;   Vector3Type* P_;  Vector6Type* S_j;
    Matrix6Type* M_;  Vector6Type* S_i;  Vector6Type* c_;  Vector6Type* p_;

    //! Method -> member function
    template<typename Derived>
    int operator()(CRTPInterface<Derived> & BaseType) const {
      BaseType.TwCbPb(zeroFlag, vi, (*R_).derived(),  (*P_).derived(), (*S_j).derived(), (*M_).derived(),
                      (*S_i).derived(), (*c_).derived(), (*p_).derived());
      return 0;
    }

  };


  //! Visitor for U, u & invD expressions
  //!------------------------------------------------------------------------------!//
  template<typename ScalarType, typename Vector6Type, typename Matrix6Type>
  class UuiD_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor
    UuiD_visitor( ) : u(nullptr), iD(nullptr), U_(nullptr), P_A_(nullptr), M_A_(nullptr) {}

    //! Members
    ScalarType* u; ScalarType* iD; ScalarType tau;
    Vector6Type* U_;  Vector6Type* P_A_;  Matrix6Type* M_A_;

    //! Method -> member function
    template<typename Derived>
    int operator()(CRTPInterface<Derived> & BaseType) const {
      ScalarType & u_ = (*u);  ScalarType & iD_ = (*iD);
      BaseType.UuiD(u_, iD_, tau, (*U_).derived(), (*P_A_).derived(), (*M_A_).derived());
      return 0;
    }

  };


  //! Visitor for preparing inertial expressions
  //!------------------------------------------------------------------------------!//
  template<typename ScalarType, typename Vector6Type, typename Matrix6Type>
  class PreIner_visitor : public boost::static_visitor<int> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor
    PreIner_visitor( ) : U_(nullptr), c_(nullptr), P_a_(nullptr),
                         M_a_(nullptr), P_A_(nullptr), M_A_(nullptr) {}

    //! Members
    ScalarType u, iD;
    Vector6Type* U_;    Vector6Type* c_;    Vector6Type* P_a_;
    Matrix6Type* M_a_;  Vector6Type* P_A_;  Matrix6Type* M_A_;

    //! Method -> member function
    template<typename Derived>
    int operator()(CRTPInterface<Derived> & BaseType) const {
      BaseType.PreIner(u, iD, (*U_).derived(), (*c_).derived(), (*P_a_).derived(),
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
    InerProj_visitor( ) : P_(nullptr), R_(nullptr), P_a_(nullptr), P_A_(nullptr), M_a_(nullptr), M_A_(nullptr) {}

    //! Members
    bool P_z;
    Vector3Type* P_;    Matrix3Type* R_;
    Vector6Type* P_a_;  Vector6Type* P_A_;
    Matrix6Type* M_a_;  Matrix6Type* M_A_;

    //! Method -> member function
    template<typename Derived>
    int operator()(CRTPInterface<Derived> & BaseType) const {
      BaseType.InerProj(P_z, (*P_).derived(), (*R_).derived(), (*P_a_).derived(), (*P_A_).derived(),
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
    Accel_visitor( ) : ddq(nullptr), P_(nullptr), R_(nullptr), c_(nullptr),
      U_(nullptr), Acc_i_(nullptr), Acc_j_(nullptr) {}

    //! Members
    bool zeroFlag;
    ScalarType u, iD;
    ScalarType* ddq;

    Vector3Type* P_;  Matrix3Type* R_;  Vector6Type* c_;
    Vector6Type* U_;  Vector6Type* Acc_i_;  Vector6Type* Acc_j_;

    //! Method -> member function
    template<typename Derived>
    int operator()(CRTPInterface<Derived> & BaseType) const {
      BaseType.Accel(zeroFlag, u, iD, ddq, (*P_).derived(), (*R_).derived(), (*c_).derived(),
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
    Accel_root_visitor( ) : ddq(nullptr), P_(nullptr), R_(nullptr), U_(nullptr), Acc_i_(nullptr) {}

    //! Members
    ScalarType u, iD;
    ScalarType* ddq;

    Vector3Type* P_;  Matrix3Type* R_;  Vector6Type* U_;  Vector6Type* Acc_i_;

    //! Method -> member function
    template<typename Derived>
    int operator()(CRTPInterface<Derived> & BaseType) const {
      BaseType.AccelRoot(u, iD, ddq, (*P_).derived(), (*R_).derived(), (*U_).derived(), (*Acc_i_).derived());
      return 0;
    }

  };


}

#endif // GEOMBD_VISITORS_HXX
