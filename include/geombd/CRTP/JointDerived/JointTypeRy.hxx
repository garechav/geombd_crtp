/**
 *    \file include/geombd/CRTP/JointDerived/JointTypeRy.hxx
 *    \author Alvaro Paz, Gustavo Arechavaleta
 *    \version 1.0
 *    \date 2021
 *
 *    Derived class for Ry joint type
 *    Copyright (c) 2021 Cinvestav
 *    This library is distributed under the MIT License.
 */

#ifndef GEOMBD_JOINT_TYPE_REVOLUTE_Y_HXX
#define GEOMBD_JOINT_TYPE_REVOLUTE_Y_HXX

#define EIGEN_NO_DEBUG
#include "Eigen/Core"
#include "Eigen/Geometry"

namespace geoCRTP{


  //! Derived -> Joint Type Revolute On Y Axis
  //!------------------------------------------------------------------------------!//
  struct JointTypeRy : public CRTPInterface<JointTypeRy>
  {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! FK Forward Declaration.
    template<typename ScalarType, typename Vector3Type, typename Matrix3Type>
    inline static void
    runFK(const ScalarType & qi,
          typename Eigen::MatrixBase<Vector3Type> & S,
          typename Eigen::MatrixBase<Matrix3Type> & R);


    //! TwCbPb at Root Forward Declaration.
    template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename Matrix6Type>
    inline static void
    runTCP_root(const ScalarType & vi,
                const Eigen::MatrixBase<Vector3Type> & S,
                Eigen::MatrixBase<Vector6Type> & S_i,
                Eigen::MatrixBase<Vector6Type> & p_,
                const Eigen::MatrixBase<Matrix6Type> & M_);


    //! TwCbPb Forward Declaration.
    template<typename ScalarType, typename Matrix3Type, typename Matrix6Type,
             typename Vector3Type, typename Vector6Type>
    inline static void
    runTwCbPb(bool zeroFlag,
              const ScalarType & vi,
              const Eigen::MatrixBase<Vector3Type> & S,
              const Eigen::MatrixBase<Matrix3Type> & R_,
              const Eigen::MatrixBase<Vector3Type> & P_,
              const Eigen::MatrixBase<Vector6Type> & S_l,
              const Eigen::MatrixBase<Matrix6Type> & M_,
              Eigen::MatrixBase<Vector6Type> & S_i,
              Eigen::MatrixBase<Vector6Type> & c_,
              Eigen::MatrixBase<Vector6Type> & p_);


    //! UuinvD Forward Declaration.
    template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename Matrix6Type>
    inline static void
    runUuiD(ScalarType & u,
            ScalarType & iD,
            const ScalarType tau,
            const Eigen::MatrixBase<Vector3Type> & S,
            Eigen::MatrixBase<Vector6Type> & U_r,
            const Eigen::MatrixBase<Vector6Type> & P_A_r,
            const Eigen::MatrixBase<Matrix6Type> & M_A_r);


    //! Preparing Inertias Forward Declaration.
    template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename Matrix6Type>
    inline static void
    runPreIner(ScalarType u,
               ScalarType iD,
               const Eigen::MatrixBase<Vector3Type> & S,
               const Eigen::MatrixBase<Vector6Type> & U_r,
               const Eigen::MatrixBase<Vector6Type> & c_r,
               Eigen::MatrixBase<Vector6Type> & P_a_r,
               Eigen::MatrixBase<Matrix6Type> & M_a_r,
               const Eigen::MatrixBase<Vector6Type> & P_A_r,
               const Eigen::MatrixBase<Matrix6Type> & M_A_r);


    //! Forward Declaration for inertial back-projection.
    template<typename Vector3Type, typename Matrix3Type, typename Vector6Type, typename Matrix6Type>
    inline static void
    runInerProj(bool P_z,
                const Eigen::MatrixBase<Vector3Type> & S,
                const Eigen::MatrixBase<Vector3Type> & P_r,
                const Eigen::MatrixBase<Matrix3Type> & R_r,
                const Eigen::MatrixBase<Vector6Type> & P_a_r,
                Eigen::MatrixBase<Vector6Type> & P_A_r,
                const Eigen::MatrixBase<Matrix6Type> & M_a_r,
                Eigen::MatrixBase<Matrix6Type> & M_A_r);


    //! Forward Declaration for Acceleration.
    template<typename ScalarType, typename Vector3Type, typename Matrix3Type, typename Vector6Type>
    inline static void
    runAccel(bool zeroFlag,
             ScalarType u,
             ScalarType iD,
             ScalarType* ddq,
             const Eigen::MatrixBase<Vector3Type> & S,
             const Eigen::MatrixBase<Vector3Type> & P_r,
             const Eigen::MatrixBase<Matrix3Type> & R_r,
             const Eigen::MatrixBase<Vector6Type> & c_r,
             const Eigen::MatrixBase<Vector6Type> & U_r,
             Eigen::MatrixBase<Vector6Type> & Acc_i_r,
             const Eigen::MatrixBase<Vector6Type> & Acc_j_r);


    //! Forward Declaration for Acceleration at Root.
    template<typename ScalarType, typename Vector3Type, typename Matrix3Type, typename Vector6Type>
    inline static void
    runAccelRoot(ScalarType u,
                 ScalarType iD,
                 ScalarType* ddq,
                 const Eigen::MatrixBase<Vector3Type> & S,
                 const Eigen::MatrixBase<Vector3Type> & P_r,
                 const Eigen::MatrixBase<Matrix3Type> & R_r,
                 const Eigen::MatrixBase<Vector6Type> & U_r,
                 Eigen::MatrixBase<Vector6Type> & Acc_i_r);

  };


  //! Member-function declaration
  //!------------------------------------------------------------------------------!//
  //!------------------------------------------------------------------------------!//
  //!------------------------------------------------------------------------------!//


  //! Forward Kinematics Declaration
  //! EIGEN_ASM_COMMENT("MyBegin");
  //!------------------------------------------------------------------------------!//
  template<typename ScalarType, typename Vector3Type, typename Matrix3Type>
  inline void
  JointTypeRy::runFK(const ScalarType & qi,
                     typename Eigen::MatrixBase<Vector3Type> & S,
                     typename Eigen::MatrixBase<Matrix3Type> & R) {
    static ScalarType sqi, cqi;
    ::geo::SINCOS<ScalarType>(qi, &sqi, &cqi);

    R.coeffRef(0,0) =  cqi;  R.coeffRef(0,2) = sqi;
    R.coeffRef(2,0) = -sqi;  R.coeffRef(2,2) = cqi;
  }


  //! TwCbPb at Root Declaration.
  template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename Matrix6Type>
  inline void
  JointTypeRy::runTCP_root(const ScalarType & vi,
                           const Eigen::MatrixBase<Vector3Type> & S,
                           Eigen::MatrixBase<Vector6Type> & S_i,
                           Eigen::MatrixBase<Vector6Type> & p_,
                           const Eigen::MatrixBase<Matrix6Type> & M_) {

    S_i.coeffRef(3) = vi;

    Matrix6Type & M = const_cast<Matrix6Type &>(M_.derived());

    ScalarType vi2 = vi * vi;

    p_.coeffRef(0) =  M.coeff(2,4)*vi2;
    p_.coeffRef(2) = -M.coeff(0,4)*vi2;
    p_.coeffRef(3) =  M.coeff(5,4)*vi2;
    p_.coeffRef(5) = -M.coeff(3,4)*vi2;
  }


  //! Apply the Ad operator on se(3) element to get the Twist and bias
  //! TODO: Twist is only used here, thus consider to store it only when multiple children
  //!------------------------------------------------------------------------------!//
  template<typename ScalarType, typename Matrix3Type, typename Matrix6Type,
           typename Vector3Type, typename Vector6Type>
  inline void
  JointTypeRy::runTwCbPb(bool zeroFlag,
                         const ScalarType & vi,
                         const Eigen::MatrixBase<Vector3Type> & S,
                         const Eigen::MatrixBase<Matrix3Type> & R_,
                         const Eigen::MatrixBase<Vector3Type> & P_,
                         const Eigen::MatrixBase<Vector6Type> & S_l,
                         const Eigen::MatrixBase<Matrix6Type> & M_,
                         Eigen::MatrixBase<Vector6Type> & S_i,
                         Eigen::MatrixBase<Vector6Type> & c_,
                         Eigen::MatrixBase<Vector6Type> & p_) {

    EIGEN_STATIC_ASSERT(Vector6Type::ColsAtCompileTime == 1,
                        YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_A_MATRIX);

    //! Twist
    //!------------------------------------------------------------------------------!//
    typedef const Eigen::Block<Vector6Type,3,1> constSegment3;
    typedef Eigen::Block<Vector6Type,3,1> Segment3;

    Vector6Type & S_ = const_cast<Vector6Type &>(S_l.derived());
    constSegment3 & Sup   = S_.template segment<3>(0);
    constSegment3 & Sdown = S_.template segment<3>(3);

    Segment3 Up_   = S_i.template segment<3>(0);
    Segment3 Down_ = S_i.template segment<3>(3);

    if(zeroFlag) {
        Down_.noalias() = Sup - P_.cross(Sdown);  //tmp
        Up_.noalias()   = R_.transpose()*Down_;
      } else {
        Up_.noalias()   = R_.transpose()*Sup;
      }

    Down_.noalias() = R_.transpose()*Sdown;

    Down_.coeffRef(1) += vi;

    //! C bias
    //!------------------------------------------------------------------------------!//
    Segment3 Cup   = c_.template segment<3>(0);
    Segment3 Cdown = c_.template segment<3>(3);    

    Cup.coeffRef(0) = -Up_.coeff(2)*vi;  // -ad(Sy) effect
    Cup.coeffRef(2) =  Up_.coeff(0)*vi;

    Cdown.coeffRef(0) = -Down_.coeff(2)*vi;
    Cdown.coeffRef(2) =  Down_.coeff(0)*vi;

    //! P bias
    //!------------------------------------------------------------------------------!//
    typedef const Eigen::Block<Matrix6Type,3,3> constBlock3;
    typedef Eigen::Block<Matrix6Type,3,3> Block3;

    Matrix6Type & M = const_cast<Matrix6Type &>(M_.derived());
    constBlock3 & UL = M.template block<3,3>(0,0);
    constBlock3 & UR = M.template block<3,3>(0,3);
    constBlock3 & LL = M.template block<3,3>(3,0);
    constBlock3 & LR = M.template block<3,3>(3,3);

    typename Vector3Type::PlainObject muUp, muDo;
//    typedef typename GEOMBD_EIGEN_PLAIN_TYPE(Vector3Type) PlainVector3Type;
//    PlainVector3Type muUp, muDo;

    Segment3 p_Up = p_.template segment<3>(0);
    Segment3 p_Do = p_.template segment<3>(3);

    //! M is symmetric and its shape is known
    muUp.noalias() = UL*Up_ + UR*Down_;
    muDo.noalias() = LL*Up_ + LR*Down_;

    //! TODO: Cross is very expensive, consider enable staticness with traits expressions
    p_Up.noalias() = Skew(Down_)*muUp;
    p_Do.noalias() = Skew(Up_)*muUp;
    p_Do.noalias() += Skew(Down_)*muDo;
  }

//  EIGEN_ALWAYS_INLINE EIGEN_DEVICE_FUNC

  //! Implementation for U, u & invD.
  template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename Matrix6Type>
  inline void
  JointTypeRy::runUuiD(ScalarType & u,
                       ScalarType & iD,
                       const ScalarType tau,
                       const Eigen::MatrixBase<Vector3Type> & S,
                       Eigen::MatrixBase<Vector6Type> & U_r,
                       const Eigen::MatrixBase<Vector6Type> & P_A_r,
                       const Eigen::MatrixBase<Matrix6Type> & M_A_r) {

    EIGEN_STATIC_ASSERT(Vector6Type::ColsAtCompileTime == 1,
                        YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_A_MATRIX);

    U_r = M_A_r.template col(4);
    iD = 1 / U_r.coeff(4);
    u = tau - P_A_r.coeff(4);
  }


  //! Implementation for preparing inertial expressions.
  template<typename ScalarType, typename Vector3Type, typename Vector6Type, typename Matrix6Type>
  inline void
  JointTypeRy::runPreIner(ScalarType u,
                          ScalarType iD,
                          const Eigen::MatrixBase<Vector3Type> & S,
                          const Eigen::MatrixBase<Vector6Type> & U_r,
                          const Eigen::MatrixBase<Vector6Type> & c_r,
                          Eigen::MatrixBase<Vector6Type> & P_a_r,
                          Eigen::MatrixBase<Matrix6Type> & M_a_r,
                          const Eigen::MatrixBase<Vector6Type> & P_A_r,
                          const Eigen::MatrixBase<Matrix6Type> & M_A_r) {

    EIGEN_STATIC_ASSERT(Vector6Type::ColsAtCompileTime == 1,
                        YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_A_MATRIX);

    //! Ma
    typename Vector6Type::PlainObject & U_ = const_cast<Vector6Type &>(U_r.derived());
    typename Matrix6Type::PlainObject & M_A_ = const_cast<Matrix6Type &>(M_A_r.derived());

    typename Vector6Type::PlainObject UD_ = U_*iD;
    M_a_r.noalias() = M_A_ - U_*UD_.transpose();

    //! Pa
    typename Vector6Type::PlainObject & P_A_ = const_cast<Vector6Type &>(P_A_r.derived());
    typename Vector6Type::PlainObject & c_ = const_cast<Vector6Type &>(c_r.derived());

    P_a_r.noalias() = P_A_ + M_a_r*c_ + u*UD_;
  }


  //! Implementation for inertial back-projection.
  template<typename Vector3Type, typename Matrix3Type, typename Vector6Type, typename Matrix6Type>
  inline void
  JointTypeRy::runInerProj(bool P_z,
                           const Eigen::MatrixBase<Vector3Type> & S,
                           const Eigen::MatrixBase<Vector3Type> & P_,
                           const Eigen::MatrixBase<Matrix3Type> & R_,
                           const Eigen::MatrixBase<Vector6Type> & P_a_r,
                           Eigen::MatrixBase<Vector6Type> & P_A_r,
                           const Eigen::MatrixBase<Matrix6Type> & M_a_r,
                           Eigen::MatrixBase<Matrix6Type> & M_A_r) {

    //! Ad*M*Ad
    typedef const Eigen::Block<Matrix6Type,3,3> constBlock3;
    typedef Eigen::Block<Matrix6Type,3,3> Block3;

    typename Matrix6Type::PlainObject & M_a_ = const_cast<Matrix6Type &>(M_a_r.derived());
    constBlock3 & Ai = M_a_.template block<3,3>(0,0);
    constBlock3 & Bi = M_a_.template block<3,3>(0,3);
    constBlock3 & Di = M_a_.template block<3,3>(3,3);

//    typename GEOMBD_EIGEN_PLAIN_TYPE(Matrix6Type) Mtmp;
    typename Matrix6Type::PlainObject Mtmp;

    Block3 Ao = Mtmp.template block<3,3>(0,0);
    Block3 Bo = Mtmp.template block<3,3>(0,3);
    Block3 Co = Mtmp.template block<3,3>(3,0);
    Block3 Do = Mtmp.template block<3,3>(3,3);

    typedef typename Matrix3Type::Scalar MyScalar;
    MyScalar sq, cq;
    sq = R_.coeff(0,2);  cq = R_.coeff(0,0);

    ::geo::Mat3ProjRy(sq, cq, Ai, Ao);

    Do.noalias() = R_*Bi; // tmp variable
    Co.noalias() = Do*R_.transpose();

    ::geo::Mat3ProjRy(sq, cq, Di, Do);

    Bo = Co;
    if(P_z){
        //! Linear
        typename Matrix3Type::PlainObject SkP = Skew(P_);
        Bo.noalias() -= Ao*SkP;

        typename Matrix3Type::PlainObject Dtmp1;
        Dtmp1.noalias() = SkP*Bo;
        Do.noalias() += Dtmp1;
        Do.noalias() -= Co.transpose()*SkP;
      }

    Co = Bo.transpose();

    M_A_r.noalias() += Mtmp;


    //! Ad*Xi*Ad
    typedef const Eigen::Block<Vector6Type,3,1> constSegment3;
    typedef Eigen::Block<Vector6Type,3,1> Segment3;

    Vector6Type & P_a_ = const_cast<Vector6Type &>(P_a_r.derived());
    constSegment3 & Pa_up   = P_a_.template segment<3>(0);
    constSegment3 & Pa_down = P_a_.template segment<3>(3);

    typename Vector6Type::PlainObject P_A_i;

    Segment3 Up_   = P_A_i.template segment<3>(0);
    Segment3 Down_ = P_A_i.template segment<3>(3);

    Up_.noalias()   = R_*Pa_up;
    Down_.noalias() = R_*Pa_down;

    if(P_z) Down_.noalias() += P_.cross(Up_);  // If P != 0

    P_A_r.noalias() += P_A_i;
  }


  //! Implementation for acceleration expressions.
  template<typename ScalarType, typename Vector3Type, typename Matrix3Type, typename Vector6Type>
  inline void
  JointTypeRy::runAccel(bool zeroFlag,
                        ScalarType u,
                        ScalarType iD,
                        ScalarType* ddq,
                        const Eigen::MatrixBase<Vector3Type> & S,
                        const Eigen::MatrixBase<Vector3Type> & P_r,
                        const Eigen::MatrixBase<Matrix3Type> & R_r,
                        const Eigen::MatrixBase<Vector6Type> & c_r,
                        const Eigen::MatrixBase<Vector6Type> & U_r,
                        Eigen::MatrixBase<Vector6Type> & Acc_i_r,
                        const Eigen::MatrixBase<Vector6Type> & Acc_j_r) {

    EIGEN_STATIC_ASSERT(Vector6Type::ColsAtCompileTime == 1,
                        YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_A_MATRIX);

    //! Acceleration bias.
    //!------------------------------------------------------------------------------!//
    typedef const Eigen::Block<Vector6Type,3,1> constSegment3;
    typedef Eigen::Block<Vector6Type,3,1> Segment3;

    Vector6Type & Acc_j_ = const_cast<Vector6Type &>(Acc_j_r.derived());
    constSegment3 & AjUp = Acc_j_.template segment<3>(0);
    constSegment3 & AjDown = Acc_j_.template segment<3>(3);

    typename Vector6Type::PlainObject Acc_a;

    Segment3 AaUp = Acc_a.template segment<3>(0);
    Segment3 AaDown = Acc_a.template segment<3>(3);

    if(zeroFlag) {
        AaDown.noalias() = AjUp - P_r.cross(AjDown);  //tmp
        AaUp.noalias()   = R_r.transpose()*AaDown;
      } else {
        AaUp.noalias()   = R_r.transpose()*AjUp;
      }

    AaDown.noalias() = R_r.transpose()*AjDown;

    Acc_a.noalias() += c_r;

    //! Joint acceleration.
    //!------------------------------------------------------------------------------!//
    (*ddq) = iD * ( u - U_r.transpose()*Acc_a );

    //! Update spatial acceleration.
    //!------------------------------------------------------------------------------!//
    Acc_i_r = Acc_a;
    Acc_i_r(4) += (*ddq);

  }


  //! Implementation for acceleration expressions at root.
  template<typename ScalarType, typename Vector3Type, typename Matrix3Type, typename Vector6Type>
  inline void
  JointTypeRy::runAccelRoot(ScalarType u,
                            ScalarType iD,
                            ScalarType* ddq,
                            const Eigen::MatrixBase<Vector3Type> & S,
                            const Eigen::MatrixBase<Vector3Type> & P_r,
                            const Eigen::MatrixBase<Matrix3Type> & R_r,
                            const Eigen::MatrixBase<Vector6Type> & U_r,
                            Eigen::MatrixBase<Vector6Type> & Acc_i_r) {

    EIGEN_STATIC_ASSERT(Vector6Type::ColsAtCompileTime == 1,
                        YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_A_MATRIX);

    //! Acceleration bias.
    //!------------------------------------------------------------------------------!//
    //! g = [ 0, 0, 9.81, 0, 0, 0 ]^T
    //! Acc_a = Ad(G[Sy])*g = [ -9.81*sq, 0, 9.81*cq, 0, 0, 0 ]^T
    ScalarType sq, cq;  sq = R_r.coeff(0,2);  cq = R_r.coeff(0,0);
    Acc_i_r.setZero();
    Acc_i_r.coeffRef(0) = -9.81*sq;  Acc_i_r.coeffRef(2) = 9.81*cq;

    //! Joint acceleration.
    //!------------------------------------------------------------------------------!//
    ScalarType ddq_;
    ddq_ = u - U_r.transpose()*Acc_i_r;
    ddq_ *= iD;
    (*ddq) = ddq_;

    //! Update spatial acceleration.
    //!------------------------------------------------------------------------------!//
    Acc_i_r.coeffRef(4) = ddq_;
  }


} // end namespace geoCRTP

#endif // GEOMBD_JOINT_TYPE_REVOLUTE_Y_HXX
