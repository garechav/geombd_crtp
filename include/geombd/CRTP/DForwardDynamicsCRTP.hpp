/**
 *    \file include/geombd/CRTP/DForwardDynamicsCRTP.hpp
 *    \author Alvaro Paz, Gustavo Arechavaleta
 *    \version 1.0
 *    \date 2021
 *
 *    Class to implement the Articulated-Body Algorithm derivative wrt state
 *    Copyright (c) 2021 Cinvestav
 *    This library is distributed under the MIT License.
 */

#ifndef GEOMBD_FORWARD_DYNAMICS_DIFFERENTIATION_CRTP_HPP
#define GEOMBD_FORWARD_DYNAMICS_DIFFERENTIATION_CRTP_HPP

#define EIGEN_NO_DEBUG
#define EIGEN_MPL2_ONLY
#define EIGEN_UNROLLING_LIMIT 30

//#include "geombd/CRTP/utils/sincos.hxx"
//#include "geombd/CRTP/utils/functors.hxx"
//#include "geombd/CRTP/traits/GeoOperator"

#include "geombd/CRTP/DJointBase.hxx"
#include "geombd/pinocchio/container/aligned-vector.hpp"
#include "geombd/io/parser.hpp"

namespace geo {

  const int maxBody = 32;

  template<typename ScalarType>
  class FwdDynDifCRTP
  {
    //!------------------------------------------------------------------------------!//
    //! Types
    //!------------------------------------------------------------------------------!//

    typedef int index_t;
    typedef ScalarType real_t;

    typedef Eigen::Matrix<index_t, 6, 1> Vector6int;

    typedef Eigen::Matrix<real_t, 6, 1> SpatialVector;
    typedef Eigen::Matrix<real_t, 6, 6> SpatialMatrix;

    typedef Eigen::Matrix<index_t, Eigen::Dynamic, 1, 0, maxBody, 1> VectorXint;
    typedef Eigen::Matrix<real_t, Eigen::Dynamic, 1, 0, 6*maxBody, 1>  VectorXr;

    typedef Eigen::Matrix<real_t , 2 , 1> Vector2r;
    typedef Eigen::Matrix<real_t , 2 , 2> Matrix2r;

    typedef Eigen::Matrix<real_t , 3 , 1> Vector3r;
    typedef Eigen::Matrix<real_t , 3 , 3> Matrix3r;

    typedef Eigen::Matrix<real_t , 4 , 1> Vector4r;
    typedef Eigen::Matrix<real_t , 4 , 4> Matrix4r;

    //!------------------------------------------------------------------------------!//

    typedef Eigen::Matrix<real_t, 6, Eigen::Dynamic, 0, 6, maxBody> D_SpatialVector;
    typedef Eigen::Matrix<real_t, Eigen::Dynamic, 6, Eigen::RowMajor, 6*maxBody, 6> D_SpatialMatrix;

    typedef Eigen::Matrix<real_t, 1, Eigen::Dynamic, Eigen::RowMajor, 1, maxBody> RowVectorXr;

    typedef Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, 0, maxBody, 2*maxBody> MatrixXr;

    //!------------------------------------------------------------------------------!//
    //! Constructors and Destructors
    //!------------------------------------------------------------------------------!//

  public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW //128-bit alignment


    //! Custom constructor
    FwdDynDifCRTP( const std::shared_ptr< Robot >& robot )
    {

      EIGEN_MAKE_ALIGNED_OPERATOR_NEW

          this->n = robot->nq;

      //! Pre-allocation for enhancement
      //!------------------------------------------------------------------------------!//
      ddq = VectorXr::Zero(n);
      PRE.resize(n);  SUC.resize(n);  PRESUC.resize(n);
      PRE_n.resize(n);  SUC_n.resize(n);  PRESUC_n.resize(n);
      this->NP = VectorXint::Zero(n,1);  this->NS = VectorXint::Zero(n,1);

      for ( index_t i = 1; i < n+1; i++ )
        {
          auto body = robot->operator [](i);
          auto joint = robot->operator ()(body->parent_joints[0]);
          j = body->parent_links[0];

          TransConst.push_back(joint->origin.xyz);
          RotConst.push_back(Matrix3r::Identity());
          Trans.push_back(Vector3r::Zero());
          Rot.push_back(Matrix3r::Identity());
          Twists.push_back(SpatialVector::Zero());
          Cbias.push_back(SpatialVector::Zero());
          Pbias.push_back(SpatialVector::Zero());
          MConst.push_back(body->inertial_properties.spatial);
          M_a.push_back(SpatialMatrix::Zero());
          M_A.push_back(SpatialMatrix::Zero());
          P_a.push_back(SpatialVector::Zero());
          P_A.push_back(SpatialVector::Zero());

          u.push_back(0);  invD.push_back(0);
          U.push_back(SpatialVector::Zero());

          Accel.push_back(SpatialVector::Zero());

          parent.push_back(j-1);

          real_t Pi_sum = TransConst.back().cwiseAbs().sum();
          if(Pi_sum == 0) {
              P_zero.push_back(false);
            } else {
              P_zero.push_back(true);
            }

          //! Robot topology
          //!------------------------------------------------------------------------------!//
          Pre = body->Pre;         Suc = body->Suc;
          nP = Pre.size();         nS = Suc.size();  //! Both consider the i-th body ID
          NP.coeffRef(i-1) = nP;  NS.coeffRef(i-1) = nS;

          Suc.erase(Suc.begin());

          PRE.at(i-1) = Pre;  SUC.at(i-1) = Suc;  PreSuc = Pre;
          PreSuc.insert(std::end(PreSuc), std::begin(Suc), std::end(Suc));
          PRESUC.at(i-1) = PreSuc;

          if(nP == 1) {
              rootFlag.push_back(false);
            } else {
              rootFlag.push_back(true);
            }

          Pre_n = VectorXint::Zero(Pre.size(), 1);
          for(int w_ = 0; w_ < Pre.size(); w_++) {
              Pre_n(w_) = Pre[w_];
            }
          PRE_n.at(i-1) = Pre_n;

          Suc_n = VectorXint::Zero(Suc.size(), 1);
          for(int w_ = 0; w_ < Suc.size(); w_++) {
              Suc_n(w_) = Suc[w_];
            }
          SUC_n.at(i-1) = Suc_n;

          PreSuc_n = VectorXint::Zero(PreSuc.size(), 1);
          for(int w_ = 0; w_ < PreSuc.size(); w_++) {
              PreSuc_n(w_) = PreSuc[w_];
            }
          PRESUC_n.at(i-1) = PreSuc_n;

          Vector6int auxInd;
          auto body_parent = robot->operator [](body->parent_links[0]);

          index_t n_child_p = body_parent->child_links.size();
          if (j)  auxInd << i-1, j, nP, nS, NP(j-1), n_child_p;
          else  auxInd << i-1, j, nP, nS, -1, n_child_p;

          indexVec.push_back( auxInd );

          if (body->child_links.size() > 1)  branchFlag.push_back(true);
          else branchFlag.push_back(false);

          if (nS-1) isLeaf.push_back(true);
          else isLeaf.push_back(false);

          //! Fill containers with proper-dimensions matrices
          //!------------------------------------------------------------------------------!//
          D_SpatialVector D_SVaux = D_SpatialVector::Zero(6, nP);
          D_q_V.push_back( D_SVaux );    D_dq_V.push_back( D_SVaux );
          D_q_c.push_back( D_SVaux );    D_dq_c.push_back( D_SVaux );
          D_q_p.push_back( D_SVaux );    D_dq_p.push_back( D_SVaux );

          D_U_h.push_back( D_SpatialVector::Zero(6, nS-1) );
          D_U_v.push_back( VectorXr::Zero(6*(nS-1), 1) );
          D_invD.push_back( RowVectorXr::Zero(1, nS-1) );
          D_M_A.push_back( D_SpatialMatrix::Zero(6*(nS-1), 6) );

          D_q_u.push_back( RowVectorXr::Zero(1, nP+nS-1) );        D_dq_u.push_back( RowVectorXr::Zero(1, nP+nS-1) );
          D_SVaux = D_SpatialVector::Zero(6, nP+nS-1);
          D_q_Pa.push_back( D_SVaux );   D_dq_Pa.push_back( D_SVaux );
          D_q_PA.push_back( D_SVaux );   D_dq_PA.push_back( D_SVaux );

          D_q_A.push_back( D_SpatialVector::Zero(6, n) );   D_dq_A.push_back( D_SpatialVector::Zero(6, n) );

          D_ddq = MatrixXr::Zero(n,2*n);

          //! Handling screw axes
          //!------------------------------------------------------------------------------!//
          Vector3r axis = joint->axis;
          Screw_w.push_back(axis);

          if (joint->type == Joint::joint_type::prismatic) {
              D_dq_V[i-1].template rightCols<1>().template segment<3>(0) = axis;
              //! code here for Joint Types Px Py and Pz
              //!------------------------------------------------------------------------------!//
            }
          if (joint->type == Joint::joint_type::revolute) {
              D_dq_V[i-1].template rightCols<1>().template segment<3>(3) = axis;
              //!------------------------------------------------------------------------------!//
              if((real_t)axis(0) == 1) {
//                  JointTypesVec.push_back( D_JointTypeRx{} );
                }
              //!------------------------------------------------------------------------------!//
              if((real_t)axis(1) == 1) {
//                  JointTypesVec.push_back( D_JointTypeRy{} );
                }
              //!------------------------------------------------------------------------------!//
              if((real_t)axis(2) == 1) {
                  JointTypesVec.push_back( D_JointTypeRz{} );
                }
              //!------------------------------------------------------------------------------!//
              if((real_t)axis(1) != 0 && (real_t)axis(2) != 0) {
//                  JointTypesVec.push_back( D_JointTypeRxyz{} );
                }

            }

        }

    }


    ~FwdDynDifCRTP(){}

    //!------------------------------------------------------------------------------!//
    //! Members
    //!------------------------------------------------------------------------------!//

  public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Dof
    index_t n;

    index_t ID, IDj, j;

    SpatialVector Screw;

    std::vector< bool > P_zero, rootFlag, branchFlag, isLeaf;

    std::vector< index_t > parent;

    //! Joint-types variant
    typedef boost::variant< D_JointTypeRz > CV;

    //! Joint-types vector
    std::vector< CV > JointTypesVec;
    typename std::vector< CV >::iterator JointTypesIter;

    D_FwdKin_visitor<ScalarType, Vector3r, Matrix3r> visitorFK;
    D_TCP01_visitor<ScalarType, Vector3r, SpatialVector> visitorTCP01;
    D_TCP02_visitor<Vector3r, D_SpatialVector> visitorTCP02;
    D_TCProot_visitor<ScalarType, Vector3r, SpatialVector, SpatialMatrix, D_SpatialVector> visitorTCProot;

    Inertia01_visitor<Vector3r, SpatialVector, SpatialMatrix, VectorXr, D_SpatialMatrix> visitorInertia01;
    Inertia02_visitor<ScalarType, Vector3r, SpatialVector, RowVectorXr, D_SpatialVector> visitorInertia02;
    Inertia03_visitor<index_t, Vector3r, Matrix3r, SpatialVector, SpatialMatrix, D_SpatialMatrix> visitorInertia03;
    Leaf01_visitor<ScalarType, Vector3r, SpatialVector, SpatialMatrix, RowVectorXr, D_SpatialVector> visitorLeaf01;
    Leaf02_visitor<Vector3r, Matrix3r, SpatialVector, SpatialMatrix, D_SpatialMatrix> visitorLeaf02;

    Accel01_visitor<Vector3r, SpatialVector> visitorAccel01;
    Accel02_visitor<ScalarType, Vector3r, SpatialVector, RowVectorXr, D_SpatialVector> visitorAccel02;
    AccelRoot_visitor<ScalarType, Vector3r, Matrix3r, SpatialVector, D_SpatialVector, RowVectorXr, MatrixXr> visitorRoot;


    PINOCCHIO_ALIGNED_STD_VECTOR( Vector3r ) Trans;
    PINOCCHIO_ALIGNED_STD_VECTOR( Matrix3r ) Rot;
    PINOCCHIO_ALIGNED_STD_VECTOR( Vector3r ) TransConst;
    PINOCCHIO_ALIGNED_STD_VECTOR( Matrix3r ) RotConst;
    PINOCCHIO_ALIGNED_STD_VECTOR( Vector3r ) Screw_w;
    PINOCCHIO_ALIGNED_STD_VECTOR( SpatialVector ) Twists;
    PINOCCHIO_ALIGNED_STD_VECTOR( SpatialVector ) Cbias;
    PINOCCHIO_ALIGNED_STD_VECTOR( SpatialVector ) Pbias;
    PINOCCHIO_ALIGNED_STD_VECTOR( SpatialMatrix ) MConst;
    PINOCCHIO_ALIGNED_STD_VECTOR( SpatialMatrix ) M_a;
    PINOCCHIO_ALIGNED_STD_VECTOR( SpatialMatrix ) M_A;
    PINOCCHIO_ALIGNED_STD_VECTOR( SpatialVector ) P_a;
    PINOCCHIO_ALIGNED_STD_VECTOR( SpatialVector ) P_A;

    PINOCCHIO_ALIGNED_STD_VECTOR( SpatialVector ) U;
    PINOCCHIO_ALIGNED_STD_VECTOR( real_t ) u;
    PINOCCHIO_ALIGNED_STD_VECTOR( real_t ) invD;

    PINOCCHIO_ALIGNED_STD_VECTOR( SpatialVector ) Accel;
    VectorXr ddq;

    //!------------------------------------------------------------------------------!//

    index_t nP, nS, PS, nPj, nSj;
    VectorXint NP, NS;

    std::vector< index_t > Pre, Suc, PreSuc;
    std::vector< std::vector< index_t > > PRE, SUC, PRESUC;

    VectorXint Pre_n, Suc_n, PreSuc_n;
    std::vector< VectorXint > PRE_n, SUC_n, PRESUC_n;

    SpatialMatrix D_p_aux;

    PINOCCHIO_ALIGNED_STD_VECTOR( D_SpatialVector ) D_q_V;   PINOCCHIO_ALIGNED_STD_VECTOR( D_SpatialVector ) D_dq_V;
    PINOCCHIO_ALIGNED_STD_VECTOR( D_SpatialVector ) D_q_c;   PINOCCHIO_ALIGNED_STD_VECTOR( D_SpatialVector ) D_dq_c;


    //    std::vector< D_RowSpatialVector, Eigen::aligned_allocator<D_SpatialVector> > D_q_c;
    //    std::vector< D_RowSpatialVector, Eigen::aligned_allocator<D_SpatialVector> > D_dq_c;


    PINOCCHIO_ALIGNED_STD_VECTOR( Vector6int ) indexVec;
    PINOCCHIO_ALIGNED_STD_VECTOR( D_SpatialVector ) D_U_h;   PINOCCHIO_ALIGNED_STD_VECTOR( VectorXr ) D_U_v;  // D_U horizontally and vertically arranged
    PINOCCHIO_ALIGNED_STD_VECTOR( RowVectorXr ) D_invD;
    PINOCCHIO_ALIGNED_STD_VECTOR( D_SpatialMatrix ) D_M_A;

    PINOCCHIO_ALIGNED_STD_VECTOR( RowVectorXr ) D_q_u;       PINOCCHIO_ALIGNED_STD_VECTOR( RowVectorXr ) D_dq_u;
    PINOCCHIO_ALIGNED_STD_VECTOR( D_SpatialVector ) D_q_p;  PINOCCHIO_ALIGNED_STD_VECTOR( D_SpatialVector ) D_dq_p;
    PINOCCHIO_ALIGNED_STD_VECTOR( D_SpatialVector ) D_q_Pa;  PINOCCHIO_ALIGNED_STD_VECTOR( D_SpatialVector ) D_dq_Pa;
    PINOCCHIO_ALIGNED_STD_VECTOR( D_SpatialVector ) D_q_PA;  PINOCCHIO_ALIGNED_STD_VECTOR( D_SpatialVector ) D_dq_PA;

    PINOCCHIO_ALIGNED_STD_VECTOR( D_SpatialVector ) D_q_A;   PINOCCHIO_ALIGNED_STD_VECTOR( D_SpatialVector ) D_dq_A;

    MatrixXr D_ddq;

    //!------------------------------------------------------------------------------!//
    //! Methods
    //!------------------------------------------------------------------------------!//

  public:

    //!------------------------------------------------------------------------------!//
    //! Articulated-Body Algorithm
    //!------------------------------------------------------------------------------!//
    template<typename ConfigVectorT, typename TangentVectorT, typename CotangentVectorT>
    const void
    D_aba(const Eigen::MatrixBase<ConfigVectorT> & q,
          const Eigen::MatrixBase<TangentVectorT> & v,
          const Eigen::MatrixBase<CotangentVectorT> & tau)
    {
      typedef typename ConfigVectorT::Scalar ScalarT;

      //! First Recursion
      ID = 0;
      for ( JointTypesIter = JointTypesVec.begin(); JointTypesIter != JointTypesVec.end(); JointTypesIter++ ){
          //!------------------------------------------------------------------------------!//
          //! Forward kinematics
          //!------------------------------------------------------------------------------!//
          //! Call CRTP visitor.
          //!------------------------------------------------------------------------------!//
          visitorFK.qi = q(ID);  visitorFK.S = &Screw_w[ID];  visitorFK.R = &Rot[ID];
          boost::apply_visitor( visitorFK, *JointTypesIter );

          IDj = parent[ID];
          //! Twist and c, p bias
          //!------------------------------------------------------------------------------!//
          if(ID > 0){
              //! Twist
              //!------------------------------------------------------------------------------!//
              typedef const Eigen::Block<SpatialVector,3,1> constSegment3;
              typedef Eigen::Block<SpatialVector,3,1> Segment3;

              SpatialVector & S_ = const_cast<SpatialVector &>(Twists[IDj].derived());
              constSegment3 & Sup   = S_.template segment<3>(0);
              constSegment3 & Sdown = S_.template segment<3>(3);

              Segment3 Up_   = Twists[ID].template segment<3>(0);
              Segment3 Down_ = Twists[ID].template segment<3>(3);

              int icol = D_q_V[ID].cols();

              //!------------------------------------------------------------------------------!//
              static D_SpatialVector D_q_V_i_, D_dq_V_i_;

              if(icol == 2) {
                  D_q_V_i_ = D_SpatialVector::Zero(6,1);
                  D_dq_V_i_ = D_dq_V[IDj];
                }
              if(icol == 7) {
                  D_q_V_i_ = D_q_V[IDj];
                  D_dq_V_i_ = D_dq_V[IDj];
                }

              //!------------------------------------------------------------------------------!//
              if(P_zero[ID]) {
                  Matrix3r SkP = Skew(TransConst[ID]);
                  Down_.noalias() = Sup - SkP*Sdown;  //tmp
                  Up_.noalias()   = Rot[ID].transpose()*Down_;
                  //!------------------------------------------------------------------------------!//
                  D_q_V_i_.template topRows<3>().noalias() -= SkP*D_q_V_i_.template bottomRows<3>();
                  D_q_V_i_.template topRows<3>().applyOnTheLeft(Rot[ID].transpose());

                  D_dq_V_i_.template topRows<3>().noalias() -= SkP*D_dq_V_i_.template bottomRows<3>();
                  D_dq_V_i_.template topRows<3>().applyOnTheLeft(Rot[ID].transpose());
                } else {
                  Up_.noalias() = Rot[ID].transpose()*Sup;
                  //!------------------------------------------------------------------------------!//
                  D_q_V_i_.template topRows<3>().applyOnTheLeft(Rot[ID].transpose());

                  D_dq_V_i_.template topRows<3>().applyOnTheLeft(Rot[ID].transpose());
                }

              Down_.noalias() = Rot[ID].transpose()*Sdown;
              //!------------------------------------------------------------------------------!//
              D_q_V_i_.template bottomRows<3>().applyOnTheLeft(Rot[ID].transpose());
              D_dq_V_i_.template bottomRows<3>().applyOnTheLeft(Rot[ID].transpose());

              //! Call CRTP visitor.
              //!------------------------------------------------------------------------------!//
              visitorTCP01.vi = v(ID);  visitorTCP01.Sw_ = &Screw_w[ID];  visitorTCP01.S_ = &Twists[ID];  visitorTCP01.C_ = &Cbias[ID];
              boost::apply_visitor( visitorTCP01, *JointTypesIter );

              //! Twist differentiation -> taking advantage of c bias computation.
              //!------------------------------------------------------------------------------!//
              D_q_V_i_.resize(Eigen::NoChange, icol);
              D_q_V_i_.template rightCols<1>() = Cbias[ID];

              D_dq_V_i_.resize(Eigen::NoChange, icol);
              D_dq_V_i_.template rightCols<1>().template segment<3>(3) = Screw_w[ID];

              //! Store
              if(icol==6) {
                  D_q_V[ID] = D_q_V_i_;
                  D_dq_V[ID] = D_dq_V_i_;
                }

              Cbias[ID] *= v(ID);

              //! C bias differentiation -> taking advantage of twist computation.
              //!------------------------------------------------------------------------------!//
              D_SpatialVector D_q_c_aux, D_dq_c_aux;
              D_q_c_aux.noalias() = D_q_V_i_*v(ID);

              D_dq_c_aux.noalias() = D_dq_V_i_*v(ID);
              D_dq_c_aux.template rightCols<1>() += Twists[ID];

              //! Call CRTP visitor.
              //!------------------------------------------------------------------------------!//
              visitorTCP02.D_q_c_ = &D_q_c[ID];       visitorTCP02.D_dq_c_ = &D_dq_c[ID];     visitorTCP02.Sw_ = &Screw_w[ID];
              visitorTCP02.D_q_c_aux_ = &D_q_c_aux;   visitorTCP02.D_dq_c_aux_ = &D_dq_c_aux;
              boost::apply_visitor( visitorTCP02, *JointTypesIter );

              //! P bias
              //!------------------------------------------------------------------------------!//
              typedef const Eigen::Block<SpatialMatrix,3,3> constBlock3;
              typedef Eigen::Block<SpatialMatrix,3,3> Block3;

              SpatialMatrix & M = const_cast<SpatialMatrix &>(MConst[ID].derived());
              constBlock3 & UL = M.template block<3,3>(0,0);
              constBlock3 & UR = M.template block<3,3>(0,3);
              constBlock3 & LL = M.template block<3,3>(3,0);
              constBlock3 & LR = M.template block<3,3>(3,3);

              typename Vector3r::PlainObject muUp, muDo;

              Segment3 p_Up = Pbias[ID].template segment<3>(0);
              Segment3 p_Do = Pbias[ID].template segment<3>(3);

              //! M is symmetric and its shape is known
              muUp.noalias() = UL*Up_;  muUp.noalias() += UR*Down_;
              muDo.noalias() = LL*Up_;  muDo.noalias() += LR*Down_;

              //! TODO: Cross is very expensive, consider enable staticness with traits expressions
              p_Up.noalias() = Down_.cross(muUp);
              p_Do.noalias() = Up_.cross(muUp);
              p_Do.noalias() += Down_.cross(muDo);

              //! P bias differentiation
              //!------------------------------------------------------------------------------!//
              SpatialMatrix D_p_, D_aux;

              Block3 UL_ = D_p_.template block<3,3>(0,0);  Block3 UR_ = D_p_.template block<3,3>(0,3);
              Block3 LL_ = D_p_.template block<3,3>(3,0);  Block3 LR_ = D_p_.template block<3,3>(3,3);

              //! adDual
              UL_.noalias() = Skew(Down_);  UR_.setZero();
              LL_.noalias() = Skew(Up_);    LR_.noalias() = UL_;

              D_aux.noalias() = D_p_*M.derived();  // tmp

              D_p_.noalias() = D_aux;

              //! adBar
              UR_.noalias() -= Skew(muUp);
              LL_.noalias() -= Skew(muUp);  LR_.noalias() -= Skew(muDo);

              D_q_p[ID].noalias()  = D_p_*D_q_V_i_;
              D_dq_p[ID].noalias() = D_p_*D_dq_V_i_.derived();
            } else {
              //! Call CRTP visitor.
              //!------------------------------------------------------------------------------!//
              visitorTCProot.vi = v(ID);       visitorTCProot.S_ = &Screw_w[ID];  visitorTCProot.T_ = &Twists[ID];
              visitorTCProot.p_ = &Pbias[ID];  visitorTCProot.M_ = &MConst[ID];   visitorTCProot.D_dq_p_ = &D_dq_p[ID];
              boost::apply_visitor( visitorTCProot, *JointTypesIter );
            }
          if (branchFlag[ID]) {
              if (isLeaf[ID]) {
                  D_q_PA[ID].noalias() = D_q_p[ID];
                  D_dq_PA[ID].noalias() = D_dq_p[ID];
                }
            } else {
              D_q_PA[ID].setZero();   D_q_PA[ID].leftCols( NP.coeffRef(ID) ) = D_q_p[ID];
              D_dq_PA[ID].setZero();  D_dq_PA[ID].leftCols( NP.coeffRef(ID) ) = D_dq_p[ID];
            }
          ID++;
        }


      //! Initialize inertial terms
      M_A = MConst;  P_A = Pbias;

      //! Second Recursion
      ID = n - 1;
      for ( JointTypesIter = JointTypesVec.end()-1; JointTypesIter != JointTypesVec.begin()-1; JointTypesIter-- ) {
          nP = NP(ID);    nS = NS(ID)-1;

          IDj = parent[ID];
          //! If the current body is not a leaf
          //!------------------------------------------------------------------------------!//
          if (nS) {

              index_t nPj, nChild;
              nPj = indexVec[ID].coeff(4);  nChild = indexVec[ID].coeff(5);
              static D_SpatialMatrix D_M_A_i_;

              //! if it's a parent's leaf or has multiple children
              if(nS == 1) {
                  D_M_A_i_ = D_M_A[ID];
                }
              if(ID == 5) {   // change condition to no. of children of i > 1
                  D_M_A_i_ = D_M_A[ID];
                }

              //! Call CRTP visitor.
              //!------------------------------------------------------------------------------!//
              visitorInertia01.Sw_ = &Screw_w[ID];    visitorInertia01.U_ = &U[ID];          visitorInertia01.M_A_ = &M_A[ID];
              visitorInertia01.D_U_v_ = &D_U_v[ID];  visitorInertia01.D_M_A_i_ = &D_M_A_i_;
              boost::apply_visitor( visitorInertia01, *JointTypesIter );

              D_U_h[ID] = Eigen::Map<D_SpatialVector>(D_U_v[ID].derived().data(), 6, nS);

              u[ID] = tau(ID);

              //! Call CRTP visitor.
              //!------------------------------------------------------------------------------!//
              visitorInertia02.invD_ = &invD[ID];      visitorInertia02.u_ = &u[ID];            visitorInertia02.Sw_ = &Screw_w[ID];
              visitorInertia02.U_ = &U[ID];            visitorInertia02.P_A_ = &P_A[ID];
              visitorInertia02.D_invD_ = &D_invD[ID];  visitorInertia02.D_q_u_ = &D_q_u[ID];    visitorInertia02.D_dq_u_ = &D_dq_u[ID];
              visitorInertia02.D_U_h_ = &D_U_h[ID];    visitorInertia02.D_q_PA_ = &D_q_PA[ID];  visitorInertia02.D_dq_PA_ = &D_dq_PA[ID];
              boost::apply_visitor( visitorInertia02, *JointTypesIter );


              //! Differentiating inertial back-projection to parent. Only if parent is not the root.
              //!------------------------------------------------------------------------------!//
              if (rootFlag[ID]) {

                  //! Prepare inertial expression M_a.
                  //!------------------------------------------------------------------------------!//
                  SpatialVector UD_ = U[ID]*invD[ID];
                  M_a[ID] = M_A[ID];
                  M_a[ID].noalias() -= UD_*U[ID].transpose();

                  //! Prepare inertial expression P_a.
                  //!------------------------------------------------------------------------------!//
                  P_a[ID] = P_A[ID];
                  P_a[ID].noalias() += M_a[ID]*Cbias[ID];
                  P_a[ID].noalias() += u[ID]*UD_;

                  //! Prepare inertial expression differentiation D_M_a.
                  //!------------------------------------------------------------------------------!//
                  SpatialMatrix Mtmp, AdjointDual;  VectorXr kronAux;  D_SpatialMatrix D_MtmpI;

                  //!------------------------------------------------------------------------------!//
                  D_MtmpI.noalias() = D_U_v[ID]*UD_.transpose();
                  D_M_A_i_.noalias() -= D_MtmpI;              //! D_M_A_i_ is D_M_a_i_ from this point
                  Mtmp.noalias() = U[ID]*U[ID].transpose();

                  AdDual(TransConst[ID].derived(), Rot[ID].derived(), AdjointDual);
                  index_t ite = 0;
                  for(index_t iter = 0; iter < nS*6; iter += 6){
                      D_M_A_i_.middleRows(iter, 6).noalias() -= D_MtmpI.middleRows(iter, 6).transpose() + D_invD[ID].coeffRef(ite)*Mtmp;
                      ite++;
                    }

                  kronAux.noalias() = D_M_A_i_*Cbias[ID];

                  //! then D_MtmpII = D_M_a and D_MtmpI is AdjointDual*D_M_a*Adjoint

                  //! Prepare inertial expression differentiation D_q_Pa.
                  //!------------------------------------------------------------------------------!//
                  D_SpatialVector D_Vec, D_VecII;
                  D_Vec = Eigen::Map<D_SpatialVector>(kronAux.data(),6,nS);

                  D_Vec.noalias() += (u[ID]*invD[ID])*D_U_h[ID] + (u[ID]*U[ID])*D_invD[ID];

                  D_q_Pa[ID].leftCols(nP).noalias() = M_a[ID]*D_q_c[ID];
                  D_q_Pa[ID].rightCols(nS) = D_Vec;

                  D_q_Pa[ID].noalias() += D_q_PA[ID] + UD_*D_q_u[ID];

                  //! Prepare inertial expression differentiation D_dq_Pa.
                  //!------------------------------------------------------------------------------!//
                  D_dq_Pa[ID].noalias() = D_dq_PA[ID] + UD_*D_dq_u[ID];
                  D_dq_Pa[ID].leftCols(nP).noalias() += M_a[ID]*D_dq_c[ID];

                  //!------------------------------------------------------------------------------!//
                  //!                         INERTIAL BACK PROJECTION                             !//
                  //!------------------------------------------------------------------------------!//

                  //! Back projection of P_a.
                  //!------------------------------------------------------------------------------!//
                  typedef Eigen::Block<SpatialVector,3,1> Segment3;

                  Segment3 Pa_up   = P_a[ID].template segment<3>(0);
                  Segment3 Pa_down = P_a[ID].template segment<3>(3);

                  SpatialVector P_A_i;

                  Segment3 Up_   = P_A_i.template segment<3>(0);
                  Segment3 Down_ = P_A_i.template segment<3>(3);

                  Up_.noalias()   = Rot[ID]*Pa_up;
                  Down_.noalias() = Rot[ID]*Pa_down;

                  if(P_zero[ID]) Down_.noalias() += TransConst[ID].cross(Up_);  // If P != 0

                  P_A[IDj].noalias() += P_A_i;

                  //! Make room for new D_M_i
                  //!------------------------------------------------------------------------------!//
                  D_M_A_i_.resize(D_M_A_i_.rows()+6, Eigen::NoChange);

                  //! Call CRTP visitor.
                  //!------------------------------------------------------------------------------!//
                  visitorInertia03.P_z_ = P_zero[ID];  visitorInertia03.nS_ = nS;         visitorInertia03.Sw_ = &Screw_w[ID];   visitorInertia03.P_ = &TransConst[ID];
                  visitorInertia03.R_ = &Rot[ID];      visitorInertia03.P_a_ = &P_a[ID];  visitorInertia03.P_A_i_ = &P_A_i;
                  visitorInertia03.M_a_ = &M_a[ID];    visitorInertia03.Mtmp_ = &Mtmp;    visitorInertia03.D_M_A_i_ = &D_M_A_i_;
                  boost::apply_visitor( visitorInertia03, *JointTypesIter );

                  M_A[IDj].noalias() += Mtmp;

                  //!------------------------------------------------------------------------------!//
                  //!                                  BRANCHING                                   !//
                  //!------------------------------------------------------------------------------!//
                  if (nChild == 1) {  //! Serial Body
                      //!-------------------------------------------------------------------------  D_MA
                      //            D_M_A_i_.template topRows<6>() = Mtmp;
                      //!-------------------------------------------------------------------------  D_q_PA
                      D_q_PA[IDj].noalias() += AdjointDual*D_q_Pa[ID];

                      D_q_PA[IDj].col(nP-1) += AdjointDual*P_A_i;
                      //!-------------------------------------------------------------------------  D_dq_PA
                      D_dq_PA[IDj].noalias() += AdjointDual*D_dq_Pa[ID];
                    } else {  //! Branched Body
                      //!-------------------------------------------------------------------------  D_MA
                      //            D_M_Aj_.template middleRows<6>((ID-j)*6) = Mtmp;
                      D_M_A[IDj].middleRows((ID-IDj+1)*6, 6*(nS+1)) = D_M_A_i_;
                      //!-------------------------------------------------------------------------  D_q_PA
                      D_Vec.noalias() = AdjointDual*D_q_Pa[ID];
                      D_Vec.leftCols(nPj) += D_q_PA[IDj].leftCols(nPj);

                      D_VecII = D_q_PA[IDj].rightCols(D_q_PA[IDj].cols()-nPj);

                      D_q_PA[IDj].derived().resize(Eigen::NoChange, D_Vec.cols()+D_VecII.cols());
                      D_q_PA[IDj] << D_Vec, D_VecII;

                      D_q_PA[IDj].col(nP-1) += AdjointDual*P_A_i;
                      //!-------------------------------------------------------------------------  D_dq_PA
                      D_Vec.noalias() = AdjointDual*D_dq_Pa[ID];
                      D_Vec.leftCols(nPj) += D_dq_PA[IDj].leftCols(nPj);

                      D_VecII = D_dq_PA[IDj].rightCols(D_dq_PA[IDj].cols()-nPj);

                      D_dq_PA[IDj].derived().resize(Eigen::NoChange, D_Vec.cols()+D_VecII.cols());
                      D_dq_PA[IDj] << D_Vec, D_VecII;
                    }

                }

            } else {

              u[ID] = tau(ID);

              //! Call CRTP visitor.
              //!------------------------------------------------------------------------------!//
              visitorLeaf01.invD_ = &invD[ID];      visitorLeaf01.u_ = &u[ID];            visitorLeaf01.U_ = &U[ID];             visitorLeaf01.Sw_ = &Screw_w[ID];
              visitorLeaf01.P_A_ = &P_A[ID];        visitorLeaf01.M_A_ = &M_A[ID];        visitorLeaf01.D_q_u_ = &D_q_u[ID];
              visitorLeaf01.D_dq_u_ = &D_dq_u[ID];  visitorLeaf01.D_q_PA_ = &D_q_PA[ID];  visitorLeaf01.D_dq_PA_ = &D_dq_PA[ID];
              boost::apply_visitor( visitorLeaf01, *JointTypesIter );

              //! Prepare inertial expressions.
              //!------------------------------------------------------------------------------!//
              //! Ma
              SpatialVector UD_ = U[ID]*invD[ID];
              M_a[ID].noalias() = M_A[ID] - UD_*U[ID].transpose();
              
              //! Pa
              P_a[ID].noalias() = P_A[ID] + M_a[ID]*Cbias[ID] + u[ID]*UD_;

              //! Ad*Xi*
              typedef Eigen::Block<SpatialVector,3,1> Segment3;
              
              Segment3 Pa_up = P_a[ID].template segment<3>(0);
              Segment3 Pa_dw = P_a[ID].template segment<3>(3);
              
              typename SpatialVector::PlainObject P_A_i;
              
              Segment3 Up_   = P_A_i.template segment<3>(0);
              Segment3 Down_ = P_A_i.template segment<3>(3);
              
              Up_.noalias()   = Rot[ID]*Pa_up;
              Down_.noalias() = Rot[ID]*Pa_dw;
              
              if(P_zero[ID]) Down_.noalias() += TransConst[ID].cross(Up_);  // If P != 0
              
              P_A[IDj].noalias() += P_A_i;

              //! Call CRTP visitor.
              //!------------------------------------------------------------------------------!//
              visitorLeaf02.P_z_ = P_zero[ID];  visitorLeaf02.Sw_ = &Screw_w[ID];  visitorLeaf02.P_ = &TransConst[ID];  visitorLeaf02.R_ = &Rot[ID];  visitorLeaf02.P_A_i_ = &P_A_i;
              visitorLeaf02.P_a_ = &P_a[ID];    visitorLeaf02.M_a_ = &M_a[ID];     visitorLeaf02.M_A_j_ = &M_A[IDj];  visitorLeaf02.D_M_A_j_ = &D_M_A[IDj];
              boost::apply_visitor( visitorLeaf02, *JointTypesIter );

              //              typename GEOMBD_EIGEN_PLAIN_TYPE(SpatialMatrix) AdjointDual;
              SpatialMatrix AdjointDual;

              //! P_a differentiation -> D_q_Pa and D_dq_Pa.
              //!------------------------------------------------------------------------------!//
              D_q_Pa[ID] = D_q_p[ID];
              D_q_Pa[ID].noalias() += M_a[ID]*D_q_c[ID];
              D_q_Pa[ID].noalias() += UD_*D_q_u[ID];

              D_dq_Pa[ID] = D_dq_p[ID];
              D_dq_Pa[ID].noalias() += M_a[ID]*D_dq_c[ID];
              D_dq_Pa[ID].noalias() += UD_*D_dq_u[ID];

              //! D_q_Pa back projection.
              //!------------------------------------------------------------------------------!//
              AdDual(TransConst[ID].derived(), Rot[ID].derived(), AdjointDual);
              D_q_PA[IDj].template rightCols<1>() = AdjointDual*P_A_i;
              D_q_PA[IDj].noalias() += AdjointDual*D_q_Pa[ID];
              
              //! D_dq_Pa back projection.
              //!------------------------------------------------------------------------------!//
              D_dq_PA[IDj].noalias() += AdjointDual*D_dq_Pa[ID];

            }
          ID--;
        }

      visitorRoot.D_ddq_ = &D_ddq;

      //! Third Recursion
      ID = 0;
      for ( JointTypesIter = JointTypesVec.begin(); JointTypesIter != JointTypesVec.end(); JointTypesIter++ ){
          //! Spatial acceleration
          IDj = parent[ID];
          if( ID ){

              //! Acceleration bias and its derivative.
              //!------------------------------------------------------------------------------!//
              typedef const Eigen::Block<SpatialVector,3,1> constSegment3;
              typedef Eigen::Block<SpatialVector,3,1> Segment3;

              constSegment3 & AjUp = Accel[IDj].template segment<3>(0);
              constSegment3 & AjDown = Accel[IDj].template segment<3>(3);

              typename SpatialVector::PlainObject Aa_, AdAj;

              Segment3 AaUp = Aa_.template segment<3>(0);
              Segment3 AaDown = Aa_.template segment<3>(3);

              if(P_zero[ID]) {
                  //! delete cross product
                  AaDown.noalias() = AjUp - TransConst[ID].cross(AjDown);  //tmp
                  AaUp.noalias()   = Rot[ID].transpose()*AaDown;
                } else {
                  AaUp.noalias()   = Rot[ID].transpose()*AjUp;
                }

              AaDown.noalias() = Rot[ID].transpose()*AjDown;

              //! Call CRTP visitor.
              //!------------------------------------------------------------------------------!//
              visitorAccel01.Sw_ = &Screw_w[ID];  visitorAccel01.AdAj_ = &AdAj;  visitorAccel01.Aa_ = &Aa_;
              boost::apply_visitor( visitorAccel01, *JointTypesIter );

              Aa_.noalias() += Cbias[ID];
              //!------------------------------------------------------------------------------!//
              Eigen::Matrix<ScalarType, 6, 6> AdjointDual;
              AdDual(TransConst[ID].derived(), Rot[ID].derived(), AdjointDual);
              typename D_SpatialVector::PlainObject D_q_Aa_, D_dq_Aa_;

              D_q_Aa_.noalias() = AdjointDual.transpose()*D_q_A[IDj];
              D_q_c[ID].template rightCols<1>() += AdAj;

              for(int j=0; j < PRE_n[ID].size(); j++)
                D_q_Aa_.col( PRE_n[ID](j) ) += D_q_c[ID].col(j);

              D_dq_Aa_.noalias() = AdjointDual.transpose()*D_dq_A[IDj];

              for(int j=0; j < PRE_n[ID].size(); j++)
                D_dq_Aa_.col( PRE_n[ID](j) ) += D_dq_c[ID].col(j);


              //! Joint acceleration and its derivatives D_q_ddq and D_dq_ddq.
              //!------------------------------------------------------------------------------!//
              ScalarType ddq__, ddq_;
              ddq__ = u[ID] - U[ID].transpose()*Aa_;  ddq_ = invD[ID]*ddq__;
              ddq(ID) = ddq_;
              //!-------------------------------------------------------
              typename RowVectorXr::PlainObject D_q_ddq, D_dq_ddq;
              D_q_ddq = -invD[ID]*U[ID].transpose()*D_q_Aa_;

              for(int j=0; j < PRESUC_n[ID].size(); j++)
                D_q_ddq( PRESUC_n[ID](j) ) += invD[ID]*D_q_u[ID](j);

              D_dq_ddq = -invD[ID]*U[ID].transpose()*D_dq_Aa_;

              for(int j=0; j < PRESUC_n[ID].size(); j++)
                D_dq_ddq( PRESUC_n[ID](j) ) += invD[ID]*D_dq_u[ID](j);


              //! Special treatment when body is not a leaf.
              //!------------------------------------------------------------------------------!//
              if (isLeaf[ID]) {

                  for(int j=0; j < SUC_n[ID].size(); j++)
                    D_q_ddq( SUC_n[ID](j) ) += ddq__*D_invD[ID](j) - invD[ID]*Aa_.transpose()*D_U_h[ID].col(j);

                  Accel[ID] = Aa_;  D_q_A[ID] = D_q_Aa_;  D_dq_A[ID] = D_dq_Aa_;

                  //! Call CRTP visitor.
                  //!------------------------------------------------------------------------------!//
                  visitorAccel02.ddq_ = &ddq_;  visitorAccel02.Sw_ = &Screw_w[ID];  visitorAccel02.A_ = &Accel[ID];  visitorAccel02.D_q_A_ = &D_q_A[ID];
                  visitorAccel02.D_dq_A_ = &D_dq_A[ID];  visitorAccel02.D_q_ddq_ = &D_q_ddq;  visitorAccel02.D_dq_ddq_ = &D_dq_ddq;
                  boost::apply_visitor( visitorAccel02, *JointTypesIter );

                }

              //! Fill up D_ddq.
              //!------------------------------------------------------------------------------!//
              D_ddq.row(ID) << D_q_ddq, D_dq_ddq;

            } else {

              visitorRoot.u = u[ID];         visitorRoot.iD = invD[ID];         visitorRoot.ddq = &ddq[ID];
              visitorRoot.S = &Screw_w[ID];  visitorRoot.P_ = &TransConst[ID];  visitorRoot.R_ = &Rot[ID];
              visitorRoot.U_ = &U[ID];       visitorRoot.Acc_i_ = &Accel[ID];
              //!-------------------------------------------------------
              visitorRoot.D_U_h_ = &D_U_h[ID];  visitorRoot.D_invD_ = &D_invD[ID];
              visitorRoot.D_q_u_ = &D_q_u[ID];  visitorRoot.D_dq_u_ = &D_dq_u[ID];
              visitorRoot.D_q_A_ = &D_q_A[ID];  visitorRoot.D_dq_A_ = &D_dq_A[ID];

              boost::apply_visitor( visitorRoot, *JointTypesIter );

            }
          ID++;
        }

      //! --> Alternative FOREACH
      // BOOST_FOREACH(CV & c, JointTypesVec)
      // {  boost::apply_visitor( visitor, c );  }

      //! --> Alternative std::for_each
      // std::for_each(JointTypesVec.begin(), JointTypesVec.end(), boost::apply_visitor(visitor));

    }

  };

} // namespace geombd

#endif // GEOMBD_FORWARD_DYNAMICS_DIFFERENTIATION_CRTP_HPP
