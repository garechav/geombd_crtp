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
    D_TCP01_visitor<ScalarType, SpatialVector> visitorTCP01;
    D_TCP02_visitor<D_SpatialVector> visitorTCP02;
    D_TCProot_visitor<ScalarType, Vector3r, SpatialVector, SpatialMatrix, D_SpatialVector> visitorTCProot;

    Inertia01_visitor<SpatialVector, SpatialMatrix, VectorXr, D_SpatialMatrix> visitorInertia01;
    Inertia02_visitor<ScalarType, SpatialVector, RowVectorXr, D_SpatialVector> visitorInertia02;


    //    D_TCP_root_visitor<ScalarType, Vector3r, SpatialVector, SpatialMatrix, D_SpatialVector> TCP_rootVis;
    //    D_TwCbPb_visitor<ScalarType, Matrix3r, SpatialMatrix, Vector3r, SpatialVector, D_SpatialVector> TwCbPbVis;

    //    D_InertiaLeaf_visitor<ScalarType, Vector3r, Matrix3r, SpatialVector, SpatialMatrix, RowVectorXr, D_SpatialVector, D_SpatialMatrix> InerLeafVis;
    //    D_Inertial_visitor<ScalarType, Vector6int, Vector3r, Matrix3r, SpatialVector, SpatialMatrix, VectorXr, RowVectorXr, D_SpatialVector, D_SpatialMatrix> InertialVis;

    //    D_Accel_root_visitor<ScalarType, Vector3r, Matrix3r, SpatialVector, D_SpatialVector, RowVectorXr, MatrixXr> Accel_rootVis;
    //    D_Accel_visitor<ScalarType, index_t, Vector3r, Matrix3r, SpatialVector, D_SpatialVector, RowVectorXr, MatrixXr, VectorXint> AccelVis;


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
              visitorTCP01.vi = v(ID);  visitorTCP01.S_ = &Twists[ID];  visitorTCP01.C_ = &Cbias[ID];
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
              visitorTCP02.D_q_c_ = &D_q_c[ID];       visitorTCP02.D_dq_c_ = &D_dq_c[ID];
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
              visitorInertia01.U_ = &U[ID];          visitorInertia01.M_A_ = &M_A[ID];
              visitorInertia01.D_U_v_ = &D_U_v[ID];  visitorInertia01.D_M_A_i_ = &D_M_A_i_;
              boost::apply_visitor( visitorInertia01, *JointTypesIter );

              D_U_h[ID] = Eigen::Map<D_SpatialVector>(D_U_v[ID].derived().data(), 6, nS);

              u[ID] = tau(ID);

              //! Call CRTP visitor.
              //!------------------------------------------------------------------------------!//
              visitorInertia02.invD_ = &invD[ID];  visitorInertia02.u_ = &u[ID];
              visitorInertia02.U_ = &U[ID];  visitorInertia02.P_A_ = &P_A[ID];
              visitorInertia02.D_invD_ = &D_invD[ID];  visitorInertia02.D_q_u_ = &D_q_u[ID];  visitorInertia02.D_dq_u_ = &D_dq_u[ID];
              visitorInertia02.D_U_h_ = &D_U_h[ID];  visitorInertia02.D_q_PA_ = &D_q_PA[ID];  visitorInertia02.D_dq_PA_ = &D_dq_PA[ID];
              boost::apply_visitor( visitorInertia02, *JointTypesIter );


              //! Differentiating inertial back-projection to parent. Only if parent is not the root.           54 seg
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
                  //!------------------------------------------------------------------------------!//           53.9 seg
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
                    }                                                                               ///            62.3 seg

                  kronAux.noalias() = D_M_A_i_*Cbias[ID];

                  //! then D_MtmpII = D_M_a and D_MtmpI is AdjointDual*D_M_a*Adjoint
//!------------------------------------------------------------------------------!//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                  D_M_A_i_.resize(D_M_A_i_.rows()+6, Eigen::NoChange);                     ///      63.0 seg
                  D_Mat6ProjRz(P_zero[ID], nS, TransConst[ID].derived(), Rot[ID].derived(), D_M_A_i_);
//!------------------------------------------------------------------------------!//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                  //! Prepare inertial expression differentiation D_q_Pa.
                  //!------------------------------------------------------------------------------!//              73.8 seg
                  D_SpatialVector D_Vec, D_VecII;
                  D_Vec = Eigen::Map<D_SpatialVector>(kronAux.data(),6,nS);

                  D_Vec.noalias() += (u[ID]*invD[ID])*D_U_h[ID] + (u[ID]*U[ID])*D_invD[ID];   //75.41

                  D_q_Pa[ID].leftCols(nP).noalias() = M_a[ID]*D_q_c[ID];
                  D_q_Pa[ID].rightCols(nS) = D_Vec;           // 79.05

                  D_q_Pa[ID].noalias() += D_q_PA[ID] + UD_*D_q_u[ID];  //  78.46

                  //! Prepare inertial expression differentiation D_dq_Pa.
                  //!------------------------------------------------------------------------------!//
                  D_dq_Pa[ID].noalias() = D_dq_PA[ID] + UD_*D_dq_u[ID];
                  D_dq_Pa[ID].leftCols(nP).noalias() += M_a[ID]*D_dq_c[ID];

                  //!------------------------------------------------------------------------------!//     83.61 seg
                  //!                         INERTIAL BACK PROJECTION                             !//
                  //!------------------------------------------------------------------------------!//
                  //! Back projection of M_a.
                  //!------------------------------------------------------------------------------!//
//!------------------------------------------------------------------------------!//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                  Mat6ProjRz(P_zero[ID], TransConst[ID].derived(), Rot[ID].derived(), M_a[ID], Mtmp);
                  M_A[IDj].noalias() += Mtmp;
//!------------------------------------------------------------------------------!//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                  //! Back projection of P_a.
                  //!------------------------------------------------------------------------------!//    83.99
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

//!------------------------------------------------------------------------------!//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                  //! Back projection of the partial derivative D_M_a.
                  //! M_a differentiation -> Single-motion-revolute screws enable skew-symmetric properties.
                  //!------------------------------------------------------------------------------!//
                  SpatialMatrix M_a_S, M_a_S_;
                  M_a_S.setZero();
                  M_a_S.template row(0) = -M_a[ID].template row(1);  // mimicking effect ad_dual(-Sz)*Ma
                  M_a_S.template row(1) =  M_a[ID].template row(0);
                  M_a_S.template row(3) = -M_a[ID].template row(4);
                  M_a_S.template row(4) =  M_a[ID].template row(3);

                  M_a_S_.noalias() = M_a_S + M_a_S.transpose();   // mimicking effect ad_dual(-Sz)*Ma - Ma*ad(Sz)

                  //! Transform M_a_S_ since it is now symmetric.
                  Mat6ProjRz(P_zero[ID], TransConst[ID].derived(), Rot[ID].derived(), M_a_S_, Mtmp);
//!------------------------------------------------------------------------------!//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                  //!------------------------------------------------------------------------------!//
                  //!                                  BRANCHING                                   !//
                  //!------------------------------------------------------------------------------!//
//!------------------------------------------------------------------------------!//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                  P_A_i << -P_a[ID].coeff(1), P_a[ID].coeff(0), 0, -P_a[ID].coeff(4), P_a[ID].coeff(3), 0;
//!------------------------------------------------------------------------------!//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                  D_M_A_i_.template topRows<6>() = Mtmp;

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








              //              InertialVis.indexVec = &indexVec[ID];
              //              InertialVis.u = &u[ID];   InertialVis.iD = &invD[ID];   InertialVis.tau = tau(ID);
              //              InertialVis.S = &Screw_w[ID];  InertialVis.U_ = &U[ID];  InertialVis.c_ = &Cbias[ID];
              //              InertialVis.P_A_ = &P_A[ID];  InertialVis.M_A_ = &M_A[ID];
              //              InertialVis.P_a_ = &P_a[ID];  InertialVis.M_a_ = &M_a[ID];
              //              InertialVis.P_ = &TransConst[ID];  InertialVis.R_ = &Rot[ID];  InertialVis.P_z = P_zero[ID];
              //              InertialVis.P_Aj_ = &P_A[IDj];  InertialVis.M_Aj_ = &M_A[IDj];
              //              //!-------------------------------------------------------
              //              InertialVis.rootFlag = rootFlag[ID];
              //              InertialVis.D_U_h_ = &D_U_h[ID];    InertialVis.D_U_v_ = &D_U_v[ID];
              //              InertialVis.D_invD_ = &D_invD[ID];  InertialVis.D_M_A_ = &D_M_A[ID];    InertialVis.D_M_Aj_ = &D_M_A[IDj];
              //              //!-------------------------------------------------------
              //              InertialVis.D_q_u_ = &D_q_u[ID];      InertialVis.D_dq_u_ = &D_dq_u[ID];
              //              InertialVis.D_q_Pa_ = &D_q_Pa[ID];    InertialVis.D_dq_Pa_ = &D_dq_Pa[ID];
              //              InertialVis.D_q_PA_ = &D_q_PA[ID];    InertialVis.D_dq_PA_ = &D_dq_PA[ID];
              //              InertialVis.D_q_PAj_ = &D_q_PA[IDj];  InertialVis.D_dq_PAj_ = &D_dq_PA[IDj];
              //              InertialVis.D_q_c_ = &D_q_c[ID];      InertialVis.D_dq_c_ = &D_dq_c[ID];
              //              //!-------------------------------------------------------
              //              boost::apply_visitor( InertialVis, *JointTypesIter );










            } else {
//              InerLeafVis.u = &u[ID];   InerLeafVis.iD = &invD[ID];    InerLeafVis.tau = tau(ID);
//              InerLeafVis.S = &Screw_w[ID];  InerLeafVis.U_ = &U[ID];  InerLeafVis.c_ = &Cbias[ID];
//              InerLeafVis.P_A_ = &P_A[ID];  InerLeafVis.M_A_ = &M_A[ID];
//              InerLeafVis.P_a_ = &P_a[ID];  InerLeafVis.M_a_ = &M_a[ID];
//              InerLeafVis.P_ = &TransConst[ID];  InerLeafVis.R_ = &Rot[ID];  InerLeafVis.P_z = P_zero[ID];
//              InerLeafVis.P_Aj_ = &P_A[IDj];  InerLeafVis.M_Aj_ = &M_A[IDj];
//              //!-------------------------------------------------------
//              InerLeafVis.D_M_Aj_ = &D_M_A[IDj];
//              InerLeafVis.D_q_u_ = &D_q_u[ID];      InerLeafVis.D_dq_u_ = &D_dq_u[ID];
//              InerLeafVis.D_q_p_ = &D_q_p[ID];      InerLeafVis.D_dq_p_ = &D_dq_p[ID];
//              InerLeafVis.D_q_Pa_ = &D_q_Pa[ID];    InerLeafVis.D_dq_Pa_ = &D_dq_Pa[ID];
//              InerLeafVis.D_q_PA_ = &D_q_PA[ID];    InerLeafVis.D_dq_PA_ = &D_dq_PA[ID];
//              InerLeafVis.D_q_PAj_ = &D_q_PA[IDj];  InerLeafVis.D_dq_PAj_ = &D_dq_PA[IDj];
//              InerLeafVis.D_q_c_ = &D_q_c[ID];      InerLeafVis.D_dq_c_ = &D_dq_c[ID];
//              //!-------------------------------------------------------
//              boost::apply_visitor( InerLeafVis, *JointTypesIter );
            }
          ID--;
        }













      //      AccelVis.D_ddq_ = &D_ddq;  Accel_rootVis.D_ddq_ = &D_ddq;

      //      //! Third Recursion
      //      ID = 0;
      //      for ( JointTypesIter = JointTypesVec.begin(); JointTypesIter != JointTypesVec.end(); JointTypesIter++ ){
      //          //! Spatial acceleration
      //          IDj = parent[ID];
      //          if( ID ){
      //              AccelVis.zeroFlag = P_zero[ID];  AccelVis.u = u[ID];  AccelVis.iD = invD[ID];
      //              AccelVis.ddq = &ddq[ID];   AccelVis.S = &Screw_w[ID];  AccelVis.P_ = &TransConst[ID];
      //              AccelVis.R_ = &Rot[ID];    AccelVis.c_ = &Cbias[ID];   AccelVis.U_ = &U[ID];
      //              AccelVis.Acc_i_ = &Accel[ID];  AccelVis.Acc_j_ = &Accel[ IDj ];
      //              //!-------------------------------------------------------
      //              AccelVis.isLeaf = isLeaf[ID];  AccelVis.ID = ID;
      //              AccelVis.Pre_ = &PRE_n[ID];  AccelVis.Suc_ = &SUC_n[ID];  AccelVis.PreSuc_ = &PRESUC_n[ID];
      //              AccelVis.D_U_h_ = &D_U_h[ID];    AccelVis.D_invD_ = &D_invD[ID];
      //              AccelVis.D_q_u_ = &D_q_u[ID];    AccelVis.D_dq_u_ = &D_dq_u[ID];
      //              AccelVis.D_q_c_ = &D_q_c[ID];    AccelVis.D_dq_c_ = &D_dq_c[ID];
      //              AccelVis.D_q_A_ = &D_q_A[ID];    AccelVis.D_dq_A_ = &D_dq_A[ID];
      //              AccelVis.D_q_Aj_ = &D_q_A[IDj];  AccelVis.D_dq_Aj_ = &D_dq_A[IDj];

      //              boost::apply_visitor( AccelVis, *JointTypesIter );
      //            } else {
      //              Accel_rootVis.u = u[ID];  Accel_rootVis.iD = invD[ID];  Accel_rootVis.ddq = &ddq[ID];
      //              Accel_rootVis.S = &Screw_w[ID];  Accel_rootVis.P_ = &TransConst[ID];  Accel_rootVis.R_ = &Rot[ID];
      //              Accel_rootVis.U_ = &U[ID];  Accel_rootVis.Acc_i_ = &Accel[ID];
      //              //!-------------------------------------------------------
      //              Accel_rootVis.D_U_h_ = &D_U_h[ID];  Accel_rootVis.D_invD_ = &D_invD[ID];
      //              Accel_rootVis.D_q_u_ = &D_q_u[ID];  Accel_rootVis.D_dq_u_ = &D_dq_u[ID];
      //              Accel_rootVis.D_q_A_ = &D_q_A[ID];  Accel_rootVis.D_dq_A_ = &D_dq_A[ID];

      //              boost::apply_visitor( Accel_rootVis, *JointTypesIter );
      //            }
      //          ID++;
      //        }

      //! --> Alternative FOREACH
      // BOOST_FOREACH(CV & c, JointTypesVec)
      // {  boost::apply_visitor( visitor, c );  }

      //! --> Alternative std::for_each
      // std::for_each(JointTypesVec.begin(), JointTypesVec.end(), boost::apply_visitor(visitor));

    }

  };

} // namespace geombd

#endif // GEOMBD_FORWARD_DYNAMICS_DIFFERENTIATION_CRTP_HPP
