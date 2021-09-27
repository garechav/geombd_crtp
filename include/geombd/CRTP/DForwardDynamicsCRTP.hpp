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

#include "geombd/CRTP/DJointBase.hxx"
#include "geombd/pinocchio/container/aligned-vector.hpp"
#include "geombd/io/parser.hpp"

namespace geo {

  const int maxBody = 32;

  template<typename ScalarType>
  class FwdDynDifCRTP
  {
    // --------------------------------------------
    // Types
    // --------------------------------------------

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

    //!-------------------------------------------------------

    typedef Eigen::Matrix<real_t, 6, Eigen::Dynamic, 0, 6, maxBody> D_SpatialVector;
    typedef Eigen::Matrix<real_t, Eigen::Dynamic, 6, 0, 6*maxBody, 6> D_SpatialMatrix;

    typedef Eigen::Matrix<real_t, 1, Eigen::Dynamic, Eigen::RowMajor, 1, maxBody> RowVectorXr;

    typedef Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, 0, maxBody, 2*maxBody> MatrixXr;

    // --------------------------------------------
    // Constructors and Destructors
    // --------------------------------------------

  public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW //128-bit alignment


    //! Custom constructor
    FwdDynDifCRTP( const std::shared_ptr< Robot >& robot )
    {

      EIGEN_MAKE_ALIGNED_OPERATOR_NEW

      this->n = robot->nq;

      //! Pre-allocation for enhancement
      //!-------------------------------------------------------
      ddq = VectorXr::Zero(n);
      PRE.resize(n);  SUC.resize(n);  PRESUC.resize(n);
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
          //!-------------------------------------------------------
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
          //!-------------------------------------------------------
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
          //!-------------------------------------------------------
          Vector3r axis = joint->axis;
          Screw_w.push_back(axis);

          if (joint->type == Joint::joint_type::prismatic) {
              D_dq_V[i-1].template rightCols<1>().template segment<3>(0) = axis;
              //!-------------------------------------------------------
              //! code here for Joint Types Px Py and Pz
            }
          if (joint->type == Joint::joint_type::revolute) {
              D_dq_V[i-1].template rightCols<1>().template segment<3>(3) = axis;
//              //!-------------------------------------------------------
//              if((real_t)axis(0) == 1) {
//                  JointTypesVec.push_back( D_JointTypeRx{} );
//                }
//              //!-------------------------------------------------------
//              if((real_t)axis(1) == 1) {
//                  JointTypesVec.push_back( D_JointTypeRy{} );
//                }
//              //!-------------------------------------------------------
//              if((real_t)axis(2) == 1) {
//                  JointTypesVec.push_back( D_JointTypeRz{} );
//                }
//              //!-------------------------------------------------------
//              if((real_t)axis(1) != 0 && (real_t)axis(2) != 0) {
                  JointTypesVec.push_back( D_JointTypeRxyz{} );
//                }

            }

        }

    }


    ~FwdDynDifCRTP(){}

    // --------------------------------------------
    // Members
    // --------------------------------------------

  public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Dof
    index_t n;

    index_t ID, IDj, j;

    SpatialVector Screw;

    std::vector< bool > P_zero, rootFlag, branchFlag, isLeaf;

    std::vector< index_t > parent;

    //! Joint-types variant
    typedef boost::variant<D_JointTypeRx, D_JointTypeRy, D_JointTypeRz, D_JointTypeRxyz> CV;

    //! Joint-types vector
    std::vector< CV > JointTypesVec;
    typename std::vector< CV >::iterator JointTypesIter;

    D_FwdKin_visitor<ScalarType, Vector3r, Matrix3r> visitorFK;
    D_TCP_root_visitor<ScalarType, Vector3r, SpatialVector, SpatialMatrix, D_SpatialVector> TCP_rootVis;
    D_TwCbPb_visitor<ScalarType, Matrix3r, SpatialMatrix, Vector3r, SpatialVector, D_SpatialVector> TwCbPbVis;

    D_InertiaLeaf_visitor<ScalarType, Vector3r, Matrix3r, SpatialVector, SpatialMatrix, RowVectorXr, D_SpatialVector, D_SpatialMatrix> InerLeafVis;
    D_Inertial_visitor<ScalarType, Vector6int, Vector3r, Matrix3r, SpatialVector, SpatialMatrix, VectorXr, RowVectorXr, D_SpatialVector, D_SpatialMatrix> InertialVis;

    D_Accel_root_visitor<ScalarType, Vector3r, Matrix3r, SpatialVector, D_SpatialVector, RowVectorXr, MatrixXr> Accel_rootVis;
    D_Accel_visitor<ScalarType, index_t, Vector3r, Matrix3r, SpatialVector, D_SpatialVector, RowVectorXr, MatrixXr> AccelVis;


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

    //!-------------------------------------------------------

    index_t nP, nS, PS, nPj, nSj;
    VectorXint NP, NS;

    std::vector< index_t > Pre, Suc, PreSuc;
    std::vector< std::vector< index_t > > PRE, SUC, PRESUC;

    SpatialMatrix D_p_aux;

    PINOCCHIO_ALIGNED_STD_VECTOR( D_SpatialVector ) D_q_V;   PINOCCHIO_ALIGNED_STD_VECTOR( D_SpatialVector ) D_dq_V;
    PINOCCHIO_ALIGNED_STD_VECTOR( D_SpatialVector ) D_q_c;   PINOCCHIO_ALIGNED_STD_VECTOR( D_SpatialVector ) D_dq_c;

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

    // --------------------------------------------
    // Methods
    // --------------------------------------------

  public:


    //! Articulated-Body Algorithm
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
          //! Forward kinematics
          //!-------------------------------------------------------
          visitorFK.qi = q(ID);  visitorFK.S = &Screw_w[ID];  visitorFK.R = &Rot[ID];
          boost::apply_visitor( visitorFK, *JointTypesIter );

          IDj = parent[ID];
          //! Twist and c, p bias
          //!-------------------------------------------------------
          if(ID > 0){
              TwCbPbVis.zeroFlag = P_zero[ID];  TwCbPbVis.vi = v(ID);  TwCbPbVis.S = &Screw_w[ID];  TwCbPbVis.R_ = &Rot[ID];
              TwCbPbVis.P_ = &TransConst[ID];   TwCbPbVis.S_j = &Twists[IDj];     TwCbPbVis.M_ = &MConst[ID];
              TwCbPbVis.S_i = &Twists[ID];      TwCbPbVis.c_ = &Cbias[ID];  TwCbPbVis.p_ = &Pbias[ID];
              //!-------------------------------------------------------
              TwCbPbVis.D_q_Vj_ = &D_q_V[IDj];  TwCbPbVis.D_dq_Vj_ = &D_dq_V[IDj];  TwCbPbVis.D_q_V_ = &D_q_V[ID];
              TwCbPbVis.D_dq_V_ = &D_dq_V[ID];  TwCbPbVis.D_q_c_ = &D_q_c[ID];      TwCbPbVis.D_dq_c_ = &D_dq_c[ID];
              TwCbPbVis.D_q_p_ = &D_q_p[ID];    TwCbPbVis.D_dq_p_ = &D_dq_p[ID];
              //!-------------------------------------------------------
              boost::apply_visitor( TwCbPbVis, *JointTypesIter );
            } else {
              TCP_rootVis.vi = v(ID);  TCP_rootVis.S = &Screw_w[ID];  TCP_rootVis.S_i = &Twists[ID];
              TCP_rootVis.p_ = &Pbias[ID];  TCP_rootVis.M_ = &MConst[ID];  TCP_rootVis.D_dq_p_ = &D_dq_p[ID];
              //!-------------------------------------------------------
              boost::apply_visitor( TCP_rootVis, *JointTypesIter );
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
              InertialVis.indexVec = &indexVec[ID];
              InertialVis.u = &u[ID];   InertialVis.iD = &invD[ID];   InertialVis.tau = tau(ID);
              InertialVis.S = &Screw_w[ID];  InertialVis.U_ = &U[ID];  InertialVis.c_ = &Cbias[ID];
              InertialVis.P_A_ = &P_A[ID];  InertialVis.M_A_ = &M_A[ID];
              InertialVis.P_a_ = &P_a[ID];  InertialVis.M_a_ = &M_a[ID];
              InertialVis.P_ = &TransConst[ID];  InertialVis.R_ = &Rot[ID];  InertialVis.P_z = P_zero[ID];
              InertialVis.P_Aj_ = &P_A[IDj];  InertialVis.M_Aj_ = &M_A[IDj];
              //!-------------------------------------------------------
              InertialVis.rootFlag = rootFlag[ID];
              InertialVis.D_U_h_ = &D_U_h[ID];    InertialVis.D_U_v_ = &D_U_v[ID];
              InertialVis.D_invD_ = &D_invD[ID];  InertialVis.D_M_A_ = &D_M_A[ID];    InertialVis.D_M_Aj_ = &D_M_A[IDj];
              //!-------------------------------------------------------
              InertialVis.D_q_u_ = &D_q_u[ID];      InertialVis.D_dq_u_ = &D_dq_u[ID];
              InertialVis.D_q_Pa_ = &D_q_Pa[ID];    InertialVis.D_dq_Pa_ = &D_dq_Pa[ID];
              InertialVis.D_q_PA_ = &D_q_PA[ID];    InertialVis.D_dq_PA_ = &D_dq_PA[ID];
              InertialVis.D_q_PAj_ = &D_q_PA[IDj];  InertialVis.D_dq_PAj_ = &D_dq_PA[IDj];
              InertialVis.D_q_c_ = &D_q_c[ID];      InertialVis.D_dq_c_ = &D_dq_c[ID];
              //!-------------------------------------------------------
              boost::apply_visitor( InertialVis, *JointTypesIter );
            } else {
              InerLeafVis.u = &u[ID];   InerLeafVis.iD = &invD[ID];    InerLeafVis.tau = tau(ID);
              InerLeafVis.S = &Screw_w[ID];  InerLeafVis.U_ = &U[ID];  InerLeafVis.c_ = &Cbias[ID];
              InerLeafVis.P_A_ = &P_A[ID];  InerLeafVis.M_A_ = &M_A[ID];
              InerLeafVis.P_a_ = &P_a[ID];  InerLeafVis.M_a_ = &M_a[ID];
              InerLeafVis.P_ = &TransConst[ID];  InerLeafVis.R_ = &Rot[ID];  InerLeafVis.P_z = P_zero[ID];
              InerLeafVis.P_Aj_ = &P_A[IDj];  InerLeafVis.M_Aj_ = &M_A[IDj];
              //!-------------------------------------------------------
              InerLeafVis.D_M_Aj_ = &D_M_A[IDj];
              InerLeafVis.D_q_u_ = &D_q_u[ID];      InerLeafVis.D_dq_u_ = &D_dq_u[ID];
              InerLeafVis.D_q_p_ = &D_q_p[ID];      InerLeafVis.D_dq_p_ = &D_dq_p[ID];
              InerLeafVis.D_q_Pa_ = &D_q_Pa[ID];    InerLeafVis.D_dq_Pa_ = &D_dq_Pa[ID];
              InerLeafVis.D_q_PA_ = &D_q_PA[ID];    InerLeafVis.D_dq_PA_ = &D_dq_PA[ID];
              InerLeafVis.D_q_PAj_ = &D_q_PA[IDj];  InerLeafVis.D_dq_PAj_ = &D_dq_PA[IDj];
              InerLeafVis.D_q_c_ = &D_q_c[ID];      InerLeafVis.D_dq_c_ = &D_dq_c[ID];
              //!-------------------------------------------------------
              boost::apply_visitor( InerLeafVis, *JointTypesIter );
            }
          ID--;
        }

      AccelVis.D_ddq_ = &D_ddq;  Accel_rootVis.D_ddq_ = &D_ddq;

      //! Third Recursion
      ID = 0;
      for ( JointTypesIter = JointTypesVec.begin(); JointTypesIter != JointTypesVec.end(); JointTypesIter++ ){
          //! Spatial acceleration
          IDj = parent[ID];
          if( ID ){
              AccelVis.zeroFlag = P_zero[ID];  AccelVis.u = u[ID];  AccelVis.iD = invD[ID];
              AccelVis.ddq = &ddq[ID];   AccelVis.S = &Screw_w[ID];  AccelVis.P_ = &TransConst[ID];
              AccelVis.R_ = &Rot[ID];    AccelVis.c_ = &Cbias[ID];   AccelVis.U_ = &U[ID];
              AccelVis.Acc_i_ = &Accel[ID];  AccelVis.Acc_j_ = &Accel[ IDj ];
              //!-------------------------------------------------------
              AccelVis.isLeaf = isLeaf[ID];  AccelVis.ID = ID;
              AccelVis.Pre_ = &PRE[ID];  AccelVis.Suc_ = &SUC[ID];  AccelVis.PreSuc_ = &PRESUC[ID];
              AccelVis.D_U_h_ = &D_U_h[ID];    AccelVis.D_invD_ = &D_invD[ID];
              AccelVis.D_q_u_ = &D_q_u[ID];    AccelVis.D_dq_u_ = &D_dq_u[ID];
              AccelVis.D_q_c_ = &D_q_c[ID];    AccelVis.D_dq_c_ = &D_dq_c[ID];
              AccelVis.D_q_A_ = &D_q_A[ID];    AccelVis.D_dq_A_ = &D_dq_A[ID];
              AccelVis.D_q_Aj_ = &D_q_A[IDj];  AccelVis.D_dq_Aj_ = &D_dq_A[IDj];

              boost::apply_visitor( AccelVis, *JointTypesIter );
            } else {
              Accel_rootVis.u = u[ID];  Accel_rootVis.iD = invD[ID];  Accel_rootVis.ddq = &ddq[ID];
              Accel_rootVis.S = &Screw_w[ID];  Accel_rootVis.P_ = &TransConst[ID];  Accel_rootVis.R_ = &Rot[ID];
              Accel_rootVis.U_ = &U[ID];  Accel_rootVis.Acc_i_ = &Accel[ID];
              //!-------------------------------------------------------
              Accel_rootVis.D_U_h_ = &D_U_h[ID];  Accel_rootVis.D_invD_ = &D_invD[ID];
              Accel_rootVis.D_q_u_ = &D_q_u[ID];  Accel_rootVis.D_dq_u_ = &D_dq_u[ID];
              Accel_rootVis.D_q_A_ = &D_q_A[ID];  Accel_rootVis.D_dq_A_ = &D_dq_A[ID];

              boost::apply_visitor( Accel_rootVis, *JointTypesIter );
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
