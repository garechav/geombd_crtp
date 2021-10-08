/**
 *    \file include/geombd/CRTP/ForwardDynamics.hpp
 *    \author Alvaro Paz, Gustavo Arechavaleta
 *    \version 1.0
 *    \date 2021
 *
 *    Class to implement the Articulated-Body Algorithm
 *    Copyright (c) 2021 Cinvestav
 *    This library is distributed under the MIT License.
 */

#ifndef GEOMBD_FORWARD_DYNAMICS_CRTP_HPP
#define GEOMBD_FORWARD_DYNAMICS_CRTP_HPP

#include "geombd/CRTP/JointBase.hxx"

#include "geombd/pinocchio/container/aligned-vector.hpp"

#define EIGEN_NO_DEBUG
#include "geombd/io/parser.hpp"


namespace geoCRTP
{
  const int maxBody = 32;

  template<typename ScalarType>
  class FwdDynCRTP
  {
    // --------------------------------------------
    // Types
    // --------------------------------------------

    typedef int index_t;
    typedef ScalarType real_t;

    typedef Eigen::Matrix<real_t, 6, 1> SpatialVector;
    typedef Eigen::Matrix<real_t, 6, 6> SpatialMatrix;

    typedef Eigen::Matrix<real_t, Eigen::Dynamic, 1, 0, maxBody, 1> VectorXr;

    typedef Eigen::Matrix<real_t , 2 , 1> Vector2r;
    typedef Eigen::Matrix<real_t , 2 , 2> Matrix2r;

    typedef Eigen::Matrix<real_t , 3 , 1> Vector3r;
    typedef Eigen::Matrix<real_t , 3 , 3> Matrix3r;

    typedef Eigen::Matrix<real_t , 4 , 1> Vector4r;
    typedef Eigen::Matrix<real_t , 4 , 4> Matrix4r;

    // --------------------------------------------

    typedef Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, 0, maxBody, maxBody> MatrixXr;
    typedef Eigen::Matrix<real_t, 1, Eigen::Dynamic, Eigen::RowMajor, 1, maxBody> RowVectorXr;
    typedef Eigen::Matrix<real_t, 6, Eigen::Dynamic, 0, 6, 2*maxBody> D_SpatialVector;

    // --------------------------------------------
    // Constructors and Destructors
    // --------------------------------------------

  public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW //128-bit alignment

    //! Custom constructor
    FwdDynCRTP( const std::shared_ptr< Robot >& robot )
    {

      this->n = robot->nq;

      this->ddq = VectorXr::Zero(n);
      this->inv_H = MatrixXr::Zero(n,n);
      this->Fcrb = D_SpatialVector::Zero(6,n);
      this->Pc = D_SpatialVector::Zero(6,n);

      int ID = 0;

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

          D_Fcrb.push_back(D_SpatialVector::Zero(6,2*n));
          D_Pc.push_back(D_SpatialVector::Zero(6,n));

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
          NP.push_back( body->Pre.size() );
          NS.push_back( body->Suc.size() );

          //! Handling screw axes
          //!-------------------------------------------------------
          Vector3r axis = joint->axis;          
          Screw_w.push_back(axis);

          if (joint->type == Joint::joint_type::prismatic) {
              //!-------------------------------------------------------
              //! code here for Joint Types Px Py and Pz
            }
          if (joint->type == Joint::joint_type::revolute) {
              //!-------------------------------------------------------
              if((real_t)axis(0) == 1) {
                  JointTypesVec.push_back( JointTypeRx{} );
                }
              //!-------------------------------------------------------
              if((real_t)axis(1) == 1) {
                  JointTypesVec.push_back( JointTypeRy{} );
                }
              //!-------------------------------------------------------
              if((real_t)axis(2) == 1) {
                  JointTypesVec.push_back( JointTypeRz{} );
                }
              //!-------------------------------------------------------
              if((real_t)axis(1) != 0 && (real_t)axis(2) != 0) {
                  JointTypesVec.push_back( JointTypeRxyz{} );
                }

            }

          ID++;
        }

    }


    ~FwdDynCRTP(){}

    // --------------------------------------------
    // Members
    // --------------------------------------------

  public:

    //! Dof
    index_t n;

    index_t ID, j;

    std::vector< bool > P_zero;

    std::vector< index_t > parent;

    //! Joint-types variant
    typedef boost::variant<JointTypeRx, JointTypeRy, JointTypeRz, JointTypeRxyz> CV;

    //! Joint-types vector
    std::vector< CV > JointTypesVec;
    typename std::vector< CV >::iterator JointTypesIter;

    //! Articulated-Body Algorithm Visitors
    //!------------------------------------------------------------------------------!//
    FwdKin_visitor<ScalarType, Vector3r, Matrix3r> visitorFK;
    TCP_root_visitor<ScalarType, Vector3r, SpatialVector, SpatialMatrix> TCP_rootVis;
    TwCbPb_visitor<ScalarType, Matrix3r, SpatialMatrix, Vector3r, SpatialVector> TwCbPbVis;

    UuiD_visitor<ScalarType, Vector3r, SpatialVector, SpatialMatrix> UuiDVis;
    PreIner_visitor<ScalarType, Vector3r, SpatialVector, SpatialMatrix> PreInerVis;
    InerProj_visitor<Vector3r, Matrix3r, SpatialVector, SpatialMatrix> InerProjVis;

    Accel_visitor<ScalarType, Vector3r, Matrix3r, SpatialVector> AccelVis;
    Accel_root_visitor<ScalarType, Vector3r, Matrix3r, SpatialVector> Accel_rootVis;

    //! Inverse Inertia Matrix Visitors
    //!------------------------------------------------------------------------------!//
    UinvD_visitor<ScalarType, Vector3r, SpatialVector, SpatialMatrix> UinvDVis;
    invInertia_visitor<index_t, ScalarType, Vector3r, D_SpatialVector, RowVectorXr> invInertiaVis;
    HSelector_visitor<index_t, Vector3r, D_SpatialVector, MatrixXr> HSelectorVis;

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

    PINOCCHIO_ALIGNED_STD_VECTOR( D_SpatialVector ) D_Fcrb;
    PINOCCHIO_ALIGNED_STD_VECTOR( D_SpatialVector ) D_Pc;

    VectorXr ddq;
    MatrixXr inv_H, inv_H_t;
    RowVectorXr iHrow;
    D_SpatialVector Fcrb, Pc;
    std::vector<int> NP, NS;
    SpatialVector Faux, Faux_;
    SpatialMatrix M_aux, M_aux_out;
    SpatialMatrix Adjoint;

    // --------------------------------------------
    // Methods
    // --------------------------------------------

  public:


    //! Articulated-Body Algorithm
    //!------------------------------------------------------------------------------!//
    template<typename ConfigVectorT, typename TangentVectorT, typename CotangentVectorT>
    const void
    aba(const Eigen::MatrixBase<ConfigVectorT> & q,
        const Eigen::MatrixBase<TangentVectorT> & v,
        const Eigen::MatrixBase<CotangentVectorT> & tau)
    {
      typedef typename ConfigVectorT::Scalar ScalarT;

      //! First Recursion
      //!------------------------------------------------------------------------------!//
      ID = 0;
      for ( JointTypesIter = JointTypesVec.begin(); JointTypesIter != JointTypesVec.end(); JointTypesIter++ ){
          //! Forward kinematics
          visitorFK.qi = q(ID);  visitorFK.S = &Screw_w[ID];  visitorFK.R = &Rot[ID];
          boost::apply_visitor( visitorFK, *JointTypesIter );

          //! Twist and c, p bias
          if(ID != 0){
              TwCbPbVis.zeroFlag = P_zero[ID];  TwCbPbVis.vi = v(ID);  TwCbPbVis.S = &Screw_w[ID];  TwCbPbVis.R_ = &Rot[ID];  TwCbPbVis.P_ = &TransConst[ID];
              TwCbPbVis.S_j = &Twists[ parent[ID] ];  TwCbPbVis.M_ = &MConst[ID];
              TwCbPbVis.S_i = &Twists[ID];  TwCbPbVis.c_ = &Cbias[ID];  TwCbPbVis.p_ = &Pbias[ID];

              boost::apply_visitor( TwCbPbVis, *JointTypesIter );
            } else {
              TCP_rootVis.vi = v(ID);  TCP_rootVis.S = &Screw_w[ID];  TCP_rootVis.S_i = &Twists[ID];
              TCP_rootVis.p_ = &Pbias[ID];  TCP_rootVis.M_ = &MConst[ID];

              boost::apply_visitor( TCP_rootVis, *JointTypesIter );
            }
          ID++;
        }

      //! Initialize inertial terms
      M_A = MConst;  P_A = Pbias;

      //! Second Recursion
      //!------------------------------------------------------------------------------!//
      ID = n - 1;
      for ( JointTypesIter = JointTypesVec.end()-1; JointTypesIter != JointTypesVec.begin()-1; JointTypesIter-- ) {
          //! U, u & invD
          UuiDVis.u = &u[ID];   UuiDVis.iD = &invD[ID];   UuiDVis.tau = tau(ID);   UuiDVis.S = &Screw_w[ID];
          UuiDVis.U_ = &U[ID];  UuiDVis.P_A_ = &P_A[ID];  UuiDVis.M_A_ = &M_A[ID];
          boost::apply_visitor( UuiDVis, *JointTypesIter );

          //! Inertial back projection
          if(ID != 0){
              //! Prepare inertias
              PreInerVis.u = u[ID];        PreInerVis.iD = invD[ID];    PreInerVis.S = &Screw_w[ID];
              PreInerVis.U_ = &U[ID];      PreInerVis.c_ = &Cbias[ID];
              PreInerVis.P_a_ = &P_a[ID];  PreInerVis.M_a_ = &M_a[ID];
              PreInerVis.P_A_ = &P_A[ID];  PreInerVis.M_A_ = &M_A[ID];

              boost::apply_visitor( PreInerVis, *JointTypesIter );

              //! Inertial back-projection
              InerProjVis.S = &Screw_w[ID];  InerProjVis.P_ = &TransConst[ID];
              InerProjVis.R_ = &Rot[ID];         InerProjVis.P_z = P_zero[ID];
              InerProjVis.M_a_ = &M_a[ID];       InerProjVis.M_A_ = &M_A[ parent[ID] ];
              InerProjVis.P_a_ = &P_a[ID];       InerProjVis.P_A_ = &P_A[ parent[ID] ];

              boost::apply_visitor( InerProjVis, *JointTypesIter );
            }
          ID--;
        }

      //! Third Recursion
      //!------------------------------------------------------------------------------!//
      ID = 0;
      for ( JointTypesIter = JointTypesVec.begin(); JointTypesIter != JointTypesVec.end(); JointTypesIter++ ){
          //! Spatial acceleration
          if( ID ){
              AccelVis.zeroFlag = P_zero[ID];  AccelVis.u = u[ID];  AccelVis.iD = invD[ID];
              AccelVis.ddq = &ddq[ID];  AccelVis.S_ = &Screw_w[ID];  AccelVis.P_ = &TransConst[ID];
              AccelVis.R_ = &Rot[ID];  AccelVis.c_ = &Cbias[ID];  AccelVis.U_ = &U[ID];
              AccelVis.Acc_i_ = &Accel[ID];  AccelVis.Acc_j_ = &Accel[ parent[ID] ];

              boost::apply_visitor( AccelVis, *JointTypesIter );
            } else {
              Accel_rootVis.u = u[ID];  Accel_rootVis.iD = invD[ID];  Accel_rootVis.ddq = &ddq[ID];
              Accel_rootVis.S_ = &Screw_w[ID];  Accel_rootVis.P_ = &TransConst[ID];  Accel_rootVis.R_ = &Rot[ID];
              Accel_rootVis.U_ = &U[ID];  Accel_rootVis.Acc_i_ = &Accel[ID];

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


    //! Inverse of Generalized Inertia Matrix
    //!------------------------------------------------------------------------------!//
    template<typename ConfigVectorT>
    const void
    inverseInertiaMatrix( const Eigen::MatrixBase<ConfigVectorT> & q )
    {
      //! Regarding sub_tree = n for highly-DoF's robots

      typedef typename ConfigVectorT::Scalar ScalarT;

      typedef const Eigen::Block<SpatialMatrix,3,3> constBlock3;
      typedef Eigen::Block<SpatialMatrix,3,3> Block3;

      constBlock3 & Ai = M_aux.template block<3,3>(0,0);
      constBlock3 & Bi = M_aux.template block<3,3>(0,3);
      constBlock3 & Di = M_aux.template block<3,3>(3,3);

      Block3 Ao = M_aux_out.template block<3,3>(0,0);
      Block3 Bo = M_aux_out.template block<3,3>(0,3);
      Block3 Co = M_aux_out.template block<3,3>(3,0);
      Block3 Do = M_aux_out.template block<3,3>(3,3);

      //! First Recursion
      //!------------------------------------------------------------------------------!//
      ID = 0;
      for ( JointTypesIter = JointTypesVec.begin(); JointTypesIter != JointTypesVec.end(); JointTypesIter++ ){
          //! Forward kinematics
          visitorFK.qi = q(ID);  visitorFK.S = &Screw_w[ID];  visitorFK.R = &Rot[ID];
          boost::apply_visitor( visitorFK, *JointTypesIter );

          D_Fcrb[ID].setZero();  //! is it really neccesary?
          D_Pc[ID].setZero();  //! is it really neccesary?

          ID++;
        }

      //! Initialize inertial terms
      inv_H.setZero();  //! is it really neccesary?
      M_A = MConst;

      //! Second Recursion
      //!------------------------------------------------------------------------------!//
      ID = n - 1;
      for ( JointTypesIter = JointTypesVec.end()-1; JointTypesIter != JointTypesVec.begin()-1; JointTypesIter-- ) {
          //! U & invD
          UinvDVis.iD = &invD[ID];  UinvDVis.S = &Screw_w[ID];  UinvDVis.U_ = &U[ID];  UinvDVis.M_A_ = &M_A[ID];
          boost::apply_visitor( UinvDVis, *JointTypesIter );

          inv_H(ID, ID) = invD[ID];

          //! If the current body is not a leaf
          //!------------------------------------------------------------------------------!//
          if (NS[ID]-1) {
              Fcrb = D_Fcrb[ID].leftCols(n);

              invInertiaVis.n_ID = n-ID;   invInertiaVis.iD = &invD[ID];  invInertiaVis.S = &Screw_w[ID];
              invInertiaVis.Fcrb_ = &Fcrb; invInertiaVis.iHrow_ = &iHrow;
              boost::apply_visitor( invInertiaVis, *JointTypesIter );

              inv_H.block(ID,ID,1,n-ID) -= iHrow;

              if (NP[ID]!=1) {
                  Fcrb.rightCols(n-ID) += U[ID]*inv_H.block(ID,ID,1,n-ID);

                  Pc.topRightCorner(3,n-ID) = Rot[ID]*Fcrb.topRightCorner(3,n-ID);  //!  tmp
                  Pc.bottomRightCorner(3,n-ID) = Rot[ID]*Fcrb.bottomRightCorner(3,n-ID);  //!  tmp

                  if(P_zero[ID]) Pc.bottomRightCorner(3,n-ID) += Skew(TransConst[ID])*Pc.topRightCorner(3,n-ID);  // If P != 0

                  Fcrb.rightCols(n-ID) = D_Fcrb[ parent[ID] ].middleCols(ID,n-ID) + Pc.rightCols(n-ID);

                  D_Fcrb[ parent[ID] ].leftCols(n) = Fcrb;
                }

            } else {
              Fcrb = D_Fcrb[ parent[ID] ].leftCols(n);

              Faux_ = U[ID]*inv_H(ID,ID);

              Faux.template segment<3>(0) = Rot[ID]*Faux_.template segment<3>(0);
              Faux.template segment<3>(3) = Rot[ID]*Faux_.template segment<3>(3);

              if(P_zero[ID]) Faux.template segment<3>(3) += TransConst[ID].cross( Faux.template segment<3>(0) );  // If P != 0

              Fcrb.col(ID) += Faux;

              D_Fcrb[ parent[ID] ].leftCols(n) = Fcrb;
            }
          inv_H(ID, ID) = invD[ID];

          if (NP[ID]!=1) {
              M_aux.noalias() = M_A[ID] - invD[ID]*U[ID]*U[ID].transpose();

              Do.noalias() = Rot[ID]*Ai; // tmp variable
              Ao.noalias() = Do*Rot[ID].transpose();

              Bo.noalias() = Rot[ID]*Di; // tmp variable
              Do.noalias() = Bo*Rot[ID].transpose();

              Bo.noalias() = Rot[ID]*Bi; // tmp variable
              Co.noalias() = Bo*Rot[ID].transpose();

              Bo = Co;
              if(P_zero[ID]){
                  //! Linear
                  typename Matrix3r::PlainObject SkP = Skew(TransConst[ID]);
                  Bo.noalias() -= Ao*SkP;

                  typename Matrix3r::PlainObject Dtmp1;
                  Dtmp1.noalias() = SkP*Bo;
                  Do.noalias() += Dtmp1;
                  Do.noalias() -= Co.transpose()*SkP;
                }
              Co = Bo.transpose();

              M_A[ parent[ID] ] += M_aux_out;
            }

//          std::cout<<"Fcrb "<<ID<<" = "<<D_Fcrb[ parent[ID] ].cwiseAbs().sum()<<std::endl;

          ID--;
        }


      //! Third Recursion
      //!------------------------------------------------------------------------------!//
      ID = 0;
      for ( JointTypesIter = JointTypesVec.begin(); JointTypesIter != JointTypesVec.end(); JointTypesIter++ ){
          typename Vector3r::PlainObject & P_ = const_cast<Vector3r &>(TransConst[ID].derived());
          typename Matrix3r::PlainObject & R_ = const_cast<Matrix3r &>(Rot[ID].derived());

          ::geo::AdDual(P_, R_, Adjoint);

          if (NP[ID]!=1) {
              Faux = invD[ID]*Adjoint*U[ID];  //! tmp (iD*U^{T}*A^{T})
              inv_H.block(ID,ID,1,n-ID) -= Faux.transpose()*D_Fcrb[ parent[ID] ].middleCols(ID,n-ID);
            }

          HSelectorVis.n = n;      HSelectorVis.ID = ID;  HSelectorVis.S_ = &Screw_w[ID];
          HSelectorVis.Pc_ = &Pc;  HSelectorVis.iH_total_ = &inv_H;
          boost::apply_visitor( HSelectorVis, *JointTypesIter );

          if (NP[ID]!=1) {
              Pc.rightCols(n-ID) += Adjoint.transpose()*D_Fcrb[ parent[ID] ].middleCols(ID,n-ID);
            }

          D_Fcrb[ID].leftCols(n) = Pc;

          ID++;
        }

      //! Forcing symmetry
      for ( short int ID = 0; ID < n; ID++ ) inv_H(ID,ID) /= 2;

      inv_H_t = inv_H.transpose();
      inv_H += inv_H_t;

    }

  };

} // namespace geombdCRTP

#endif // GEOMBD_FORWARD_DYNAMICS_CRTP_HPP
