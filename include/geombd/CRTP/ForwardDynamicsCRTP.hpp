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

#define EIGEN_NO_DEBUG
#include "geombd/io/parser.hpp"


namespace geoCRTP
{

  template<typename ScalarType>
  class FwdDynCRTP
  {
    // --------------------------------------------
    // Types
    // --------------------------------------------

    typedef short int index_t;
    typedef ScalarType real_t;

    typedef Eigen::Matrix<real_t, 6, 1> SpatialVector;
    typedef Eigen::Matrix<real_t, 6, 6> SpatialMatrix;

    typedef  Eigen::Matrix<real_t, Eigen::Dynamic, 1, 0, 30, 1> VectorXr;

    typedef Eigen::Matrix<real_t , 2 , 1> Vector2r;
    typedef Eigen::Matrix<real_t , 2 , 2> Matrix2r;

    typedef Eigen::Matrix<real_t , 3 , 1> Vector3r;
    typedef Eigen::Matrix<real_t , 3 , 3> Matrix3r;

    typedef Eigen::Matrix<real_t , 4 , 1> Vector4r;
    typedef Eigen::Matrix<real_t , 4 , 4> Matrix4r;

    // --------------------------------------------
    // Constructors and Destructors
    // --------------------------------------------

  public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW //128-bit alignment

    //! Custom constructor
    FwdDynCRTP( const std::shared_ptr< Robot >& robot )
    {

      this->n = robot->nq;

      ddq = VectorXr::Zero(n);

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

    FwdKin_visitor<ScalarType, Vector3r, Matrix3r> visitorFK;
    TCP_root_visitor<ScalarType, Vector3r, SpatialVector, SpatialMatrix> TCP_rootVis;
    TwCbPb_visitor<ScalarType, Matrix3r, SpatialMatrix, Vector3r, SpatialVector> TwCbPbVis;

    UuiD_visitor<ScalarType, Vector3r, SpatialVector, SpatialMatrix> UuiDVis;
    PreIner_visitor<ScalarType, Vector3r, SpatialVector, SpatialMatrix> PreInerVis;
    InerProj_visitor<Vector3r, Matrix3r, SpatialVector, SpatialMatrix> InerProjVis;

    Accel_visitor<ScalarType, Vector3r, Matrix3r, SpatialVector> AccelVis;
    Accel_root_visitor<ScalarType, Vector3r, Matrix3r, SpatialVector> Accel_rootVis;


    std::vector< Vector3r, Eigen::aligned_allocator<Vector3r> > Trans, TransConst, Screw_w;
    std::vector< Matrix3r, Eigen::aligned_allocator<Matrix3r> > Rot, RotConst;
    std::vector< SpatialVector, Eigen::aligned_allocator<SpatialVector> > Twists, Cbias, Pbias, P_a, P_A, U, Accel;
    std::vector< SpatialMatrix, Eigen::aligned_allocator<SpatialMatrix> > MConst, M_a, M_A;
    std::vector< real_t, Eigen::aligned_allocator<real_t> > u, invD;


    VectorXr ddq;

    // --------------------------------------------
    // Methods
    // --------------------------------------------

  public:


    //! Articulated-Body Algorithm
    template<typename ConfigVectorT, typename TangentVectorT, typename CotangentVectorT>
    const void
    aba(const Eigen::MatrixBase<ConfigVectorT> & q,
        const Eigen::MatrixBase<TangentVectorT> & v,
        const Eigen::MatrixBase<CotangentVectorT> & tau)
    {
      typedef typename ConfigVectorT::Scalar ScalarT;

      //! First Recursion
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

  };

} // namespace geombdCRTP

#endif // GEOMBD_FORWARD_DYNAMICS_CRTP_HPP
