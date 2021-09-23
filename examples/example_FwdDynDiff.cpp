/**
 *    \file examples/example_FwdDynDiff.cpp
 *    \author Alvaro Paz, Gustavo Arechavaleta
 *    \version 1.0
 *    \date 2021
 *
 *    Example to test the ENHANCED ABA Differentiation wrt state
 *    Copyright (c) 2021 Cinvestav
 *    This library is distributed under the MIT License.
 */

//#define EIGEN_RUNTIME_NO_MALLOC

#define EIGEN_NO_DEBUG
#define EIGEN_MPL2_ONLY
#define EIGEN_UNROLLING_LIMIT 30

#include "geombd/CRTP/DForwardDynamicsCRTP.hpp"

#include <iostream>
#include <chrono>
#include "geombd/io/parser.hpp"


#define __FU_PATH_PREFIX__ "../../geombd_crtp/data/TROmodels/"
//std::string urdf_dir = __FU_PATH_PREFIX__ "lwr.urdf"; //7
std::string urdf_dir = __FU_PATH_PREFIX__ "nao_inertial_python.urdf"; //24
//std::string urdf_dir = __FU_PATH_PREFIX__ "HRP2.urdf"; //28
//std::string urdf_dir = __FU_PATH_PREFIX__ "atlas.urdf"; //30


//! Set time variables
const int M = 100000;    // sample size;
typedef double DataType;
int n;

typedef Eigen::Matrix<DataType, Eigen::Dynamic, 1> VectorXr;
VectorXr q, dq, tau;


void loop_ABA ( std::shared_ptr< geo::FwdDynDifCRTP< DataType > > robotDynamics ) {
  //! Time loop for forward dynamics
  for (int k = 0 ; k < M ; k++) {
      robotDynamics->D_aba(q.derived(), dq.derived(), tau.derived());
    }
}


int main(){
  //! Light-weight parser
  //!------------------------------------------------------------------------------!//
  auto robot = Robot::build_model(urdf_dir);
//  robot.value()->dump(std::cout);
//  robot.value()->print_tree(std::cout);
//  robot.value()->print_featherstone(std::cout);

  //! Forward dynamics object pointer
  auto robotDynamics = std::make_shared< geo::FwdDynDifCRTP< DataType > >( robot.value() );

  auto t1 = std::chrono::high_resolution_clock::now();
  auto t2 = std::chrono::high_resolution_clock::now();

  n = robot.value()->nq;
  q   = VectorXr::LinSpaced(n, 1, 2*3.1416);
  dq  = VectorXr::LinSpaced(n, 1, 2*3.1416);
  tau = VectorXr::LinSpaced(n, 1, 2*3.1416);

  //  Eigen::internal::set_is_malloc_allowed(false); //! false to enable it

  //! Perform the ABA loop
  t1 = std::chrono::high_resolution_clock::now();
  loop_ABA( robotDynamics );
  t2 = std::chrono::high_resolution_clock::now();

  auto t_total = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

  //  Eigen::internal::set_is_malloc_allowed(true);

  std::cout<<"Forward dynamics + D = "<<t_total/M<<std::endl;


//  auto ddq_crtp = robotDynamics->ddq;
//  std::cout<<"ddq = "<<ddq_crtp.cwiseAbs().sum()<<std::endl;
//  std::cout << std::scientific << std::setprecision(20) << ddq_crtp.transpose() << std::endl;
//  //!-------------------------------------------------------
//  auto D_ddq_crtp = robotDynamics->D_ddq;
//  std::cout << "D_ddq_crtp = "<< std::endl<< std::scientific << std::setprecision(20) << D_ddq_crtp << std::endl;


//  //!------------------------------------------------------------------------------!//
//  //!--------------Numeric verification through finite differences-----------------!//
//  //!------------------------------------------------------------------------------!//
//  auto ddq_ = robotDynamics->ddq;
//  auto D_ddq_ = robotDynamics->D_ddq;

//  Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic> numericD_ddq(D_ddq_.rows(), D_ddq_.cols());

//  double inc_s = pow(2,-24);  //induced variation

//  for(short int iter = 0 ; iter < n ; iter++){
//      //! wrt q
//      auto q_aux1 = q;
//      q_aux1(iter) += inc_s;

//      robotDynamics->D_aba(q_aux1, dq, tau);

//      auto temporal_ddq = robotDynamics->ddq;

//      numericD_ddq.col(iter) = (temporal_ddq-ddq_)/inc_s;

//      //! wrt dq
//      auto dq_aux1 = dq;
//      dq_aux1(iter) += inc_s;

//      robotDynamics->D_aba(q, dq_aux1, tau);

//      temporal_ddq = robotDynamics->ddq;

//      numericD_ddq.col(iter+n) = (temporal_ddq-ddq_)/inc_s;
//    }

//  std::cout << "-------------------------- D_ddq" << std::endl;
//  std::cout << "Analytic: " << D_ddq_.cwiseAbs().sum() << std::endl;
//  std::cout << "Numeric: " << numericD_ddq.cwiseAbs().sum() << std::endl;
//  auto error_in_D_ddq = D_ddq_ - numericD_ddq;
//  std::cout << "Error: " << error_in_D_ddq.eval().cwiseAbs().sum() << std::endl;
//  std::cout << "--------------------------" << std::endl;

  return 0;
}
