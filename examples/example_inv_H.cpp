/**
 *    \file examples/example_FwdDyn.cpp
 *    \author Alvaro Paz, Gustavo Arechavaleta
 *    \version 1.0
 *    \date 2021
 *
 *    Example to test the ABA
 *    Copyright (c) 2021 Cinvestav
 *    This library is distributed under the MIT License.
 */

//#define EIGEN_RUNTIME_NO_MALLOC
#define EIGEN_NO_DEBUG

//#define EIGEN_DONT_VECTORIZE

#include "geombd/CRTP/ForwardDynamicsCRTP.hpp"

#include <iostream>
#include <chrono>

#include "geombd/io/parser.hpp"

#define __FU_PATH_PREFIX__ "../../data/TROmodels/"
//std::string urdf_dir = __FU_PATH_PREFIX__ "lwr.urdf"; //7
std::string urdf_dir = __FU_PATH_PREFIX__ "nao_inertial_python.urdf"; //24
//std::string urdf_dir = __FU_PATH_PREFIX__ "HRP2.urdf"; //28
//std::string urdf_dir = __FU_PATH_PREFIX__ "atlas.urdf"; //30

//! Set time variables
const int M = 100000;    // sample size;
typedef double DataType;
int n;

typedef Eigen::Matrix<DataType, Eigen::Dynamic, 1> VectorXr;
VectorXr q;


void loop_iH ( std::shared_ptr< geoCRTP::FwdDynCRTP< DataType > > robotDynamics ) {
  //! Time loop for forward dynamics
  for (int k = 0 ; k < M ; k++) {
      robotDynamics->inverseInertiaMatrix( q.derived() );
    }
}

int main(){
  //! Light-weight parser
  //!------------------------------------------------------------------------------!//
  auto robot = Robot::build_model(urdf_dir);

  //! Forward dynamics object pointer
  auto robotDynamics = std::make_shared< geoCRTP::FwdDynCRTP< DataType > >( robot.value() );

  auto t1 = std::chrono::high_resolution_clock::now();
  auto t2 = std::chrono::high_resolution_clock::now();

  n = robot.value()->nq;
  q   = VectorXr::LinSpaced(n, 1, 2*3.1416);

  //  Eigen::internal::set_is_malloc_allowed(false); //! false to enable it

  //! Perform the inv H loop
  t1 = std::chrono::high_resolution_clock::now();
  loop_iH( robotDynamics );
  t2 = std::chrono::high_resolution_clock::now();

  auto t_total = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

  //  Eigen::internal::set_is_malloc_allowed(true);

  std::cout<<"Inverse of Inertia Matrix = "<<t_total/M<<std::endl;
//  std::cout<<"inv(H) = "<<std::endl<< std::scientific << std::setprecision(20) <<robotDynamics->inv_H<<std::endl;
  std::cout<<"sum(sum(abs( inv(H) ))) = "<< std::scientific << std::setprecision(20) <<robotDynamics->inv_H.cwiseAbs().sum()<<std::endl;
  std::cout<<"error = "<<robotDynamics->inv_H.cwiseAbs().sum()-1.91743564889628964011e+05<<std::endl;

  return 0;
}

