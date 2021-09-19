/**
 *    \file examples/example_FwdDyn.cc
 *    \author Alvaro Paz, Gustavo Arechavaleta
 *    \version 1.0
 *    \date 2020
 *
 *    Example to test the ABA
 */

//#define EIGEN_RUNTIME_NO_MALLOC
#define EIGEN_NO_DEBUG

//#define EIGEN_DONT_VECTORIZE

#include "geombd/CRTP/ForwardDynamicsCRTP.hpp"

#include <iostream>
#include <chrono>

#include "geombd/io/parser.hpp"

#define __FU_PATH_PREFIX__ "../data/TROmodels/"
std::string urdf_dir = __FU_PATH_PREFIX__ "nao_inertial_python.urdf";

//! Set time variables
const int M = 100000;    // sample size;
typedef double DataType;
int n;

typedef Eigen::Matrix<DataType, Eigen::Dynamic, 1> VectorXr;
VectorXr q, dq, tau;


void loop_ABA ( std::shared_ptr< geoCRTP::FwdDynCRTP< DataType > > robotDynamics ) {
  //! Time loop for forward dynamics
  for (int k = 0 ; k < M ; k++) {
      robotDynamics->aba(q.derived(), dq.derived(), tau.derived());
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
  dq  = VectorXr::LinSpaced(n, 1, 2*3.1416);
  tau = VectorXr::LinSpaced(n, 1, 2*3.1416);

  //  Eigen::internal::set_is_malloc_allowed(false); //! false to enable it

  //! Perform the ABA loop
  t1 = std::chrono::high_resolution_clock::now();
  loop_ABA( robotDynamics );
  t2 = std::chrono::high_resolution_clock::now();

  auto t_total = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

  //  Eigen::internal::set_is_malloc_allowed(true);

  std::cout<<"Forward dynamics = "<<t_total/M<<std::endl;

  std::cout<<"FD = "<<robotDynamics->ddq.transpose()<<std::endl;














//  Eigen::internal::set_is_malloc_allowed(false);
//  Eigen::MatrixXd R1 = Eigen::MatrixXd::Random(23,23);
//  Eigen::MatrixXd R2 = Eigen::MatrixXd::Random(23,23);
//  R1 += R1*R2;
//  Eigen::internal::set_is_malloc_allowed(true);

//  Eigen::internal::set_is_malloc_allowed(true);

//  std::cout<<"Forward dynamics = "<<t_total_00/M<<std::endl;




  return 0;
}

