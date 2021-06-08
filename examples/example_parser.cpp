/**
 *    \file examples/example_parser.cc
 *    \author Alvaro Paz, Gustavo Arechavaleta
 *    \version 1.0
 *    \date 2020
 *
 *    Example to verify the parser and multibody
 */

#include <iostream>
#include <vector>
#include <chrono>

#include "geombd/core.h"
#include "geombd/io/parser.hpp"

#define __FU_PATH_PREFIX__ "../GeoMBD/"
// #define __FU_PATH_PREFIX__ "../../../"

//std::string urdf_nao = __FU_PATH_PREFIX__ "data/TROmodels/nao_inertial_XYZ_python.urdf";
std::string urdf_nao = __FU_PATH_PREFIX__ "data/TROmodels/atlas.urdf";

std::string naoFile = urdf_nao;



void test_comparison(const std::shared_ptr<geo::MultiBody>& old_robot ,const std::shared_ptr<Robot>& new_robot)
{
    using namespace std;

    cout << old_robot->getName() << endl;
    cout << new_robot->getName() << endl << endl;

    cout << "Degrees of freedom = " << old_robot->getDoF() << " and " << new_robot->nq << endl << endl;


    int ID = 1;  int j;
    typename std::vector< geo::Body* >::iterator bodyIterator;
    geo::DynamicBody* tempBody;
    geo::DynamicBody* parentBody;
    std::vector< geo::Body* > bodies;
    bodies = old_robot->getPubBodies();

    bool ifPrint = true;
    double total_e = 0.0;
    for ( bodyIterator = bodies.begin()+1; bodyIterator != bodies.end(); bodyIterator++ )
      {
        tempBody = (geo::DynamicBody *)(*bodyIterator);

        j = old_robot->getParentId(tempBody->getId());
        parentBody = (geo::DynamicBody *)(bodies.at(j));

//        //! Link name
//        //!-------------------------------------------------------
//        auto old_name = tempBody->getName();
//        auto new_name = new_robot->operator [](ID)->name;

////        cout << ID << " -> " << old_name << " -> " << new_name << endl;

//        //! XYZ position
//        //!-------------------------------------------------------
//        Vector3d old_P = tempBody->getHomeConfig().block<3,1>(0,3);
//        auto new_P = new_robot->operator ()(ID-1)->origin.xyz;
//        auto error_P = old_P - new_P.transpose();
//        auto abs_P_error = error_P.cwiseAbs().sum();
//        if (abs_P_error == 0.0) {
//            if (ifPrint)  cout << "P_ in " << ID-1 << " joint is OK." << endl;
//          } else {
//            cout << "WARNING in P_ " << ID-1 << " -> " << "old = " << old_P.transpose()
//                 << " and new = " << new_P << " ---> error = " << abs_P_error << endl;
//          }


//        //! Rotation
//        //!-------------------------------------------------------
//        Matrix3d R_ = tempBody->getHomeConfig().block<3,3>(0,0);
//        Matrix3d R_Id = Matrix3d::Identity();
//        Matrix3d R_error = R_ - R_Id;
//        double e_R = R_error.cwiseAbs().sum();
//        if (e_R == 0.0) {
//            if (ifPrint)  cout << "R_ in " << ID-1 << " joint is OK." << endl;
//          } else {
//            cout << "WARNING in R_ " << ID-1 << " -> " << "old = "<< R_ << " ---> error = " << e_R << endl;
//          }


//        //! Mass
//        //!-------------------------------------------------------
//        auto old_m = tempBody->getMass();
//        auto new_m = new_robot->operator [](ID)->inertial_properties.mass;
//        auto error_m = old_m - new_m;
//        auto abs_error_m = error_m;
//        if (abs_error_m == 0.0) {
//            if (ifPrint)  cout << "Mass in " << ID << " link is OK." << endl;
//          } else {
//            cout << "WARNING in Mass at " << ID << " -> " << "old = " << old_m
//                 << " and new = " << new_m << " ---> error = " << abs_error_m << endl;
//          }


//        //! Center of mass
//        //!-------------------------------------------------------
//        auto old_com = tempBody->getCoM();
//        auto new_com = new_robot->operator [](ID)->inertial_properties.origin.xyz;
//        auto error_com = old_com - new_com.transpose();
//        auto abs_error = error_com.cwiseAbs().sum();
//        if (abs_error == 0.0) {
//            if (ifPrint)  cout << "COM in " << ID << " link is OK." << endl;
//          } else {
//            cout << "WARNING in COM " << ID << " -> " << "old = " << old_com.transpose()
//                 << " and new = " << new_com << " ---> error = " << abs_error << endl;
//          }


//        //! Inertia tensor
//        //!-------------------------------------------------------
//        Matrix3d old_I = tempBody->getInertiaTensor();
//        Matrix3d new_I = new_robot->operator [](ID)->inertial_properties.inertia.value;
//        Matrix3d error_I = old_I - new_I;
//        double abs_error_I = error_I.cwiseAbs().sum();
//        if (abs_error_I == 0.0) {
//            if (ifPrint)  cout << "Inertia in " << ID << " link is OK." << endl;
//          } else {
//            cout << "WARNING in Inertia " << ID << " -> " << "old = " << old_I
//                 << " and new = " << new_I << " ---> error = " << abs_error_I << endl;
//          }


//        //! Spatial inertia
//        //!-------------------------------------------------------
//        typedef Matrix<double, 6, 6> SpatialMatrix;
//        SpatialMatrix old_spatial = tempBody->getSpatialInertia();
//        SpatialMatrix new_spatial = new_robot->operator [](ID)->inertial_properties.spatial;
//        SpatialMatrix error_spatial = old_spatial - new_spatial;
//        double abs_error_spatial = error_spatial.cwiseAbs().sum();
//        if (abs_error_spatial == 0.0) {
//            if (ifPrint)  cout << "Spatial inertia in " << ID << " link is OK." << endl;
//          } else {
//            cout << "WARNING in spatial inertia " << ID << " -> " << "old = " << endl << old_spatial
//                 << " and new = " << endl << new_spatial << " ---> error = " << abs_error_spatial << endl;
//          }


//        //! Cummulative error in parameters
//        //!-------------------------------------------------------
//        total_e +=  abs_P_error + e_R + abs_error_m + abs_error + abs_error_I + abs_error_spatial;


        //!-------------------------------------------------------!//
        //!                   Robot Topology                      !//
        //!-------------------------------------------------------!//

//        auto pred = old_robot->getPredecessorJoints(tempBody->getId());
//        cout << ID-1 << " -> predOld => ";
//        for (int iter = 0; iter < pred.size(); iter++)  cout << pred[iter]-1 << ", ";
//        cout << endl;
        //!-------------------------------------------------------
        auto preNew = new_robot->operator [](ID)->Pre;
        cout << ID << " -> predNew => ";
        for (int iter = 0; iter < preNew.size(); iter++)  cout << preNew[iter]+1 << ", ";
        cout << endl;
        //!-------------------------------------------------------
//        auto suce = old_robot->getSubTreeBody(tempBody->getId());
//        cout << ID-1 << " -> suceOld => ";
//        for (int iter = 0; iter < suce.size(); iter++)  cout << suce[iter]-1 << ", ";
//        cout << endl;
        //!-------------------------------------------------------
        auto sucNew = new_robot->operator [](ID)->Suc;
        cout << ID << " -> suceNew => ";
        for (int iter = 0; iter < sucNew.size(); iter++)  cout << sucNew[iter]+1 << ", ";
        cout << endl << endl << endl;

        ID++;
      }

    if (total_e == 0.0) {
        cout << "Kinematic and dynamic parameters were succesfully verifyed." << endl;
      } else {
        cout << "WARNING -> Parameter mismatchs have been found." << endl;
      }

}


int main( ){


  //! Multibody
  //!------------------------------------------------------------------------------!//
  geo::World world = geo::World();
  std::string sFile = naoFile;
  int robotId = world.getRobotsVector()->size();

  world.loadMultiBodyURDF(sFile,robotId, geo::kNAO);

  std::shared_ptr<geo::MultiBody> old_robot = world.getRobot(0);


  //! Light-weight parser
  //!------------------------------------------------------------------------------!//
  std::cout << "Loading filename [" << urdf_nao << "]...\n";
  auto new_robot = Robot::build_model(urdf_nao);
  new_robot.value()->dump(std::cout);
  new_robot.value()->print_tree(std::cout);
  new_robot.value()->print_featherstone(std::cout);


  //! Light-weight parser
  //!------------------------------------------------------------------------------!//
  test_comparison(old_robot, new_robot.value());


  return 0;
}

