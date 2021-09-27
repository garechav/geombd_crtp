/**
 *    \file include/geombd/io/parser.hpp
 *    \author Keny Ordaz
 *    \version 1.0
 *    \date 2021
 *
 *    Class to parse a multibody system from URDF
 *    Copyright (c) 2021 Cinvestav
 *    This library is distributed under the MIT License.
 */

#ifndef HEADER_IO_PARSER_HPP
#define HEADER_IO_PARSER_HPP
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <sstream>
#include <optional>
#include <memory>
#include <unordered_map>
//#include <chrono>
#include <Eigen/Core>
//#include "pugixml.hpp"

namespace pugi {
    class xml_node;
}
using namespace Eigen;

// support {{{1
RowVectorXd string_to_vector(const std::string& vec, int size=3);

class Pose
{
public:
    Pose(): xyz(0, 0, 0), rpy(0, 0, 0) {}

    void init_rpy(double roll, double pitch, double yaw)
    {
        auto phi = roll / 2.0;
        auto the = pitch / 2.0;
        auto psi = yaw / 2.0;
        auto sphi = sin(phi);  auto cphi = cos(phi);
        auto sthe = sin(the);  auto cthe = cos(the);
        auto spsi = sin(psi);  auto cpsi = cos(psi);

        this->x = sphi * cthe * cpsi - cphi * sthe * spsi;
        this->y = cphi * sthe * cpsi + sphi * cthe * spsi;
        this->z = cphi * cthe * spsi - sphi * sthe * cpsi;
        this->w = cphi * cthe * cpsi + sphi * sthe * spsi;

        this->normalize();
    };
    void normalize()
    {
        double s = sqrt(this->x * this->x +
                this->y * this->y +
                this->z * this->z +
                this->w * this->w);
        if (s == 0.0) {
            this->x = 0.0;
            this->y = 0.0;
            this->z = 0.0;
            this->w = 1.0;
        } else {
            this->x /= s;
            this->y /= s;
            this->z /= s;
            this->w /= s;
        }
    };
    double x, y, z, w;
    RowVector3d xyz;
    RowVector3d rpy;
};

class Axis
{
public:
    Axis(): xyz(1.0, 0., 0.) {}

    RowVector3d xyz;
};

class Inertia { // optional zero inertia
public:
    Inertia(double ixx=0.0, double ixy=0.0, double ixz=0.0, double iyy=0.0, double iyz=0.0, double izz=0.0)
        : ixx(ixx), ixy(ixy), ixz(ixz), iyy(iyy), iyz(iyz), izz(izz) 
        { setMatrix(ixx, ixy, ixz, iyy, iyz, izz); }
    void setMatrix(double ixx=0.0, double ixy=0.0, double ixz=0.0, double iyy=0.0, double iyz=0.0, double izz=0.0)
        { value << ixx, ixy, ixz, ixy, iyy, iyz, ixz, iyz, izz; }

    double ixx;
    double ixy;
    double ixz;
    double iyy;
    double iyz;
    double izz;
    Matrix3d value;
};

class Inertial {
public:
    Inertial(): origin(), mass(0.0), inertia(), spatial() {}
    Inertial(const pugi::xml_node & node);
    std::string to_string() const {
        static IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", " << ", ";");
        static IOFormat OctaveFmt(StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
        std::ostringstream out;
        out << "mass (" << mass << "),\n\t origin[xyz] (" << origin.xyz.format(CommaInitFmt);
        out << "), origin[rpy] (" << origin.rpy.format(CommaInitFmt) << "); inertia tensor:\n";
        out << inertia.value.format(OctaveFmt);
        return out.str();
 
    }
    Pose origin;
    double mass;
    Inertia inertia;
    Matrix<double, 6, 6> spatial;
};
// support }}}1

class Element // {{{2
{
public:
    Element(const std::string& name, int id): name(name), id(id), tree_id(-1) {}
    std::string name;
    int id;
    int tree_id;
    //~Element() { std::cout << __PRETTY_FUNCTION__ << "\n"; }
};  // }}}2

class Body : public Element  // {{{1
{
public:
    Body(const std::string& name="unassigned", int id=-1);
    std::string to_string() const;
    void addChild(int child_id);
    void addParent(int parent_id);
    void addJoint(int joint_id);
    void setInertialProperties(const Inertial& value);
    bool has_inertial;
    Inertial inertial_properties;
    std::vector<int> child_links;
    std::vector<int> parent_links;
    std::vector<int> parent_joints;
    int parent;    
    std::vector<int> Pre;
    std::vector<int> Suc;
};  // }}}1

class Joint : public Element // {{{1
{
    // MAke enum for joint types
public:
    enum joint_type {
        revolute = 'r',  /*!< The type of the joint is revolute. */
        prismatic = 'p', /*!< The type of the joint is prismatic. */
        fixed = 'f',      /*!< The type of the joint is glue. */
        notype = 'n'      /*!< There is no information about the type. */ 
    };

    Joint(const std::string& name="unassigned", int id=-1): nq(0), Element(name, id) {}
    Joint(const std::string& name, unsigned int id,
            const std::string& parent, const std::string& child, const std::string& type)
        : /*name(name), id(id),*/ predecessor(parent), successor(child)
        , axis(1.0, 0, 0), Element(name, id) {
            if (type == "revolute") {
                this->type = Joint::joint_type::revolute;
            } else if (type == "prismatic") {
                this->type = Joint::joint_type::prismatic;
            } else if (type == "fixed") {
                this->type = Joint::joint_type::fixed;
            } else {
                this->type = Joint::joint_type::notype;

            }

            switch(this->type) {
            case joint_type::revolute:
            case joint_type::prismatic:
                // TODO: check limits
                setLimit();
                nq = 1;
                break;
            default:
                nq = 0;// Do nothing?;
            }
        }
    //~Joint() { std::cout << __PRETTY_FUNCTION__ << "\n"; }

    std::string to_string() {
        std::ostringstream out;
        out << "joint " << name << ", id " << id << ", type " << std::string(1, type);
        out << ", axis " << axis;
        if (type == joint_type::revolute || type == joint_type::prismatic) {
            out << ", limits(" << upper << ", " << lower << ", " << effort << ", " << velocity << ") ";
        }
        out << joint_axis_indicator();
        out << "\t/" << successor << "/" << predecessor << "/ ";
        return out.str();
    }

    std::string joint_axis_indicator() {
        std::ostringstream out;
        out << "[";
        switch(type) {
        case Joint::joint_type::revolute:
            out << "⚙";
            break;
        case Joint::joint_type::prismatic:
            out << "↕";
            break;
        case Joint::joint_type::fixed:
            out << "⚓";
            break;
        case Joint::joint_type::notype:
            out << "⦻";
        }
        if (type == joint_type::revolute || type == joint_type::prismatic) {
            if (axis == RowVector3d(1.0, 0.0, 0.0))
                out << "+X";
            else if (axis == RowVector3d(0.0, 1.0, 0.0))
                out << "+Y";
            else if (axis == RowVector3d(-1.0, 0.0, 0.0))
                out << "-X";
            else if (axis == RowVector3d(0.0, -1.0, 0.0))
                out << "-Y";
            else if (axis == RowVector3d(0.0, 0.0, 1.0))
                out << "+Z";
            else if (axis == RowVector3d(0.0, 0.0, -1.0))
                out << "-Z";
            else
                out << "⦻";
        }

        out << "]";
        return out.str();
    }

    void setLimit(double upper=0.0, double lower=0.0, double effort=0.0, double velocity=0.0) {
        this->upper = upper; this->lower = lower; this->effort = effort; this->velocity = velocity;
    }
    
    void setOrigin(const RowVector3d& xyz, const RowVector3d& rpy) {
        origin.xyz = xyz;
        origin.rpy = rpy;
    }

    void setDampingAndFriction(double damping=0.0, double friction=0.0) {} //FIXME: If required
    std::string predecessor;
    std::string successor;
    joint_type type;
    RowVector3d axis;
    Pose origin;
    double upper, lower, effort, velocity;
    double damping, friction;
    int nq;
};  // }}}1

class Robot { // {{{1
public:
    Robot(const std::string& name="unknown");
    void addJoint(const std::string& name, unsigned int id,
            const std::string& parent, const std::string& child, const std::string& type);
    void addBody(const std::string& name, unsigned int id);
    void dump(std::ostream & os);
    void print_featherstone(std::ostream & os);
    void print_tree(std::ostream & os, const std::string & prefix, int id, int parent_id=-1);
    void print_tree(std::ostream & os);
    int getBodyID(const std::string& frame_name) const;
    int getJointID(const std::string& name) const;
    std::shared_ptr<Joint> & operator()(int i) { // robot(j)
        return joints[i];
    }
    std::shared_ptr<Body> & operator[](int i) { // robot[i]
        return bodies[i];
    }
    int getRootBodyID() const;
    const std::string& getName() const { return name; } 
    static std::optional<std::shared_ptr<Robot>> build_model(const std::string& urdf_file);
    void add_pre(std::shared_ptr<Body> body, std::shared_ptr<Body> body_s);
    void add_suc(std::shared_ptr<Body> body, std::shared_ptr<Body> body_p);
    void build_PSvector( );
    void build_spatial_inertias( );
protected:
    bool load_from_urdf(const pugi::xml_node & node);
    void visit_spanning_tree(int root_id);
    void build_graph();
protected:
    std::string name;
    std::vector<std::shared_ptr<Joint>> joints;
    std::vector<std::shared_ptr<Body>> bodies;
    std::unordered_map<std::string,int> joint_dico;
    std::unordered_map<std::string,int> joint_pc_dico;
    std::unordered_map<std::string,int> link_dico;
public:
    std::vector<int> body_l;
    std::vector<int> body_m;
    std::vector<int> joint_p;
    std::vector<int> joint_s;
    Eigen::RowVectorXi tree_ids;
    Eigen::RowVectorXi j_tree_ids;
public:
    int nq;
private:
    int _count;
};  // }}}1


#endif
