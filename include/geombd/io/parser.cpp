/*#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <sstream>
#include <optional>
#include <memory>
#include <unordered_map>
#include <chrono>
#include <Eigen/Dense>*/
#include "pugixml.hpp"

#include "parser.hpp"

using namespace Eigen;

// support {{{1
RowVectorXd string_to_vector(const std::string& vec, int size)
{
   RowVectorXd res(size);
   std::istringstream ss(vec);
   for (auto i=0; i < size; i++) {
       ss >> res[i];
   }
   return res;
}

/*
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
*/

Inertial::Inertial(const pugi::xml_node & node): Inertial() {

    auto node_o = node.child("origin"); //(optional: defaults to identity if not specified)
    if (!node_o.empty()) { // TODO: Should default to identity if it is not indicated in the URDF
        auto xyz = string_to_vector(node_o.attribute("xyz").value());
        auto rpy = string_to_vector(node_o.attribute("rpy").value());
        origin.xyz = xyz;
        origin.rpy = rpy;
    }
    mass = node.child("mass").attribute("value").as_double(); //default 0.0
    auto node_i = node.child("inertia");
    if (!node_i.empty()) {
        inertia.setMatrix(node_i.attribute("ixx").as_double(),
                node_i.attribute("ixy").as_double(),
                node_i.attribute("ixz").as_double(),
                node_i.attribute("iyy").as_double(),
                node_i.attribute("iyz").as_double(),
                node_i.attribute("izz").as_double());
    }

}

/*

    Inertial(): origin(), mass(0.0), inertia() {}

    std::string Inertial::to_string() const {
        static IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", " << ", ";");
        static IOFormat OctaveFmt(StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
        std::ostringstream out;
        out << "mass (" << mass << "),\n\t origin[xyz] (" << origin.xyz.format(CommaInitFmt);
        out << "), origin[rpy] (" << origin.rpy.format(CommaInitFmt) << "); inertia tensor:\n";
        out << inertia.value.format(OctaveFmt);
        return out.str();
 
    }

*/


Body::Body(const std::string& name, int id)
    : parent(-1), has_inertial(false), Element(name, id)
{ 
}

std::string Body::to_string() const
{
    std::ostringstream out;
    out << "body " << name << ", id " << id;
    if (parent_links.size() > 0)
        out << "; parent (";
    for (auto && parent : parent_links) {
        out << parent << " ";
    }
    out << ")";
    if (child_links.size() > 0) {
        out << "; children ( ";
        for (auto && child : child_links) {
            out << child << " ";
        }
        out << ")";
    }
    if (has_inertial) {
        out << " has inertial properties: ";
        out << inertial_properties.to_string();
    }
    //out << "\n";
    return out.str();
}

void Body::addChild(int child_id) { child_links.push_back(child_id); }
void Body::addParent(int parent_id) { parent_links.push_back(parent_id); }
void Body::addJoint(int joint_id) { parent_joints.push_back(joint_id); }
void Body::setInertialProperties(const Inertial& value) 
{
    inertial_properties = value; has_inertial = true;
}


/*
class Joint : public Element // {{{1
{
    // MAke enum for joint types
public:
    enum joint_type {
        revolute = 'r',  // The type of the joint is revolute. 
        prismatic = 'p', // The type of the joint is prismatic. 
        fixed = 'f',      // The type of the joint is glue. 
        notype = 'n'      // There is no information about the type.  
    };

    Joint(const std::string& name="unassigned", int id=-1): nq(0), Element(name, id) {}
    Joint(const std::string& name, unsigned int id,
            const std::string& parent, const std::string& child, const std::string& type)
        : predecessor(parent), successor(child)
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
*/

Robot::Robot(const std::string& name): name(name), nq(0), joint_p(), joint_s(), body_l(), body_m(), _count(0) {}

void Robot::addJoint(const std::string& name, unsigned int id,
        const std::string& parent, const std::string& child, const std::string& type) {
    auto joint = std::make_shared<Joint>(name, id, parent, child, type);
    joints.push_back(joint);
    joint_dico[name] = id;
    joint_pc_dico[parent + ":" + child] = id;
    //auto parent_id = getBodyID(parent);
    //auto child_id = getBodyID(child);
    auto parent_id = link_dico[parent];
    auto child_id = link_dico[child];

    bodies[parent_id]->addChild(child_id);
    bodies[child_id]->addParent(parent_id);
    bodies[child_id]->addJoint(id);
    this->nq += joints[id]->nq;
    //std::cout << __FUNCTION__ << " dofs: " << this->nq << "\n";
}

void Robot::addBody(const std::string& name, unsigned int id){

    bodies.push_back(std::make_shared<Body>(name, id));
    link_dico[name] = id;
}

void Robot::dump(std::ostream & os) {
    os << "robot name is: " << name << " with " << nq << "dof\n";
    os << "root link: " << bodies[getRootBodyID()]->name << "\n";
    for (auto && body : bodies) {
        os << body->to_string() << "\n";
    }
    for (auto && joint : joints) {
        os << joint->to_string() << "\n";
    }
    // λμ 11 daniel. 12 alvaro, 1 carla y viernes con Yang 12
    os << "\n";
}

void Robot::print_featherstone(std::ostream & os) {
    auto n = joint_p.size();
    auto ii = Eigen::RowVectorXi::LinSpaced(n,0,n-1);
    Eigen::MatrixXi feather(4, n);
    //feather.row(0) = ii;

    Map<RowVectorXi> p_(joint_p.data(), n);
    Map<RowVectorXi> s_(joint_s.data(), n);

    feather << ii ,  p_ , s_, ii;

    feather.row(3) = feather.block(1,0, 2, n).colwise().minCoeff();
    os << feather<< std::endl;
    os << "p(i): ";
    for (auto && i : joint_p) {
        os << " " << i;
    }
    os << "\ns(i): ";
    for (auto && i : joint_s) {
        os << " " << i;
    }
    os << "\nλ(i): ";
    for (auto && i : body_l) {
        os << " " << i;
    }
    os << "\nlambda: " << this->j_tree_ids;
    os << "\nmap: " << this->tree_ids;
    os << "\n";
}
void Robot::print_tree(std::ostream & os, const std::string & prefix, int id, int parent_id) {
    static const std::string spc_pref{"   "};
    auto root = bodies[id];
    auto nb_children = root->child_links.size();
    os << prefix;
    if (parent_id >= 0) {
        auto pre = bodies[parent_id]->name;
        auto curr = root->name;
        auto j_id = this->joint_pc_dico[pre + ":" + curr];
        os << "[" << j_id << "|" << joints[j_id]->tree_id << "|" << joints[j_id]->name << "] ";
    }
    os << "/" << root->tree_id << "/" << root->name << "/" ;
    if (nb_children == 1) {
        os << " has " <<  nb_children << " child\n";
    } else if (nb_children > 1) {
        os << " has " << nb_children << " children\n";
    } else 
        os << " \n";

    for (auto && child : root->child_links) {
        print_tree(os, spc_pref + prefix, child, id);
    }

}
void Robot::print_tree(std::ostream & os) {
    auto root_id = getRootBodyID();
    print_tree(os, "", root_id);
}

int Robot::getBodyID(const std::string& frame_name) const
{
    auto id = this->link_dico.find(frame_name);
    if (id == this->link_dico.end()) {
        return -1;
    } else {
        auto [key, value] = (*id);
        //std::cout << "Found at " << key << " wth id: " << value << std::endl;
        return value;
    }
}

int Robot::getJointID(const std::string& name) const
{
    auto id = this->joint_dico.find(name);
    if (id == this->joint_dico.end()) {
        return -1;
    } else {
        auto [key, value] = (*id);
        //std::cout << "Found at " << key << " with id: " << value << std::endl;
        return value;
    }
}

int Robot::getRootBodyID() const {
    auto root = -1;
    for (auto && body : bodies) {
        if (body->parent_links.size() == 0) {
            root = body->id;
            break;
        }
    }
    return root;
}

        void Robot::add_pre(std::shared_ptr<Body> body, std::shared_ptr<Body> body_p) {
          if (body_p->id != 0) body->Pre.insert(body->Pre.begin(), body_p->id-1);  // -1 for zero init and 0 for one init

          auto parent_v = body_p->parent_links;
          for (auto && parent_i : parent_v) {
              add_pre(body, bodies[parent_i]);
            }
        }

        void Robot::add_suc(std::shared_ptr<Body> body, std::shared_ptr<Body> body_p) {
          body->Suc.push_back(body_p->id-1);  // -1 for zero init and 0 for one init
          auto child_v = body_p->child_links;
          for (auto && child_i : child_v) {
              add_suc(body, bodies[child_i]);
            }
        }

        void Robot::build_PSvector( ) {
          for (auto && body = bodies.begin()+1; body != bodies.end(); body++) {
              add_pre(*body, *body);
              add_suc(*body, *body);
            }
        }

        void Robot::build_spatial_inertias( ) {
          for (auto && body = bodies.begin()+1; body != bodies.end(); body++) {
              if ((*body)->has_inertial) {
                  double mass = (*body)->inertial_properties.mass;
                  auto com = (*body)->inertial_properties.origin.xyz;
                  Matrix3d skew_com;
                  skew_com <<       0, -com(2),  com(1),
                               com(2),       0, -com(0),
                              -com(1),  com(0),       0;
                  (*body)->inertial_properties.spatial << mass*Matrix3d::Identity(),
                      -mass*skew_com,
                      mass*skew_com,
                      (*body)->inertial_properties.inertia.value - mass*skew_com*skew_com;
                }
            }
        }

bool Robot::load_from_urdf(const pugi::xml_node & node) {
    auto id = 0;
    for (pugi::xml_node link: node.children("link"))  // http://wiki.ros.org/urdf/XML/link
    {
        auto name =  link.attribute("name").value();
        this->addBody(name, id);
        auto inertial = link.child("inertial");
        if (!inertial.empty()) { // optional
            Inertial inertial_properties{inertial};
            bodies[id]->setInertialProperties(inertial_properties);  // modify body [id]
        }
        id++;
    }
    id = 0;
    for (pugi::xml_node joint: node.children("joint")) // http://wiki.ros.org/urdf/XML/joint 
    {
        auto name =  joint.attribute("name").value();
        auto type =  joint.attribute("type").value();
        auto predecessor = joint.child("parent").attribute("link").value();
        auto successor = joint.child("child").attribute("link").value();
        this->addJoint(name, id, predecessor, successor, type);
        auto axis = joint.child("axis");
        if (!axis.empty()) joints[id]->axis = string_to_vector(axis.attribute("xyz").value()); //optional: defaults to (1,0,0
        auto origin = joint.child("origin");
        if (!origin.empty()) { // TODO: Should default to identity if it is not indicated in the URDF
            auto xyz = string_to_vector(origin.attribute("xyz").value());
            auto rpy = string_to_vector(origin.attribute("rpy").value());

            joints[id]->setOrigin(xyz, rpy);  // modify joint (id)
        }
        auto limit = joint.child("limit");
        if (!limit.empty()) { // required for revolute and prismatic
            auto upper = limit.attribute("upper").as_double(); // optional, 0.0
            auto lower = limit.attribute("lower").as_double(); // optional, defaults to 0
            auto effort = limit.attribute("effort").as_double();
            auto velocity = limit.attribute("velocity").as_double();

            joints[id]->setLimit(upper, lower, effort, velocity);
        }
        // TODO: Comment if no dynamics properties are considered. 
        /*auto dynamics = joint.child("dynamics");
          if (!dynamics.empty()) {
          auto damping = dynamics.attribute("damping").as_double(); // optional, defaults to 0
          auto friction = dynamics.attribute("friction").as_double(); // optional, defaults to 0

          joints[id]->setDampingAndFriction(damping, friction);
          }*/
        id++;
    }

    this->build_graph();
    return true;
}

void Robot::visit_spanning_tree(int root_id) {
    auto root_name = bodies[root_id]->name;
    for (auto && i : bodies[root_id]->child_links){
        auto curr_name = bodies[i]->name;
        auto & joint_ref = joints[joint_pc_dico[root_name + ":" + curr_name]]->tree_id;
        auto & ref = bodies[i]->tree_id;
        if (ref == -1) {
            ref = ++_count;
            //if (joint_ref == -1)
            //    joint_ref = _count;
        }
        visit_spanning_tree(i);
    }   
}

void Robot::build_graph() {
    auto root_id = getRootBodyID();
    if (root_id != 0) {
        std::cerr << "Warning: first body is not root.\n";
    }
    auto nb_joints = joints.size();
    auto nb_links = bodies.size();

    if (nb_joints > nb_links - 1) {
        std::cerr << "Warning: there maybe cycles (closed kinematic chains).\n"; // kinematic trees and closed-loop systems.
    }

    this->_count = 0;
    bodies[root_id]->tree_id = 0;
    visit_spanning_tree(root_id);

    //j_tree_ids.resize(nb_joints + 1);
    //j_tree_ids[0] = -1;
    j_tree_ids.setConstant(nb_joints + 1, -1);
    joint_p.resize(nb_joints);
    joint_s.resize(nb_joints);
    auto i = 0;
    for (auto && joint : joints) {
        joint_p[i] = link_dico[joint->predecessor];
        joint_s[i] = link_dico[joint->successor];
        i++;
    }

    tree_ids.resize(nb_links);
    for (auto && link : bodies) {
        tree_ids[link->tree_id] = link->id;
    }
    /*for (auto i = 1; i < nb_links; i++) {
      auto idx = bodies[tree_ids[i]]->parent_joints[0];
      j_tree_ids[i] = idx;
      }*/

    for (auto && joint : joints) {
        //std::cout << "j id: (" << joint->id << ";" << joint->tree_id << ")  ";
        //j_tree_ids[joint->tree_id] = joint->id;
    }
    //std::set<int> joint_list{};

    //std::cout << "Regular numbering:\n" << tree_ids << "\n";

    body_l.resize(nb_links);
    body_m.resize(nb_links); //TODO: validate
    tree_ids.resize(nb_links);
    i = 0;
    auto j=0;
    for (auto && link : bodies) {
        if (link->parent_links.size() == 0) {
            body_l[i] = -1;
        } else {
            body_l[i] = std::min(joint_p[j], joint_s[j]);
            j++;
        }
        //std::cout << i << " < λ(i)= " << body_l[i];
        i++; 

    }

}


/**
 *
 *
 * https://github.com/ros/urdfdom/blob/master/xsd/urdf.xsd
 */
std::optional<std::shared_ptr<Robot>> Robot::build_model(const std::string& urdf_file)
{
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(urdf_file.c_str());
    if (!result)
        return {};
    auto node = doc.child("robot");
    auto name = node.attribute("name").value(); // http://wiki.ros.org/urdf/XML/robot
    auto robot = std::make_shared<Robot>(name);
    robot->load_from_urdf(node); // TODO: modify to send false for unknown joint types
    robot->build_PSvector( );
    robot->build_spatial_inertias( );    
    return {robot};

}
