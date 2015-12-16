//
// Copyright (c) 2015 CNRS
//
// This file is part of Pinocchio
// Pinocchio is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// Pinocchio is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Lesser Public License for more details. You should have
// received a copy of the GNU Lesser General Public License along with
// Pinocchio If not, see
// <http://www.gnu.org/licenses/>.

#ifndef __se3_urdf_geom_hpp__
#define __se3_urdf_geom_hpp__

#include <urdf_model/model.h>
#include <urdf_parser/urdf_parser.h>

#include <iostream>
#include <boost/foreach.hpp>
#include "pinocchio/multibody/model.hpp"

#include <hpp/fcl/collision_object.h>
#include <hpp/fcl/collision.h>
#include <hpp/fcl/shape/geometric_shapes.h>
#include "pinocchio/multibody/parser/from-collada-to-fcl.hpp"

#include <exception>

namespace urdf
{
  typedef boost::shared_ptr<ModelInterface> ModelInterfacePtr;
  typedef boost::shared_ptr<const Joint> JointConstPtr;
  typedef boost::shared_ptr<const Link> LinkConstPtr;
  typedef boost::shared_ptr<Link> LinkPtr;
  typedef boost::shared_ptr<const Inertial> InertialConstPtr;
}

namespace se3
{
  namespace urdf
  {


fcl::CollisionObject retrieveCollisionGeometry (const ::urdf::LinkConstPtr & link, std::string meshRootDir)
  {
    boost::shared_ptr < ::urdf::Collision> collision = link->collision;
    boost::shared_ptr < fcl::CollisionGeometry > geometry;

    // Handle the case where collision geometry is a mesh
    if (collision->geometry->type == ::urdf::Geometry::MESH)
    {
      boost::shared_ptr < ::urdf::Mesh> collisionGeometry = boost::dynamic_pointer_cast< ::urdf::Mesh> (collision->geometry);
      std::string collisionFilename = collisionGeometry->filename;

      std::string full_path = fromURDFMeshPathToAbsolutePath(collisionFilename, meshRootDir);

      ::urdf::Vector3 scale = collisionGeometry->scale;

      // Create FCL mesh by parsing Collada file.
      PolyhedronPtrType  polyhedron (new PolyhedronType);

      loadPolyhedronFromResource (full_path, scale, polyhedron);
      geometry = polyhedron;
    }
    

    // Compute body position in world frame.
    // MatrixHomogeneousType position = computeBodyAbsolutePosition (link, collision->origin);
    if (!geometry)
    {
      throw std::runtime_error(std::string("The polyhedron retrived is empty"));
    }
    fcl::CollisionObject collisionObject (geometry, fcl::Transform3f());

    return collisionObject; // TO CHECK: what happens if geometry is empty ?
    
  }

    void parseTreeWithGeom( ::urdf::LinkConstPtr link, Model & model,GeometryModel & model_geom, std::string meshRootDir, const SE3 & placementOffset = SE3::Identity()) throw (std::invalid_argument)
{


  ::urdf::JointConstPtr joint = link->parent_joint;
  SE3 nextPlacementOffset = SE3::Identity(); // OffSet of the next link. In case we encounter a fixed joint, we need to propagate the length of its attached body to next joint.

  // std::cout << " *** " << link->name << "    < attached by joint ";
  // if(joint)
  //   std::cout << "#" << link->parent_joint->name << std::endl;
  // else std::cout << "###ROOT" << std::endl;

  // std::cout << "placementOffset: " << placementOffset << std::endl;

  if(joint!=NULL)
  {
    assert(link->getParent()!=NULL);

    if (!link->inertial && joint->type != ::urdf::Joint::FIXED)
    {
      const std::string exception_message (link->name + " - spatial inertia information missing.");
      throw std::invalid_argument(exception_message);
    }

    Model::Index parent = (link->getParent()->parent_joint==NULL) ? (model.existJointName("root_joint") ? model.getJointId("root_joint") : 0) :
                                                                    model.getJointId( link->getParent()->parent_joint->name );
    //std::cout << joint->name << " === " << parent << std::endl;

    const SE3 & jointPlacement = placementOffset*convertFromUrdf(joint->parent_to_joint_origin_transform);

    const Inertia & Y = (link->inertial) ?  convertFromUrdf(*link->inertial) :
                                          Inertia::Identity();

    bool visual = (link->visual) ? true : false;


    //std::cout << "Parent = " << parent << std::endl;
    //std::cout << "Placement = " << (Matrix4)jointPlacement << std::endl;

    switch(joint->type)
    {
      case ::urdf::Joint::REVOLUTE:
      case ::urdf::Joint::CONTINUOUS: // Revolute with no joint limits
      {

        Eigen::VectorXd maxEffort;
        Eigen::VectorXd velocity;
        Eigen::VectorXd lowerPosition;
        Eigen::VectorXd upperPosition;

        if (joint->limits)
        {
          maxEffort.resize(1);      maxEffort     << joint->limits->effort;
          velocity.resize(1);       velocity      << joint->limits->velocity;
          lowerPosition.resize(1);  lowerPosition << joint->limits->lower;
          upperPosition.resize(1);  upperPosition << joint->limits->upper;
        }

        Eigen::Vector3d jointAxis(Eigen::Vector3d::Zero());
        AxisCartesian axis = extractCartesianAxis(joint->axis);
        switch(axis)
        {
          case AXIS_X:
            model.addBody(  parent, JointModelRX(), jointPlacement, Y,
                            maxEffort, velocity, lowerPosition, upperPosition,
                            joint->name,link->name, visual );
            break;
          case AXIS_Y:
            model.addBody(  parent, JointModelRY(), jointPlacement, Y,
                            maxEffort, velocity, lowerPosition, upperPosition,
                            joint->name,link->name, visual );
            break;
          case AXIS_Z:
            model.addBody(  parent, JointModelRZ(), jointPlacement, Y,
                            maxEffort, velocity, lowerPosition, upperPosition,
                            joint->name,link->name, visual );
            break;
          case AXIS_UNALIGNED:
            jointAxis= Eigen::Vector3d( joint->axis.x,joint->axis.y,joint->axis.z );
            jointAxis.normalize();
            model.addBody(  parent, JointModelRevoluteUnaligned(jointAxis), 
                            jointPlacement, Y,
                            maxEffort, velocity, lowerPosition, upperPosition,
                            joint->name,link->name, visual );
            break;
          default:
            assert( false && "Fatal Error while extracting revolute joint axis");
            break;
        }
        break;
      }
      case ::urdf::Joint::PRISMATIC:
      {
        Eigen::VectorXd maxEffort = Eigen::VectorXd(0.);
        Eigen::VectorXd velocity = Eigen::VectorXd(0.);
        Eigen::VectorXd lowerPosition = Eigen::VectorXd(0.);
        Eigen::VectorXd upperPosition = Eigen::VectorXd(0.);

        if (joint->limits)
        {
          maxEffort.resize(1);      maxEffort     << joint->limits->effort;
          velocity.resize(1);       velocity      << joint->limits->velocity;
          lowerPosition.resize(1);  lowerPosition << joint->limits->lower;
          upperPosition.resize(1);  upperPosition << joint->limits->upper;
        }
        AxisCartesian axis = extractCartesianAxis(joint->axis);   
        switch(axis)
        {
          case AXIS_X:
            model.addBody(  parent, JointModelPX(), jointPlacement, Y,
                            maxEffort, velocity, lowerPosition, upperPosition,
                            joint->name,link->name, visual );
            break;
          case AXIS_Y:
            model.addBody(  parent, JointModelPY(), jointPlacement, Y,
                            maxEffort, velocity, lowerPosition, upperPosition,
                            joint->name,link->name, visual );
            break;
          case AXIS_Z:
            model.addBody(  parent, JointModelPZ(), jointPlacement, Y,
                            maxEffort, velocity, lowerPosition, upperPosition,
                            joint->name,link->name, visual );
            break;
          case AXIS_UNALIGNED:
            std::cerr << "Bad axis = (" << joint->axis.x <<"," << joint->axis.y << "," << joint->axis.z << ")" << std::endl;
            assert(false && "Only X, Y or Z axis are accepted." );
            break;
          default:
            assert( false && "Fatal Error while extracting prismatic joint axis");
            break;
        }
        break;
      }
      case ::urdf::Joint::FIXED:
      {
        // In case of fixed joint, if link has inertial tag:
        //    -add the inertia of the link to his parent in the model
        // Otherwise do nothing.
        // In all cases:
        //    -let all the children become children of parent
        //    -inform the parser of the offset to apply
        //    -add fixed body in model to display it in gepetto-viewer
        if (link->inertial)
        {
          model.mergeFixedBody(parent, jointPlacement, Y); //Modify the parent inertia in the model
        }

        SE3 ptjot_se3 = convertFromUrdf(link->parent_joint->parent_to_joint_origin_transform);

        //transformation of the current placement offset
        nextPlacementOffset = placementOffset*ptjot_se3;

        //add the fixed Body in the model for the viewer
        model.addFixedBody(parent,nextPlacementOffset,link->name,visual);

        BOOST_FOREACH(::urdf::LinkPtr child_link,link->child_links) 
        {
          child_link->setParent(link->getParent() );  //skip the fixed generation
        }
        break;
      }
      default:
      {
        std::cerr << "The joint type " << joint->type << " is not supported." << std::endl;
        assert(false && "Only revolute, prismatic and fixed joints are accepted." );
        break;
      }
    }
    if (link->collision)
    {
      fcl::CollisionObject collision_object = retrieveCollisionGeometry(link, meshRootDir);
      SE3 geomPlacement = convertFromUrdf(link->collision->origin);
      std::string collision_object_name = link->name ;
      model_geom.addGeomObject(model.getJointId(joint->name), collision_object, geomPlacement, collision_object_name);
    }      
  }
  else if (link->getParent() != NULL)
  {
    const std::string exception_message (link->name + " - joint information missing.");
    throw std::invalid_argument(exception_message);
  }

  BOOST_FOREACH(::urdf::LinkConstPtr child,link->child_links)
  {
    parseTreeWithGeom(child, model, model_geom, meshRootDir, nextPlacementOffset);
  }
}


    template <typename D>
    void parseTreeWithGeom( ::urdf::LinkConstPtr link, Model & model, GeometryModel & model_geom, std::string meshRootDir, const SE3 & placementOffset , const JointModelBase<D> &  root_joint  )
    {
      const Inertia & Y = (link->inertial) ?
      convertFromUrdf(*link->inertial)
      : Inertia::Identity();
      model.addBody( 0, root_joint, placementOffset, Y , "root_joint", link->name, true );
      BOOST_FOREACH(::urdf::LinkConstPtr child,link->child_links)
      {
        parseTreeWithGeom(child, model, model_geom, meshRootDir, SE3::Identity());
      }
    }


    template <typename D>
    std::pair<Model, GeometryModel > buildModelAndGeom( const std::string & filename, const std::string & meshRootDir, const JointModelBase<D> &  root_joint )
    {
      Model model; GeometryModel model_geom;

      ::urdf::ModelInterfacePtr urdfTree = ::urdf::parseURDFFile (filename);
      parseTreeWithGeom(urdfTree->getRoot(), model, model_geom, meshRootDir, SE3::Identity(), root_joint);
      return std::make_pair(model, model_geom);
    }

    std::pair<Model, GeometryModel > buildModelAndGeom( const std::string & filename, const std::string & meshRootDir)
    {
      Model model; GeometryModel model_geom;

      ::urdf::ModelInterfacePtr urdfTree = ::urdf::parseURDFFile (filename);
      parseTreeWithGeom(urdfTree->getRoot(), model, model_geom, meshRootDir, SE3::Identity());
      return std::make_pair(model, model_geom);
    }

  } // namespace urdf
} // namespace se3

#endif // ifndef __se3_urdf_geom_hpp__

