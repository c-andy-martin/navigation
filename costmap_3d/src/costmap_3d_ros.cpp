/*********************************************************************
 *
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2019, Badger Technologies LLC
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of Willow Garage, Inc. nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * Author: C. Andy Martin
 *********************************************************************/
// Uncomment to compile in auto profiling support for optimizing queries
//#define COSTMAP_3D_ROS_AUTO_PROFILE_QUERY 1
#if (COSTMAP_3D_ROS_AUTO_PROFILE_QUERY) > 0
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>
#endif
#include <cmath>
#include <costmap_3d/costmap_3d_ros.h>
#include <costmap_3d/layered_costmap_3d.h>
#include <tf/transform_datatypes.h>
#include <visualization_msgs/Marker.h>

namespace costmap_3d
{

Costmap3DROS::Costmap3DROS(std::string name, tf::TransformListener& tf)
  : super(name, tf),
    private_nh_("~/" + name + "/"),
    layered_costmap_3d_(layered_costmap_),
    initialized_(false),
    plugin_loader_("costmap_3d", "costmap_3d::Layer3D"),
    publisher_(private_nh_, &layered_costmap_3d_, "costmap_3d"),
    dsrv_(ros::NodeHandle("~/" + name + "/costmap_3d")),
    get_plan_cost_action_srv_(
        private_nh_,
        "get_plan_cost_3d",
        std::bind(&Costmap3DROS::getPlanCost3DActionCallback, this, std::placeholders::_1),
        false)
{
  {
    std::lock_guard<std::mutex> initialize_lock(initialized_mutex_);

    if (!private_nh_.getParam("costmap_3d/footprint_mesh_resource", footprint_mesh_resource_))
    {
      ROS_ERROR("Unable to find footprint_mesh_resource parameter");
    }

    XmlRpc::XmlRpcValue my_list;
    if (private_nh_.getParam("costmap_3d/plugins", my_list))
    {
      for (int32_t i = 0; i < my_list.size(); ++i)
      {
        std::string pname = static_cast<std::string>(my_list[i]["name"]);
        std::string type = static_cast<std::string>(my_list[i]["type"]);

        try
        {
          boost::shared_ptr<Layer3D> plugin = plugin_loader_.createInstance(type);
          ROS_INFO_STREAM("3D costmap \"" << name << "\" "
                          << "Using 3D plugin \"" << pname << "\" of type \"" << type << "\"");
          plugin->initialize(&layered_costmap_3d_, name + "/costmap_3d/" + pname, &tf_);
          layered_costmap_3d_.addPlugin(plugin);
        }
        catch (class_loader::ClassLoaderException e)
        {
          ROS_ERROR_STREAM("Unable to load plugin \"" << pname << "\" of type \"" << type << "\". "
                           << "Error: " << e.what());
        }
      }
    }
    if (layered_costmap_3d_.getPlugins().size() == 0)
    {
      ROS_WARN_STREAM("3D costmap \"" << name << "\" has no 3D plugin layers, will be 2D only");
    }
    footprint_pub_ = private_nh_.advertise<visualization_msgs::Marker>("footprint_3d", 1, true);
    initialized_ = true;
  }

  // Now that the object is initialized, advertise dynamic reconfigure and services
  dsrv_.setCallback(
      std::bind(
          &Costmap3DROS::reconfigureCB,
          this,
          std::placeholders::_1,
          std::placeholders::_2));

  get_plan_cost_action_srv_.start();
  get_plan_cost_srv_ = private_nh_.advertiseService(
      "get_plan_cost_3d",
      &Costmap3DROS::getPlanCost3DServiceCallback,
      this);
}

Costmap3DROS::~Costmap3DROS()
{
}

void Costmap3DROS::reconfigureCB(Costmap3DConfig &config, uint32_t level)
{
  geometry_msgs::Point min, max;

  layered_costmap_3d_.setBounds(config.map_z_min, config.map_z_max);

  footprint_3d_padding_ = config.footprint_3d_padding;

  publishFootprint();
}

void Costmap3DROS::updateMap()
{
  std::lock_guard<std::mutex> initialize_lock(initialized_mutex_);

  if (initialized_ && !isPaused())
  {
    // First update the 3D map.
    // We update 3D first, in case any 3D layers need to affect 2D layers.
    // Get global pose of robot
    tf::Stamped < tf::Pose > pose;
    if (getRobotPose (pose))
    {
      geometry_msgs::Pose pose_msg;
      tf::poseTFToMsg(pose, pose_msg);
      layered_costmap_3d_.updateMap(pose_msg);
    }

    // Now update the 2D map
    super::updateMap();
  }
}

void Costmap3DROS::start()
{
  // check if we're stopped or just paused
  if (isStopped())
  {
    layered_costmap_3d_.activate();
  }

  super::start();
}

void Costmap3DROS::stop()
{
  super::stop();
  layered_costmap_3d_.deactivate();
}

void Costmap3DROS::resetLayers()
{
  layered_costmap_3d_.reset();
  super::resetLayers();
}

std::shared_ptr<Costmap3DQuery> Costmap3DROS::getBufferedQuery(
        const std::string& footprint_mesh_resource,
        double padding)
{
  const std::string& query_mesh(getFootprintMeshResource(footprint_mesh_resource));
  padding = getFootprintPadding(padding);
  std::lock_guard<LayeredCostmap3D> lock(layered_costmap_3d_);
  return std::shared_ptr<Costmap3DQuery>(new Costmap3DQuery(layered_costmap_3d_.getCostmap3D(), query_mesh, padding));
}

std::shared_ptr<Costmap3DQuery> Costmap3DROS::getAssociatedQuery(
        const std::string& footprint_mesh_resource,
        double padding)
{
  const std::string& query_mesh(getFootprintMeshResource(footprint_mesh_resource));
  padding = getFootprintPadding(padding);
  std::lock_guard<LayeredCostmap3D> lock(layered_costmap_3d_);
  return std::shared_ptr<Costmap3DQuery>(new Costmap3DQuery(&layered_costmap_3d_, query_mesh, padding));
}

std::set<std::string> Costmap3DROS::getLayerNames()
{
  // Begin with any 2D layers
  std::set<std::string> plugin_names(super::getLayerNames());

  std::lock_guard<LayeredCostmap3D> lock(layered_costmap_3d_);

  // Add the 3D layers to the list
  for (auto plugin : layered_costmap_3d_.getPlugins())
  {
    plugin_names.insert(plugin->getName());
  }

  return plugin_names;
}

void Costmap3DROS::resetBoundingBox(geometry_msgs::Point min, geometry_msgs::Point max, const std::set<std::string>& layers)
{
  {
    // Its not OK to hold the lock while calling super::resetBoundingBox even
    // though the lock is recursive because it will break the condition variable
    // wait.
    std::lock_guard<LayeredCostmap3D> lock(layered_costmap_3d_);

    // Reset all 3D layers first, then 2D layers, because Costmap2DROS
    // resetBoundingBox waits for the next update
    for (auto plugin : layered_costmap_3d_.getPlugins())
    {
      // Only reset layers that are in the layer set
      // Match either the whole layer name, or just the final name after the
      // final '/'.
      const std::string& plugin_full_name(plugin->getName());
      std::string plugin_last_name_only;
      int slash = plugin_full_name.rfind('/');
      if( slash != std::string::npos )
      {
        plugin_last_name_only = plugin_full_name.substr(slash+1);
      }
      ROS_INFO_STREAM("resetBoundingBox consider 3D layer: " << plugin_full_name <<
                      " last name: " << plugin_last_name_only);
      if (layers.find(plugin_full_name) != layers.end()
          || plugin_last_name_only.size() > 0 && layers.find(plugin_last_name_only) != layers.end())
      {
        ROS_INFO_STREAM("resetBoundingBox 3D layer " << plugin->getName());
        plugin->resetBoundingBox(min, max);
      }
    }
  }

  super::resetBoundingBox(min, max, layers);
}

const std::string& Costmap3DROS::getFootprintMeshResource(const std::string& alt_mesh)
{
  // TODO: need to verify that the alt file actually exists
  if (alt_mesh.empty())
  {
    return footprint_mesh_resource_;
  }
  return alt_mesh;
}

double Costmap3DROS::getFootprintPadding(double alt_padding)
{
  if (!std::isfinite(alt_padding))
  {
    return footprint_3d_padding_;
  }
  return alt_padding;
}

void Costmap3DROS::publishFootprint()
{
  // TODO: publish a footprint with the padding applied
  visualization_msgs::Marker marker;
  marker.header.frame_id = getBaseFrameID();
  marker.header.stamp = ros::Time();
  marker.type = visualization_msgs::Marker::MESH_RESOURCE;
  marker.action = visualization_msgs::Marker::ADD;
  marker.pose.position.x = 0.0;
  marker.pose.position.y = 0.0;
  marker.pose.position.z = 0.0;
  marker.pose.orientation.x = 0.0;
  marker.pose.orientation.y = 0.0;
  marker.pose.orientation.z = 0.0;
  marker.pose.orientation.w = 1.0;
  marker.scale.x = 1;
  marker.scale.y = 1;
  marker.scale.z = 1;
  marker.color.a = 0.5;
  marker.color.r = 0.5;
  marker.color.g = 0.5;
  marker.color.b = 0.5;
  marker.frame_locked = true;
  marker.mesh_resource = footprint_mesh_resource_;
  footprint_pub_.publish(marker);
}

std::shared_ptr<Costmap3DQuery> Costmap3DROS::getQuery(const std::string& footprint_mesh_resource, double padding)
{
  const std::string& query_mesh(getFootprintMeshResource(footprint_mesh_resource));
  padding = getFootprintPadding(padding);
  auto query_pair = std::make_pair(query_mesh, padding);
  auto query_it = query_map_.find(query_pair);
  if (query_it == query_map_.end())
  {
    // Query object does not exist, create it and add it to the map
    std::shared_ptr<Costmap3DQuery> query;
    query.reset(new Costmap3DQuery(&layered_costmap_3d_, query_mesh, padding));
    query_it = query_map_.insert(std::make_pair(query_pair, query)).first;
  }
  return query_it->second;
}

double Costmap3DROS::footprintCost(geometry_msgs::Pose pose,
                                   const std::string& footprint_mesh_resource,
                                   double padding)
{
  return getQuery(footprint_mesh_resource, padding)->footprintCost(pose);
}

bool Costmap3DROS::footprintCollision(geometry_msgs::Pose pose,
                                      const std::string& footprint_mesh_resource,
                                      double padding)
{
  return getQuery(footprint_mesh_resource, padding)->footprintCollision(pose);
}

double Costmap3DROS::footprintDistance(geometry_msgs::Pose pose,
                                       const std::string& footprint_mesh_resource,
                                       double padding)
{
  return getQuery(footprint_mesh_resource, padding)->footprintDistance(pose);
}

double Costmap3DROS::footprintSignedDistance(geometry_msgs::Pose pose,
                                             const std::string& footprint_mesh_resource,
                                             double padding)
{
  return getQuery(footprint_mesh_resource, padding)->footprintSignedDistance(pose);
}

void Costmap3DROS::getPlanCost3DActionCallback(
    const actionlib::SimpleActionServer<GetPlanCost3DAction>::GoalConstPtr& goal)
{
  GetPlanCost3DResult result;
  processPlanCost3D(*goal, result);
  get_plan_cost_action_srv_.setSucceeded(result);
}

bool Costmap3DROS::getPlanCost3DServiceCallback(
    GetPlanCost3DService::Request& request,
    GetPlanCost3DService::Response& response)
{
  processPlanCost3D(request, response);
  return true;
}

template <typename RequestType, typename ResponseType>
void Costmap3DROS::processPlanCost3D(RequestType& request, ResponseType& response)
{
  // Be sure the costmap is locked while querying, if necessary.
  std::unique_lock <LayeredCostmap3D> lock(layered_costmap_3d_, std::defer_lock);
  std::shared_ptr<Costmap3DQuery> query;

#if (COSTMAP_3D_ROS_AUTO_PROFILE_QUERY) > 0
  ROS_INFO("Starting getPlanCost3DServiceCallback");
  ros::Time start_time = ros::Time::now();
  kill(getpid(), 12);
#endif

  if (request.buffered)
  {
    // Buffered request, getBufferedQuery gets the lock directly
    query = getBufferedQuery(request.footprint_mesh_resource, request.padding);
  }
  else
  {
    // Unbuffered request, must hold the lock for the duration.
    lock.lock();
    query = getQuery(request.footprint_mesh_resource, request.padding);
  }

  if (request.cost_query_mode == GetPlanCost3DService::Request::COST_QUERY_MODE_DISTANCE ||
      request.cost_query_mode == GetPlanCost3DService::Request::COST_QUERY_MODE_EXACT_SIGNED_DISTANCE)
  {
    response.cost = std::numeric_limits<double>::max();
  }
  else
  {
    response.cost = 0.0;
  }

  response.pose_costs.reserve(request.poses.size());
  response.lethal_indices.reserve(request.poses.size());
  for (int i = 0; i < request.poses.size(); ++i)
  {
    bool collision = false;
    // NOTE: Costmap3D does not support time (its not a 4D costmap!) so ignore the stamp.
    const auto pose = request.poses[i].pose;
    double pose_cost = -1.0;
    // Warn if the frame doesn't match. We currently don't transform the poses.
    if (request.poses[i].header.frame_id != layered_costmap_3d_.getGlobalFrameID())
    {
      ROS_WARN_STREAM_THROTTLE(1.0, "Costmap3DROS::getPlanCost3DServiceCallback expects poses in frame " <<
                               layered_costmap_3d_.getGlobalFrameID() << " but pose was in frame " <<
                               request.poses[i].header.frame_id);
      // Leave the pose_cost lethal
    }
    else if (request.cost_query_mode == GetPlanCost3DService::Request::COST_QUERY_MODE_COLLISION_ONLY)
    {
      pose_cost = query->footprintCollision(pose) ? -1.0 : 0.0;
    }
    else if (request.cost_query_mode == GetPlanCost3DService::Request::COST_QUERY_MODE_DISTANCE)
    {
      pose_cost = query->footprintDistance(pose);
    }
    else if (request.cost_query_mode == GetPlanCost3DService::Request::COST_QUERY_MODE_EXACT_SIGNED_DISTANCE)
    {
      pose_cost = query->footprintSignedDistance(pose);
    }
    else
    {
      pose_cost = query->footprintCost(pose);
    }
    // negative is a collision
    if (pose_cost < 0.0)
    {
      collision = true;
    }
    if (request.cost_query_mode == GetPlanCost3DService::Request::COST_QUERY_MODE_DISTANCE ||
        request.cost_query_mode == GetPlanCost3DService::Request::COST_QUERY_MODE_EXACT_SIGNED_DISTANCE)
    {
      // in distance mode, the cost is the minimum distance across all poses
      response.cost = std::min(response.cost, pose_cost);
    }
    else
    {
      // in collision or cost mode, plan cost will either be negative if there is a collision
      // or the aggregate non-lethal cost.
      if (collision)
      {
        if (response.cost >= 0.0)
        {
          // this pose is in collision, but the plan hasn't seen a collision yet.
          // reset the cost to this pose's cost
          response.cost = pose_cost;
        }
        else
        {
          // otherwise, add the lethal pose_cost to the cost (making it more negative)
          response.cost += pose_cost;
        }
      }
      else
      {
        if (response.cost >= 0.0)
        {
          // not in collision, plan not in collision, add the cost
          response.cost += pose_cost;
        }
        // don't add the non-lethal cost to a lethal plan
      }
    }
    response.pose_costs.push_back(pose_cost);
    if (collision)
    {
      response.lethal_indices.push_back(i);
      if (request.lazy)
      {
        break;
      }
    }
  }

#if (COSTMAP_3D_ROS_AUTO_PROFILE_QUERY) > 0
  kill(getpid(), 12);
  ROS_INFO_STREAM("Finished getting " << request.poses.size() << " poses in " <<
                  (ros::Time::now() - start_time).toSec() << " seconds.");
#endif
}

}  // namespace costmap_3d