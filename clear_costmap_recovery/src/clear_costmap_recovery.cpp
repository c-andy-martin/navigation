/*********************************************************************
*
* Software License Agreement (BSD License)
*
*  Copyright (c) 2009, Willow Garage, Inc.
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
* Author: Eitan Marder-Eppstein
*********************************************************************/
#include <clear_costmap_recovery/clear_costmap_recovery.h>
#include <pluginlib/class_list_macros.h>
#include <boost/pointer_cast.hpp>
#include <vector>

//register this planner as a RecoveryBehavior plugin
PLUGINLIB_EXPORT_CLASS(clear_costmap_recovery::ClearCostmapRecovery, nav_core::RecoveryBehavior)

using costmap_2d::NO_INFORMATION;

namespace clear_costmap_recovery {
ClearCostmapRecovery::ClearCostmapRecovery(): global_costmap_(NULL), local_costmap_(NULL),
  tf_(NULL), initialized_(false) {}

void ClearCostmapRecovery::initialize(std::string name, tf2_ros::Buffer* tf,
    costmap_2d::Costmap2DROS* global_costmap, costmap_2d::Costmap2DROS* local_costmap){
  if(!initialized_){
    name_ = name;
    tf_ = tf;
    global_costmap_ = global_costmap;
    local_costmap_ = local_costmap;

    //get some parameters from the parameter server
    ros::NodeHandle private_nh("~/" + name_);

    private_nh.param("reset_distance", reset_distance_, 3.0);
    private_nh.param("invert_area_to_clear", invert_area_to_clear_, false);
    private_nh.param("force_updating", force_updating_, false);
    private_nh.param("affected_maps", affected_maps_, std::string("both"));
    if (affected_maps_ != "local" && affected_maps_ != "global" && affected_maps_ != "both")
    {
      ROS_WARN("Wrong value for affected_maps parameter: '%s'; valid values are 'local', 'global' or 'both'; " \
               "defaulting to 'both'", affected_maps_.c_str());
      affected_maps_ = "both";
    }

    std::vector<std::string> clearable_layers_default, clearable_layers;
    clearable_layers_default.push_back( std::string("obstacles") );
    private_nh.param("layer_names", clearable_layers, clearable_layers_default);

    for(unsigned i=0; i < clearable_layers.size(); i++) {
        ROS_INFO("Recovery behavior will clear layer '%s'", clearable_layers[i].c_str());
        clearable_layers_.insert(clearable_layers[i]);
    }

    initialized_ = true;
  }
  else{
    ROS_ERROR("You should not call initialize twice on this object, doing nothing");
  }
}

void ClearCostmapRecovery::runBehavior(){
  if(!initialized_){
    ROS_ERROR("This object must be initialized before runBehavior is called");
    return;
  }

  if(global_costmap_ == NULL || local_costmap_ == NULL){
    ROS_ERROR("The costmaps passed to the ClearCostmapRecovery object cannot be NULL. Doing nothing.");
    return;
  }

  if (!invert_area_to_clear_){
    ROS_WARN("Clearing %s costmap%s outside a square (%.2fm) large centered on the robot.", affected_maps_.c_str(),
           affected_maps_ == "both" ? "s" : "", reset_distance_);
  }else {
    ROS_WARN("Clearing %s costmap%s inside a square (%.2fm) large centered on the robot.", affected_maps_.c_str(),
           affected_maps_ == "both" ? "s" : "", reset_distance_);
  }

  ros::WallTime t0 = ros::WallTime::now();
  if (affected_maps_ == "global" || affected_maps_ == "both")
  {
    clear(global_costmap_);

    if (force_updating_)
      global_costmap_->updateMap();

    ROS_DEBUG("Global costmap cleared in %fs", (ros::WallTime::now() - t0).toSec());
  }

  t0 = ros::WallTime::now();
  if (affected_maps_ == "local" || affected_maps_ == "both")
  {
    clear(local_costmap_);

    if (force_updating_)
      local_costmap_->updateMap();

    ROS_DEBUG("Local costmap cleared in %fs", (ros::WallTime::now() - t0).toSec());
  }
}

void ClearCostmapRecovery::clear(costmap_2d::Costmap2DROS* costmap){
  std::vector<boost::shared_ptr<costmap_2d::Layer> >* plugins = costmap->getLayeredCostmap()->getPlugins();

  geometry_msgs::PoseStamped pose;

  if(!costmap->getRobotPose(pose)){
    ROS_ERROR("Cannot clear map because pose cannot be retrieved");
    return;
  }

  double x = pose.pose.position.x;
  double y = pose.pose.position.y;

  geometry_msgs::Point min, max;
  min.x = x - reset_distance_ / 2;
  min.y = y - reset_distance_ / 2;
  min.z = -std::numeric_limits<double>::max();
  max.x = x + reset_distance_ / 2;
  max.y = y + reset_distance_ / 2;
  max.z = std::numeric_limits<double>::max();

  costmap->resetBoundingBox(min, max, clearable_layers_);
}

};
