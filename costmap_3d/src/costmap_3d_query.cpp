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
#include <costmap_3d/costmap_3d_query.h>
#include <chrono>
#include <sstream>
#include <fcl/geometry/octree/octree.h>
#include <fcl/narrowphase/collision.h>
#include <fcl/narrowphase/distance_result.h>
#include <fcl/geometry/shape/sphere.h>
#include <ros/ros.h>
#include <ros/package.h>
#include <tf/transform_datatypes.h>
#include <costmap_3d/octree_solver.h>

// In some environments, pcl::io isn't available because it depends on VTK,
// with its enormous dependencies, so allow for an alternative implementation:
#if(USE_STL_LOADER)
  #include <stl_loader/stl_loader.h>
#else
  #include <pcl/io/vtk_lib_io.h>
#endif

namespace costmap_3d
{

// Define and initialize static thread_local variables.
// statics.
thread_local unsigned int Costmap3DQuery::tls_last_layered_costmap_update_number_ = 0;
thread_local Costmap3DQuery* Costmap3DQuery::tls_last_instance_ = nullptr;
thread_local Costmap3DQuery::DistanceCacheEntry Costmap3DQuery::tls_last_cache_entries_[Costmap3DQuery::MAX][Costmap3DQuery::OBSTACLES_MAX] = {};

Costmap3DQuery::Costmap3DQuery(
    const LayeredCostmap3D* layered_costmap_3d,
    const std::string& mesh_resource,
    double padding,
    Costmap3DQuery::QueryRegionScale query_region_scale,
    unsigned int pose_bins_per_meter,
    unsigned int pose_bins_per_radian,
    unsigned int pose_milli_bins_per_meter,
    unsigned int pose_milli_bins_per_radian,
    unsigned int pose_micro_bins_per_meter,
    unsigned int pose_micro_bins_per_radian)
  : layered_costmap_3d_(layered_costmap_3d),
    query_region_scale_(query_region_scale),
    pose_bins_per_meter_(pose_bins_per_meter),
    pose_bins_per_radian_(pose_bins_per_radian),
    pose_milli_bins_per_meter_(pose_milli_bins_per_meter),
    pose_milli_bins_per_radian_(pose_milli_bins_per_radian),
    pose_micro_bins_per_meter_(pose_micro_bins_per_meter),
    pose_micro_bins_per_radian_(pose_micro_bins_per_radian)
{
  init();
  updateMeshResource(mesh_resource, padding);
}

Costmap3DQuery::Costmap3DQuery(const Costmap3DConstPtr& costmap_3d,
    const std::string& mesh_resource,
    double padding,
    Costmap3DQuery::QueryRegionScale query_region_scale,
    unsigned int pose_bins_per_meter,
    unsigned int pose_bins_per_radian,
    unsigned int pose_milli_bins_per_meter,
    unsigned int pose_milli_bins_per_radian,
    unsigned int pose_micro_bins_per_meter,
    unsigned int pose_micro_bins_per_radian)
  : layered_costmap_3d_(nullptr),
    query_region_scale_(query_region_scale),
    pose_bins_per_meter_(pose_bins_per_meter),
    pose_bins_per_radian_(pose_bins_per_radian),
    pose_milli_bins_per_meter_(pose_milli_bins_per_meter),
    pose_milli_bins_per_radian_(pose_milli_bins_per_radian),
    pose_micro_bins_per_meter_(pose_micro_bins_per_meter),
    pose_micro_bins_per_radian_(pose_micro_bins_per_radian)
{
  init();
  // Make a local copy of the costmap in question
  // It would be awesome if the Costmap3D had a way to snapshot
  // or copy-on-write. As it stands, for many scenarios involving
  // space-limited local costmaps, copying a 3D costmap will only take a
  // couple milliseconds and is better than leaving the costmap locked for
  // an entire planning cycle.
  octree_ptr_.reset(new Costmap3D(*costmap_3d));
  updateMeshResource(mesh_resource, padding);
  // For a buffered query, go ahead and setup the interior collision LUT now,
  // as the resolution can not change.
  checkInteriorCollisionLUT();
}

Costmap3DQuery::~Costmap3DQuery()
{
}

void Costmap3DQuery::init()
{
  clearStatistics();
  robot_model_halfspaces_.reset(new std::vector<fcl::Halfspace<FCLFloat>>);
  last_octomap_resolution_ = 0.0;
}

void Costmap3DQuery::setCacheBinSize(
      unsigned int pose_bins_per_meter,
      unsigned int pose_bins_per_radian,
      unsigned int pose_milli_bins_per_meter,
      unsigned int pose_milli_bins_per_radian,
      unsigned int pose_micro_bins_per_meter,
      unsigned int pose_micro_bins_per_radian)
{
  upgrade_lock upgradeable_lock(upgrade_mutex_);
  if (pose_bins_per_meter != pose_bins_per_meter_ ||
      pose_bins_per_radian != pose_bins_per_radian_ ||
      pose_milli_bins_per_meter != pose_milli_bins_per_meter_ ||
      pose_milli_bins_per_radian != pose_milli_bins_per_radian_ ||
      pose_micro_bins_per_meter != pose_micro_bins_per_meter_ ||
      pose_micro_bins_per_radian != pose_micro_bins_per_radian_)
  {
    upgrade_to_unique_lock write_lock(upgradeable_lock);
    // There is no need to invalidate the existing cache entries. If the bins
    // per meter change, either some or all of the new calls will miss the
    // existing cache entries. The existing entries are all still correct, so
    // keeping them is fine. Also, if only the thresholds change all the cache
    // entries will still work.
    pose_bins_per_meter_ = pose_bins_per_meter;
    pose_bins_per_radian_ = pose_bins_per_radian;
    pose_milli_bins_per_meter_ = pose_milli_bins_per_meter;
    pose_milli_bins_per_radian_ = pose_milli_bins_per_radian;
    pose_micro_bins_per_meter_ = pose_micro_bins_per_meter;
    pose_micro_bins_per_radian_ = pose_micro_bins_per_radian;
    ROS_INFO_STREAM("Set cache bins per meter: " << pose_bins_per_meter_);
    ROS_INFO_STREAM("Set cache bins per radian: " << pose_bins_per_radian_);
    ROS_INFO_STREAM("Set milli-cache bins per meter: " << pose_milli_bins_per_meter_);
    ROS_INFO_STREAM("Set milli-cache bins per radian: " << pose_milli_bins_per_radian_);
    ROS_INFO_STREAM("Set micro-cache bins per meter: " << pose_micro_bins_per_meter_);
    ROS_INFO_STREAM("Set micro-cache bins per radian: " << pose_micro_bins_per_radian_);
    calculateCacheDistanceThresholds();
  }
}

void Costmap3DQuery::setCacheThresholdParameters(
    bool threshold_two_d_mode,
    double threshold_factor)
{
  upgrade_lock upgradeable_lock(upgrade_mutex_);
  if (threshold_two_d_mode != threshold_two_d_mode_ ||
      threshold_factor != threshold_factor_)
  {
    upgrade_to_unique_lock write_lock(upgradeable_lock);
    threshold_two_d_mode_ = threshold_two_d_mode;
    threshold_factor_ = threshold_factor;
    ROS_INFO_STREAM("Set threshold 2D mode: " << threshold_two_d_mode_ ? "true" : "false");
    ROS_INFO_STREAM("Set threshold factor: " << threshold_factor_);
    calculateCacheDistanceThresholds();
  }
}

void Costmap3DQuery::calculateCacheDistanceThresholds()
{
  // The distance thresholds must be greater than the maximum translation error
  // plus the maximum rotation error, otherwise the query may return
  // no-collision when one is actually present. In a robot with full
  // three-dimensional movement this would equate to the length of the diagonal
  // of a bin plus the chord length corresponding to the size of a rotational
  // bin given the maximum radius of the robot mesh. For two-dimensional
  // movement of a robot on a plane aligned with the costmap, the threshold
  // is lower at only the diagonal of the bin projected into the plane as
  // a square and the chord length given the radius of the robot mesh projected
  // into the plane (the robot footprint radius).
  Eigen::Vector4f origin(0.0, 0.0, 0.0, 0.0), max_pt(0.0, 0.0, 0.0, 0.0);
  double diagonal_factor;
  if (threshold_two_d_mode_)
  {
    // Create a 2D version of the mesh cloud by truncating the Z value.
    pcl::PointCloud<pcl::PointXYZ> mesh_points_2d(*robot_mesh_points_);
    for (auto& point : mesh_points_2d)
    {
      point.z = 0.0;
    }
    pcl::getMaxDistance(mesh_points_2d, origin, max_pt);
    diagonal_factor = std::sqrt(2.0);
  }
  else
  {
    pcl::getMaxDistance(*robot_mesh_points_, origin, max_pt);
    diagonal_factor = std::sqrt(3.0);
  }
  const Eigen::Vector3f max_pt3 = max_pt.head<3>();
  double mesh_radius = max_pt3.norm();
  milli_cache_threshold_ = diagonal_factor / pose_milli_bins_per_meter_ +
    2.0 * mesh_radius * std::sin(0.5 / pose_milli_bins_per_radian_);
  micro_cache_threshold_ = diagonal_factor / pose_micro_bins_per_meter_ +
    2.0 * mesh_radius * std::sin(0.5 / pose_micro_bins_per_radian_);
  milli_cache_threshold_ *= threshold_factor_;
  micro_cache_threshold_ *= threshold_factor_;
  ROS_INFO_STREAM("Calculated mesh radius: " << mesh_radius << "m");
  ROS_INFO_STREAM("Calculated milli-distance cache threshold: " << milli_cache_threshold_ << "m");
  ROS_INFO_STREAM("Calculated micro-distance cache threshold: " << micro_cache_threshold_ << "m");
}

void Costmap3DQuery::updateCostmap(const Costmap3DConstPtr& costmap_3d)
{
  // Copy the given costmap
  std::shared_ptr<const octomap::OcTree> new_octree(new Costmap3D(*costmap_3d));
  upgrade_lock upgrade_lock(upgrade_mutex_);
  checkCostmap(upgrade_lock, new_octree);
}

// Caller must hold at least a shared lock
bool Costmap3DQuery::needCheckCostmap(std::shared_ptr<const octomap::OcTree> new_octree)
{
  bool need_update = false;
  {
    if (layered_costmap_3d_)
    {
      if (layered_costmap_3d_->getCostmap3D() != octree_ptr_ ||
         (last_octomap_resolution_ != layered_costmap_3d_->getCostmap3D()->getResolution()) ||
         (last_layered_costmap_update_number_ != layered_costmap_3d_->getNumberOfUpdates()))
      {
        need_update = true;
      }
    }
    else if (new_octree)
    {
      if (new_octree != octree_ptr_)
      {
        need_update = true;
      }
    }
  }
  return need_update;
}

// Borrowed from https://artificial-mind.net/blog/2021/10/09/unordered-map-badness
template <class Map>
static double unordered_map_badness(Map const& map)
{
  auto const lambda = map.size() / double(map.bucket_count());

  auto cost = 0.;
  for (auto const& entry : map)
    cost += map.bucket_size(map.bucket(entry.first));
  cost /= map.size();

  return std::max(0., cost / (1 + lambda) - 1);
}

template <class Map>
static void printDistanceCacheDebug(
    const Map& cache,
    const std::string& prefix)
{
  ROS_DEBUG_STREAM_NAMED(
      "query_distance_cache", prefix <<
      " size: " << cache.size() <<
      " bucket count: " << cache.bucket_count() <<
      " badness: " << unordered_map_badness(cache));
}

// Caller must already hold an upgradable shared lock via the passed
// upgrade_lock, which we can use to upgrade to exclusive access if necessary.
void Costmap3DQuery::checkCostmap(Costmap3DQuery::upgrade_lock& upgrade_lock,
                                  std::shared_ptr<const octomap::OcTree> new_octree)
{
  // Check if we need to update w/ just the upgrade lock held. This way only
  // one thread has to pay to upgrade to an exclusive lock (which forces no new
  // readers and waits for all existing readers to finish).
  bool need_update = needCheckCostmap(new_octree);
  if (need_update)
  {
    // get write access
    upgrade_to_unique_lock write_lock(upgrade_lock);
    // While it isn't possible for the cache locks to be held once the overall
    // upgrade mutex is locked (since they are only ever grabbed while the
    // shared read lock is held), go ahead and lock them too for completeness.
    unique_lock distance_cache_lock(distance_cache_mutex_);
    unique_lock milli_distance_cache_lock(milli_distance_cache_mutex_);
    unique_lock micro_distance_cache_lock(micro_distance_cache_mutex_);
    unique_lock exact_distance_cache_lock(exact_distance_cache_mutex_);

    if (layered_costmap_3d_)
    {
      if (layered_costmap_3d_->getCostmap3D() != octree_ptr_)
      {
        // The octomap has been reallocated, change where we are pointing.
        octree_ptr_ = layered_costmap_3d_->getCostmap3D();
        nonlethal_octree_ptr_.reset();
      }
    }
    else if (new_octree)
    {
      if (new_octree != octree_ptr_)
      {
        // There is a new octomap allocated, change where we are pointing.
        octree_ptr_ = new_octree;
        nonlethal_octree_ptr_.reset();
      }
    }
    checkInteriorCollisionLUT();
    // The costmap has been updated since the last query, reset our caches

    // Delete any distance cache entries that have had their corresponding
    // octomap cells removed. It is fine to keep entries in the presence of
    // additions, as the entry only defines an upper bound. Because the size of
    // the tree is limited, the size of the cache has a natural limit. If this
    // limit is ever too large, a separate cache size may need to be set.
    {
      printDistanceCacheDebug(distance_cache_, "distance cache: ");
      size_t start_size = distance_cache_.size();
      auto it = distance_cache_.begin();
      while (it != distance_cache_.end())
      {
        bool erase = true;
        Costmap3DIndex index;
        unsigned int depth;
        if (it->second.getCostmapIndexAndDepth(*octree_ptr_, &index, &depth))
        {
          auto* node = octree_ptr_->search(index, depth);
          if (node && !octree_ptr_->nodeHasChildren(node))
          {
            // The node exists and has no children (so its a leaf and not an inner node)
            // Keep this entry.
            erase = false;
          }
        }

        if (erase)
          it = distance_cache_.erase(it);
        else
          ++it;
      }
      ROS_DEBUG_STREAM_NAMED("query_distance_cache",
          "distance cache removed " << start_size - distance_cache_.size() << " entries");
    }

    // Delete any milli distance cache entries that have had their
    // corresponding octomap cells removed. If the cell has not been removed,
    // invalidate the milli cache entry from being used in the fast-path, but
    // let it still be used as an upper-bound for the slow path, as it greatly
    // reduces subsequent query times. Also note that the caches are always
    // updated after a successful query, so the staleness of the bound will be
    // fixed the next query that hits the cache.
    {
      printDistanceCacheDebug(milli_distance_cache_, "milli-distance cache: ");
      size_t start_size = milli_distance_cache_.size();
      auto it = milli_distance_cache_.begin();
      while (it != milli_distance_cache_.end())
      {
        bool erase = true;
        Costmap3DIndex index;
        unsigned int depth;
        if (it->second.getCostmapIndexAndDepth(*octree_ptr_, &index, &depth))
        {
          auto* node = octree_ptr_->search(index, depth);
          if (node && !octree_ptr_->nodeHasChildren(node))
          {
            // The node exists and has no children (so its a leaf and not an inner node)
            // Keep this entry.
            erase = false;
          }
        }

        if (erase)
        {
          it = milli_distance_cache_.erase(it);
        }
        else
        {
          // Fast-path only used for finite recorded distance
          it->second.distance = std::numeric_limits<FCLFloat>::infinity();
          ++it;
        }
      }
      ROS_DEBUG_STREAM_NAMED("query_distance_cache",
          "milli-distance cache removed " << start_size - milli_distance_cache_.size() << " entries");
    }

    // Drop the micro cache. It is not worth iterating over the (relatively
    // large) micro cache.  New entries in the costmap make the upper bound
    // fairly worthless and we can not use the fast path which is the micro
    // cache's strong point.
    printDistanceCacheDebug(micro_distance_cache_, "micro-distance cache: ");
    micro_distance_cache_.clear();
    // We must drop the exact cache after the costmap has changed
    printDistanceCacheDebug(exact_distance_cache_, "exact distance cache: ");
    exact_distance_cache_.clear();
    printStatistics();
    clearStatistics();
    if (layered_costmap_3d_)
    {
      last_layered_costmap_update_number_ = layered_costmap_3d_->getNumberOfUpdates();
    }
  }
  // Must re-check TLS as the layred costmap update number has changed.
  checkTLS();
}

void Costmap3DQuery::checkTLS()
{
  // Check if any thread local state must be invalidated.
  // We must handle if either the costmap has been updated since this thread's
  // storage was invalidated, or if a different instance of the object is
  // being used in the same thread.
  if (tls_last_layered_costmap_update_number_ != last_layered_costmap_update_number_ ||
      tls_last_instance_ != this)
  {
    for (unsigned int i=0; i<MAX; ++i)
    {
      for (unsigned int j=0; j<OBSTACLES_MAX; ++j)
      {
        tls_last_cache_entries_[i][j].clear();
      }
    }
    tls_last_layered_costmap_update_number_ = last_layered_costmap_update_number_;
    tls_last_instance_ = this;
  }
}

void Costmap3DQuery::clearStatistics()
{
  queries_since_clear_ = 0;
  empties_since_clear_ = 0;
  hits_since_clear_ = 0;
  fast_milli_hits_since_clear_ = 0;
  slow_milli_hits_since_clear_ = 0;
  fast_micro_hits_since_clear_ = 0;
  slow_micro_hits_since_clear_ = 0;
  exact_hits_since_clear_ = 0;
  reuse_results_since_clear_ = 0;
  misses_since_clear_us_ = 0;
  hits_since_clear_us_ = 0;
  fast_milli_hits_since_clear_us_ = 0;
  slow_milli_hits_since_clear_us_ = 0;
  fast_micro_hits_since_clear_us_ = 0;
  slow_micro_hits_since_clear_us_ = 0;
  exact_hits_since_clear_us_ = 0;
  reuse_results_since_clear_us_ = 0;
  hit_fcl_bv_distance_calculations_ = 0;
  hit_fcl_primative_distance_calculations_ = 0;
  miss_fcl_bv_distance_calculations_ = 0;
  miss_fcl_primative_distance_calculations_ = 0;
}

void Costmap3DQuery::printStatistics()
{
  unsigned int cache_misses = queries_since_clear_ -
      empties_since_clear_ -
      hits_since_clear_ -
      fast_milli_hits_since_clear_ -
      slow_milli_hits_since_clear_ -
      fast_micro_hits_since_clear_ -
      slow_micro_hits_since_clear_ -
      reuse_results_since_clear_ -
      exact_hits_since_clear_;
  double hit_ratio = (double)hits_since_clear_ / queries_since_clear_;
  double empty_ratio = (double)empties_since_clear_ / queries_since_clear_;
  double fast_milli_hit_ratio = (double)fast_milli_hits_since_clear_ / queries_since_clear_;
  double slow_milli_hit_ratio = (double)slow_milli_hits_since_clear_ / queries_since_clear_;
  double fast_micro_hit_ratio = (double)fast_micro_hits_since_clear_ / queries_since_clear_;
  double slow_micro_hit_ratio = (double)slow_micro_hits_since_clear_ / queries_since_clear_;
  double exact_hit_ratio = (double)exact_hits_since_clear_ / queries_since_clear_;
  uint64_t total_us = misses_since_clear_us_ +
      hits_since_clear_us_ +
      fast_milli_hits_since_clear_us_ +
      slow_milli_hits_since_clear_us_ +
      fast_micro_hits_since_clear_us_ +
      slow_micro_hits_since_clear_us_ +
      reuse_results_since_clear_us_ +
      exact_hits_since_clear_us_;
  std::ostringstream ss;
  ss << "Costmap3DQuery statistics:"
      "\n\tqueries this cycle: " << queries_since_clear_ <<
      "\n\tcache misses: " << cache_misses <<
      "\n\tcache hits: " << hits_since_clear_ <<
      "\n\tempty query ratio: " << empty_ratio <<
      "\n\tcache hit ratio: " << hit_ratio <<
      "\n\tslow milli cache hits: " << slow_milli_hits_since_clear_ <<
      "\n\tslow milli cache hit ratio: " << slow_milli_hit_ratio <<
      "\n\tfast milli cache hits: " << fast_milli_hits_since_clear_ <<
      "\n\tfast milli cache hit ratio: " << fast_milli_hit_ratio <<
      "\n\tslow micro cache hits: " << slow_micro_hits_since_clear_ <<
      "\n\tslow micro cache hit ratio: " << slow_micro_hit_ratio <<
      "\n\tfast micro cache hits: " << fast_micro_hits_since_clear_ <<
      "\n\tfast micro cache hit ratio: " << fast_micro_hit_ratio <<
      "\n\texact cache hits: " << exact_hits_since_clear_ <<
      "\n\texact cache hit ratio: " << exact_hit_ratio <<
      "\n\treuse past results: " << reuse_results_since_clear_ <<
      "\n\ttotal usecs: " << total_us <<
      "\n\tmiss usecs/query: " << (double)misses_since_clear_us_ / cache_misses <<
      "\n\thit usecs/query: " << (double)hits_since_clear_us_ / hits_since_clear_ <<
      "\n\tslow milli hit usecs/query: " << (double)slow_milli_hits_since_clear_us_ / slow_milli_hits_since_clear_ <<
      "\n\tfast milli hit usecs/query: " << (double)fast_milli_hits_since_clear_us_ / fast_milli_hits_since_clear_ <<
      "\n\tslow micro hit usecs/query: " << (double)slow_micro_hits_since_clear_us_ / slow_micro_hits_since_clear_ <<
      "\n\tfast micro hit usecs/query: " << (double)fast_micro_hits_since_clear_us_ / fast_micro_hits_since_clear_ <<
      "\n\texact hit usecs/query: " << (double)exact_hits_since_clear_us_ / exact_hits_since_clear_ <<
      "\n\treuse results usecs/query: " << (double)reuse_results_since_clear_us_ / reuse_results_since_clear_ <<
      "\n\tmiss FCL BV distance calculations: " << miss_fcl_bv_distance_calculations_ <<
      "\n\tmiss FCL primative distance calculations: " << miss_fcl_primative_distance_calculations_ <<
      "\n\thit FCL BV distance calculations: " << hit_fcl_bv_distance_calculations_ <<
      "\n\thit FCL primative distance calculations: " << hit_fcl_primative_distance_calculations_;
  // Leverage the fact that the print arguments are only run if the log is
  // enabled to make it simple to track if statistics tracking should be
  // enabled. There is no need to pay for the atomic operations for statistics
  // tracking if it is disabled.
  bool should_track_statistics = false;
  ROS_DEBUG_STREAM_NAMED(
      "query_statistics",
      // Weird construction here to get the side effect of setting
      // should_track_statistics when the named debug is enabled. Also skips
      // printing the first cycle as it will have invalid (zero) statistics
      (((should_track_statistics = true) && track_statistics_) ? ss.str() :
      "Costmap3DQuery: begin tracking statistics"));
  track_statistics_ = should_track_statistics;
}

void Costmap3DQuery::addPCLPolygonToFCLTriangles(
    const pcl::Vertices& polygon,
    std::vector<fcl::Triangle>* fcl_triangles)
{
  // Assume the polygons are convex. Break them into triangles.
  const std::size_t zero_index = polygon.vertices[0];
  for (int i=1; i < polygon.vertices.size() - 1; ++i)
  {
    fcl_triangles->push_back(fcl::Triangle(zero_index, polygon.vertices[i], polygon.vertices[i+1]));
  }
}


void Costmap3DQuery::addPCLPolygonMeshToRobotModel(
    const pcl::PolygonMesh& pcl_mesh,
    double padding,
    FCLRobotModel* robot_model)
{
  robot_mesh_points_.reset(new pcl::PointCloud<pcl::PointXYZ>);
  pcl::fromPCLPointCloud2(pcl_mesh.cloud, *robot_mesh_points_);

  padPoints(robot_mesh_points_, padding);

  std::vector<fcl::Vector3<FCLFloat>> fcl_points;
  std::vector<fcl::Triangle> fcl_triangles;

  for (auto pcl_point : *robot_mesh_points_)
  {
    fcl_points.push_back(convertPCLPointToFCL<FCLFloat>(pcl_point));
  }

  for (auto polygon : pcl_mesh.polygons)
  {
    addPCLPolygonToFCLTriangles(polygon, &fcl_triangles);
  }

  robot_model->addSubModel(fcl_points, fcl_triangles);
}

void Costmap3DQuery::updateMeshResource(const std::string& mesh_resource, double padding)
{
  unique_lock write_lock(upgrade_mutex_);
  std::string filename = getFileNameFromPackageURL(mesh_resource);
  if (filename.size() == 0)
  {
    return;
  }
  int pcl_rv = pcl::io::loadPolygonFileSTL(filename, robot_mesh_);
  if (pcl_rv < 0)
  {
    ROS_ERROR_STREAM("Costmap3DQuery: unable to load STL mesh file " << filename
                     << " query object will always return collision!");
    return;
  }
  robot_model_.reset(new FCLRobotModel());
  robot_model_->beginModel();
  addPCLPolygonMeshToRobotModel(robot_mesh_, padding, robot_model_.get());
  robot_model_->endModel();
  robot_model_->computeLocalAABB();

  crop_hull_.setHullCloud(robot_mesh_points_);
  crop_hull_.setHullIndices(robot_mesh_.polygons);
  crop_hull_.setCropOutside(true);

  // Calculate halfspaces for each mesh triangle
  robot_model_halfspaces_->clear();
  robot_model_halfspaces_->resize(robot_model_->num_tris);
  for (unsigned int i = 0; i < robot_model_->num_tris; ++i)
  {
    fcl::Triangle tri = robot_model_->tri_indices[i];
    (*robot_model_halfspaces_)[i] = convertTriangleToHalfspace<FCLFloat>(
        fcl::TriangleP<FCLFloat>(
            robot_model_->vertices[tri[0]],
            robot_model_->vertices[tri[1]],
            robot_model_->vertices[tri[2]]));
  }
  calculateCacheDistanceThresholds();
}

std::string Costmap3DQuery::getFileNameFromPackageURL(const std::string& url)
{
  /* Unfortunately the resource retriever does not have a way to get a path from a URL
   * (it only returns the contents of a URL in memory), and equally unfortunate, PCL does not
   * have a way to parse an STL from memory. Therefore we have to duplicate some of the
   * (slightly-modified) resource retriever code here. */
  // Reference: https://github.com/ros/resource_retriever/blob/kinetic-devel/src/retriever.cpp
  std::string mod_url = url;
  if (url.find("package://") == 0)
  {
    mod_url.erase(0, strlen("package://"));
    size_t pos = mod_url.find("/");
    if (pos == std::string::npos)
    {
      ROS_ERROR_STREAM("Costmap3DQuery: Could not parse package:// format URL "
                       << url
                       << " query object will always return collision!");
      return "";
    }

    std::string package = mod_url.substr(0, pos);
    mod_url.erase(0, pos);
    std::string package_path = ros::package::getPath(package.c_str());

    if (package_path.empty())
    {
      ROS_ERROR_STREAM("Costmap3DQuery: Package [" << package << "] from URL "
                       << url << " does not exist, "
                       << " query object will always return collision!");
      return "";
    }

    mod_url = package_path + mod_url;
  }
  else if (url.find("file://") == 0)
  {
    mod_url.erase(0, strlen("file://"));
  }

  return mod_url;
}

double Costmap3DQuery::footprintCost(const geometry_msgs::Pose& pose, Costmap3DQuery::QueryRegion query_region)
{
  // TODO: implement as cost query. For now, just translate a collision to cost
  return footprintCollision(pose, query_region) ? -1.0 : 0.0;
}

bool Costmap3DQuery::footprintCollision(const geometry_msgs::Pose& pose, Costmap3DQuery::QueryRegion query_region)
{
  // It is more correct and even more efficient to query the distance to find
  // collisions than it is to use FCL to directly find octomap collisions.
  // This is because our distance query correctly handles interior collisions,
  // which requires finding the nearest octomap box, which an FCL collision
  // will not do.
  DistanceOptions opts;
  opts.query_region = query_region;
  // For collision only queries, we only need to check when bounding volumes overlap.
  // Make cache-misses in such cases very fast by limiting the bound distance.
  // Bound to just above zero, as zero would be considered a collision.
  // Remember that min for double is the smallest positive (normalized) double.
  opts.distance_limit = std::numeric_limits<double>::min();
  // We want exact answers for collision only checks.
  opts.relative_error = 0.0;
  return calculateDistance(pose, opts) <= 0.0;
}

// Discern if the given octomap box is an interior collision and adjust
// distance or signed distance appropriately.
// This is done by using PCL's point-in-a-mesh call which works on concave
// meshes that represent closed polyhedra. FCL handles concave meshes as
// surfaces, not volumes.
double Costmap3DQuery::handleDistanceInteriorCollisions(
      const DistanceCacheEntry& cache_entry,
      const geometry_msgs::Pose& pose)
{
  FCLFloat distance;

  // Turn pose into tf
  const fcl::Transform3<FCLFloat> pose_tf(costmap_3d::poseToFCLTransform<FCLFloat>(pose));

  // Start with the interior collision check as it is very fast.
  // Find out if the center of the box is inside the given footprint volume mesh.
  distance = interior_collision_lut_.distance(
      *cache_entry.octomap_box,
      cache_entry.octomap_box_tf,
      pose_tf,
      pose_tf.inverse());

  if (distance < 0.0)
  {
    return distance;
  }

  // We need to calculate the distance between the mesh triangle at the new
  // pose and the box.
  //
  // As of the time this code was written, the normal FCL API does not
  // allow box/triangle distance or signed distance queries.
  // Yet FCL internally does such checks all the time, so use the
  // internal mechanism for now.
  fcl::detail::GJKSolver_libccd<FCLFloat> solver;
  solver.shapeTriangleDistance(
      *cache_entry.octomap_box,
      cache_entry.octomap_box_tf,
      cache_entry.mesh_triangle->a,
      cache_entry.mesh_triangle->b,
      cache_entry.mesh_triangle->c,
      pose_tf,
      &distance);

  // Box/triangle intersect, use penetration depth with box/halfspace model.
  if (distance < 0.0)
  {
    distance = boxHalfspaceSignedDistance(
        *cache_entry.octomap_box,
        cache_entry.octomap_box_tf,
        cache_entry.mesh_triangle_id,
        pose_tf);
  }

  return distance;
}

Costmap3DQuery::FCLFloat Costmap3DQuery::boxHalfspaceSignedDistance(
    const fcl::Box<FCLFloat>& box,
    const fcl::Transform3<FCLFloat>& box_tf,
    int mesh_triangle_id,
    const fcl::Transform3<FCLFloat>& mesh_tf) const
{
  fcl::Halfspace<FCLFloat> halfspace(fcl::transform(
      (*robot_model_halfspaces_)[mesh_triangle_id], mesh_tf));
  return costmap_3d::boxHalfspaceSignedDistance<FCLFloat>(box, box_tf, halfspace);
}

double Costmap3DQuery::calculateDistance(const geometry_msgs::Pose& pose,
                                         const Costmap3DQuery::DistanceOptions& opts)
{
  // Make convenience aliases of some of the options
  const QueryRegion query_region = opts.query_region;
  const QueryObstacles query_obstacles = opts.query_obstacles;
  // Let compiler optimize track_statistics by adding it to the stack.
  const bool track_statistics = track_statistics_;

  std::chrono::high_resolution_clock::time_point start_time;
  if (track_statistics)
  {
    start_time = std::chrono::high_resolution_clock::now();
    queries_since_clear_.fetch_add(1, std::memory_order_relaxed);
  }

  // Use the passed in distance limit up to the "min" float (close to zero)
  // We do not want to use a zero or negative distance limit from the options
  // because that will prevent us from correctly searching the octomap/mesh
  FCLFloat pose_distance = std::max(opts.distance_limit, std::numeric_limits<FCLFloat>::min());

  // Handle reuse results first before acquiring any locks.
  // It is already up to the caller using reuse results to ensure the costmap
  // hasn't updated in-between calls, so it is pointless to acquire the lock to
  // check.
  if (opts.reuse_past_result)
  {
    // Reuse the past result (if there was one) and directly return the new distance.
    // In cases where the caller knows this is the correct behavior (such as
    // when making minor perturbations to estimate derivatives), this is
    // faster than having to calculate the hash and find the cache entry.
    // Do not check if the entry is in the query region still, it is assumed
    // the caller wants to ignore the potential error this would cause, as this
    // interface is mainly useful for very small perturbations.
    if (tls_last_cache_entries_[query_region][query_obstacles])
    {
      double distance = handleDistanceInteriorCollisions(
          tls_last_cache_entries_[query_region][query_obstacles],
          pose);
      if (track_statistics)
      {
        reuse_results_since_clear_.fetch_add(1, std::memory_order_relaxed);
        reuse_results_since_clear_us_.fetch_add(
            std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start_time).count(),
            std::memory_order_relaxed);
      }
      return distance;
    }
    else
    {
      // It is possible for the caller to have just previously called on this
      // thread and got no results because of the way the query was limited.
      // Assume that this query is limited the same way, and return the current
      // pose_distance, which will be set from the passed in distance limit.
      // This will be the exact same distance returned by the previous call,
      // which is the desired result when the previous query resulted in
      // nothing found.
      return pose_distance;
    }
  }

  shared_lock read_lock(upgrade_mutex_);

  if (!robot_model_)
  {
    // We failed to create a robot model.
    // The failure would have been logged, so simply return collision.
    return -1.0;
  }

  // Check if an upgrade lock will be necessary. Do not do this with an upgrade
  // lock held, as only one thread can have an upgrade lock at a time.
  bool need_upgrade_lock = false;
  if (needCheckCostmap())
  {
    need_upgrade_lock = true;
  }
  if (query_obstacles == NONLETHAL_ONLY && !nonlethal_octree_ptr_)
  {
    need_upgrade_lock = true;
  }

  if (need_upgrade_lock)
  {
    // Drop the shared lock to get an upgrade lock. The only safe way to do
    // this and avoid deadlock/livelock is to drop the shared lock and then get
    // an upgrade lock. Because the lock is dropped, we must re-check the
    // conditions after acquiring the upgrade mutex (as it may no longer be
    // necessary). Because this very slow path is only encountered at most once
    // or twice per update cycle (depending on if non-lethal queries are
    // issued and if this is a buffered query) this is OK.
    read_lock.unlock();
    read_lock.release();
    upgrade_lock upgradeable_lock(upgrade_mutex_);
    checkCostmap(upgradeable_lock);
    if (query_obstacles == NONLETHAL_ONLY && !nonlethal_octree_ptr_)
    {
      // We do not yet have a nonlethal copy of the tree, create one now and save it until
      // the next costmap update. In order to do this safely we must grab the write lock.
      // This is OK, it will only happen once per update cycle.
      upgrade_to_unique_lock write_lock(upgradeable_lock);
      nonlethal_octree_ptr_ = std::static_pointer_cast<const Costmap3D>(octree_ptr_)->nonlethalOnly();
    }
    // Re-acquire the shared lock. This can be done atomically by downgrading
    // the upgrade lock to a shared lock.
    read_lock = shared_lock(std::move(upgradeable_lock));
  }

  std::shared_ptr<const octomap::OcTree> octree_to_query;
  if (query_obstacles == NONLETHAL_ONLY)
  {
    octree_to_query = nonlethal_octree_ptr_;
  }
  else
  {
    octree_to_query = octree_ptr_;
  }
  assert(octree_to_query);

  // FCL does not correctly handle an empty octomap.
  if (octree_to_query->size() == 0 || !octree_to_query->isNodeOccupied(octree_to_query->getRoot()))
  {
    if (track_statistics)
    {
      empties_since_clear_.fetch_add(1, std::memory_order_relaxed);
    }
    return std::numeric_limits<double>::max();
  }

  DistanceCacheKey exact_cache_key(pose, query_region, query_obstacles);
  {
    shared_lock cache_read_lock(exact_distance_cache_mutex_);
    auto exact_cache_entry = exact_distance_cache_.find(exact_cache_key);
    if (exact_cache_entry != exact_distance_cache_.end())
    {
      double distance = exact_cache_entry->second.distance;
      // Be sure to update the TLS last cache entry.
      // We do not need the write lock to update thread local storage.
      tls_last_cache_entries_[query_region][query_obstacles] = exact_cache_entry->second;
      if (track_statistics)
      {
        exact_hits_since_clear_.fetch_add(1, std::memory_order_relaxed);
        exact_hits_since_clear_us_.fetch_add(
            std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start_time).count(),
            std::memory_order_relaxed);
      }
      return distance;
    }
  }

  fcl::DistanceRequest<FCLFloat> request;
  fcl::DistanceResult<FCLFloat> result;

  // Setup the regions of interest corresponding to the query region.
  // This has to be done prior to checking the distance caches, so the cache entry
  // can be checked to ensure it is still inside the region.
  // Store the region on the stack to avoid dynamic memory allocation.
  RegionsOfInterestAtPose rois(query_region, query_region_scale_, pose);

  DistanceCacheKey milli_cache_key(pose, query_region, query_obstacles, pose_milli_bins_per_meter_, pose_milli_bins_per_radian_);
  bool milli_hit = false;
  DistanceCacheEntry cache_entry;
  {
    shared_lock cache_read_lock(milli_distance_cache_mutex_);
    auto milli_cache_entry = milli_distance_cache_.find(milli_cache_key);
    if (milli_cache_entry != milli_distance_cache_.end())
    {
      cache_entry = milli_cache_entry->second;
    }
  }
  if (cache_entry &&
      rois.distanceCacheEntryInside(cache_entry))
  {
    double distance = handleDistanceInteriorCollisions(
        cache_entry,
        pose);
    if ((!opts.exact_signed_distance && distance <= 0.0) ||
        std::isfinite(cache_entry.distance) && (
          cache_entry.distance > milli_cache_threshold_ ||
          distance > milli_cache_threshold_))
    {
      // Be sure to update the TLS last cache entry.
      // We do not need the write lock to update thread local storage.
      tls_last_cache_entries_[query_region][query_obstacles] = cache_entry;
      if (track_statistics)
      {
        fast_milli_hits_since_clear_.fetch_add(1, std::memory_order_relaxed);
        fast_milli_hits_since_clear_us_.fetch_add(
            std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start_time).count(),
            std::memory_order_relaxed);
      }
      return distance;
    }
    else
    {
      // we are too close to directly use the milli-cache, but use the
      // calculated distance as the pose_distance upper bound.
      if (distance < pose_distance)
      {
        milli_hit = true;
        pose_distance = distance;
        cache_entry.setupResult(&result);
      }
    }
  }

  DistanceCacheKey micro_cache_key(pose, query_region, query_obstacles, pose_micro_bins_per_meter_, pose_micro_bins_per_radian_);
  bool micro_hit = false;
  cache_entry.clear();
  {
    shared_lock cache_read_lock(micro_distance_cache_mutex_);
    auto micro_cache_entry = micro_distance_cache_.find(micro_cache_key);
    if (micro_cache_entry != micro_distance_cache_.end())
    {
      cache_entry = micro_cache_entry->second;
    }
  }
  if (cache_entry &&
      rois.distanceCacheEntryInside(cache_entry))
  {
    double distance = handleDistanceInteriorCollisions(
        cache_entry,
        pose);
    if ((!opts.exact_signed_distance && distance <= 0.0) ||
        cache_entry.distance > micro_cache_threshold_ ||
        distance > micro_cache_threshold_)
    {
      // Be sure to update the TLS last cache entry.
      // We do not need the write lock to update thread local storage.
      tls_last_cache_entries_[query_region][query_obstacles] = cache_entry;
      if (track_statistics)
      {
        fast_micro_hits_since_clear_.fetch_add(1, std::memory_order_relaxed);
        fast_micro_hits_since_clear_us_.fetch_add(
            std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start_time).count(),
            std::memory_order_relaxed);
      }
      return distance;
    }
    else
    {
      // we are too close to directly use the micro-cache, but use the
      // calculated distance as the pose_distance upper bound.
      if (distance < pose_distance)
      {
        micro_hit = true;
        pose_distance = distance;
        cache_entry.setupResult(&result);
      }
    }
  }

  DistanceCacheKey cache_key(pose, query_region, query_obstacles, pose_bins_per_meter_, pose_bins_per_radian_);
  bool cache_hit = false;
  // Assume a milli hit or micro hit is much better than the normal distance
  // cache at setting an upper bound, skip the normal cache lookup in such
  // cases.
  if (!milli_hit && !micro_hit)
  {
    cache_entry.clear();
    {
      shared_lock cache_read_lock(distance_cache_mutex_);
      auto distance_cache_entry = distance_cache_.find(cache_key);
      if (distance_cache_entry != distance_cache_.end())
      {
        cache_entry = distance_cache_entry->second;
      }
    }
    if (cache_entry &&
        rois.distanceCacheEntryInside(cache_entry))
    {
      // Cache hit, find the distance between the mesh triangle at the new pose
      // and the octomap box, and use this as our initial guess in the result.
      // This greatly prunes the search tree, yielding a big increase in runtime
      // performance.
      double distance = handleDistanceInteriorCollisions(
          cache_entry,
          pose);
      if (distance < pose_distance)
      {
        cache_hit = true;
        pose_distance = distance;
        cache_entry.setupResult(&result);
      }
    }
  }

  // Check the distance from the last cache entry too, and use it if it is a
  // better match. This is especialy useful if we miss every cache, but have a
  // last entry, and the queries are being done on a path. First check the TLS
  // which will invalidate it if it is no longer valid.
  checkTLS();
  if (tls_last_cache_entries_[query_region][query_obstacles] &&
      rois.distanceCacheEntryInside(tls_last_cache_entries_[query_region][query_obstacles]))
  {
    double last_entry_distance = handleDistanceInteriorCollisions(
        tls_last_cache_entries_[query_region][query_obstacles],
        pose);
    if (last_entry_distance < pose_distance)
    {
      pose_distance = last_entry_distance;
      tls_last_cache_entries_[query_region][query_obstacles].setupResult(&result);
    }
  }

  request.rel_err = opts.relative_error;
  request.enable_nearest_points = true;
  request.enable_signed_distance = true;

  result.min_distance = pose_distance;

  double distance;

  // Because FCL's OcTree/Mesh distance treats the Mesh as hollow, we must
  // use our own distance code which treats the Mesh as a closed mesh
  // defining a filled volume. Otherwise, there are cases where an octomap
  // box is inside the mesh, but the distance is positive.
  // The solver and octree_solver need to be on the stack as they are not
  // thread-safe.
  FCLSolver solver;
  // Use our interior collision LUT to model the robot as a volume.
  // Use box-halfspace distance to model box/mesh penetrations for signed distance.
  OcTreeMeshSolver<FCLSolver> octree_solver(
      &solver,
      std::bind(&InteriorCollisionLUT<FCLFloat>::distance,
                &interior_collision_lut_,
                std::placeholders::_1,
                std::placeholders::_2,
                std::placeholders::_3,
                std::placeholders::_4,
                std::placeholders::_5),
      std::bind(&Costmap3DQuery::boxHalfspaceSignedDistance,
                this,
                std::placeholders::_1,
                std::placeholders::_2,
                std::placeholders::_3,
                std::placeholders::_4));
  if (query_obstacles == NONLETHAL_ONLY)
  {
    octree_solver.setUncertainOnly(true);
  }

  if (opts.exact_signed_distance)
  {
    octree_solver.setExactSignedDistance(true);
  }

  if (result.min_distance > 0.0 || opts.exact_signed_distance)
  {
    fcl::OcTree<FCLFloat> fcl_octree(octree_to_query);
    // Always setup the correct occupancy limits
    fcl_octree.setFreeThres(octomap::probability(FREE));
    fcl_octree.setOccupancyThres(octomap::probability(LETHAL));
    rois.setupFCLOctree(&fcl_octree);

    octree_solver.distance(
        &fcl_octree,
        robot_model_.get(),
        fcl::Transform3<FCLFloat>::Identity(),
        poseToFCLTransform<FCLFloat>(pose),
        request,
        &result);
  }
  distance = result.min_distance;

  // Note that it is possible for the result to be empty. The octomap might
  // only contain non-lethal leaves and we may have missed every cache.
  // Or the query region may be set to something other than ALL and there are
  // no map entries in the queried region.
  // If we get no result primitives, do not add null pointers to the cache!
  if (result.primitive1 && result.primitive2)
  {
    DistanceCacheEntry new_entry(result);
    new_entry.distance = distance;

    // Update distance caches.
    // While it may seem expensive to copy the cache entries into the
    // caches, it prevents cache aliasing and avoids dynamic memory.
    {
      unique_lock cache_write_lock(distance_cache_mutex_);
      distance_cache_[cache_key] = new_entry;
    }
    {
      unique_lock cache_write_lock(milli_distance_cache_mutex_);
      milli_distance_cache_[milli_cache_key] = new_entry;
    }
    {
      unique_lock cache_write_lock(micro_distance_cache_mutex_);
      micro_distance_cache_[micro_cache_key] = new_entry;
    }
    {
      unique_lock cache_write_lock(exact_distance_cache_mutex_);
      exact_distance_cache_[exact_cache_key] = new_entry;
    }
    tls_last_cache_entries_[query_region][query_obstacles] = new_entry;
  }
  else
  {
    // Need to erase the last cache since there was no result.
    // This should only happen if the query was restricted such that no octomap
    // cells were considered. One way this could happen is if the distance was
    // limited by the caller. If the next call is a reuse_past_result query, we
    // must set the last cache to nullptr to prevent erroneously using some
    // other result.
    // No need to lock TLS
    tls_last_cache_entries_[query_region][query_obstacles].clear();
  }

  if (track_statistics)
  {
    if (micro_hit)
    {
      slow_micro_hits_since_clear_.fetch_add(1, std::memory_order_relaxed);
      slow_micro_hits_since_clear_us_.fetch_add(
          std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start_time).count(),
          std::memory_order_relaxed);
      hit_fcl_bv_distance_calculations_.fetch_add(result.bv_distance_calculations);
      hit_fcl_primative_distance_calculations_.fetch_add(result.primative_distance_calculations);
    }
    else if (milli_hit)
    {
      slow_milli_hits_since_clear_.fetch_add(1, std::memory_order_relaxed);
      slow_milli_hits_since_clear_us_.fetch_add(
          std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start_time).count(),
          std::memory_order_relaxed);
      hit_fcl_bv_distance_calculations_.fetch_add(result.bv_distance_calculations);
      hit_fcl_primative_distance_calculations_.fetch_add(result.primative_distance_calculations);
    }
    else if (cache_hit)
    {
      hits_since_clear_.fetch_add(1, std::memory_order_relaxed);
      hits_since_clear_us_.fetch_add(
          std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start_time).count(),
          std::memory_order_relaxed);
      hit_fcl_bv_distance_calculations_.fetch_add(result.bv_distance_calculations);
      hit_fcl_primative_distance_calculations_.fetch_add(result.primative_distance_calculations);
    }
    else
    {
      misses_since_clear_us_.fetch_add(
          std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start_time).count(),
          std::memory_order_relaxed);
      miss_fcl_bv_distance_calculations_.fetch_add(result.bv_distance_calculations);
      miss_fcl_primative_distance_calculations_.fetch_add(result.primative_distance_calculations);
    }
  }

  return distance;
}

double Costmap3DQuery::footprintDistance(const geometry_msgs::Pose& pose,
                                         Costmap3DQuery::QueryRegion query_region,
                                         bool reuse_past_result,
                                         double relative_error)
{
  DistanceOptions opts;
  opts.query_region = query_region;
  opts.reuse_past_result = reuse_past_result;
  opts.relative_error = relative_error;
  return calculateDistance(pose, opts);
}

double Costmap3DQuery::footprintSignedDistance(const geometry_msgs::Pose& pose,
                                               Costmap3DQuery::QueryRegion query_region,
                                               bool reuse_past_result,
                                               double relative_error,
                                               bool exact_signed_distance)
{
  DistanceOptions opts;
  opts.query_region = query_region;
  opts.reuse_past_result = reuse_past_result;
  opts.relative_error = relative_error;
  opts.signed_distance = true;
  opts.exact_signed_distance = exact_signed_distance;
  return calculateDistance(pose, opts);
}

double Costmap3DQuery::footprintDistance(const geometry_msgs::Pose& pose,
                                         const Costmap3DQuery::DistanceOptions& opts)
{
  return calculateDistance(pose, opts);
}

void Costmap3DQuery::checkInteriorCollisionLUT()
{
  if (last_octomap_resolution_ != octree_ptr_->getResolution())
  {
    // Resolution changed, need to setup our interior collision LUT
    last_octomap_resolution_ = octree_ptr_->getResolution();
    const double box_size = last_octomap_resolution_;
    // Instead of adding the orientation to the LUT, just be more than double the
    // spatial resolution. This gives acceptable results without using much
    // memory.
    const double lut_res = box_size / 2.5;
    interior_collision_lut_.setup(box_size, lut_res, *robot_model_, crop_hull_, robot_model_halfspaces_);
  }
}

}  // namespace costmap_3d
