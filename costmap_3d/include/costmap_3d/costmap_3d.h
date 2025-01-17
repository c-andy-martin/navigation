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
#ifndef COSTMAP_3D_COSTMAP_3D_H_
#define COSTMAP_3D_COSTMAP_3D_H_

#include <memory>
#include <string>

#include <geometry_msgs/Point.h>
#include <octomap/octomap.h>

namespace costmap_3d
{

/* OccupancyOcTree Log odds values for various costmap states. */
typedef float Cost;
// Sentinel. Don't choose NAN as it can't be pruned in octomap
const Cost UNKNOWN = -200.0;
const Cost FREE = -10.0;
// A "nonlethal" Cost is any Cost strictly between FREE and LETHAL
// Specifically, UNKNOWN is *not* nonlethal. Nonlethal is meant to
// represent a sensed object that is desirable to avoid but would
// likely be OK to collide with.
const Cost LETHAL = 10.0;
/* A log-odds between FREE and LETHAL represents a non-lethal,
 * but greater than zero cost for that space */

// Default depth for costmap octrees. This can be overriden in the constructor.
// The default is 20 to be able to represent a large space (~5km in any
// direction) at 1cm resolution. Increasing the depth allows larger spaces to
// be reprsented at a slight increased CPU cost per costmap query.
constexpr unsigned int DEFAULT_DEPTH = 20;

using Costmap3DIndex = octomap::OcTreeKey;
using Costmap3DIndexEntryType = octomap::key_type;

class Costmap3D : public octomap::OcTree
{
public:
  explicit Costmap3D(double resolution, unsigned int depth = DEFAULT_DEPTH);
  // Load a Costmap3D from a path to a binary octomap file
  explicit Costmap3D(std::string path);
  Costmap3D(const Costmap3D& rhs);

  // Return a new Costmap3D contatinaing only the nonlethal cells.
  // This is useful for speeding up nonlethal only queries, as the tree only
  // stores the maximum cost in the inner nodes it is impossible to cut off the
  // searches early in the presence of lethal obstacles. This nonlethal only
  // tree can be constructed when the costmap updates and stored and reused for
  // all nonlethal queries.
  std::shared_ptr<Costmap3D> nonlethalOnly() const;
protected:
  virtual void init();
};

typedef std::shared_ptr<Costmap3D> Costmap3DPtr;
typedef std::shared_ptr<const Costmap3D> Costmap3DConstPtr;

octomap::point3d toOctomapPoint(const geometry_msgs::Point& pt);
geometry_msgs::Point fromOctomapPoint(const octomap::point3d& point);

}  // namespace costmap_3d

#endif  // COSTMAP_3D_COSTMAP_3D_H_
