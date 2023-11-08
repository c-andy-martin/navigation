/*********************************************************************
 *
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2024, Badger Technologies LLC
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
#ifndef COSTMAP_3D_BIN_POSE_H_
#define COSTMAP_3D_BIN_POSE_H_

#include <limits>
#include <cmath>

#include <Eigen/Dense>
#include <geometry_msgs/Pose.h>

namespace costmap_3d
{

inline geometry_msgs::Pose binPose(
    const geometry_msgs::Pose& pose,
    int bins_per_meter,
    int bins_per_rotation)
{
  geometry_msgs::Pose rv;
  // std::round is slow on AVX2, but floor is fast, and where the
  // quantization happens is not important
  rv.position.x = std::floor(pose.position.x * bins_per_meter) / bins_per_meter;
  rv.position.y = std::floor(pose.position.y * bins_per_meter) / bins_per_meter;
  rv.position.z = std::floor(pose.position.z * bins_per_meter) / bins_per_meter;

  // Decompose uniformly from SO(3) to the volume of a cube by solving this
  // uniform distribution from "Uniform Random Rotations" by K. Shoemake in
  // reverse:
  //
  // x = std::sqrt(u0) * cos(2*M_PI*u1)
  // y = std::sqrt(u0) * sin(2*M_PI*u1)
  // z = std::sqrt(1.0 - u0) * cos(2*M_PI*u2)
  // w = std::sqrt(1.0 - u0) * sin(2*M_PI*u2)
  //
  // Note x/y/z/w can be in any order and the distribution is still evenly
  // spaced in SO(3), so the choice of which way to assign x/y/z/w is arbitrary
  // but needs to remain consistent. The order above is chosen so rotations
  // around the z axis are binned evenly on rotation angle (but x/y axes are
  // not). This is a good choice for navigation, where 2D navigation will
  // always use rotations only about the z-axis, and 3D navigation may use any
  // orientation, but would often have rotation axes near the z-axis as well.

  // Go ahead and map any rotation in the "lower" hyper-hemisphere to its
  // equivalent antipode rotation. This prevents us from having to map antipode
  // rotations to the same bins later on, which would be complicated due to
  // rounding leaving the numbers of matching rotations slightly different.
  bool negate =
    pose.orientation.w < 0.0 || (pose.orientation.w == 0.0 && (
        pose.orientation.z < 0.0 || (pose.orientation.z == 0.0 && (
            pose.orientation.y < 0.0 || (pose.orientation.y == 0.0 && (
                pose.orientation.x < 0.0))))));
  // Keep the quaternion in the same order as the geometry message.
  Eigen::Array4d q(
      negate ? -pose.orientation.x : pose.orientation.x,
      negate ? -pose.orientation.y : pose.orientation.y,
      negate ? -pose.orientation.z : pose.orientation.z,
      negate ? -pose.orientation.w : pose.orientation.w);

  Eigen::Array3d u;
  // There are several ways to solve for u0. The simplest is:
  //
  // x = std::sqrt(u0) * cos(2*M_PI*u1)
  // y = std::sqrt(u0) * sin(2*M_PI*u1)
  //
  // So
  //
  // x^2 = u0 * cos^2(2*M_PI*u1)
  // y^2 = u0 * sin^2(2*M_PI*u1)
  //
  // So
  //
  // x^2 + y^2 = u0 * cos^2(2*M_PI*u1) + u0 * sin^2(2*M_PI*u1)
  //
  // Which simplifies to:
  //
  // x^2 + y^2 = u0
  //
  u(0) = q(0) * q(0) + q(1) * q(1);
  // Check against solving for u0 with z/w
  assert(std::abs(u(0) - (1.0 - (q(2) * q(2) + q(3) * q(3)))) < 1e-6);
  double angle1;
  double angle2;
  // Solve for u1 and u2 using atan2 to avoid singularities.
  // This works because:
  // y / x = tan(2*M_PI*u1)
  // w / z = tan(2*M_PI*u2)
  //
  // However, in the event that both coordinates are some form of zero, atan2
  // will return +/-pi or +/-0 depending on the signs of the zeros.
  // Yet those should all bin to the same value, so specifically handle that.
  // The only way that both coordinates for angle1 are zero is if u0 is 0.
  // The only way that both coordinates for angle2 are zero is if u0 is 1.
  // Also handle clamping u0 just below 1 if it is at 1 to avoid placing it
  // in its own bin when performing the floor operation down below.
  if (u(0) == 0.0)
  {
    angle1 = 0.0;
    angle2 = std::atan2(q(3), q(2));
  }
  else if (u(0) >= 1.0)
  {
    u(0) = 1.0 - std::numeric_limits<double>::epsilon();
    angle1 = std::atan2(q(1), q(0));
    angle2 = 0.0;
  }
  else
  {
    angle1 = std::atan2(q(1), q(0));
    angle2 = std::atan2(q(3), q(2));
  }
  u(1) = angle1 / (2.0 * M_PI);
  u(2) = angle2 / (2.0 * M_PI);
  // If the angle1 is at the end of the possible range, wrap it around back to
  // a half negative rotation. This keeps it in the proper bin, and also keeps
  // the floor below from generating a too-high value (the end of the range is
  // not inclusive).
  if (u(1) >= 0.5)
  {
    u(1) = -0.5;
  }
  // For angle2, because w is always >= 0.0, clamp any value at the ends to the
  // very last bin. 0.5 or -0.5 may happen when w is +/-0.
  // If this were not done quaternions with w at +/-0 would fall into their own
  // tiny bins.
  if (std::abs(u(2)) >= 0.5)
  {
    u(2) = 0.5 - std::numeric_limits<double>::epsilon();
  }
  // Verify the proper ranges when assertions are enabled (during testing).
  assert(u(0) >= 0.0);
  assert(u(0) < 1.0);
  assert(u(1) >= -0.5);
  assert(u(1) < 0.5);
  assert(u(2) >= 0.0);
  assert(u(2) < 0.5);

  // Need to double the number of buckets to get the correct size in rotations.
  // This is because a value of one of the u values of .5 represents an entire
  // rotation.
  const double bucket_factor = bins_per_rotation * 2;
  u *= bucket_factor;
  // Avoid zero values for the binned u values, as that would create huge
  // buckets at zero (since multiplying by zero destroys the significance of
  // the other variable). Instead, represent each bin as its midpoint.
  u = (u.array().floor() + Eigen::Array3d(0.5, 0.5, 0.5)) / bucket_factor;
  assert(u(0) > 0.0);
  assert(u(0) < 1.0);
  assert(u(1) > -0.5);
  assert(u(1) < 0.5);
  assert(u(2) > 0.0);
  assert(u(2) < 0.5);
  // Stash away the angles in locals to make it easy for the compiler to
  // optimize the std::sin and std::cos to a single call to sincos.
  double angle_u1 = 2 * M_PI * u(1);
  double angle_u2 = 2 * M_PI * u(2);
  double sin_angle_u1 = std::sin(angle_u1);
  double cos_angle_u1 = std::cos(angle_u1);
  double sin_angle_u2 = std::sin(angle_u2);
  double cos_angle_u2 = std::cos(angle_u2);
  q(0) = std::sqrt(u(0)) * cos_angle_u1;
  q(1) = std::sqrt(u(0)) * sin_angle_u1;
  q(2) = std::sqrt(1.0 - u(0)) * cos_angle_u2;
  q(3) = std::sqrt(1.0 - u(0)) * sin_angle_u2;

  assert(std::isfinite(q(0)));
  assert(std::isfinite(q(1)));
  assert(std::isfinite(q(2)));
  assert(std::isfinite(q(3)));
  // The middle of bins should never be on a zero, which is important to keep
  // from mapping too many bins together at a zero. Also, this avoids returning
  // negative zero which aids in hashing using the raw bits of the floats.
  assert(q(0) != 0.0);
  assert(q(1) != 0.0);
  assert(q(2) != 0.0);
  assert(q(3) > 0.0);

  rv.orientation.x = q(0);
  rv.orientation.y = q(1);
  rv.orientation.z = q(2);
  rv.orientation.w = q(3);

  return rv;
}

// Returns maximum angular distance for two rotations that hit the same bin
// (from furthest corner to furthest corner). If the desired error is to the
// center of the bin (the binned pose value), simply call with double the
// bins_per_rotation value (to get the max bin corner to the center bin
// position).
inline double binPoseAngularDistanceLimit(int bins_per_rotation, bool two_d_mode = false)
{
  if (two_d_mode)
  {
    // If running in 2D mode, the pose bins are evenly spaced when the
    // rotation axis is the z-axis. In such cases, the true maximum angular error
    // is just 2 * M_PI / bins_per_rotation.
    return 2 * M_PI / bins_per_rotation;
  }
  // The largest bins are those at u0 = 1 or u0 = 0. To find the maximum
  // angular bin size (error) simply measure the angular distance along the
  // diagonal of the bin.
  //
  // So compare the corresponding quaternions at u0 = u1 = u2 = 0.0, and the
  // corresponding corner u0 = u1 = u2 = 1 / (2 * bins_per_rotation)
  //
  // So q0 (at the corner) would be:
  //
  // w = std::sqrt(1.0 - u0) * sin(2*M_PI*u1)
  // x = std::sqrt(u0) * cos(2*M_PI*u2)
  // y = std::sqrt(u0) * sin(2*M_PI*u2)
  // z = std::sqrt(1.0 - u0) * cos(2*M_PI*u1)
  // w = 0
  // x = 0
  // y = 0
  // z = 1
  //
  // Because the angular bin size is found using the dot product of the
  // corners, and the corner being examined has zero coefficients at w/x/y,
  // only the z coefficient is relevant.
  //
  // So z at the corresponding corner:
  // z_corner = std::sqrt(1 - .5 / bins_per_rotation) * cos(M_PI / bins_per_rotation)
  //
  // The angular distance from the original quaternion to the binned quaternion is therefore:
  // 2 * acos(q dot binned_q)
  // 2 * acos(z_corner)
  // so:
  // 2 * acos(std::sqrt(1 - .5 / bins_per_rotation) * cos(M_PI / bins_per_rotation))
  //
  return 2 * acos(std::sqrt(1.0 - .5 / bins_per_rotation) * cos(M_PI / bins_per_rotation));
}

}  // namespace costmap_3d

#endif  // COSTMAP_3D_BIN_POSE_H_
