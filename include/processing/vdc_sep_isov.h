//! @file vdc_sep_isov.h
//! @brief Isosurface vertex separation algorithms
//! @details Provides different methods to separate active cubes to ensure
//! non-adjacent isosurface vertices for improved triangulation quality.

#ifndef VDC_SEP_ISOV_H
#define VDC_SEP_ISOV_H

#include "core/vdc_type.h"
#include "processing/vdc_grid.h"
#include <vector>

//! @brief Unified separation of active cubes parameterized by distance and splits.
/*!
 * Applies greedy selection in a refined grid of factor (K+1) with clearance
 * measured in subcubes. Cube centers are snapped to power-of-two fractions
 * based on K to reduce floating-point error.
 *
 * @param activeCubes Input vector of active cubes (modified by sorting)
 * @param grid Scalar grid metadata
 * @param isovalue Target isovalue
 * @param sep_dist Separation distance D in subcubes (D<=1 disables filtering)
 * @param sep_split Number of splits K per axis (refine factor = K+1)
 * @return Separated subset of cubes
 */
std::vector<Cube> separate_active_cubes(
    std::vector<Cube> &activeCubes,
    const UnifiedGrid &grid,
    float isovalue,
    int sep_dist,
    int sep_split);

//! @brief Compute accurate iso-crossing point using edge-intersection centroids
/*!
 * Finds all edge-isovalue intersections and returns their centroid.
 * Used for separation logic, not for Delaunay triangulation.
 *
 * @param grid The scalar grid
 * @param i,j,k Cube indices
 * @param isovalue Target isovalue
 * @return Centroid of edge-crossing points
 */
Point compute_iso_crossing_point_accurate(
    const UnifiedGrid &grid,
    int i, int j, int k,
    float isovalue);

// ============================================================================
// DEBUG UTILITIES (temporary)
// ============================================================================

//! @brief Check all triples of iso-crossing points for near-collinearity
void check_collinear_isocrossings(const std::vector<Cube>& activeCubes,
                                  const UnifiedGrid& grid,
                                  float isovalue,
                                  double angle_threshold = 5.0);

#endif // VDC_SEP_ISOV_H
