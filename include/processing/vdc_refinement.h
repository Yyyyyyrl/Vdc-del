//! @file vdc_refinement.h
//! @brief Facet-centric Delaunay refinement to improve small angles near the isosurface.

#ifndef VDC_REFINEMENT_H
#define VDC_REFINEMENT_H

#include "core/vdc_commandline.h"
#include "processing/vdc_grid.h"

//! @brief Statistics reported by the refinement pass.
struct RefinementStats
{
    int iterations_run = 0;            //!< Number of refinement iterations executed.
    std::size_t candidate_facets = 0;  //!< Facets matching angle criteria.
    std::size_t bipolar_facets = 0;    //!< Those whose dual Voronoi edge is bipolar.
    std::size_t inserted_points = 0;   //!< Points actually inserted into the Delaunay.
};

//! @brief Refine poorly shaped facets whose dual Voronoi edge is bipolar.
/*!
 * Scans finite Delaunay facets, checks their corner angles, and if the
 * dual Voronoi edge is bipolar with respect to the target isovalue, inserts a
 * new vertex snapped to the nearest active cube (or subcell). Refinement
 * triggers when the minimum angle drops below the configured threshold or the
 * maximum angle exceeds the user-provided limit.
 *
 * @param dt Delaunay triangulation to refine (modified in place).
 * @param grid Scalar grid used for value interpolation.
 * @param active_cubes Active cubes detected for the current isovalue.
 * @param params Global parameters controlling thresholds and resolution.
 * @return Collected refinement statistics.
 */
RefinementStats refine_delaunay(Delaunay &dt,
                                const UnifiedGrid &grid,
                                const std::vector<Cube> &active_cubes,
                                const VdcParam &params);

#endif // VDC_REFINEMENT_H
