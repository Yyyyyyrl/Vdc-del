//! @file vdc_func.h
//! @brief Delaunay triangulation construction functions for vdc-del.
//!
//! This is a simplified version that only includes functions needed for
//! the Delaunay-based isosurface extraction algorithm.

#ifndef VDC_FUNC_H
#define VDC_FUNC_H

#include "core/vdc_type.h"
#include "core/vdc_utilities.h"
#include "core/vdc_commandline.h"
#include "core/vdc_debug.h"
#include "processing/vdc_grid.h"

//! @brief Constructs a Delaunay triangulation from active cube centers and boundary dummy points.
/*!
 * This function builds the Delaunay triangulation that forms the basis for
 * the dual Voronoi diagram. It adds:
 * - Active cube centers as regular vertices
 * - Dummy boundary points from grid facets
 *
 * @param dt The Delaunay triangulation to populate.
 * @param grid The unified grid containing the scalar field.
 * @param grid_facets Grid facets for generating boundary dummy points.
 * @param vdc_param Parameters controlling the construction.
 * @param activeCubeCenters Centers of active cubes to be inserted as Delaunay vertices.
 */
void construct_delaunay_triangulation(
    Delaunay &dt,
    UnifiedGrid &grid,
    const std::vector<std::vector<GridFacets>> &grid_facets,
    VdcParam &vdc_param,
    std::vector<Point> &activeCubeCenters
);

//! @brief Adds dummy points from a facet for Voronoi diagram bounding.
/*!
 * This function generates dummy points on the boundary of the domain
 * to ensure proper Voronoi cell bounding.
 *
 * @param facet The grid facet to generate dummy points from.
 * @param grid The unified grid.
 * @param supersample_multiplier Factor to scale dummy point spacing.
 * @return Vector of dummy points for this facet.
 */
std::vector<Point> add_dummy_from_facet(
    const GridFacets &facet,
    const UnifiedGrid &grid,
    double supersample_multiplier
);

#endif // VDC_FUNC_H
