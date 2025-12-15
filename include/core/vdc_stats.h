//! @file vdc_stats.h
//! @brief Header file for collecting and reporting summary statistics.
//!
//! Simplified version for vdc-del that doesn't depend on Voronoi structures.

#ifndef VDC_STATS_H
#define VDC_STATS_H

#include "core/vdc_type.h"
#include "processing/vdc_grid.h"
#include <array>
#include <cstddef>
#include <vector>

// Forward declaration
struct DelaunayIsosurface;

//! @brief Structure to hold summary statistics for the Delaunay-based pipeline.
struct SummaryStats
{
    std::size_t active_cubes = 0;                   //!< Number of active cubes
    std::size_t delaunay_vertices = 0;              //!< Number of vertices in the Delaunay triangulation
    std::size_t delaunay_cells = 0;                 //!< Number of cells in the Delaunay triangulation
    std::size_t active_vertices = 0;                //!< Number of active Delaunay vertices
    std::size_t isosurface_facets = 0;              //!< Number of Delaunay facets marked as isosurface
    std::size_t total_cycles = 0;                   //!< Total number of facet cycles
    std::size_t multi_cycle_vertices = 0;           //!< Number of vertices with multiple cycles
    std::size_t iso_vertices = 0;                   //!< Number of vertices in the extracted isosurface
    std::size_t iso_triangles = 0;                  //!< Number of triangles in the extracted isosurface
};

//! @brief Collects summary statistics from the Delaunay-based pipeline.
/*!
 * @param activeCubes Vector of active cubes in the grid
 * @param dt The Delaunay triangulation
 * @param iso_surface The extracted isosurface
 * @return SummaryStats structure containing all collected statistics
 */
SummaryStats collect_summary_stats_delaunay(
    const std::vector<Cube> &activeCubes,
    const Delaunay &dt,
    const DelaunayIsosurface &iso_surface
);

//! @brief Prints a formatted summary report of collected statistics.
/*!
 * @param stats The SummaryStats structure containing statistics to report
 */
void print_summary_report(const SummaryStats &stats);

#endif // VDC_STATS_H
