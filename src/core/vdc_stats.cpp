//! @file vdc_stats.cpp
//! @brief Implementation of statistics collection and reporting for vdc-del.

#include "core/vdc_stats.h"
#include "processing/vdc_del_isosurface.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>

SummaryStats collect_summary_stats_delaunay(
    const std::vector<Cube> &activeCubes,
    const Delaunay &dt,
    const DelaunayIsosurface &iso_surface)
{
    SummaryStats stats;

    stats.active_cubes = activeCubes.size();
    stats.delaunay_vertices = static_cast<std::size_t>(
        std::distance(dt.finite_vertices_begin(), dt.finite_vertices_end()));
    stats.delaunay_cells = static_cast<std::size_t>(
        std::distance(dt.finite_cells_begin(), dt.finite_cells_end()));

    // Count active vertices, isosurface facets, and cycles
    std::size_t active_count = 0;
    std::size_t total_cycles = 0;
    std::size_t multi_cycle_count = 0;

    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
        if (vit->info().active) {
            active_count++;
            std::size_t num_cycles = vit->info().facet_cycles.size();
            total_cycles += num_cycles;
            if (num_cycles > 1) {
                multi_cycle_count++;
            }
        }
    }

    // Count isosurface facets
    std::size_t isosurface_facet_count = 0;
    for (auto cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit) {
        for (int i = 0; i < 4; ++i) {
            if (cit->info().facet_is_isosurface[i]) {
                isosurface_facet_count++;
            }
        }
    }

    stats.active_vertices = active_count;
    stats.isosurface_facets = isosurface_facet_count;
    stats.total_cycles = total_cycles;
    stats.multi_cycle_vertices = multi_cycle_count;
    stats.iso_vertices = iso_surface.isovertices.size();
    stats.iso_triangles = iso_surface.triangles.size();

    return stats;
}

void print_summary_report(const SummaryStats &stats)
{
    std::cout << "\n[SUMMARY] Run statistics (Delaunay-based)\n";
    std::cout << "  Active cubes: " << stats.active_cubes << "\n";
    std::cout << "  Delaunay finite vertices: " << stats.delaunay_vertices
              << ", finite cells: " << stats.delaunay_cells << "\n";
    std::cout << "  Active Delaunay vertices: " << stats.active_vertices << "\n";
    std::cout << "  Isosurface facets: " << stats.isosurface_facets << "\n";
    std::cout << "  Total facet cycles: " << stats.total_cycles
              << " (multi-cycle vertices: " << stats.multi_cycle_vertices << ")\n";
    std::cout << "  Iso-surface vertices: " << stats.iso_vertices
              << ", triangles: " << stats.iso_triangles << "\n";
}
