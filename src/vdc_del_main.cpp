//! @file vdc_del_main.cpp
//! @brief Main entry point for the Delaunay-based VDC isosurface extraction.
//!
//! This program implements the Delaunay-based algorithm for isosurface extraction
//! as described in vdc-DelaunayBased.txt. It works directly in the Delaunay domain
//! rather than constructing a full Voronoi diagram.

#include "core/vdc_commandline.h"
#include "core/vdc_timing.h"
#include "core/vdc_debug.h"
#include "processing/vdc_grid.h"
#include "processing/vdc_func.h"
#include "processing/vdc_del_isosurface.h"
#include "processing/vdc_sep_isov.h"
#include "vdc_io.h"

#include <iostream>
#include <string>

//! @brief Print program header
void print_header() {
    std::cout << "========================================\n";
    std::cout << " vdc-del: Delaunay-based Dual Contouring\n";
    std::cout << "========================================\n\n";
}

//! @brief Main entry point
int main(int argc, char* argv[]) {
    print_header();

    // ========================================================================
    // Step 1: Parse command-line arguments
    // ========================================================================
    VdcParam param;
    parse_arguments(argc, argv, param);

    TimingStats::getInstance().startTimer("Total", "");

    std::cout << "Input file: " << param.file_path << "\n";
    std::cout << "Isovalue: " << param.isovalue << "\n";
    std::cout << "Output format: " << param.output_format << "\n\n";

    // ========================================================================
    // Step 2: Load grid data
    // ========================================================================
    TimingStats::getInstance().startTimer("LoadGrid", "Total");
    std::cout << "Loading grid data...\n";

    UnifiedGrid grid = load_nrrd_data(param.file_path);

    std::cout << "  Grid dimensions: "
              << grid.num_cells[0] << " x "
              << grid.num_cells[1] << " x "
              << grid.num_cells[2] << "\n";
    std::cout << "  Grid spacing: "
              << grid.spacing[0] << " x "
              << grid.spacing[1] << " x "
              << grid.spacing[2] << "\n";

    TimingStats::getInstance().stopTimer("LoadGrid", "Total");

    // ========================================================================
    // Step 3: Optional supersampling
    // ========================================================================
    if (param.supersample && param.supersample_r > 1) {
        TimingStats::getInstance().startTimer("Supersample", "Total");
        std::cout << "\nSupersampling by factor " << param.supersample_r << "...\n";

        grid = supersample_grid(grid, param.supersample_r);

        std::cout << "  New grid dimensions: "
                  << grid.num_cells[0] << " x "
                  << grid.num_cells[1] << " x "
                  << grid.num_cells[2] << "\n";

        TimingStats::getInstance().stopTimer("Supersample", "Total");
    }

    if (grid.boundary_crosses_isovalue(param.isovalue)) {
        grid.zero_boundary_shell();
        std::cout << "[INFO] Clamped boundary voxels to minimum scalar value.\n";
    }

    // ========================================================================
    // Step 4: Find active cubes
    // ========================================================================
    TimingStats::getInstance().startTimer("FindActiveCubes", "Total");
    std::cout << "\nFinding active cubes...\n";

    std::vector<Cube> activeCubes;
    find_active_cubes(grid, param.isovalue, activeCubes);

    std::cout << "  Found " << activeCubes.size() << " active cubes\n";

    TimingStats::getInstance().stopTimer("FindActiveCubes", "Total");

    if (activeCubes.empty()) {
        std::cerr << "Error: No active cubes found at isovalue " << param.isovalue << "\n";
        return 1;
    }

    // ========================================================================
    // Step 4b: Optional separation (if -sep_dist/-sep_split specified)
    // ========================================================================
    if (param.sep) {
        TimingStats::getInstance().startTimer("Separation", "Total");
        std::cout << "\nSeparating active cubes (sep_dist=" << param.sep_dist
                  << ", sep_split=" << param.sep_split << ")...\n";

        std::vector<Cube> separated = separate_active_cubes(
            activeCubes, grid, param.isovalue, param.sep_dist, param.sep_split);
        activeCubes.swap(separated);

        std::cout << "  Remaining cubes after separation: " << activeCubes.size() << "\n";
        TimingStats::getInstance().stopTimer("Separation", "Total");

        if (activeCubes.empty()) {
            std::cerr << "Error: Separation removed all active cubes.\n";
            return 1;
        }
    }

    // ========================================================================
    // Step 5: Construct Delaunay triangulation
    // ========================================================================
    TimingStats::getInstance().startTimer("Delaunay", "Total");
    std::cout << "\nConstructing Delaunay triangulation...\n";

    // Create grid facets for dummy points at boundaries
    std::vector<std::vector<GridFacets>> grid_facets = create_grid_facets(activeCubes);

    // Get cube centers for Delaunay vertices
    std::vector<Point> cubeCenters = get_cube_centers(activeCubes);

    // Use accurate iso-crossing points when explicitly requested or when
    // separation uses a refined subgrid (keeps site selection consistent).
    const bool use_isov_sites = param.position_delv_on_isov || (param.sep && param.sep_split > 0);
    if (use_isov_sites) {
        cubeCenters = get_cube_accurate_iso_crossing_points(activeCubes);
    }

    // Construct Delaunay triangulation
    Delaunay dt;
    construct_delaunay_triangulation(dt, grid, grid_facets, param, cubeCenters);

    std::cout << "  Delaunay vertices: " << dt.number_of_vertices() << "\n";
    std::cout << "  Delaunay cells: " << dt.number_of_finite_cells() << "\n";

    // Store per-site isosurface sample points in vertex info.
    // `info().index` is assigned by `construct_delaunay_triangulation` and matches
    // the active-cube list order for non-dummy vertices.
    const std::vector<Point> isoCrossings = get_cube_accurate_iso_crossing_points(activeCubes);
    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
        if (vit->info().is_dummy) {
            continue;
        }
        const int site_idx = vit->info().index;
        if (site_idx >= 0 && static_cast<size_t>(site_idx) < isoCrossings.size()) {
            vit->info().isov = isoCrossings[static_cast<size_t>(site_idx)];
        } else {
            // Fallback (should not happen): use the Delaunay site itself.
            vit->info().isov = vit->point();
        }
    }

    TimingStats::getInstance().stopTimer("Delaunay", "Total");

    // ========================================================================
    // Step 6: Compute cell circumcenters and scalar values
    // ========================================================================
    TimingStats::getInstance().startTimer("CellProcessing", "Total");
    std::cout << "\nComputing cell circumcenters and scalars...\n";

    compute_cell_circumcenters_and_scalars(dt, grid, param.isovalue);

    TimingStats::getInstance().stopTimer("CellProcessing", "Total");

    // ========================================================================
    // Step 7: Mark isosurface facets
    // ========================================================================
    TimingStats::getInstance().startTimer("MarkFacets", "Total");
    std::cout << "Marking isosurface facets...\n";

    mark_isosurface_facets(dt);

    TimingStats::getInstance().stopTimer("MarkFacets", "Total");

    // ========================================================================
    // Step 8: Mark active vertices
    // ========================================================================
    TimingStats::getInstance().startTimer("MarkVertices", "Total");
    std::cout << "Marking active vertices...\n";

    mark_active_vertices(dt);

    TimingStats::getInstance().stopTimer("MarkVertices", "Total");

    // ========================================================================
    // Step 9: Compute facet cycles around active vertices
    // ========================================================================
    TimingStats::getInstance().startTimer("CycleDetection", "Total");
    std::cout << "\nComputing facet cycles...\n";

    compute_facet_cycles(dt);

    // Optional mod_cyc pass to fix non-manifolds from problematic matchings
    if (param.mod_cyc) {
        modify_cycles_pass(dt);
    }

    TimingStats::getInstance().stopTimer("CycleDetection", "Total");

    // ========================================================================
    // Step 10: Compute cycle isovertices
    // ========================================================================
    TimingStats::getInstance().startTimer("IsovertexComputation", "Total");
    std::cout << "Computing isovertex positions...\n";

    compute_cycle_isovertices(dt, grid, param.isovalue, param.position_delv_on_isov);

    TimingStats::getInstance().stopTimer("IsovertexComputation", "Total");

    // ========================================================================
    // Step 11: Generate triangles
    // ========================================================================
    TimingStats::getInstance().startTimer("TriangleGeneration", "Total");
    std::cout << "\nGenerating isosurface triangles...\n";

    DelaunayIsosurface iso_surface;
    generate_isosurface_triangles(dt, iso_surface);

    std::cout << "  Generated " << iso_surface.num_vertices() << " vertices\n";
    std::cout << "  Generated " << iso_surface.num_triangles() << " triangles\n";

    TimingStats::getInstance().stopTimer("TriangleGeneration", "Total");

    // ========================================================================
    // Step 12: Output mesh
    // ========================================================================
    TimingStats::getInstance().startTimer("Output", "Total");

    // Determine output filename
    std::string output_filename = param.output_filename;
    if (output_filename.empty()) {
        // Generate default filename from input
        size_t last_slash = param.file_path.find_last_of("/\\");
        size_t last_dot = param.file_path.find_last_of(".");
        std::string basename = param.file_path.substr(
            last_slash == std::string::npos ? 0 : last_slash + 1,
            last_dot == std::string::npos ? std::string::npos : last_dot - (last_slash == std::string::npos ? 0 : last_slash + 1)
        );
        output_filename = basename + "_del." + param.output_format;
    }

    std::cout << "\nWriting output to " << output_filename << "...\n";

    if (param.output_format == "ply") {
        write_ply_delaunay(output_filename, iso_surface);
    } else {
        write_off_delaunay(output_filename, iso_surface);
    }

    TimingStats::getInstance().stopTimer("Output", "Total");

    // ========================================================================
    // Print timing statistics if requested
    // ========================================================================
    TimingStats::getInstance().stopTimer("Total", "");

    if (param.timing_stats) {
        std::cout << "\n";
        TimingStats::getInstance().printReport();
    }

    std::cout << "\nDone!\n";
    return 0;
}
