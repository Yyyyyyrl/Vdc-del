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
#include "processing/vdc_refinement.h"
#include "processing/vdc_sep_isov.h"
#include "vdc_io.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

//! @brief Print program header
void print_header() {
    std::cout << "========================================\n";
    std::cout << " vdc-del: Delaunay-based Dual Contouring\n";
    std::cout << "========================================\n\n";
}

int main(int argc, char* argv[]) {
    // ========================================================================
    // Step 1: Parse command-line arguments
    // ========================================================================
    VdcParam param;
    parse_arguments(argc, argv, param);

    indicator = !param.terse;

    TimingStats::getInstance().startTimer("Total", "");

    if (!param.terse) {
        print_header();
        std::cout << "Input file: " << param.file_path << "\n";
        std::cout << "Isovalue: " << param.isovalue << "\n";
        std::cout << "Output format: " << param.output_format << "\n\n";
    }

    // ========================================================================
    // Step 2: Load grid data
    // ========================================================================
    TimingStats::getInstance().startTimer("LoadGrid", "Total");
    if (!param.terse) {
        std::cout << "Loading grid data...\n";
    }

    UnifiedGrid grid = load_nrrd_data(param.file_path);

    if (!param.terse) {
        std::cout << "  Grid dimensions: "
                  << grid.num_cells[0] << " x "
                  << grid.num_cells[1] << " x "
                  << grid.num_cells[2] << "\n";
        std::cout << "  Grid spacing: "
                  << grid.spacing[0] << " x "
                  << grid.spacing[1] << " x "
                  << grid.spacing[2] << "\n";
    }

    TimingStats::getInstance().stopTimer("LoadGrid", "Total");

    // ========================================================================
    // Step 3: Optional supersampling
    // ========================================================================
    if (param.supersample && param.supersample_r > 1) {
        TimingStats::getInstance().startTimer("Supersample", "Total");
        if (!param.terse) {
            std::cout << "\nSupersampling by factor " << param.supersample_r << "...\n";
        }

        grid = supersample_grid(grid, param.supersample_r);

        if (!param.terse) {
            std::cout << "  New grid dimensions: "
                      << grid.num_cells[0] << " x "
                      << grid.num_cells[1] << " x "
                      << grid.num_cells[2] << "\n";
        }

        TimingStats::getInstance().stopTimer("Supersample", "Total");
    }

    if (grid.boundary_crosses_isovalue(param.isovalue)) {
        grid.zero_boundary_shell();
        if (!param.terse) {
            std::cout << "[INFO] Clamped boundary voxels to minimum scalar value.\n";
        }
    }

    // ========================================================================
    // Step 4: Find active cubes
    // ========================================================================
    TimingStats::getInstance().startTimer("FindActiveCubes", "Total");
    if (!param.terse) {
        std::cout << "\nFinding active cubes...\n";
    }

    std::vector<Cube> activeCubes;
    find_active_cubes(grid, param.isovalue, activeCubes);
    std::vector<Cube> activeCubes_pre_sep;
    if (!param.dump_site_selection_json.empty()) {
        activeCubes_pre_sep = activeCubes;
    }

    if (!param.terse) {
        std::cout << "  Found " << activeCubes.size() << " active cubes\n";
    }

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
        if (!param.terse) {
            std::cout << "\nSeparating active cubes (sep_dist=" << param.sep_dist
                      << ", sep_split=" << param.sep_split << ")...\n";
        }

        std::vector<Cube> separated = separate_active_cubes(
            activeCubes, grid, param.isovalue, param.sep_dist, param.sep_split);
        activeCubes.swap(separated);

        if (!param.terse) {
            std::cout << "  Remaining cubes after separation: " << activeCubes.size() << "\n";
        }
        TimingStats::getInstance().stopTimer("Separation", "Total");

        if (activeCubes.empty()) {
            std::cerr << "Error: Separation removed all active cubes.\n";
            return 1;
        }
    }
    if (!param.dump_site_selection_json.empty()) {
        const std::vector<Cube>& pre = activeCubes_pre_sep.empty() ? activeCubes : activeCubes_pre_sep;
        vdc_debug::dump_site_selection_json(
            param.dump_site_selection_json,
            grid,
            param.isovalue,
            param.sep_dist,
            param.sep_split,
            pre,
            activeCubes);
    }

    // ========================================================================
    // Step 5: Construct Delaunay triangulation
    // ========================================================================
    TimingStats::getInstance().startTimer("Delaunay", "Total");
    if (!param.terse) {
        std::cout << "\nConstructing Delaunay triangulation...\n";
    }

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

    TimingStats::getInstance().stopTimer("Delaunay", "Total");

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
            vit->info().has_isov_sample = true;
        }
    }

    if (param.refine_small_angles) {
        TimingStats::getInstance().startTimer("Refinement", "Total");
        RefinementStats refine_stats = refine_delaunay(dt, grid, activeCubes, param);
        TimingStats::getInstance().stopTimer("Refinement", "Total");

        if (!param.terse) {
            std::cout << "  Refinement iterations: " << refine_stats.iterations_run
                      << ", candidates: " << refine_stats.candidate_facets
                      << ", bipolar facets: " << refine_stats.bipolar_facets
                      << ", inserted sites: " << refine_stats.inserted_points
                      << "\n";
        }
    }

    if (!param.terse) {
        std::cout << "  Delaunay vertices: " << dt.number_of_vertices() << "\n";
        std::cout << "  Delaunay cells: " << dt.number_of_finite_cells() << "\n";
    }

    // Optional: Output Delaunay triangulation to OFF file
    if (param.out_delv) {
        // Derive filename from mesh output: delv_<basename>.off or delv_<basename>_crop.off
        std::string delv_filename;
        {
            std::filesystem::path mesh_path(param.output_filename);
            std::string stem = mesh_path.stem().string();
            if (stem.empty()) stem = "mesh";
            if (param.out_delv_has_bbox) {
                delv_filename = "delv_" + stem + "_crop.off";
            } else {
                delv_filename = "delv_" + stem + ".off";
            }
        }
        vdc_debug::write_delv_off(
            dt,
            delv_filename,
            param.out_delv_has_bbox,
            param.out_delv_bbox_min,
            param.out_delv_bbox_max);
    }

    // ========================================================================
    // Step 6: Compute cell circumcenters and scalar values
    // ========================================================================
    TimingStats::getInstance().startTimer("CellProcessing", "Total");
    if (!param.terse) {
        std::cout << "\nComputing cell circumcenters and scalars...\n";
    }

    compute_cell_circumcenters_and_scalars(dt, grid, param.isovalue);

    TimingStats::getInstance().stopTimer("CellProcessing", "Total");

    //Flip cell signs to remove very small dihedral angles between isosurface facets.
    if (param.flip_small_dihedral_cells) {
        TimingStats::getInstance().startTimer("CellProcessing", "SmallDihedralFlip");
        if (!param.terse) {
            std::cout << "Flipping small-dihedral cells...\n";
        }
        // Convert angle threshold from degrees to cosine (done once here for efficiency)
        const double cos_threshold = std::cos(param.dihedral_angle_threshold_deg * M_PI / 180.0);
        const int flipped = flip_cell_signs_for_small_isofacet_dihedral(
            dt, cos_threshold);
        if (!param.terse) {
            std::cout << "  Flipped " << flipped << " cells.\n";
        }
        TimingStats::getInstance().stopTimer("CellProcessing", "SmallDihedralFlip");
    }

    // ========================================================================
    // Step 7: Mark isosurface facets
    // ========================================================================
    TimingStats::getInstance().startTimer("MarkFacets", "Total");
    if (!param.terse) {
        std::cout << "Marking isosurface facets...\n";
    }

    mark_isosurface_facets(dt);

    TimingStats::getInstance().stopTimer("MarkFacets", "Total");

    if (!param.dump_facet_example_json.empty()) {
        vdc_debug::dump_isosurface_facet_example_json(param.dump_facet_example_json, param, dt);
    }

    // ========================================================================
    // Step 8: Mark active vertices
    // ========================================================================
    TimingStats::getInstance().startTimer("MarkVertices", "Total");
    if (!param.terse) {
        std::cout << "Marking active vertices...\n";
    }

    mark_active_vertices(dt);

    TimingStats::getInstance().stopTimer("MarkVertices", "Total");

    // ========================================================================
    // Step 9: Compute facet cycles around active vertices
    // ========================================================================
    TimingStats::getInstance().startTimer("CycleDetection", "Total");
    if (!param.terse) {
        std::cout << "\nComputing facet cycles...\n";
    }

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
    if (!param.terse) {
        std::cout << "Computing isovertex positions...\n";
    }

    std::string multi_isov_trace_dir;
    if (param.multi_isov_trace) {
        multi_isov_trace_dir = vdc_debug::prepare_multi_isov_trace_dir(argv[0], param);
    }

    const CycleIsovertexOptions isovertex_options{
        param.position_multi_isov_on_delv,
        param.multi_isov_trace,
        param.foldover,
        param.use_sep_dir,
        multi_isov_trace_dir,
    };

    compute_cycle_isovertices(
        dt, grid, param.isovalue, param.position_delv_on_isov, isovertex_options);

    TimingStats::getInstance().stopTimer("IsovertexComputation", "Total");

    if (!param.dump_multicycle_json.empty()) {
        vdc_debug::dump_multicycle_vertex_example_json(param.dump_multicycle_json, param, dt);
    }

    // ========================================================================
    // Step 11: Generate triangles
    // ========================================================================
    TimingStats::getInstance().startTimer("TriangleGeneration", "Total");
    if (!param.terse) {
        std::cout << "\nGenerating isosurface triangles...\n";
    }

    DelaunayIsosurface iso_surface;
    generate_isosurface_triangles(dt, iso_surface);

    // Set vertex scale from physical spacing for correct output with non-uniform grids
    iso_surface.vertex_scale = {grid.physical_spacing[0], grid.physical_spacing[1], grid.physical_spacing[2]};

    if (!param.dump_isovertex_map.empty()) {
        vdc_debug::dump_isovertex_map_txt(param.dump_isovertex_map, iso_surface);
    }

    if (param.multi_isov_trace && !multi_isov_trace_dir.empty()) {
        vdc_debug::dump_duplicate_isovertex_positions(multi_isov_trace_dir, iso_surface);
    }

    // maybe_trace_isovertex_stars(argv[0], param, dt, iso_surface); // TODO: implement or remove

    const size_t num_vertices = iso_surface.num_vertices();
    const size_t num_triangles = iso_surface.num_triangles();

    if (!param.terse) {
        std::cout << "  Generated " << num_vertices << " vertices\n";
        std::cout << "  Generated " << num_triangles << " triangles\n";
    }

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

    if (!param.terse) {
        std::cout << "\nWriting output to " << output_filename << "...\n";
    }

    bool wrote_ok = false;
    if (param.output_format == "ply") {
        wrote_ok = write_ply_delaunay(output_filename, iso_surface);
    } else {
        wrote_ok = write_off_delaunay(output_filename, iso_surface);
    }

    TimingStats::getInstance().stopTimer("Output", "Total");

    // ========================================================================
    // Print timing statistics if requested
    // ========================================================================
    TimingStats::getInstance().stopTimer("Total", "");

    if (!wrote_ok) {
        return 1;
    }

    if (param.terse) {
        std::cout << num_vertices << " isosurface vertices. "
                  << num_triangles << " isosurface triangles.\n";
        std::cout << "Wrote output to file: " << output_filename << "\n";
        return 0;
    }

    std::cout << "Wrote output to file: " << output_filename << "\n";

    if (param.timing_stats) {
        std::cout << "\n";
        TimingStats::getInstance().printReport();
    }

    std::cout << "\nDone!\n";
    return 0;
}
