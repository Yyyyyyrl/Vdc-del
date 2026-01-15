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
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
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

static std::string format_float_for_path(float v) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(6) << static_cast<double>(v);
    std::string s = ss.str();
    while (!s.empty() && s.back() == '0') {
        s.pop_back();
    }
    if (!s.empty() && s.back() == '.') {
        s.pop_back();
    }
    if (s.empty()) {
        s = "0";
    }
    return s;
}

static std::string trace_config_tag(const VdcParam& param) {
    std::vector<std::string> parts;

    if (param.supersample) {
        parts.push_back("sup" + std::to_string(param.supersample_r));
    }

    if (param.sep) {
        std::ostringstream ss;
        ss << "sep" << param.sep_split << "_dist" << param.sep_dist;
        parts.push_back(ss.str());
    }

    if (!param.multi_isov) {
        parts.push_back("single_isov");
    }
    if (param.position_delv_on_isov) {
        parts.push_back("delv_on_isov");
    }
    if (!param.mod_cyc) {
        parts.push_back("no_modcyc");
    }

    if (parts.empty()) {
        return "default";
    }

    std::ostringstream out;
    for (size_t i = 0; i < parts.size(); ++i) {
        if (i > 0) {
            out << "_";
        }
        out << parts[i];
    }
    return out.str();
}

/**
 * @brief Write Delaunay triangulation boundary facets to an OFF file.
 *
 * Outputs all finite triangular facets (boundary of the Delaunay tetrahedra).
 * This lets users compare the isosurface mesh against the Delaunay structure.
 *
 * @param dt The Delaunay triangulation.
 * @param filename Output filename.
 * @param has_bbox If true, filter to only include facets within the bounding box.
 * @param bbox_min Minimum corner of bounding box (x, y, z).
 * @param bbox_max Maximum corner of bounding box (x, y, z).
 */
static void write_delv_off(
    const Delaunay& dt,
    const std::string& filename,
    bool has_bbox = false,
    const double* bbox_min = nullptr,
    const double* bbox_max = nullptr
) {
    // Helper to check if a point is within the bounding box
    auto in_bbox = [&](const Point& p) -> bool {
        if (!has_bbox || !bbox_min || !bbox_max) return true;
        return p.x() >= bbox_min[0] && p.x() <= bbox_max[0] &&
               p.y() >= bbox_min[1] && p.y() <= bbox_max[1] &&
               p.z() >= bbox_min[2] && p.z() <= bbox_max[2];
    };

    // Collect finite vertices that are in bbox and assign indices
    std::unordered_map<Vertex_handle, int> vertex_index;
    std::vector<Point> vertices;
    int idx = 0;
    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
        if (in_bbox(vit->point())) {
            vertex_index[vit] = idx++;
            vertices.push_back(vit->point());
        }
    }

    // Collect finite facets where all 3 vertices are in the bbox
    std::vector<std::array<int, 3>> facets;
    for (auto fit = dt.finite_facets_begin(); fit != dt.finite_facets_end(); ++fit) {
        Cell_handle c = fit->first;
        int i = fit->second;

        std::array<int, 3> tri;
        int ti = 0;
        bool all_in_bbox = true;
        for (int j = 0; j < 4; ++j) {
            if (j == i) continue;
            Vertex_handle vh = c->vertex(j);
            if (dt.is_infinite(vh)) goto skip_facet;
            auto it = vertex_index.find(vh);
            if (it == vertex_index.end()) {
                all_in_bbox = false;
                break;
            }
            tri[ti++] = it->second;
        }
        if (!all_in_bbox) goto skip_facet;
        // Sort triangle indices for consistent orientation
        if (i % 2 == 0) {
            std::swap(tri[0], tri[1]);
        }
        facets.push_back(tri);
        skip_facet:;
    }

    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: Failed to open " << filename << " for writing.\n";
        return;
    }

    out << "OFF\n";
    out << vertices.size() << " " << facets.size() << " 0\n";
    out << std::setprecision(17);
    for (const auto& p : vertices) {
        out << p.x() << " " << p.y() << " " << p.z() << "\n";
    }
    for (const auto& f : facets) {
        out << "3 " << f[0] << " " << f[1] << " " << f[2] << "\n";
    }

    std::cout << "Wrote Delaunay triangulation to: " << filename << " ("
              << vertices.size() << " vertices, " << facets.size() << " facets)";
    if (has_bbox) {
        std::cout << " [cropped to bbox]";
    }
    std::cout << "\n";
}

static void dump_duplicate_isovertex_positions(
    const std::filesystem::path& trace_root,
    const DelaunayIsosurface& iso_surface
) {
    if (trace_root.empty()) {
        return;
    }

    struct Key {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;

        bool operator==(const Key& other) const noexcept {
            return x == other.x && y == other.y && z == other.z;
        }
    };

    struct KeyHash {
        size_t operator()(const Key& k) const noexcept {
            const size_t h1 = std::hash<double>{}(k.x);
            const size_t h2 = std::hash<double>{}(k.y);
            const size_t h3 = std::hash<double>{}(k.z);
            size_t h = h1;
            h ^= (h2 + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
            h ^= (h3 + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
            return h;
        }
    };

    std::unordered_map<Key, std::vector<int>, KeyHash> groups;
    groups.reserve(iso_surface.isovertices.size());

    // Quantize to ~8 decimal places (matching OFF file output precision)
    // This ensures positions that appear identical after rounding are grouped
    constexpr double QUANT = 1e8;

    const auto& s = iso_surface.vertex_scale;
    for (size_t i = 0; i < iso_surface.isovertices.size(); ++i) {
        const Point& p = iso_surface.isovertices[i];
        const Key k{
            std::round(p.x() * s[0] * QUANT) / QUANT,
            std::round(p.y() * s[1] * QUANT) / QUANT,
            std::round(p.z() * s[2] * QUANT) / QUANT,
        };
        groups[k].push_back(static_cast<int>(i));
    }

    std::vector<std::pair<Key, std::vector<int>>> dups;
    dups.reserve(32);
    for (auto& kv : groups) {
        if (kv.second.size() > 1) {
            dups.emplace_back(kv.first, std::move(kv.second));
        }
    }

    std::sort(dups.begin(), dups.end(), [](const auto& a, const auto& b) {
        return a.second.size() > b.second.size();
    });

    std::error_code ec;
    std::filesystem::create_directories(trace_root, ec);

    std::ofstream out(trace_root / "global_duplicate_isovertex_positions.txt");
    if (!out) {
        return;
    }

    out << std::setprecision(17);
    out << "count_groups " << dups.size() << "\n";

    if (dups.empty()) {
        return;
    }

    for (const auto& g : dups) {
        out << "pos " << g.first.x << " " << g.first.y << " " << g.first.z
            << " count " << g.second.size() << "\n";
        for (const int isov_idx : g.second) {
            const int delv = (isov_idx >= 0 && static_cast<size_t>(isov_idx) < iso_surface.isovertex_delaunay_vertex.size())
                                 ? iso_surface.isovertex_delaunay_vertex[static_cast<size_t>(isov_idx)]
                                 : -1;
            const int cyc = (isov_idx >= 0 && static_cast<size_t>(isov_idx) < iso_surface.isovertex_cycle_index.size())
                                 ? iso_surface.isovertex_cycle_index[static_cast<size_t>(isov_idx)]
                                 : -1;
            out << "  isov " << isov_idx << " delv " << delv << " cycle " << cyc << "\n";
        }
    }
}

//! @brief Main entry point
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
        write_delv_off(dt, delv_filename,
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

    // ========================================================================
    // Step 7: Mark isosurface facets
    // ========================================================================
    TimingStats::getInstance().startTimer("MarkFacets", "Total");
    if (!param.terse) {
        std::cout << "Marking isosurface facets...\n";
    }

    mark_isosurface_facets(dt);

    TimingStats::getInstance().stopTimer("MarkFacets", "Total");

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
        // Root trace outputs at repo root, independent of the current working directory.
        std::filesystem::path repo_root;
        {
            std::error_code ec;
            std::filesystem::path exe_path = std::filesystem::weakly_canonical(argv[0], ec);
            if (ec || exe_path.empty()) {
                ec.clear();
                exe_path = std::filesystem::absolute(argv[0], ec);
            }

            // Typical dev layout: <repo>/{build,build-local}/vdc-del
            if (!exe_path.empty()) {
                repo_root = exe_path.parent_path().parent_path();
                if (!std::filesystem::exists(repo_root / "CMakeLists.txt")) {
                    repo_root.clear();
                }
            }
            if (repo_root.empty()) {
                repo_root = std::filesystem::current_path();
            }
        }

        const std::filesystem::path base = "simple_multi_failures_trace";
        const std::filesystem::path input_path(param.file_path);
        const std::string dataset = input_path.stem().string().empty() ? "dataset" : input_path.stem().string();
        const std::string iso = format_float_for_path(param.isovalue);
        const std::string cfg = trace_config_tag(param);
        const std::filesystem::path dir = base / dataset / ("iso" + iso + "_" + cfg);

        // Clear any existing dumps so repeated runs don't accumulate stale cases.
        {
            std::error_code rm_ec;
            std::filesystem::remove_all(dir / "local", rm_ec);
            rm_ec.clear();
            std::filesystem::remove_all(dir / "final", rm_ec);
        }

        std::error_code ec;
        std::filesystem::create_directories(dir, ec);
        if (ec) {
            std::cerr << "[DEL-SELFI-TRACE] Warning: failed to create directory: " << dir
                      << " (" << ec.message() << ")\n";
        } else {
            multi_isov_trace_dir = dir.string();
            std::cerr << "[DEL-SELFI-TRACE] Dumping multi-isov trace under: " << multi_isov_trace_dir
                      << " (subdirs: local/, final/)\n";
        }
    }

    const CycleIsovertexOptions isovertex_options{
        param.position_multi_isov_on_delv,
        param.multi_isov_trace,
        param.foldover,
        multi_isov_trace_dir,
    };

    compute_cycle_isovertices(
        dt, grid, param.isovalue, param.position_delv_on_isov, isovertex_options);

    TimingStats::getInstance().stopTimer("IsovertexComputation", "Total");

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

    if (param.multi_isov_trace && !multi_isov_trace_dir.empty()) {
        dump_duplicate_isovertex_positions(multi_isov_trace_dir, iso_surface);
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
