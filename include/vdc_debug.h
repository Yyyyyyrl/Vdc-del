//! @file vdc_debug.h
//! @brief Debugging utilities, logging, and dump helpers.

#ifndef VDC_DEBUG_H
#define VDC_DEBUG_H

#include "core/vdc_type.h"
#include "processing/vdc_grid.h"

#include <array>
#include <filesystem>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_set>
#include <vector>

// Global debug flags (kept for compatibility).
extern bool debug;
extern bool indicator;

#define DEBUG_PRINT(msg) do { if (debug) { std::cerr << msg << std::endl; } } while(0)

void print_facet(Facet f);
void print_cell(Delaunay::Cell c);
void write_dummy_points(UnifiedGrid& grid, std::vector<Point> dummy_points);

struct DelaunayIsosurface;
struct VdcParam;
struct CycleIsovertexOptions;

namespace vdc_debug {

void dump_site_selection_json(
    const std::filesystem::path& path,
    const UnifiedGrid& grid,
    float isovalue,
    int sep_dist,
    int sep_split,
    const std::vector<Cube>& pre_sep,
    const std::vector<Cube>& post_sep);

void dump_isosurface_facet_example_json(
    const std::filesystem::path& path,
    const VdcParam& param,
    const Delaunay& dt);

void dump_multicycle_vertex_example_json(
    const std::filesystem::path& path,
    const VdcParam& param,
    const Delaunay& dt);

void write_delv_off(
    const Delaunay& dt,
    const std::string& filename,
    bool has_bbox = false,
    const double* bbox_min = nullptr,
    const double* bbox_max = nullptr,
    const std::array<double, 3>& vertex_scale = std::array<double, 3>{1.0, 1.0, 1.0});

void dump_duplicate_isovertex_positions(
    const std::filesystem::path& trace_root,
    const DelaunayIsosurface& iso_surface);

void dump_isovertex_map_txt(
    const std::filesystem::path& path,
    const DelaunayIsosurface& iso_surface);

std::string prepare_multi_isov_trace_dir(const char* argv0, const VdcParam& param);

// Self-intersection trace toggles (env-driven).
bool trace_selfi_enabled(int vertex_index);
const std::string& trace_selfi_dump_dir();
bool log_unresolved_A_vertices();
bool log_selfi_within_cycle_vertices();
bool trace_selfi_intersection_details_enabled();

void maybe_dump_selfi_stage(
    const Delaunay& dt,
    Vertex_handle v,
    const std::vector<Point>& cycle_isovertices,
    const std::vector<Cell_handle>& cell_by_index,
    const char* stage);

void maybe_dump_selfi_cycle_metadata(
    const Delaunay& dt,
    Vertex_handle v,
    const std::vector<Cell_handle>& cell_by_index);

// simple_multi_failures-style dumps for multi-cycle failures.
void dump_simple_multi_failure_stage(
    const std::filesystem::path& out_dir,
    Vertex_handle v,
    const std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    const char* stage);

void dump_multi_isov_trace_case(
    const std::filesystem::path& out_dir,
    Vertex_handle v,
    const std::vector<Point>& baseline_positions,
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    bool use_sep_dir);

// -------------------------------------------------------------------------
// Move-cap trace infrastructure (moved from vdc_del_cycles.cpp)
// -------------------------------------------------------------------------

struct MoveCapTraceCandidate {
    int source_cell_index = -1;
    int source_facet_index = -1;
    int incident_slot = -1; // 0: positive-side cell, 1: mirror-side cell
    int delta_cell_index = -1;
    int t_facet_idx = -1;
    int tprime_facet_idx = -1;
    int v_local = -1;
    bool acute_dihedral = false;
    double dot_nt_ntp = std::numeric_limits<double>::quiet_NaN();
    bool strict_mode = false;
    int flag_preserve_orientation = -1; // 1 true, 0 false, -1 not evaluated
    double dot_dir_ntp_away = std::numeric_limits<double>::quiet_NaN();
    int direction_pass = -1; // 1 true, 0 false, -1 not evaluated
    bool accepted_for_bound = false;
    std::string reason = "uninitialized";
    double dist = std::numeric_limits<double>::quiet_NaN();
    bool has_tet_geometry = false;
    std::array<int, 4> tet_vertex_indices = {-1, -1, -1, -1};
    std::array<Point, 4> tet_points = {
        Point(0.0, 0.0, 0.0),
        Point(0.0, 0.0, 0.0),
        Point(0.0, 0.0, 0.0),
        Point(0.0, 0.0, 0.0)};
};

struct MoveCapTraceWriteResult {
    bool wrote = false;
    std::filesystem::path path;
    int preserve_true = 0;
    int candidates = 0;
};

struct CycleMoveCapBound;  // forward declaration (defined in vdc_del_cycles.cpp)

void fill_move_cap_trace_tet_geometry(
    const Delaunay& dt,
    const Cell_handle& Delta,
    int t_facet_idx,
    Vertex_handle w0,
    MoveCapTraceCandidate* out);

MoveCapTraceWriteResult write_move_cap_trace_file(
    Vertex_handle v,
    int cycle_idx,
    bool strict_enabled,
    bool directional_gate_enabled,
    const Vector3* move_dir,
    double direction_gate_eps,
    double maxdist,
    int acute_dihedral_constraints,
    int strict_preserve_hits,
    int strict_preserve_rejects,
    int direction_keep,
    int direction_skip,
    const std::vector<MoveCapTraceCandidate>& candidates,
    const std::string& trace_dir);

// -------------------------------------------------------------------------
// Sep-dir normals trace infrastructure (moved from vdc_del_cycles.cpp)
// -------------------------------------------------------------------------

struct FacetNormalDumpEntry {
    int cell_index = -1;
    int facet_index = -1;
    int cell_sign = 0; // +1 POS, -1 NEG
    int neighbor_cell_index = -1;
    int neighbor_facet_index = -1;
    int neighbor_sign = 0;
    int normal_points_to_cell_index = -1;
    int normal_points_to_sign = 0;
    int normal_points_away_cell_index = -1;
    int normal_points_away_sign = 0;
    bool normal_flip_applied = false;
    int facet_vidx[3] = {-1, -1, -1};
    Point facet_p[3];
    int opposite_vidx = -1;
    Point opposite_p;
    double nx = 0.0;
    double ny = 0.0;
    double nz = 0.0;
    double dot_to_opposite = 0.0;
    double dot_to_circumcenter = 0.0;
};

struct SepDirDumpContext {
    Vertex_handle v0;
    int cycle_index = -1;
    int num_cycles = 0;
    size_t cycle_facets_size = 0;
    size_t normals_used = 0;
    int circumcenter_across_count = 0;
    int normals_point_to_pos_count = 0;
    int normals_point_to_neg_count = 0;
    bool flip_normals_for_cycle = false;
    const std::vector<std::vector<char>>* is_separated = nullptr;
    const std::vector<int>* separation_fail_targets = nullptr;
    double cx = 0.0, cy = 0.0, cz = 0.0;
    double center_len = 0.0;
    int foldover_edge_candidates_ge4 = 0;
    // Crossover edge info fields
    bool crossover_found = false;
    int crossover_v_low = -1;
    int crossover_v_high = -1;
    int crossover_count = 0;
    bool crossover_incident_on_v0 = false;
    int crossover_other_vidx = -1;
    double crossover_dir_x = 0.0;
    double crossover_dir_y = 0.0;
    double crossover_dir_z = 0.0;
    const std::vector<Point>* normals_as_points = nullptr;
};

void dump_sep_dir_normals_trace(
    const std::string& trace_dir,
    const SepDirDumpContext& ctx,
    const std::vector<FacetNormalDumpEntry>& entries,
    const std::string& trace_suffix = "");

void finalize_multi_isov_trace(
    const CycleIsovertexOptions& options,
    const std::unordered_set<int>& local_unresolved_dumped_vertices,
    const std::vector<Vertex_handle>& multi_cycle_vertices,
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index);

} // namespace vdc_debug

#endif // VDC_DEBUG_H
