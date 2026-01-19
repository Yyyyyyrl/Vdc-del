//! @file vdc_debug.h
//! @brief Debugging utilities, logging, and dump helpers.

#ifndef VDC_DEBUG_H
#define VDC_DEBUG_H

#include "core/vdc_type.h"
#include "processing/vdc_grid.h"

#include <filesystem>
#include <iostream>
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
    const double* bbox_max = nullptr);

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
    const std::vector<Cell_handle>& cell_by_index);

void finalize_multi_isov_trace(
    const CycleIsovertexOptions& options,
    const std::unordered_set<int>& local_unresolved_dumped_vertices,
    const std::vector<Vertex_handle>& multi_cycle_vertices,
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index);

} // namespace vdc_debug

#endif // VDC_DEBUG_H
