//! @file vdc_sep_isov.cpp
//! @brief Implementation of parameterized isosurface vertex separation

#include "processing/vdc_sep_isov.h"
#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <vector>

// ============================================================================
// COMMON UTILITIES
// ============================================================================

//! @brief Compute minimum distance to bounding box boundary in grid space
static inline int min_distance_to_boundary(int i, int j, int k, const UnifiedGrid &grid) {
    const int nx = grid.num_cells[0];
    const int ny = grid.num_cells[1];
    const int nz = grid.num_cells[2];
    int dist_i = std::min(i, (nx - 2) - i);
    int dist_j = std::min(j, (ny - 2) - j);
    int dist_k = std::min(k, (nz - 2) - k);
    return std::min({dist_i, dist_j, dist_k});
}

//! @brief Compute linear index in a refined grid of factor `factor`
static inline int linear_cell_index_refined(int i, int j, int k, const UnifiedGrid &grid, int factor) {
    const int rx = factor * (grid.num_cells[0] - 1);
    const int ry = factor * (grid.num_cells[1] - 1);
    return k * rx * ry + j * rx + i;
}

//! @brief Decode subgrid index to local coordinates
static inline void decode_subgrid_loc(int subgrid_index, int factor, int loc[3]) {
    int index = subgrid_index;
    loc[0] = index % factor;
    index /= factor;
    loc[1] = index % factor;
    index /= factor;
    loc[2] = index;
}

// ============================================================================
// ACCURATE ISO-CROSSING POINT COMPUTATION
// ============================================================================

Point compute_iso_crossing_point_accurate(
    const UnifiedGrid &grid,
    int i, int j, int k,
    float isovalue)
{
    static const int cubeVertices[8][3] = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
        {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}
    };
    static const int cubeEdges[12][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0},
        {4, 5}, {5, 6}, {6, 7}, {7, 4},
        {0, 4}, {1, 5}, {2, 6}, {3, 7}
    };

    const float dx = grid.spacing[0];
    const float dy = grid.spacing[1];
    const float dz = grid.spacing[2];
    const float min_x = grid.min_coord[0];
    const float min_y = grid.min_coord[1];
    const float min_z = grid.min_coord[2];

    float base_x = i * dx + min_x;
    float base_y = j * dy + min_y;
    float base_z = k * dz + min_z;

    std::array<Point, 8> vertices;
    std::array<float, 8> scalarValues;
    for (int v = 0; v < 8; ++v) {
        vertices[v] = Point(
            base_x + cubeVertices[v][0] * dx,
            base_y + cubeVertices[v][1] * dy,
            base_z + cubeVertices[v][2] * dz);
        scalarValues[v] = grid.get_value(
            i + cubeVertices[v][0],
            j + cubeVertices[v][1],
            k + cubeVertices[v][2]);
    }

    std::vector<Point> intersectionPoints;
    for (const auto &edge : cubeEdges) {
        int idx1 = edge[0];
        int idx2 = edge[1];
        float val1 = scalarValues[idx1];
        float val2 = scalarValues[idx2];

        if ((val1 < isovalue && val2 >= isovalue) || (val1 >= isovalue && val2 < isovalue)) {
            Point intersect = interpolate(vertices[idx1], vertices[idx2], val1, val2, isovalue, grid);
            intersectionPoints.push_back(intersect);
        }
    }

    if (!intersectionPoints.empty()) {
        Point centroid(0, 0, 0);
        for (const auto &pt : intersectionPoints) {
            centroid = Point(
                centroid.x() + pt.x(),
                centroid.y() + pt.y(),
                centroid.z() + pt.z());
        }
        float n = static_cast<float>(intersectionPoints.size());
        return Point(centroid.x() / n, centroid.y() / n, centroid.z() / n);
    }

    return Point(
        base_x + 0.5f * dx,
        base_y + 0.5f * dy,
        base_z + 0.5f * dz);
}

// ============================================================================
// SUBGRID HELPERS
// ============================================================================

//! @brief Map iso-crossing location to subgrid index using refine factor (K+1)
static int determine_subgrid_index(
    const Point &isoCrossingPoint,
    const Cube &cube,
    const UnifiedGrid &grid,
    int factor,
    int loc_out[3])
{
    const float dx = grid.spacing[0];
    const float dy = grid.spacing[1];
    const float dz = grid.spacing[2];
    const float min_x = grid.min_coord[0];
    const float min_y = grid.min_coord[1];
    const float min_z = grid.min_coord[2];

    float rel_x = (isoCrossingPoint.x() - (cube.indices[0] * dx + min_x)) / dx;
    float rel_y = (isoCrossingPoint.y() - (cube.indices[1] * dy + min_y)) / dy;
    float rel_z = (isoCrossingPoint.z() - (cube.indices[2] * dz + min_z)) / dz;

    rel_x = std::max(0.0f, std::min(1.0f, static_cast<float>(rel_x)));
    rel_y = std::max(0.0f, std::min(1.0f, static_cast<float>(rel_y)));
    rel_z = std::max(0.0f, std::min(1.0f, static_cast<float>(rel_z)));

    loc_out[0] = static_cast<int>(rel_x * factor);
    loc_out[1] = static_cast<int>(rel_y * factor);
    loc_out[2] = static_cast<int>(rel_z * factor);

    if (loc_out[0] >= factor) loc_out[0] = factor - 1;
    if (loc_out[1] >= factor) loc_out[1] = factor - 1;
    if (loc_out[2] >= factor) loc_out[2] = factor - 1;

    return loc_out[0] + factor * loc_out[1] + factor * factor * loc_out[2];
}

//! @brief Allowed offsets per split count (binary-friendly fractions)
static std::vector<float> allowed_offsets(int sep_split) {
    if (sep_split <= 0) {
        return {0.5f};
    }
    if (sep_split == 1) {
        return {0.25f, 0.5f, 0.75f};
    }
    if (sep_split == 2 || sep_split == 3) {
        return {0.125f, 0.5f, 0.875f};
    }
    // K >= 4
    return {1.0f / 16.0f, 5.0f / 16.0f, 0.5f, 11.0f / 16.0f, 15.0f / 16.0f};
}

//! @brief Snap a relative offset to nearest allowed fraction
static float quantize_offset(float raw_offset, int sep_split) {
    const std::vector<float> offsets = allowed_offsets(sep_split);
    float snapped = offsets.front();
    float best_dist = std::abs(raw_offset - snapped);
    for (size_t idx = 1; idx < offsets.size(); ++idx) {
        float dist = std::abs(raw_offset - offsets[idx]);
        if (dist < best_dist || (std::abs(dist - best_dist) < 1e-6f && offsets[idx] < snapped)) {
            best_dist = dist;
            snapped = offsets[idx];
        }
    }
    return std::max(0.0f, std::min(1.0f, snapped));
}

//! @brief Compute cube center snapped to power-of-two fractions
static Point compute_snapped_center(
    const Cube &cube,
    const UnifiedGrid &grid,
    int sep_split,
    int factor)
{
    int loc[3];
    decode_subgrid_loc(cube.isov_subgrid_index, factor, loc);

    float base_x = cube.indices[0];
    float base_y = cube.indices[1];
    float base_z = cube.indices[2];

    float rel_x = (static_cast<float>(loc[0]) + 0.5f) / static_cast<float>(factor);
    float rel_y = (static_cast<float>(loc[1]) + 0.5f) / static_cast<float>(factor);
    float rel_z = (static_cast<float>(loc[2]) + 0.5f) / static_cast<float>(factor);

    rel_x = quantize_offset(rel_x, sep_split);
    rel_y = quantize_offset(rel_y, sep_split);
    rel_z = quantize_offset(rel_z, sep_split);

    const float dx = grid.spacing[0];
    const float dy = grid.spacing[1];
    const float dz = grid.spacing[2];
    const float min_x = grid.min_coord[0];
    const float min_y = grid.min_coord[1];
    const float min_z = grid.min_coord[2];

    return Point(
        (base_x + rel_x) * dx + min_x,
        (base_y + rel_y) * dy + min_y,
        (base_z + rel_z) * dz + min_z);
}

// ============================================================================
// UNIFIED SEPARATION
// ============================================================================

std::vector<Cube> separate_active_cubes(
    std::vector<Cube> &activeCubes,
    const UnifiedGrid &grid,
    float isovalue,
    int sep_dist,
    int sep_split)
{
    const int factor = std::max(1, sep_split + 1);
    const int clearance = std::max(0, sep_dist - 1);

    // Compute accurate iso-crossings and subgrid indices
    for (Cube &cube : activeCubes) {
        cube.accurateIsoCrossing = compute_iso_crossing_point_accurate(
            grid, cube.indices[0], cube.indices[1], cube.indices[2], isovalue);

        int loc[3];
        int sub_idx = determine_subgrid_index(cube.accurateIsoCrossing, cube, grid, factor, loc);
        cube.isov_subgrid_index = sub_idx;
    }

    // Sort by distance to boundary with deterministic tie-break
    std::sort(activeCubes.begin(), activeCubes.end(),
              [&grid](const Cube &a, const Cube &b) {
                  int dist_a = min_distance_to_boundary(a.indices[0], a.indices[1], a.indices[2], grid);
                  int dist_b = min_distance_to_boundary(b.indices[0], b.indices[1], b.indices[2], grid);
                  if (dist_a == dist_b) {
                      if (a.indices[2] != b.indices[2]) return a.indices[2] < b.indices[2];
                      if (a.indices[1] != b.indices[1]) return a.indices[1] < b.indices[1];
                      return a.indices[0] < b.indices[0];
                  }
                  return dist_a < dist_b;
              });

    const int refine_nx = factor * (grid.num_cells[0] - 1);
    const int refine_ny = factor * (grid.num_cells[1] - 1);
    const int refine_nz = factor * (grid.num_cells[2] - 1);

    std::unordered_set<int> selected_indices;
    selected_indices.reserve(activeCubes.size());
    std::vector<Cube> out;
    out.reserve(activeCubes.size());

    for (Cube &cube : activeCubes) {
        int loc[3];
        decode_subgrid_loc(cube.isov_subgrid_index, factor, loc);

        int refined_loc[3] = {
            factor * cube.indices[0] + loc[0],
            factor * cube.indices[1] + loc[1],
            factor * cube.indices[2] + loc[2]
        };

        bool conflict = false;
        for (int dk = -clearance; dk <= clearance && !conflict; ++dk) {
            for (int dj = -clearance; dj <= clearance && !conflict; ++dj) {
                for (int di = -clearance; di <= clearance && !conflict; ++di) {
                    if (di == 0 && dj == 0 && dk == 0) continue;

                    int neighbor_loc[3] = {
                        refined_loc[0] + di,
                        refined_loc[1] + dj,
                        refined_loc[2] + dk
                    };

                    if (neighbor_loc[0] < 0 || neighbor_loc[1] < 0 || neighbor_loc[2] < 0) continue;
                    if (neighbor_loc[0] >= refine_nx || neighbor_loc[1] >= refine_ny || neighbor_loc[2] >= refine_nz) continue;

                    int neighbor_idx = linear_cell_index_refined(
                        neighbor_loc[0], neighbor_loc[1], neighbor_loc[2], grid, factor);
                    if (selected_indices.find(neighbor_idx) != selected_indices.end()) {
                        conflict = true;
                    }
                }
            }
        }

        if (!conflict) {
            int my_idx = linear_cell_index_refined(refined_loc[0], refined_loc[1], refined_loc[2], grid, factor);
            cube.cubeCenter = compute_snapped_center(cube, grid, sep_split, factor);
            selected_indices.insert(my_idx);
            out.push_back(cube);
        }
    }

    return out;
}
