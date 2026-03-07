//! @file vdc_del_cycles.cpp
//! @brief Implementation of cycle detection for Delaunay-based isosurface extraction.
//!
//! This file implements Steps 8-9 of the Delaunay-based VDC algorithm:
//! - Step 8: Compute cycles of isosurface facets around active vertices
//! - Step 9: Compute isosurface vertex positions for each cycle

#include "processing/vdc_del_isosurface.h"
#include "processing/vdc_del_cycles.h"
#include "core/vdc_debug.h"
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <optional>
#include <array>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <limits>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <CGAL/Triangle_3.h>
#include <CGAL/intersections.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_points_d_traits_3.h>

// ============================================================================
// Helper: Build cell index lookup
// ============================================================================

std::vector<Cell_handle> build_cell_index_vector(const Delaunay& dt) {
    std::vector<Cell_handle> cells(static_cast<size_t>(dt.number_of_finite_cells()));
    for (auto cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit) {
        const int idx = cit->info().index;
        if (idx >= 0 && static_cast<size_t>(idx) < cells.size()) {
            cells[static_cast<size_t>(idx)] = cit;
        }
    }
    return cells;
}

static inline Cell_handle lookup_cell(const std::vector<Cell_handle>& cells, int idx) {
    if (idx < 0 || static_cast<size_t>(idx) >= cells.size()) {
        return Cell_handle();
    }
    return cells[static_cast<size_t>(idx)];
}

// ============================================================================
// Helpers: Matching + adjacency on a Delaunay edge
// ============================================================================
// Note: FacetKey, FacetKeyHash, EdgeKey, EdgeKeyHash are now declared in vdc_del_cycles.h

static int find_vertex_index_in_cell(const Cell_handle& cell, const Vertex_handle& vertex) {
    for (int i = 0; i < 4; ++i) {
        if (cell->vertex(i) == vertex) {
            return i;
        }
    }
    return -1;
}

static void write_facet_cycle_indices_for_vertex(
    Vertex_handle v_handle,
    const std::vector<Cell_handle>& cell_by_index
) {
    if (v_handle == Vertex_handle()) {
        return;
    }
    const auto& cycles = v_handle->info().facet_cycles;
    for (size_t c = 0; c < cycles.size(); ++c) {
        for (const auto& [cell_index, facet_index] : cycles[c]) {
            Cell_handle cell = lookup_cell(cell_by_index, cell_index);
            if (cell == Cell_handle()) {
                continue;
            }
            if (facet_index < 0 || facet_index >= 4) {
                continue;
            }
            const int v_local = find_vertex_index_in_cell(cell, v_handle);
            if (v_local < 0 || v_local >= 4) {
                continue;
            }
            const int anchor = (facet_index + 1) % 4;
            const int slot = (v_local - anchor + 4) % 4;
            if (slot < 0 || slot >= 3) {
                continue;
            }
            // Reuse the per-facet per-vertex-slot storage (originally dualCellEdgeIndex)
            // to cache cycle indices for the Delaunay-based pipeline.
            cell->info().facet_info[facet_index].dualCellEdgeIndex[slot] = static_cast<int>(c);
        }
    }
}

static Cell_handle cell_between_facet_and_next_facet(
    const Delaunay& dt,
    const Facet& facet0,
    const Facet& facet1
) {
    // Common cell between consecutive facets is either facet0.first or facet1.first.
    const Facet facet1_mirror = dt.mirror_facet(facet1);
    if (facet0.first == facet1_mirror.first) {
        return facet0.first;
    }
    return facet1.first;
}

static std::optional<FacetKey> oriented_isosurface_facet_key(const Delaunay& dt, const Facet& facet) {
    Cell_handle cell0 = facet.first;
    if (dt.is_infinite(cell0)) {
        return std::nullopt;
    }

    const int facet_index0 = facet.second;
    Cell_handle cell1 = cell0->neighbor(facet_index0);
    if (dt.is_infinite(cell1)) {
        return std::nullopt;
    }

    const bool cell0_pos = cell0->info().flag_positive;
    const bool cell1_pos = cell1->info().flag_positive;
    if (cell0_pos == cell1_pos) {
        return std::nullopt;
    }

    if (cell0_pos) {
        // represent an isosurface facet by (positive_cell, facet_index).
        return FacetKey{cell0->info().index, facet_index0};
    }

    // Mirror into the positive cell to match the convention.
    const Facet mirror = dt.mirror_facet(facet);
    if (dt.is_infinite(mirror.first) || !mirror.first->info().flag_positive) {
        return std::nullopt;
    }
    return FacetKey{mirror.first->info().index, mirror.second};
}

static bool facet_is_valid_for_cycle_matching(
    const Delaunay& dt,
    const FacetKey& key,
    const std::vector<Cell_handle>& cell_by_index
) {
    Cell_handle cell = lookup_cell(cell_by_index, key.cell_index);
    if (cell == Cell_handle()) {
        return false;
    }
    if (dt.is_infinite(cell)) {
        return false;
    }

    if (key.facet_index < 0 || key.facet_index >= 4) {
        return false;
    }

    if (!cell->info().facet_is_isosurface[key.facet_index]) {
        return false;
    }

    // Ignore any isosurface facets involving dummy Delaunay sites.
    for (int i = 0; i < 4; ++i) {
        if (i == key.facet_index) {
            continue;
        }
        Vertex_handle vh = cell->vertex(i);
        if (dt.is_infinite(vh) || vh->info().is_dummy) {
            return false;
        }
    }

    return true;
}

static std::vector<std::pair<FacetKey, FacetKey>> compute_edge_bipolar_matching(
    const Delaunay& dt,
    const Edge& edge,
    BIPOLAR_MATCH_METHOD method,
    const std::vector<Cell_handle>& cell_by_index
) {

    Delaunay::Facet_circulator fc_start = dt.incident_facets(edge);
    if (fc_start == Delaunay::Facet_circulator()) {
        return {};
    }

    std::vector<Facet> ring_facets;
    auto fc = fc_start;
    do {
        ring_facets.push_back(*fc);
        ++fc;
    } while (fc != fc_start);

    if (ring_facets.empty()) {
        return {};
    }

    std::vector<Cell_handle> between_cells(ring_facets.size());
    for (size_t i = 0; i < ring_facets.size(); ++i) {
        const Facet& curr = ring_facets[i];
        const Facet& next = ring_facets[(i + 1) % ring_facets.size()];
        between_cells[i] = cell_between_facet_and_next_facet(dt, curr, next);

        // Require a fully finite ring; otherwise matching is ambiguous/incomplete.
        if (dt.is_infinite(between_cells[i])) {
            return {};
        }
    }

    struct BipolarFacetSlot {
        FacetKey key;
        bool between_is_positive = false; // sign of the cell sector immediately after this facet
    };

    std::vector<BipolarFacetSlot> bipolar;
    bipolar.reserve(ring_facets.size());

    for (size_t i = 0; i < ring_facets.size(); ++i) {
        const auto key_opt = oriented_isosurface_facet_key(dt, ring_facets[i]);
        if (!key_opt.has_value()) {
            continue;
        }

        const FacetKey& key = key_opt.value();
        if (!facet_is_valid_for_cycle_matching(dt, key, cell_by_index)) {
            continue;
        }

        bipolar.push_back(BipolarFacetSlot{key, between_cells[i]->info().flag_positive});
    }

    if (bipolar.size() < 2 || (bipolar.size() % 2 != 0)) {
        return {};
    }



    bool desired_between_positive = false;
    switch (method) {
        case BIPOLAR_MATCH_METHOD::SEP_NEG:
            // Separate negative regions: between-cell should be POSITIVE
            // (Algorithm: cell0.flag_positive != flag_separate_negative → skip if not matching)
            desired_between_positive = true;
            break;
        case BIPOLAR_MATCH_METHOD::SEP_POS:
            // Separate positive regions: between-cell should be NEGATIVE
            desired_between_positive = false;
            break;
        case BIPOLAR_MATCH_METHOD::UNCONSTRAINED_MATCH:
            // Default to the first consecutive pairing.
            desired_between_positive = bipolar[0].between_is_positive;
            break;
        default:
            return {};
    }

    size_t start = 0;
    if (method != BIPOLAR_MATCH_METHOD::UNCONSTRAINED_MATCH) {
        bool found = false;
        for (size_t i = 0; i < bipolar.size(); ++i) {
            if (bipolar[i].between_is_positive == desired_between_positive) {
                start = i;
                found = true;
                break;
            }
        }
        if (!found) {
            return {};
        }
    }

    std::vector<std::pair<FacetKey, FacetKey>> pairs;
    pairs.reserve(bipolar.size() / 2);
    for (size_t k = 0; k < bipolar.size(); k += 2) {
        const FacetKey& a = bipolar[(start + k) % bipolar.size()].key;
        const FacetKey& b = bipolar[(start + k + 1) % bipolar.size()].key;
        pairs.push_back({a, b});
    }
    return pairs;
}

static bool edge_has_at_least_n_valid_isosurface_facets(
    const Delaunay& dt,
    const Edge& edge,
    const std::vector<Cell_handle>& cell_by_index,
    size_t min_count
) {
    Delaunay::Facet_circulator fc_start = dt.incident_facets(edge);
    if (fc_start == Delaunay::Facet_circulator()) {
        return false;
    }

    size_t count = 0;
    auto fc = fc_start;
    do {
        const auto key_opt = oriented_isosurface_facet_key(dt, *fc);
        if (key_opt.has_value() && facet_is_valid_for_cycle_matching(dt, key_opt.value(), cell_by_index)) {
            ++count;
            if (count >= min_count) {
                return true;
            }
        }
        ++fc;
    } while (fc != fc_start);

    return false;
}

// ============================================================================
// Step 8: Compute Facet Cycles (Algorithm-aligned implementation)
// ============================================================================
// Note: FacetMatchingData is now declared in vdc_del_cycles.h


// Helper: Get index (0, 1, or 2) of vertex v in the facet (cell, facet_idx)
static int get_vertex_index_in_facet(
    Cell_handle cell, 
    int facet_idx, 
    Vertex_handle v
) {
    // Facet facet_idx is opposite vertex facet_idx
    // The three vertices of the facet are at positions (facet_idx+1)%4, (facet_idx+2)%4, (facet_idx+3)%4
    for (int i = 0; i < 3; ++i) {
        int cell_vertex_idx = (facet_idx + 1 + i) % 4;
        if (cell->vertex(cell_vertex_idx) == v) {
            return i;  // Return 0, 1, or 2 (position within the 3 vertices of the facet)
        }
    }
    return -1;  // Not found
}

void compute_facet_cycles(Delaunay& dt) {


    const std::vector<Cell_handle> cell_by_index = build_cell_index_vector(dt);

    // =========================================================================
    // Step 1: Build global matching data for all facets
    // =========================================================================
    
    std::unordered_map<FacetKey, FacetMatchingData, FacetKeyHash> facet_matching;
    
    // Track which edges we've already processed (keyed by min,max vertex indices)
    std::unordered_set<EdgeKey, EdgeKeyHash> processed_edges;
    std::vector<Edge> incident_edges;
    
    // For each active Delaunay vertex v0, match facets on its incident edges
    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
        if (!vit->info().active) continue;
        if (vit->info().is_dummy) continue;
        
        const Vertex_handle v0 = vit;
        const int v0_idx = v0->info().index;
        
        // Get all incident edges
        incident_edges.clear();
        dt.incident_edges(v0, std::back_inserter(incident_edges));
        
        for (const Edge& e : incident_edges) {
            Vertex_handle v_a = e.first->vertex(e.second);
            Vertex_handle v_b = e.first->vertex(e.third);
            
            // One of v_a or v_b should be v0
            Vertex_handle v1 = (v_a == v0) ? v_b : v_a;
            if (dt.is_infinite(v1)) continue;
            if (v1->info().is_dummy) continue;
            
            const int v1_idx = v1->info().index;
            EdgeKey ekey = EdgeKey::FromVertexIndices(v0_idx, v1_idx);
            
            // Skip if already processed
            if (processed_edges.count(ekey)) continue;
            processed_edges.insert(ekey);
            
            // Get bipolar matching for this edge
            auto pairs = compute_edge_bipolar_matching(dt, e, BIPOLAR_MATCH_METHOD::SEP_NEG, cell_by_index);
            if (pairs.empty()) continue;
            
            // For each matched pair (f0, f1), store the matching per vertex slot
            // Per algorithm: f0.cycleMatchingFacet[j0] = f1, f1.cycleMatchingFacet[j1] = f0
            for (const auto& [fkey0, fkey1] : pairs) {
                Cell_handle cell0 = lookup_cell(cell_by_index, fkey0.cell_index);
                Cell_handle cell1 = lookup_cell(cell_by_index, fkey1.cell_index);
                if (cell0 == Cell_handle() || cell1 == Cell_handle()) {
                    continue;
                }
                
                // Get vertex indices within each facet (0, 1, or 2)
                int j0 = get_vertex_index_in_facet(cell0, fkey0.facet_index, v0);
                int j1 = get_vertex_index_in_facet(cell1, fkey1.facet_index, v1);
                
                if (j0 < 0 || j1 < 0) {
                    continue;  // Facet doesn't contain the expected vertex
                }
                
                // Store the matching per algorithm
                facet_matching[fkey0].cycleMatchingFacet[j0] = fkey1;
                facet_matching[fkey0].hasMatch[j0] = true;
                
                facet_matching[fkey1].cycleMatchingFacet[j1] = fkey0;
                facet_matching[fkey1].hasMatch[j1] = true;
            }
        }
    }

    // =========================================================================
    // Step 2: For each active vertex, traverse matchings to form cycles
    // =========================================================================
    
    std::vector<Cell_handle> incident_cells;
    std::vector<FacetKey> isosurface_facets;
    std::vector<int8_t> isosurface_facet_vslots;
    std::vector<char> facet_visited;
    std::unordered_map<FacetKey, int, FacetKeyHash> facet_to_local;

    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
        if (!vit->info().active) continue;
        if (vit->info().is_dummy) continue;
        const Vertex_handle v_handle = vit;

        // Clear any existing cycles
        vit->info().facet_cycles.clear();

        // Collect all isosurface facets incident on this vertex
        isosurface_facets.clear();
        isosurface_facet_vslots.clear();
        incident_cells.clear();
        dt.incident_cells(v_handle, std::back_inserter(incident_cells));

        for (Cell_handle ch : incident_cells) {
            if (dt.is_infinite(ch)) continue;

            const int v_local = find_vertex_index_in_cell(ch, v_handle);
            if (v_local < 0) {
                continue;
            }

            for (int i = 0; i < 4; ++i) {
                if (!ch->info().facet_is_isosurface[i]) continue;
                if (v_local == i) continue;
                
                // Check for dummy vertices
                bool has_dummy = false;
                for (int j = 0; j < 4; ++j) {
                    if (j == i) continue;
                    if (ch->vertex(j)->info().is_dummy) {
                        has_dummy = true;
                        break;
                    }
                }
                if (has_dummy) continue;

                int facet_vertex_slot = -1;
                for (int t = 0; t < 3; ++t) {
                    const int cell_vertex_idx = (i + 1 + t) % 4;
                    if (cell_vertex_idx == v_local) {
                        facet_vertex_slot = t;
                        break;
                    }
                }
                if (facet_vertex_slot < 0) {
                    facet_vertex_slot = get_vertex_index_in_facet(ch, i, v_handle);
                }
                if (facet_vertex_slot < 0 || facet_vertex_slot >= 3) {
                    continue;
                }

                isosurface_facets.push_back(FacetKey{ch->info().index, i});
                isosurface_facet_vslots.push_back(static_cast<int8_t>(facet_vertex_slot));
            }
        }

        const size_t facet_count = isosurface_facets.size();
        if (facet_count == 0) {
            continue;
        }
        facet_visited.assign(facet_count, 0);

        facet_to_local.clear();
        if (facet_count > 64) {
            facet_to_local.reserve(facet_count);
            for (size_t i = 0; i < facet_count; ++i) {
                facet_to_local.emplace(isosurface_facets[i], static_cast<int>(i));
            }
        }

        const auto find_local_index = [&](const FacetKey& key) -> int {
            if (facet_count > 64) {
                auto it = facet_to_local.find(key);
                return (it == facet_to_local.end()) ? -1 : it->second;
            }
            for (size_t i = 0; i < facet_count; ++i) {
                if (isosurface_facets[i] == key) {
                    return static_cast<int>(i);
                }
            }
            return -1;
        };

        std::vector<std::vector<std::pair<int, int>>> cycles;

        for (size_t start_idx = 0; start_idx < facet_count; ++start_idx) {
            if (facet_visited[start_idx]) {
                continue;
            }

            std::vector<std::pair<int, int>> cycle;
            size_t current_idx = start_idx;

            while (true) {
                if (facet_visited[current_idx]) {
                    break;  // Completed cycle or hit already-visited facet
                }
                facet_visited[current_idx] = 1;
                const FacetKey& current_fkey = isosurface_facets[current_idx];
                cycle.emplace_back(current_fkey.cell_index, current_fkey.facet_index);

                // Find the next facet in this cycle using cycleMatchingFacet[j]
                const int j = static_cast<int>(isosurface_facet_vslots[current_idx]);
                if (j < 0 || j >= 3) {
                    break;
                }

                auto match_it = facet_matching.find(current_fkey);
                if (match_it == facet_matching.end() || !match_it->second.hasMatch[j]) {
                    break;  // No matching for this vertex slot
                }

                FacetKey next_fkey = match_it->second.cycleMatchingFacet[j];
                
                // Verify next facet is incident on v (should be in our facet list)
                const int next_idx = find_local_index(next_fkey);
                if (next_idx < 0) {
                    break;  // Next facet is not incident on this vertex
                }

                current_idx = static_cast<size_t>(next_idx);
            }

            if (!cycle.empty()) {
                cycles.push_back(cycle);
            }
        }

        // Store cycles in vertex info
        vit->info().facet_cycles = cycles;
        write_facet_cycle_indices_for_vertex(v_handle, cell_by_index);


    }

}

// ============================================================================
// Step 9: Compute Cycle Isovertices
// ============================================================================

Point compute_centroid_of_voronoi_edge_and_isosurface(
    const Delaunay& dt,
    const UnifiedGrid& grid,
    float isovalue,
    Vertex_handle v_handle,
    const std::vector<std::pair<int, int>>& cycle_facets,
    const std::vector<Cell_handle>& cell_by_index
) {

    double cx = 0, cy = 0, cz = 0;
    int numI = 0;

    for (const auto& [cell_idx, facet_idx] : cycle_facets) {
        if (facet_idx < 0 || facet_idx >= 4) {
            throw std::runtime_error("Invalid facet index in cycle facets.");
        }

        Cell_handle cellA = lookup_cell(cell_by_index, cell_idx);
        if (cellA == Cell_handle()) {
            throw std::runtime_error("Cycle facet references invalid cell index.");
        }
        Cell_handle cellB = cellA->neighbor(facet_idx);

        if (dt.is_infinite(cellB)) {
            throw std::runtime_error("Cycle facet references an infinite neighbor cell.");
        }

        // Get circumcenters (dual Voronoi vertices)
        Point wA = cellA->info().circumcenter;
        Point wB = cellB->info().circumcenter;

        // Get scalar values at circumcenters
        float sA = cellA->info().circumcenter_scalar;
        float sB = cellB->info().circumcenter_scalar;

        // Linear interpolation to find an isovalue crossing point on the Voronoi edge.
        //
        // NOTE: In the standard pipeline, facets are marked only from the positive-cell
        // side and (sA - sB) is strictly positive and the isovalue is bracketed.
        // Some experimental pipelines may flip cell signs (without changing sA/sB),
        // in which case the scalar ordering can be reversed or the isovalue may not be
        // bracketed by sA and sB. Handle these cases deterministically by:
        // - using the symmetric interpolation formula,
        // - clamping t to [0,1] so the point lies on the segment between circumcenters.
        const float denom = sB - sA;
        float t = 0.5f;
        if (std::fabs(denom) > 1e-20f) {
            t = (isovalue - sA) / denom;
            if (t < 0.0f) t = 0.0f;
            if (t > 1.0f) t = 1.0f;
        }

        // Interpolate position
        double px = wA.x() + t * (wB.x() - wA.x());
        double py = wA.y() + t * (wB.y() - wA.y());
        double pz = wA.z() + t * (wB.z() - wA.z());

        cx += px;
        cy += py;
        cz += pz;
        numI++;
    }

    if (numI == 0) {
        // According to vdc-DelaunayBased.txt this should never happen.
        throw std::runtime_error("Compute_centroid_of_Voronoi_edge_and_isosurface(): numI == 0.");
    }

    return Point(cx / numI, cy / numI, cz / numI);
}

Point project_to_sphere(const Point& point, const Point& center, double radius) {
    double dx = point.x() - center.x();
    double dy = point.y() - center.y();
    double dz = point.z() - center.z();

    double dist = std::sqrt(dx*dx + dy*dy + dz*dz);

    if (dist < 1e-10) {
        // Point is at center, return an arbitrary point on sphere
        return Point(center.x() + radius, center.y(), center.z());
    }

    double scale = radius / dist;
    return Point(
        center.x() + dx * scale,
        center.y() + dy * scale,
        center.z() + dz * scale
    );
}

//! @brief Clips an isovertex to the circumscribed sphere of a cube.
/*!
 * If the isovertex lies outside the circumscribed sphere centered at the cube center,
 * it is projected onto the sphere surface along the direction from center to isovertex.
 * The circumscribed radius is half the cube side length (approximation for stability).
 *
 * @param isovertex The isosurface vertex (centroid) to potentially clip.
 * @param cube_center The center of the active cube (Delaunay vertex).
 * @param cube_side_length The side length of the active cube.
 * @return The clipped isovertex (same as input if already inside sphere).
 */
Point clip_isovertex_to_circumscribed_sphere(
    const Point& isovertex,
    const Point& cube_center,
    double cube_side_length
) {
    // Use half the cube side length as the circumscribed radius
    // (conservative approximation; true circumscribed radius is sqrt(3)/2 * side)
    const double circumscribed_radius = 0.5 * cube_side_length;

    // Vector from cube center to isovertex
    double dx = isovertex.x() - cube_center.x();
    double dy = isovertex.y() - cube_center.y();
    double dz = isovertex.z() - cube_center.z();
    double distance = std::sqrt(dx*dx + dy*dy + dz*dz);

    // If isovertex is outside the circumscribed sphere, project it onto the sphere
    if (distance > circumscribed_radius) {
        // Normalize direction and scale to sphere radius
        double scale = circumscribed_radius / distance;
        return Point(
            cube_center.x() + dx * scale,
            cube_center.y() + dy * scale,
            cube_center.z() + dz * scale
        );
    }

    return isovertex;
}

// ============================================================================
// Step 9.b-c: Self-Intersection Detection and Resolution
// ============================================================================

//! @brief CGAL Triangle_3 type for intersection tests
typedef K::Triangle_3 Triangle_3;

enum class TriIntersectionKind {
    NONE = 0,
    DISJOINT_TRIANGLE_TRIANGLE,
    SHARED_VERTEX_T1_OPPOSITE_SEGMENT,
    SHARED_VERTEX_T2_OPPOSITE_SEGMENT,
    SHARED_VERTEX_INTERIOR_SEGMENT_T1_0,
    SHARED_VERTEX_INTERIOR_SEGMENT_T1_1,
    SHARED_VERTEX_INTERIOR_SEGMENT_T2_0,
    SHARED_VERTEX_INTERIOR_SEGMENT_T2_1,
};

struct TriIntersectionTrace {
    int shared_count = 0;
    Point shared_point = Point(0, 0, 0);
    TriIntersectionKind kind = TriIntersectionKind::NONE;
};

//! @brief Collect all triangles that would be generated for a specific cycle
/*!
 * For each facet in the cycle, determines the triangle vertices using the
 * cycle's isovertex and the isovertices of the other two Delaunay vertices.
 *
 * @param dt The Delaunay triangulation
 * @param v The vertex whose cycle we're examining
 * @param cycle_idx Index of the cycle within v's facet_cycles
 * @param cycle_isovertex The isovertex position for this cycle
 * @param cell_by_index Cell handles indexed by `cell->info().index`
 * @return Vector of Triangle_3 objects for this cycle
 */
static std::vector<Triangle_3> collect_cycle_triangles(
    const Delaunay& dt,
    Vertex_handle v,
    int cycle_idx,
    const Point& cycle_isovertex,
    const std::vector<Cell_handle>& cell_by_index
) {
    std::vector<Triangle_3> triangles;
    const auto& cycles = v->info().facet_cycles;

    if (cycle_idx < 0 || cycle_idx >= static_cast<int>(cycles.size())) {
        return triangles;
    }

    const auto& cycle_facets = cycles[cycle_idx];
    triangles.reserve(cycle_facets.size());

    for (const auto& [cell_idx, facet_idx] : cycle_facets) {
        Cell_handle ch = lookup_cell(cell_by_index, cell_idx);
        if (ch == Cell_handle()) {
            continue;
        }

        // Facet facet_idx is opposite vertex facet_idx, so its 3 vertices are
        // at indices (facet_idx+1)%4, (facet_idx+2)%4, (facet_idx+3)%4.
        Vertex_handle other_verts[2] = {Vertex_handle(), Vertex_handle()};
        int other_slots[2] = {-1, -1};
        int other_count = 0;
        for (int t = 0; t < 3; ++t) {
            const int cell_vertex_idx = (facet_idx + 1 + t) % 4;
            Vertex_handle vh = ch->vertex(cell_vertex_idx);
            if (vh == v) {
                continue;
            }
            if (other_count < 2) {
                other_verts[other_count] = vh;
                other_slots[other_count] = t;
            }
            ++other_count;
        }

        // The facet should contain v (exactly one of the 3 vertices), leaving exactly 2 others.
        if (other_count != 2) {
            continue;
        }

        // Get isovertex positions for the other two vertices
        // These vertices may also have cycles - need to find which cycle
        // contains this facet for each vertex
        Point p1 = cycle_isovertex; // This cycle's isovertex

        // For the other vertices, find their isovertex for the cycle containing this facet
        Point p2, p3;
        bool found_p2 = false, found_p3 = false;

        Vertex_handle v1 = other_verts[0];
        Vertex_handle v2 = other_verts[1];

        // Find cycle for v1 containing this facet
        int c1 = -1;
        if (facet_idx >= 0 && facet_idx < 4 &&
            other_slots[0] >= 0 && other_slots[0] < 3) {
            c1 = ch->info().facet_info[facet_idx].dualCellEdgeIndex[other_slots[0]];
        }
        if (c1 < 0) {
            c1 = find_cycle_containing_facet(v1, cell_idx, facet_idx);
        }
        if (c1 >= 0 && c1 < static_cast<int>(v1->info().cycle_isovertices.size())) {
            p2 = v1->info().cycle_isovertices[c1];
            found_p2 = true;
        }

        // Find cycle for v2 containing this facet
        int c2 = -1;
        if (facet_idx >= 0 && facet_idx < 4 &&
            other_slots[1] >= 0 && other_slots[1] < 3) {
            c2 = ch->info().facet_info[facet_idx].dualCellEdgeIndex[other_slots[1]];
        }
        if (c2 < 0) {
            c2 = find_cycle_containing_facet(v2, cell_idx, facet_idx);
        }
        if (c2 >= 0 && c2 < static_cast<int>(v2->info().cycle_isovertices.size())) {
            p3 = v2->info().cycle_isovertices[c2];
            found_p3 = true;
        }

        // Only create triangle if found all three vertices
        if (found_p2 && found_p3) {
            triangles.push_back(Triangle_3(p1, p2, p3));
        }
    }

    return triangles;
}

// ============================================================================
// Helper: Collect cycle triangles with position overrides for joint resolution
// ============================================================================
// This variant allows specifying overridden isovertex positions for specific
// neighbor vertices, enabling joint resolution of coupled multi-cycle sites.

// Forward declaration
static bool triangles_have_nontrivial_intersection(const Triangle_3& t1, const Triangle_3& t2);

using PositionOverrideMap = std::unordered_map<int, std::vector<Point>>;

static std::vector<Triangle_3> collect_cycle_triangles_with_overrides(
    const Delaunay& dt,
    Vertex_handle v,
    int cycle_idx,
    const Point& cycle_isovertex,
    const std::vector<Cell_handle>& cell_by_index,
    const PositionOverrideMap& position_overrides
) {
    std::vector<Triangle_3> triangles;
    const auto& cycles = v->info().facet_cycles;

    if (cycle_idx < 0 || cycle_idx >= static_cast<int>(cycles.size())) {
        return triangles;
    }

    const auto& cycle_facets = cycles[cycle_idx];
    triangles.reserve(cycle_facets.size());

    for (const auto& [cell_idx, facet_idx] : cycle_facets) {
        Cell_handle ch = lookup_cell(cell_by_index, cell_idx);
        if (ch == Cell_handle()) {
            continue;
        }

        Vertex_handle other_verts[2] = {Vertex_handle(), Vertex_handle()};
        int other_slots[2] = {-1, -1};
        int other_count = 0;
        for (int t = 0; t < 3; ++t) {
            const int cell_vertex_idx = (facet_idx + 1 + t) % 4;
            Vertex_handle vh = ch->vertex(cell_vertex_idx);
            if (vh == v) {
                continue;
            }
            if (other_count < 2) {
                other_verts[other_count] = vh;
                other_slots[other_count] = t;
            }
            ++other_count;
        }

        if (other_count != 2) {
            continue;
        }

        Point p1 = cycle_isovertex;
        Point p2, p3;
        bool found_p2 = false, found_p3 = false;

        Vertex_handle v1 = other_verts[0];
        Vertex_handle v2 = other_verts[1];

        // Helper to get isovertex position, checking overrides first
        auto get_isovertex = [&](Vertex_handle vh, int cycle) -> std::optional<Point> {
            const int vidx = vh->info().index;
            auto it = position_overrides.find(vidx);
            if (it != position_overrides.end()) {
                if (cycle >= 0 && cycle < static_cast<int>(it->second.size())) {
                    return it->second[cycle];
                }
            }
            // Fallback to stored positions
            if (cycle >= 0 && cycle < static_cast<int>(vh->info().cycle_isovertices.size())) {
                return vh->info().cycle_isovertices[cycle];
            }
            return std::nullopt;
        };

        int c1 = -1;
        if (facet_idx >= 0 && facet_idx < 4 &&
            other_slots[0] >= 0 && other_slots[0] < 3) {
            c1 = ch->info().facet_info[facet_idx].dualCellEdgeIndex[other_slots[0]];
        }
        if (c1 < 0) {
            c1 = find_cycle_containing_facet(v1, cell_idx, facet_idx);
        }
        if (auto pos = get_isovertex(v1, c1)) {
            p2 = *pos;
            found_p2 = true;
        }

        int c2 = -1;
        if (facet_idx >= 0 && facet_idx < 4 &&
            other_slots[1] >= 0 && other_slots[1] < 3) {
            c2 = ch->info().facet_info[facet_idx].dualCellEdgeIndex[other_slots[1]];
        }
        if (c2 < 0) {
            c2 = find_cycle_containing_facet(v2, cell_idx, facet_idx);
        }
        if (auto pos = get_isovertex(v2, c2)) {
            p3 = *pos;
            found_p3 = true;
        }

        if (found_p2 && found_p3) {
            triangles.push_back(Triangle_3(p1, p2, p3));
        }
    }

    return triangles;
}

// Helper: Count self-intersection pairs with position overrides
static int count_self_intersection_pairs_with_overrides(
    Vertex_handle v,
    const std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    const PositionOverrideMap& position_overrides,
    int max_count = std::numeric_limits<int>::max()
) {
    if (max_count <= 0) {
        return 0;
    }

    const auto& cycles = v->info().facet_cycles;
    const size_t num_cycles = cycles.size();
    if (num_cycles < 1) {
        return 0;
    }

    int count = 0;

    // Check for coincident isovertices (with tolerance to catch near-identical positions)
    constexpr double COINCIDENCE_TOL2 = 1e-16; // ~1e-8 distance squared
    for (size_t a = 0; a < num_cycles && a < cycle_isovertices.size(); ++a) {
        for (size_t b = a + 1; b < num_cycles && b < cycle_isovertices.size(); ++b) {
            const double dx = cycle_isovertices[a].x() - cycle_isovertices[b].x();
            const double dy = cycle_isovertices[a].y() - cycle_isovertices[b].y();
            const double dz = cycle_isovertices[a].z() - cycle_isovertices[b].z();
            if (dx*dx + dy*dy + dz*dz < COINCIDENCE_TOL2) {
                if (++count >= max_count) return count;
            }
        }
    }

    // Collect triangles with overrides
    std::vector<std::vector<Triangle_3>> cycle_triangles(num_cycles);
    for (size_t c = 0; c < num_cycles && c < cycle_isovertices.size(); ++c) {
        cycle_triangles[c] = collect_cycle_triangles_with_overrides(
            dt, v, static_cast<int>(c), cycle_isovertices[c], cell_by_index, position_overrides);
    }

    // Check within-cycle intersections
    for (size_t c = 0; c < num_cycles; ++c) {
        const auto& tris = cycle_triangles[c];
        for (size_t i = 0; i < tris.size(); ++i) {
            for (size_t j = i + 1; j < tris.size(); ++j) {
                if (triangles_have_nontrivial_intersection(tris[i], tris[j])) {
                    if (++count >= max_count) return count;
                }
            }
        }
    }

    // Check between-cycle intersections
    for (size_t c0 = 0; c0 < num_cycles; ++c0) {
        for (size_t c1 = c0 + 1; c1 < num_cycles; ++c1) {
            for (const auto& t0 : cycle_triangles[c0]) {
                for (const auto& t1 : cycle_triangles[c1]) {
                    if (triangles_have_nontrivial_intersection(t0, t1)) {
                        if (++count >= max_count) return count;
                    }
                }
            }
        }
    }

    return count;
}

//! @brief Check if two triangles intersect (excluding shared edges/vertices)
/*!
 * Uses CGAL's do_intersect for Triangle_3 objects.
 * Returns true if the triangles have a non-trivial intersection
 * (i.e., more than just touching at shared vertices).
 */
static bool triangles_have_nontrivial_intersection(
    const Triangle_3& t1,
    const Triangle_3& t2,
    TriIntersectionTrace* trace_out
);

static int count_self_intersection_pairs(
    Vertex_handle v,
    const std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    int max_count = std::numeric_limits<int>::max()
);

static bool triangles_have_nontrivial_intersection(
    const Triangle_3& t1,
    const Triangle_3& t2
) {
    return triangles_have_nontrivial_intersection(t1, t2, nullptr);
}

static bool triangles_have_nontrivial_intersection(
    const Triangle_3& t1,
    const Triangle_3& t2,
    TriIntersectionTrace* trace_out
) {
    TriIntersectionTrace trace;

    // Find shared vertices (by exact position, which is stable for shared mesh vertices).
    const std::array<Point, 3> t1v = {t1.vertex(0), t1.vertex(1), t1.vertex(2)};
    const std::array<Point, 3> t2v = {t2.vertex(0), t2.vertex(1), t2.vertex(2)};
    int shared_count = 0;
    Point shared_point = t1v[0];
    for (const Point& a : t1v) {
        for (const Point& b : t2v) {
            if (a == b) {
                if (shared_count == 0) {
                    shared_point = a;
                }
                ++shared_count;
                break;
            }
        }
        if (shared_count >= 2) {
            break;
        }
    }
    trace.shared_count = shared_count;
    trace.shared_point = shared_point;

    // If the triangles are disjoint in vertex set, any intersection is non-trivial.
    if (shared_count == 0) {
        if (CGAL::do_intersect(t1, t2)) {
            trace.kind = TriIntersectionKind::DISJOINT_TRIANGLE_TRIANGLE;
            if (trace_out) {
                *trace_out = trace;
            }
            return true;
        }
        if (trace_out) {
            *trace_out = trace;
        }
        return false;
    }

    // Shared edge/triangle: adjacent triangles are expected to intersect along their shared simplex.
    // Overlap beyond that would indicate degeneracy, which we ignore here for performance.
    if (shared_count >= 2) {
        if (trace_out) {
            *trace_out = trace;
        }
        return false;
    }

    // -------------------------------------------------------------------------
    // Shared single vertex: orientation-based intersection test.
    // -------------------------------------------------------------------------
    // Given triangles (a,b,e) and (c,d,e) sharing vertex e, we determine if they
    // intersect at any point other than e using CGAL orientation predicates.
    //
    // The key insight: if c and d are on the same side of plane(e,a,b), the segment
    // cd cannot cross triangle (e,a,b) (and vice versa for a,b vs plane(e,c,d)).
    //
    // Algorithm:
    //   orient_c = Orientation(e, a, b, c)
    //   orient_d = Orientation(e, a, b, d)
    //   orient_a = Orientation(e, c, d, a)
    //   orient_b = Orientation(e, c, d, b)
    //
    // The triangles intersect non-trivially iff:
    //   isNotPPorNN(orient_c, orient_d) AND
    //   isNotPPorNN(orient_a, orient_b) AND
    //   isNotPPorNN(orient_a, orient_c) AND
    //   isNotPPorNN(orient_b, orient_d)
    //
    // where isNotPPorNN returns true if the two orientations are not both positive
    // and not both negative (i.e., different signs, or at least one is coplanar).
    // -------------------------------------------------------------------------
    const Point& e = shared_point;

    std::array<Point, 2> t1_other;
    std::array<Point, 2> t2_other;
    int t1_count = 0;
    int t2_count = 0;
    for (const Point& p : t1v) {
        if (p != e && t1_count < 2) {
            t1_other[t1_count++] = p;
        }
    }
    for (const Point& p : t2v) {
        if (p != e && t2_count < 2) {
            t2_other[t2_count++] = p;
        }
    }

    if (t1_count != 2 || t2_count != 2) {
        if (trace_out) {
            *trace_out = trace;
        }
        return false;
    }

    // Name the non-shared vertices: t1 = (a, b, e), t2 = (c, d, e)
    const Point& a = t1_other[0];
    const Point& b = t1_other[1];
    const Point& c = t2_other[0];
    const Point& d = t2_other[1];

    // Helper: returns true if (o1, o2) is NOT (positive, positive) and NOT (negative, negative).
    // This means: different signs, or at least one is COPLANAR.
    auto isNotPPorNN = [](CGAL::Orientation o1, CGAL::Orientation o2) -> bool {
        if (o1 != o2) {
            return true;  // Different orientations
        }
        // Both are the same; check if both are COPLANAR (degenerate case can intersect)
        if (o1 == CGAL::COPLANAR) {
            return true;
        }
        // Both are POSITIVE or both are NEGATIVE → no intersection across the plane
        return false;
    };

    const CGAL::Orientation orient_c = CGAL::orientation(e, a, b, c);
    const CGAL::Orientation orient_d = CGAL::orientation(e, a, b, d);
    const CGAL::Orientation orient_a = CGAL::orientation(e, c, d, a);
    const CGAL::Orientation orient_b = CGAL::orientation(e, c, d, b);

    // Special case: all four points are coplanar with e.
    // The orientation-based test is too conservative here (returns "may intersect").
    // Fall back to explicit segment-triangle intersection tests.
    const bool all_coplanar =
        (orient_c == CGAL::COPLANAR) &&
        (orient_d == CGAL::COPLANAR) &&
        (orient_a == CGAL::COPLANAR) &&
        (orient_b == CGAL::COPLANAR);

    if (all_coplanar) {
        // All points are in the same plane. Check for actual overlap:
        // Does the opposite edge of one triangle intersect the other triangle?
        const Segment3 t1_opposite(a, b);
        const Segment3 t2_opposite(c, d);
        const Triangle_3 tri1(a, b, e);
        const Triangle_3 tri2(c, d, e);

        if (CGAL::do_intersect(t1_opposite, tri2) || CGAL::do_intersect(t2_opposite, tri1)) {
            trace.kind = TriIntersectionKind::SHARED_VERTEX_T1_OPPOSITE_SEGMENT;
            if (trace_out) {
                *trace_out = trace;
            }
            return true;
        }
        if (trace_out) {
            *trace_out = trace;
        }
        return false;
    }

    // The triangles intersect at a point other than e iff all four conditions hold:
    //   - c and d are not strictly on the same side of plane(e,a,b)
    //   - a and b are not strictly on the same side of plane(e,c,d)
    //   - a and c are not strictly on the same side (cross-check)
    //   - b and d are not strictly on the same side (redundant safety)
    const bool intersects =
        isNotPPorNN(orient_c, orient_d) &&
        isNotPPorNN(orient_a, orient_b) &&
        isNotPPorNN(orient_a, orient_c) &&
        isNotPPorNN(orient_b, orient_d);

    if (intersects) {
        trace.kind = TriIntersectionKind::SHARED_VERTEX_T1_OPPOSITE_SEGMENT;  // Generic shared-vertex intersection
        if (trace_out) {
            *trace_out = trace;
        }
        return true;
    }

    if (trace_out) {
        *trace_out = trace;
    }
    return false;
}

static int count_self_intersection_pairs(
    Vertex_handle v,
    const std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    int max_count
) {
    if (max_count <= 0) {
        return 0;
    }

    const auto& cycles = v->info().facet_cycles;
    const size_t num_cycles = cycles.size();
    if (num_cycles < 1) {
        return 0;
    }

    int count = 0;

    // Treat coincident cycle isovertices as an intersection condition.
    // Use tolerance to catch near-identical positions that appear the same after rounding.
    constexpr double COINCIDENCE_TOL2 = 1e-16; // ~1e-8 distance squared
    for (size_t a = 0; a < num_cycles && a < cycle_isovertices.size(); ++a) {
        for (size_t b = a + 1; b < num_cycles && b < cycle_isovertices.size(); ++b) {
            const double dx = cycle_isovertices[a].x() - cycle_isovertices[b].x();
            const double dy = cycle_isovertices[a].y() - cycle_isovertices[b].y();
            const double dz = cycle_isovertices[a].z() - cycle_isovertices[b].z();
            if (dx*dx + dy*dy + dz*dz < COINCIDENCE_TOL2) {
                ++count;
                if (count >= max_count) {
                    return count;
                }
            }
        }
    }

    std::vector<std::vector<Triangle_3>> cycle_triangles(num_cycles);
    for (size_t c = 0; c < num_cycles; ++c) {
        if (c >= cycle_isovertices.size()) {
            break;
        }
        cycle_triangles[c] = collect_cycle_triangles(
            dt, v, static_cast<int>(c), cycle_isovertices[c], cell_by_index);
    }

    for (size_t c = 0; c < num_cycles; ++c) {
        const auto& tris = cycle_triangles[c];
        for (size_t i = 0; i < tris.size(); ++i) {
            for (size_t j = i + 1; j < tris.size(); ++j) {
                if (triangles_have_nontrivial_intersection(tris[i], tris[j])) {
                    ++count;
                    if (count >= max_count) {
                        return count;
                    }
                }
            }
        }
    }

    for (size_t c0 = 0; c0 < num_cycles; ++c0) {
        for (size_t c1 = c0 + 1; c1 < num_cycles; ++c1) {
            for (const auto& t0 : cycle_triangles[c0]) {
                for (const auto& t1 : cycle_triangles[c1]) {
                    if (triangles_have_nontrivial_intersection(t0, t1)) {
                        ++count;
                        if (count >= max_count) {
                            return count;
                        }
                    }
                }
            }
        }
    }

    return count;
}

bool check_self_intersection(
    Vertex_handle v,
    const std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index
) {
    const int v_idx = v->info().index;
    const bool trace_details = vdc_debug::trace_selfi_enabled(v_idx) && vdc_debug::trace_selfi_intersection_details_enabled();
    const bool log_within = vdc_debug::log_selfi_within_cycle_vertices();

    const auto& cycles = v->info().facet_cycles;
    size_t num_cycles = cycles.size();

    if (num_cycles < 1) {
        return false;
    }

    // Between-cycle degeneracy: if two cycle isovertices are coincident (or nearly so), then
    // triangles from those cycles can touch at a single point without sharing mesh vertex indices,
    // which downstream tools (e.g. ijkmeshinfo -selfI) count as self-intersection. Treat this as
    // a self-intersection so the A-only repositioning logic can separate the cycles deterministically.
    constexpr double COINCIDENCE_TOL2 = 1e-16; // ~1e-8 distance squared
    if (num_cycles > 1) {
        for (size_t a = 0; a < num_cycles && a < cycle_isovertices.size(); ++a) {
            for (size_t b = a + 1; b < num_cycles && b < cycle_isovertices.size(); ++b) {
                const double dx = cycle_isovertices[a].x() - cycle_isovertices[b].x();
                const double dy = cycle_isovertices[a].y() - cycle_isovertices[b].y();
                const double dz = cycle_isovertices[a].z() - cycle_isovertices[b].z();
                if (dx*dx + dy*dy + dz*dz < COINCIDENCE_TOL2) {
                    if (trace_details) {
                        DEBUG_PRINT("[DEL-SELFI-INT] Vertex " << v_idx
                                    << " coincident cycle isovertices: (" << a << "," << b << ") at "
                                    << cycle_isovertices[a]);
                    }
                    return true;
                }
            }
        }
    }

    // Fast path: for a single-cycle vertex, only check for intersections within the fan.
    // With fewer than 4 triangles in the fan, all triangle pairs share an edge, and the
    // non-trivial intersection test ignores those adjacent overlaps.
    if (num_cycles == 1 && cycles[0].size() < 4) {
        return false;
    }

    // Collect triangles for each cycle
    std::vector<std::vector<Triangle_3>> cycle_triangles(num_cycles);
    for (size_t c = 0; c < num_cycles; ++c) {
        cycle_triangles[c] = collect_cycle_triangles(
            dt, v, static_cast<int>(c), cycle_isovertices[c], cell_by_index);
    }

    // Check for intersections within each cycle fan.
    for (size_t c = 0; c < num_cycles; ++c) {
        const auto& tris = cycle_triangles[c];
        for (size_t i = 0; i < tris.size(); ++i) {
            for (size_t j = i + 1; j < tris.size(); ++j) {
                TriIntersectionTrace info;
                TriIntersectionTrace* info_out = (trace_details || log_within) ? &info : nullptr;
                if (triangles_have_nontrivial_intersection(tris[i], tris[j], info_out)) {
                    return true;
                }
            }
        }
    }

    // Check all pairs of cycles for intersection.
    for (size_t i = 0; i < num_cycles; ++i) {
        for (size_t j = i + 1; j < num_cycles; ++j) {
            // Check if any triangle from cycle i intersects any triangle from cycle j
            for (const auto& t1 : cycle_triangles[i]) {
                for (const auto& t2 : cycle_triangles[j]) {
                    if (triangles_have_nontrivial_intersection(t1, t2)) {
                        return true;
                    }
                }
            }
        }
    }

    return false;
}

struct SelfIntersectionWitness {
    int cycle0 = -1;
    int cycle1 = -1;
    size_t tri0 = 0;
    size_t tri1 = 0;
};

static std::optional<SelfIntersectionWitness> find_self_intersection_witness(
    Vertex_handle v,
    const std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index
) {
    const auto& cycles = v->info().facet_cycles;
    const size_t num_cycles = cycles.size();
    if (num_cycles < 1) {
        return std::nullopt;
    }

    std::vector<std::vector<Triangle_3>> cycle_triangles(num_cycles);
    for (size_t c = 0; c < num_cycles; ++c) {
        cycle_triangles[c] = collect_cycle_triangles(
            dt, v, static_cast<int>(c), cycle_isovertices[c], cell_by_index);
    }

    // Within-cycle.
    for (size_t c = 0; c < num_cycles; ++c) {
        const auto& tris = cycle_triangles[c];
        for (size_t i = 0; i < tris.size(); ++i) {
            for (size_t j = i + 1; j < tris.size(); ++j) {
                if (triangles_have_nontrivial_intersection(tris[i], tris[j])) {
                    return SelfIntersectionWitness{
                        static_cast<int>(c), static_cast<int>(c), i, j};
                }
            }
        }
    }

    // Between cycles.
    for (size_t i = 0; i < num_cycles; ++i) {
        for (size_t j = i + 1; j < num_cycles; ++j) {
            for (size_t ti = 0; ti < cycle_triangles[i].size(); ++ti) {
                for (size_t tj = 0; tj < cycle_triangles[j].size(); ++tj) {
                    if (triangles_have_nontrivial_intersection(
                            cycle_triangles[i][ti], cycle_triangles[j][tj])) {
                        return SelfIntersectionWitness{
                            static_cast<int>(i), static_cast<int>(j), ti, tj};
                    }
                }
            }
        }
    }

    return std::nullopt;
}

double compute_sphere_radius(
    Vertex_handle v,
    const Delaunay& dt,
    const UnifiedGrid& grid,
    int sep_split,
    int supersample_r
) {
    // Suppress unused parameter warning - v is kept for potential future use.
    (void)v;
    (void)dt;

    // Sphere radius for multi-cycle isovertex projection: inscribed sphere
    // of the effective grid cell that the Delaunay vertex is located in.
    //
    // The effective cell side length accounts for:
    // - sep_split: Number of cube splits per axis (K splits -> factor K+1)
    //   e.g., -sep_split 2 divides each axis by 3
    // - supersample_r: Supersample factor
    //   e.g., -supersample 2 divides each axis by 2
    //
    // Formula: effective_side_length = grid_spacing / (subgrid_scale * supersample_scale)
    // Inscribed sphere radius = 0.5 * effective_side_length

    const double grid_spacing = std::min({grid.spacing[0], grid.spacing[1], grid.spacing[2]});
    
    // Subgrid scale from sep_split: K splits means divide by (K+1)
    const double subgrid_scale = static_cast<double>(std::max(1, sep_split + 1));
    
    // Supersample scale: factor R means divide by R
    const double supersample_scale = static_cast<double>(std::max(1, supersample_r));
    
    const double effective_side_length = grid_spacing / (subgrid_scale * supersample_scale);
    
    // Inscribed sphere radius is half the side length of the effective cube
    return 0.5 * effective_side_length;
}

// ============================================================================
// Self-intersection resolution helpers
// ============================================================================

static double vec_dot(const Vector3& a, const Vector3& b) {
    return a.x() * b.x() + a.y() * b.y() + a.z() * b.z();
}

static Vector3 vec_cross(const Vector3& a, const Vector3& b) {
    return Vector3(
        a.y() * b.z() - a.z() * b.y(),
        a.z() * b.x() - a.x() * b.z(),
        a.x() * b.y() - a.y() * b.x());
}

static double vec_norm(const Vector3& v0) {
    return std::sqrt(v0.squared_length());
}

static Vector3 vec_normalize_or(const Vector3& v0, const Vector3& fallback) {
    const double len = vec_norm(v0);
    if (len < 1e-12) {
        return fallback;
    }
    return v0 / len;
}

static Vector3 unit_dir_to_or(const Point& center, const Point& p, const Vector3& fallback) {
    return vec_normalize_or(Vector3(center, p), fallback);
}

static double squared_distance(const Point& a, const Point& b) {
    const double dx = a.x() - b.x();
    const double dy = a.y() - b.y();
    const double dz = a.z() - b.z();
    return dx * dx + dy * dy + dz * dz;
}

static bool have_min_separation(const std::vector<Point>& candidate, double min_sep2) {
    const int n = static_cast<int>(candidate.size());
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (squared_distance(candidate[i], candidate[j]) <= min_sep2) {
                return false;
            }
        }
    }
    return true;
}

// ============================================================================
//  Move-Cap: Bound multi-cycle isovertex movement by opposite facet planes
// ============================================================================

// - For a Delaunay vertex v and a specific incident cycle C, consider each cycle facet t=(v,*,*)
//   (t is an isosurface facet).
// - For each incident tetrahedron Δ that contains t (both adjacent cells),
//   let t' be the face of Δ opposite v (i.e., the face NOT containing v).
// - If the dihedral angle between t and t' is acute (< 90°), then the distance from v to
//   the plane of t' is an upper bound on how far the cycle's isovertex can move away from v
//   without crossing that plane.
// - We take maxdist = min over all such (acute) constraints, and cap the projection radius by
//   cap_r = min(r_inscribed, 0.5 * maxdist).
static int find_facet_index_in_neighbor_across(
    const Cell_handle& cell,
    const Cell_handle& neighbor
) {
    if (cell == Cell_handle() || neighbor == Cell_handle()) {
        return -1;
    }
    for (int i = 0; i < 4; ++i) {
        if (neighbor->neighbor(i) == cell) {
            return i;
        }
    }
    return -1;
}

static bool compute_unit_normal_toward_cell(
    const Delaunay& dt,
    const Cell_handle& cell,
    int facet_idx,
    Vector3* out_unit_normal
) {
    if (out_unit_normal == nullptr) {
        return false;
    }
    if (cell == Cell_handle() || dt.is_infinite(cell)) {
        return false;
    }
    if (facet_idx < 0 || facet_idx >= 4) {
        return false;
    }

    // Facet facet_idx is opposite vertex facet_idx; use CGAL's cyclic ordering.
    const Point p0 = cell->vertex((facet_idx + 1) % 4)->point();
    const Point p1 = cell->vertex((facet_idx + 2) % 4)->point();
    const Point p2 = cell->vertex((facet_idx + 3) % 4)->point();

    Vector3 n = vec_cross(Vector3(p0, p1), Vector3(p0, p2));
    const double len = vec_norm(n);
    if (len < 1e-12) {
        return false;
    }
    n = n / len;

    // Orient toward the tetrahedron interior using the vertex opposite this facet.
    // (This matches the robust orientation used in compute_cycle_separating_direction().)
    const Vertex_handle opposite_vh = cell->vertex(facet_idx);
    if (opposite_vh == Vertex_handle() || dt.is_infinite(opposite_vh)) {
        return false;
    }
    const Vector3 to_opposite(p0, opposite_vh->point());
    if (vec_dot(n, to_opposite) < 0.0) {
        n = -n;
    }

    *out_unit_normal = n;
    return true;
}

using IsosurfaceEdgeSet = std::unordered_set<EdgeKey, EdgeKeyHash>;

struct StrictTetPreservationContext {
    bool strict_mode = false;
    const IsosurfaceEdgeSet* isosurface_edges = nullptr;
};

struct CycleMoveCapBound {
    double maxdist = std::numeric_limits<double>::infinity();
    int acute_dihedral_constraints = 0;
    int strict_preserve_hits = 0;
    int strict_preserve_rejects = 0;
    int direction_keep = 0;
    int direction_skip = 0;
};

// MoveCapTraceCandidate, MoveCapTraceWriteResult, fill_move_cap_trace_tet_geometry,
// and write_move_cap_trace_file are now in vdc_debug.h / vdc_debug.cpp.
using vdc_debug::MoveCapTraceCandidate;
using vdc_debug::MoveCapTraceWriteResult;
using vdc_debug::fill_move_cap_trace_tet_geometry;
using vdc_debug::write_move_cap_trace_file;

static bool strict_tet_preservation_enabled(const StrictTetPreservationContext& context) {
    return context.strict_mode && context.isosurface_edges != nullptr;
}

static IsosurfaceEdgeSet build_isosurface_edge_membership_set(const Delaunay& dt) {
    IsosurfaceEdgeSet edges;
    edges.reserve(static_cast<size_t>(std::max<int64_t>(64, dt.number_of_finite_cells() * 3)));

    auto add_edge = [&](Vertex_handle a, Vertex_handle b) {
        if (a == Vertex_handle() || b == Vertex_handle()) {
            return;
        }
        if (dt.is_infinite(a) || dt.is_infinite(b)) {
            return;
        }
        if (!a->info().active || !b->info().active) {
            return;
        }
        if (a->info().is_dummy || b->info().is_dummy) {
            return;
        }
        edges.insert(EdgeKey::FromVertexIndices(a->info().index, b->info().index));
    };

    for (auto ch = dt.finite_cells_begin(); ch != dt.finite_cells_end(); ++ch) {
        for (int facet_idx = 0; facet_idx < 4; ++facet_idx) {
            if (!ch->info().facet_is_isosurface[facet_idx]) {
                continue;
            }
            Vertex_handle f0 = ch->vertex((facet_idx + 1) % 4);
            Vertex_handle f1 = ch->vertex((facet_idx + 2) % 4);
            Vertex_handle f2 = ch->vertex((facet_idx + 3) % 4);
            add_edge(f0, f1);
            add_edge(f1, f2);
            add_edge(f2, f0);
        }
    }

    return edges;
}

static bool is_isosurface_edge_member(
    const Delaunay& dt,
    const IsosurfaceEdgeSet& isosurface_edges,
    const Vertex_handle& a,
    const Vertex_handle& b
) {
    if (a == Vertex_handle() || b == Vertex_handle()) {
        return false;
    }
    if (dt.is_infinite(a) || dt.is_infinite(b)) {
        return false;
    }
    if (a->info().is_dummy || b->info().is_dummy) {
        return false;
    }
    return isosurface_edges.find(
               EdgeKey::FromVertexIndices(a->info().index, b->info().index)) !=
           isosurface_edges.end();
}

static bool strict_preserve_tet_orientation_for_cycle_facet(
    const Delaunay& dt,
    const Cell_handle& Delta,
    int t_facet_idx,
    Vertex_handle w0,
    const IsosurfaceEdgeSet& isosurface_edges
) {
    if (Delta == Cell_handle() || dt.is_infinite(Delta)) {
        return false;
    }
    if (t_facet_idx < 0 || t_facet_idx >= 4) {
        return false;
    }
    if (w0 == Vertex_handle() || dt.is_infinite(w0) || w0->info().is_dummy) {
        return false;
    }

    const int w0_local = find_vertex_index_in_cell(Delta, w0);
    if (w0_local < 0 || w0_local >= 4 || w0_local == t_facet_idx) {
        return false;
    }

    auto edge_dihedral_is_ge_90 = [&](int edge_a_local, int edge_b_local) -> bool {
        int opposite_locals[2] = {-1, -1};
        int opposite_count = 0;
        for (int i = 0; i < 4; ++i) {
            if (i == edge_a_local || i == edge_b_local) {
                continue;
            }
            if (opposite_count < 2) {
                opposite_locals[opposite_count++] = i;
            }
        }
        if (opposite_count != 2) {
            return false;
        }

        Vector3 n0, n1;
        if (!compute_unit_normal_toward_cell(dt, Delta, opposite_locals[0], &n0)) {
            return false;
        }
        if (!compute_unit_normal_toward_cell(dt, Delta, opposite_locals[1], &n1)) {
            return false;
        }

        // Internal dihedral >= 90deg iff dot(inward_normals) >= 0.
        // Use a tiny tolerance for near-right-angle numerical jitter.
        constexpr double kDihedralEps = 1e-12;
        return vec_dot(n0, n1) >= -kDihedralEps;
    };

    struct OpposingEdgePair {
        int e0_a;
        int e0_b;
        int e1_a;
        int e1_b;
    };
    constexpr OpposingEdgePair kOpposingPairs[3] = {
        {0, 1, 2, 3},
        {0, 2, 1, 3},
        {0, 3, 1, 2}};

    for (const OpposingEdgePair& pair : kOpposingPairs) {
        Vertex_handle e0_a = Delta->vertex(pair.e0_a);
        Vertex_handle e0_b = Delta->vertex(pair.e0_b);
        Vertex_handle e1_a = Delta->vertex(pair.e1_a);
        Vertex_handle e1_b = Delta->vertex(pair.e1_b);
        if (e0_a == Vertex_handle() || e0_b == Vertex_handle() ||
            e1_a == Vertex_handle() || e1_b == Vertex_handle()) {
            continue;
        }
        if (dt.is_infinite(e0_a) || dt.is_infinite(e0_b) ||
            dt.is_infinite(e1_a) || dt.is_infinite(e1_b)) {
            continue;
        }
        if (e0_a->info().is_dummy || e0_b->info().is_dummy ||
            e1_a->info().is_dummy || e1_b->info().is_dummy) {
            continue;
        }

        if (!is_isosurface_edge_member(dt, isosurface_edges, e0_a, e0_b)) {
            continue;
        }
        if (!is_isosurface_edge_member(dt, isosurface_edges, e1_a, e1_b)) {
            continue;
        }
        if (!edge_dihedral_is_ge_90(pair.e0_a, pair.e0_b)) {
            continue;
        }
        if (!edge_dihedral_is_ge_90(pair.e1_a, pair.e1_b)) {
            continue;
        }
        return true;
    }

    return false;
}

static CycleMoveCapBound compute_cycle_maxdist_to_opposite_face_planes(
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    Vertex_handle v,
    int cycle_idx,
    const StrictTetPreservationContext& strict_context,
    const Vector3* move_dir
) {
    CycleMoveCapBound out;
    if (v == Vertex_handle() || !v->info().active || v->info().is_dummy) {
        return out;
    }

    const auto& cycles = v->info().facet_cycles;
    if (cycle_idx < 0 || cycle_idx >= static_cast<int>(cycles.size())) {
        return out;
    }
    const auto& cycle_facets = cycles[static_cast<size_t>(cycle_idx)];
    if (cycle_facets.empty()) {
        return out;
    }

    const Point w = v->point(); // Delaunay site position (movement center)
    const bool strict_enabled = strict_tet_preservation_enabled(strict_context);
    constexpr double kDirectionGateEps = 1e-12;
    const bool directional_gate_enabled = strict_enabled && move_dir != nullptr;
    const Vector3 dir = (move_dir != nullptr ? *move_dir : Vector3(0.0, 0.0, 0.0));
    const bool trace = vdc_debug::trace_selfi_enabled(v->info().index);
    const std::string trace_dir = trace ? vdc_debug::trace_selfi_dump_dir() : std::string();
    const bool dump_trace = trace && !trace_dir.empty();
    std::vector<MoveCapTraceCandidate> trace_candidates;
    if (dump_trace) {
        trace_candidates.reserve(cycle_facets.size() * 2 + 4);
    }

    for (const auto& [cell_idx, facet_idx_pos] : cycle_facets) {
        Cell_handle cell_pos = lookup_cell(cell_by_index, cell_idx);
        if (cell_pos == Cell_handle() || dt.is_infinite(cell_pos)) {
            if (dump_trace) {
                MoveCapTraceCandidate row;
                row.source_cell_index = cell_idx;
                row.source_facet_index = facet_idx_pos;
                row.strict_mode = strict_enabled;
                row.reason = "source_cell_missing";
                trace_candidates.push_back(std::move(row));
            }
            continue;
        }
        if (facet_idx_pos < 0 || facet_idx_pos >= 4) {
            if (dump_trace) {
                MoveCapTraceCandidate row;
                row.source_cell_index = cell_idx;
                row.source_facet_index = facet_idx_pos;
                row.strict_mode = strict_enabled;
                row.reason = "source_facet_invalid";
                trace_candidates.push_back(std::move(row));
            }
            continue;
        }

        // The facet t is shared by two cells (positive + negative). Evaluate both tetrahedra Δ.
        struct IncidentCell {
            Cell_handle ch;
            int facet_idx = -1; // facet index of t within this cell
        };
        IncidentCell incident[2];
        incident[0] = {cell_pos, facet_idx_pos};
        incident[1] = {cell_pos->neighbor(facet_idx_pos), -1};
        if (incident[1].ch != Cell_handle() && !dt.is_infinite(incident[1].ch)) {
            incident[1].facet_idx = find_facet_index_in_neighbor_across(cell_pos, incident[1].ch);
            if (incident[1].facet_idx < 0) {
                incident[1].ch = Cell_handle();
            }
        } else {
            incident[1].ch = Cell_handle();
        }

        for (int incident_slot = 0; incident_slot < 2; ++incident_slot) {
            const IncidentCell& ic = incident[incident_slot];
            const Cell_handle Delta = ic.ch;
            const int t_facet_idx = ic.facet_idx;
            MoveCapTraceCandidate row;
            row.source_cell_index = cell_idx;
            row.source_facet_index = facet_idx_pos;
            row.incident_slot = incident_slot;
            row.strict_mode = strict_enabled;

            if (Delta == Cell_handle()) {
                if (dump_trace) {
                    row.reason = "missing_incident_cell";
                    trace_candidates.push_back(std::move(row));
                }
                continue;
            }
            row.delta_cell_index = Delta->info().index;
            if (t_facet_idx < 0 || t_facet_idx >= 4) {
                if (dump_trace) {
                    row.reason = "invalid_t_facet";
                    trace_candidates.push_back(std::move(row));
                }
                continue;
            }
            row.t_facet_idx = t_facet_idx;

            const int v_local = find_vertex_index_in_cell(Delta, v);
            row.v_local = v_local;
            if (v_local < 0 || v_local >= 4) {
                if (dump_trace) {
                    row.reason = "vertex_not_in_tet";
                    trace_candidates.push_back(std::move(row));
                }
                continue;
            }
            if (v_local == t_facet_idx) {
                // v must be a vertex of facet t.
                if (dump_trace) {
                    row.reason = "vertex_opposite_cycle_facet";
                    trace_candidates.push_back(std::move(row));
                }
                continue;
            }

            // t' is the face of Δ opposite v.
            const int tprime_facet_idx = v_local;
            row.tprime_facet_idx = tprime_facet_idx;
            fill_move_cap_trace_tet_geometry(dt, Delta, t_facet_idx, v, &row);

            Vector3 nt, ntp;
            if (!compute_unit_normal_toward_cell(dt, Delta, t_facet_idx, &nt)) {
                if (dump_trace) {
                    row.reason = "invalid_t_normal";
                    trace_candidates.push_back(std::move(row));
                }
                continue;
            }
            if (!compute_unit_normal_toward_cell(dt, Delta, tprime_facet_idx, &ntp)) {
                if (dump_trace) {
                    row.reason = "invalid_tprime_normal";
                    trace_candidates.push_back(std::move(row));
                }
                continue;
            }

            // Acute dihedral (< 90°) between the two faces, using consistent (cell-toward) normals.
            const double dot_nt_ntp = vec_dot(nt, ntp);
            row.dot_nt_ntp = dot_nt_ntp;
            row.acute_dihedral = (dot_nt_ntp < 0.0);
            if (!row.acute_dihedral) {
                if (dump_trace) {
                    row.reason = "non_acute_dihedral";
                    trace_candidates.push_back(std::move(row));
                }
                continue;
            }
            ++out.acute_dihedral_constraints;

            if (strict_enabled) {
                const bool preserve = strict_preserve_tet_orientation_for_cycle_facet(
                    dt, Delta, t_facet_idx, v, *strict_context.isosurface_edges);
                row.flag_preserve_orientation = preserve ? 1 : 0;
                if (preserve) {
                    ++out.strict_preserve_hits;
                } else {
                    ++out.strict_preserve_rejects;
                    if (dump_trace) {
                        row.reason = "strict_reject";
                        trace_candidates.push_back(std::move(row));
                    }
                    continue;
                }
            }

            const Vector3 ntp_away = -ntp;
            if (directional_gate_enabled) {
                const double dot_dir_ntp_away = vec_dot(dir, ntp_away);
                row.dot_dir_ntp_away = dot_dir_ntp_away;
                row.direction_pass = (dot_dir_ntp_away > kDirectionGateEps) ? 1 : 0;
                if (row.direction_pass == 0) {
                    ++out.direction_skip;
                    if (dump_trace) {
                        row.reason = "direction_gate_reject";
                        trace_candidates.push_back(std::move(row));
                    }
                    continue;
                }
                ++out.direction_keep;
            }

            // w' is the vertex of t' not shared by t: it is exactly the vertex opposite t in Δ.
            const Point wprime = Delta->vertex(t_facet_idx)->point();
            const double dist = std::abs(vec_dot(Vector3(w, wprime), ntp_away));
            row.dist = dist;
            row.accepted_for_bound = true;
            row.reason = "accepted";
            if (dump_trace) {
                trace_candidates.push_back(std::move(row));
            }
            if (dist < out.maxdist) {
                out.maxdist = dist;
            }
        }
    }

    if (dump_trace) {
        const MoveCapTraceWriteResult trace_write = write_move_cap_trace_file(
            v,
            cycle_idx,
            strict_enabled,
            directional_gate_enabled,
            (directional_gate_enabled ? &dir : nullptr),
            kDirectionGateEps,
            out.maxdist,
            out.acute_dihedral_constraints,
            out.strict_preserve_hits,
            out.strict_preserve_rejects,
            out.direction_keep,
            out.direction_skip,
            trace_candidates,
            trace_dir);
        if (trace_write.wrote) {
            DEBUG_PRINT("[DEL-MOVE-CAP-TRACE] v=" << v->info().index
                        << " cycle=" << cycle_idx
                        << " file=" << trace_write.path.string()
                        << " candidates=" << trace_write.candidates
                        << " preserve_true=" << trace_write.preserve_true
                        << " direction_keep=" << out.direction_keep
                        << " direction_skip=" << out.direction_skip);
        }
    }

    return out;
}

// ---------------------------------------------------------------------------
// Dihedral-based move gate:
// For multi-cycle vertices, treat the Delaunay site as the baseline position.
// If moving a cycle isovertex would reverse tetra orientation across an acute
// (< 90 deg) dihedral between two incident isosurface facets, reject the move
// and keep that cycle isovertex at the Delaunay site.
// ---------------------------------------------------------------------------
static bool cycle_move_reverses_acute_dihedral_tet_orientation(
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    Vertex_handle v,
    int cycle_idx,
    const Point& moved_position,
    const StrictTetPreservationContext& strict_context
) {
    if (v == Vertex_handle() || !v->info().active || v->info().is_dummy) {
        return false;
    }
    if (v->info().facet_cycles.size() < 2) {
        return false;
    }

    const auto& cycles = v->info().facet_cycles;
    if (cycle_idx < 0 || cycle_idx >= static_cast<int>(cycles.size())) {
        return false;
    }

    const Point center = v->point();
    if (squared_distance(moved_position, center) <= 1e-24) {
        return false;
    }
    const bool strict_enabled = strict_tet_preservation_enabled(strict_context);

    const auto& cycle_facets = cycles[static_cast<size_t>(cycle_idx)];
    for (const auto& [cell_idx, facet_idx_pos] : cycle_facets) {
        Cell_handle cell_pos = lookup_cell(cell_by_index, cell_idx);
        if (cell_pos == Cell_handle() || dt.is_infinite(cell_pos)) {
            continue;
        }
        if (facet_idx_pos < 0 || facet_idx_pos >= 4) {
            continue;
        }

        struct IncidentCell {
            Cell_handle ch;
            int facet_idx = -1; // facet index of cycle facet in this cell
        };
        IncidentCell incident[2];
        incident[0] = {cell_pos, facet_idx_pos};
        incident[1] = {cell_pos->neighbor(facet_idx_pos), -1};
        if (incident[1].ch != Cell_handle() && !dt.is_infinite(incident[1].ch)) {
            incident[1].facet_idx = find_facet_index_in_neighbor_across(cell_pos, incident[1].ch);
            if (incident[1].facet_idx < 0) {
                incident[1].ch = Cell_handle();
            }
        } else {
            incident[1].ch = Cell_handle();
        }

        for (const IncidentCell& ic : incident) {
            const Cell_handle Delta = ic.ch;
            const int t_facet_idx = ic.facet_idx;
            if (Delta == Cell_handle() || t_facet_idx < 0 || t_facet_idx >= 4) {
                continue;
            }

            const int v_local = find_vertex_index_in_cell(Delta, v);
            if (v_local < 0 || v_local >= 4 || v_local == t_facet_idx) {
                continue;
            }

            // t' = face opposite moved vertex v in Delta.
            const int tprime_facet_idx = v_local;

            if (strict_enabled) {
                if (!strict_preserve_tet_orientation_for_cycle_facet(
                        dt, Delta, t_facet_idx, v, *strict_context.isosurface_edges)) {
                    continue;
                }
            } else {
                // Legacy relax-mode behavior: only gate on t' faces that cross the isovalue.
                bool tprime_is_isosurface = false;
                const Cell_handle tprime_neigh = Delta->neighbor(tprime_facet_idx);
                if (tprime_neigh != Cell_handle() && !dt.is_infinite(tprime_neigh)) {
                    tprime_is_isosurface =
                        (Delta->info().flag_positive != tprime_neigh->info().flag_positive);
                }
                if (!tprime_is_isosurface) {
                    continue;
                }
            }

            Vector3 nt, ntp;
            if (!compute_unit_normal_toward_cell(dt, Delta, t_facet_idx, &nt)) {
                continue;
            }
            if (!compute_unit_normal_toward_cell(dt, Delta, tprime_facet_idx, &ntp)) {
                continue;
            }

            // Acute dihedral (< 90 deg) <=> angle between inward normals > 90 deg.
            if (vec_dot(nt, ntp) >= 0.0) {
                continue;
            }

            const Point a = Delta->vertex((tprime_facet_idx + 1) % 4)->point();
            const Point b = Delta->vertex((tprime_facet_idx + 2) % 4)->point();
            const Point c = Delta->vertex((tprime_facet_idx + 3) % 4)->point();

            const CGAL::Orientation base_orient = CGAL::orientation(center, a, b, c);
            if (base_orient == CGAL::COPLANAR) {
                continue;
            }
            const CGAL::Orientation moved_orient = CGAL::orientation(moved_position, a, b, c);
            if (moved_orient != base_orient) {
                return true;
            }
        }
    }

    return false;
}

static bool apply_dihedral_orientation_gate_to_cycle_position(
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    Vertex_handle v,
    int cycle_idx,
    Point& position,
    const StrictTetPreservationContext& strict_context
) {
    if (!cycle_move_reverses_acute_dihedral_tet_orientation(
            dt, cell_by_index, v, cycle_idx, position, strict_context)) {
        return false;
    }
    position = v->point();
    return true;
}

static int apply_dihedral_orientation_gate_to_positions(
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    Vertex_handle v,
    std::vector<Point>& positions,
    const StrictTetPreservationContext& strict_context
) {
    if (v == Vertex_handle()) {
        return 0;
    }
    int clamped = 0;
    const int n = std::min(
        static_cast<int>(positions.size()),
        static_cast<int>(v->info().facet_cycles.size()));
    for (int c = 0; c < n; ++c) {
        if (apply_dihedral_orientation_gate_to_cycle_position(
                dt, cell_by_index, v, c, positions[static_cast<size_t>(c)], strict_context)) {
            ++clamped;
        }
    }
    if (clamped > 0 && vdc_debug::trace_selfi_enabled(v->info().index)) {
        DEBUG_PRINT("[DEL-MOVE-GATE] v=" << v->info().index
                    << " clamped " << clamped
                    << " cycle move(s) to Delaunay site.");
    }
    return clamped;
}

// ============================================================================
//  Orientation/Half-space cycle separation
// ============================================================================

struct CycleBoundaryInfo {
    std::vector<Vertex_handle> boundary_verts;  // Distinct vertices (excluding v0) on cycle boundary
};

static CycleBoundaryInfo gather_cycle_boundary_info(
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    Vertex_handle v0,
    const std::vector<std::pair<int, int>>& cycle_facets
) {
    CycleBoundaryInfo out;

    std::unordered_set<int> boundary_indices;
    boundary_indices.reserve(cycle_facets.size() * 2);

    for (const auto& [cell_idx, facet_idx] : cycle_facets) {
        Cell_handle cell = lookup_cell(cell_by_index, cell_idx);
        if (cell == Cell_handle()) {
            continue;
        }
        if (dt.is_infinite(cell)) {
            continue;
        }
        if (facet_idx < 0 || facet_idx >= 4) {
            continue;
        }

        // Boundary vertices are the facet vertices other than v0.
        for (int j = 0; j < 4; ++j) {
            if (j == facet_idx) {
                continue;
            }
            Vertex_handle vh = cell->vertex(j);
            if (vh == v0) {
                continue;
            }
            if (dt.is_infinite(vh) || vh->info().is_dummy) {
                continue;
            }
            const int idx = vh->info().index;
            if (boundary_indices.insert(idx).second) {
                out.boundary_verts.push_back(vh);
            }
        }
    }

    return out;
}

// ============================================================================
//  Cycle Separating Direction via Min-Sphere of Normals
// ============================================================================

//! @brief Type definitions for CGAL Min_sphere on 3D points
typedef CGAL::Min_sphere_of_points_d_traits_3<K, double> MinSphereTraits;
typedef CGAL::Min_sphere_of_spheres_d<MinSphereTraits> MinSphere;

// Collect all finite cells incident on v0 (by index) for fast membership tests.
static std::unordered_set<int> collect_incident_finite_cell_indices_for_vertex(
    const Delaunay& dt,
    Vertex_handle v0
) {
    std::unordered_set<int> out;
    if (v0 == Vertex_handle()) {
        return out;
    }

    std::vector<Cell_handle> incident_cells;
    dt.incident_cells(v0, std::back_inserter(incident_cells));
    out.reserve(incident_cells.size());
    for (Cell_handle ch : incident_cells) {
        if (ch == Cell_handle() || dt.is_infinite(ch)) {
            continue;
        }
        out.insert(ch->info().index);
    }
    return out;
}

// Pseudocode source: SeparateMultiIsov.docx (user provided).
// Returns true iff the cells on the cycle side (positive cycle cells) do NOT reach cellA
// when traversal is constrained to finite cells incident on v0 and with all cycle facets
// (both sides) treated as barriers.
static bool are_cycle_cells_separated_from_cellA(
    Vertex_handle v0,
    const std::vector<std::pair<int, int>>& cycle_facets,
    Cell_handle cellA,
    const std::unordered_set<int>& incident_cell_indices,
    const std::vector<Cell_handle>& cell_by_index,
    const Delaunay& dt
) {
    if (v0 == Vertex_handle() || cellA == Cell_handle() || dt.is_infinite(cellA)) {
        return true;
    }
    if (cycle_facets.empty()) {
        return true;
    }

    const int cellA_idx = cellA->info().index;
    if (!incident_cell_indices.count(cellA_idx)) {
        return true;
    }

    std::unordered_set<int> visited;
    visited.reserve(incident_cell_indices.size());
    std::vector<Cell_handle> stackB;
    stackB.reserve(cycle_facets.size());

    // Initialize barrier + seed stack from cycle facets:
    // - push positive-side cycle cells;
    // - mark both positive and mirror (negative) cells as visited.
    for (const auto& [cell_idx, facet_idx] : cycle_facets) {
        Cell_handle cell_pos = lookup_cell(cell_by_index, cell_idx);
        if (cell_pos == Cell_handle() || dt.is_infinite(cell_pos)) {
            continue;
        }
        if (facet_idx < 0 || facet_idx >= 4) {
            continue;
        }
        if (!incident_cell_indices.count(cell_pos->info().index)) {
            continue;
        }

        if (visited.insert(cell_pos->info().index).second) {
            stackB.push_back(cell_pos);
        }

        Cell_handle cell_neg = cell_pos->neighbor(facet_idx);
        if (cell_neg != Cell_handle() && !dt.is_infinite(cell_neg) &&
            incident_cell_indices.count(cell_neg->info().index)) {
            visited.insert(cell_neg->info().index);
        }
    }

    // DFS/BFS over the side seeded by cycle positive cells.
    while (!stackB.empty()) {
        Cell_handle cellB = stackB.back();
        stackB.pop_back();
        if (cellB == Cell_handle() || dt.is_infinite(cellB)) {
            continue;
        }

        if (cellB->info().index == cellA_idx) {
            return false;
        }

        // Traverse to adjacent cells that are also incident on v0.
        for (int k = 0; k < 4; ++k) {
            Cell_handle cellC = cellB->neighbor(k);
            if (cellC == Cell_handle() || dt.is_infinite(cellC)) {
                continue;
            }
            const int idxC = cellC->info().index;
            if (!incident_cell_indices.count(idxC)) {
                continue;
            }
            if (visited.insert(idxC).second) {
                stackB.push_back(cellC);
            }
        }
    }

    return true;
}

// Build the non-symmetric separation matrix:
// is_separated[i][j] == 1 iff cycle i separates its cycle-side cells from a reference
// cell chosen from cycle j. (i==j is left as 0 by definition/use.)
static std::vector<std::vector<char>> determine_cycle_separation_matrix(
    Vertex_handle v0,
    const std::vector<std::vector<std::pair<int, int>>>& cycles,
    const std::vector<Cell_handle>& cell_by_index,
    const Delaunay& dt
) {
    const int num_cycles = static_cast<int>(cycles.size());
    std::vector<std::vector<char>> is_separated(
        static_cast<size_t>(num_cycles),
        std::vector<char>(static_cast<size_t>(num_cycles), 0));
    if (num_cycles <= 1 || v0 == Vertex_handle()) {
        return is_separated;
    }

    const std::unordered_set<int> incident_cell_indices =
        collect_incident_finite_cell_indices_for_vertex(dt, v0);
    if (incident_cell_indices.empty()) {
        return is_separated;
    }

    std::vector<Cell_handle> reference_cell_by_cycle(static_cast<size_t>(num_cycles), Cell_handle());
    for (int j = 0; j < num_cycles; ++j) {
        for (const auto& [cell_idx, facet_idx] : cycles[static_cast<size_t>(j)]) {
            (void)facet_idx;
            Cell_handle cj = lookup_cell(cell_by_index, cell_idx);
            if (cj == Cell_handle() || dt.is_infinite(cj)) {
                continue;
            }
            if (!incident_cell_indices.count(cj->info().index)) {
                continue;
            }
            reference_cell_by_cycle[static_cast<size_t>(j)] = cj;
            break;
        }
    }

    for (int i = 0; i < num_cycles; ++i) {
        for (int j = 0; j < num_cycles; ++j) {
            if (i == j) {
                continue;
            }
            Cell_handle cellA = reference_cell_by_cycle[static_cast<size_t>(j)];
            if (cellA == Cell_handle()) {
                // No valid reference cell for cycle j -> do not force flipping on this basis.
                is_separated[static_cast<size_t>(i)][static_cast<size_t>(j)] = 1;
                continue;
            }
            const bool sep = are_cycle_cells_separated_from_cellA(
                v0, cycles[static_cast<size_t>(i)], cellA,
                incident_cell_indices, cell_by_index, dt);
            is_separated[static_cast<size_t>(i)][static_cast<size_t>(j)] = sep ? 1 : 0;
        }
    }

    return is_separated;
}

struct CycleFoldoverEdgeInfo {
    bool found = false;
    int count = 0;
    int v_low = -1;
    int v_high = -1;
    bool incident_on_v0 = false;
    int other_vidx = -1;
    Vector3 unit_dir = Vector3(0.0, 0.0, 0.0);
};

// Detect a within-cycle crossover edge pattern:
// one (or more) cycle edges shared by >=4 cycle facets.
// When possible, return a direction along an edge incident on v0.
static CycleFoldoverEdgeInfo detect_cycle_crossover_edge_direction(
    Vertex_handle v0,
    const std::vector<std::pair<int, int>>& cycle_facets,
    const std::vector<Cell_handle>& cell_by_index,
    const Delaunay& dt,
    int* out_candidates_ge4 = nullptr
) {
    CycleFoldoverEdgeInfo out;
    if (out_candidates_ge4) {
        *out_candidates_ge4 = 0;
    }
    if (v0 == Vertex_handle() || cycle_facets.empty()) {
        return out;
    }

    struct CycleEdgeMultiplicityEntry {
        int count = 0;
        Point p_low;
        Point p_high;
        bool has_endpoints = false;
    };

    std::unordered_map<EdgeKey, CycleEdgeMultiplicityEntry, EdgeKeyHash> edge_multiplicity;
    edge_multiplicity.reserve(cycle_facets.size() * 3);

    auto count_cycle_edge = [&](const Vertex_handle& va, const Vertex_handle& vb) {
        if (va == Vertex_handle() || vb == Vertex_handle()) {
            return;
        }
        if (dt.is_infinite(va) || dt.is_infinite(vb)) {
            return;
        }

        const int ia = va->info().index;
        const int ib = vb->info().index;
        const EdgeKey key = EdgeKey::FromVertexIndices(ia, ib);
        auto& rec = edge_multiplicity[key];
        ++rec.count;
        if (!rec.has_endpoints) {
            if (key.v0 == ia) {
                rec.p_low = va->point();
                rec.p_high = vb->point();
            } else {
                rec.p_low = vb->point();
                rec.p_high = va->point();
            }
            rec.has_endpoints = true;
        }
    };

    for (const auto& [cell_idx, facet_idx] : cycle_facets) {
        Cell_handle cell = lookup_cell(cell_by_index, cell_idx);
        if (cell == Cell_handle() || dt.is_infinite(cell)) {
            continue;
        }
        if (facet_idx < 0 || facet_idx >= 4) {
            continue;
        }

        std::array<Vertex_handle, 3> facet_vhs;
        bool ok = true;
        for (int k = 0; k < 3; ++k) {
            const int cell_vertex_idx = (facet_idx + 1 + k) % 4;
            Vertex_handle vh = cell->vertex(cell_vertex_idx);
            if (vh == Vertex_handle() || dt.is_infinite(vh)) {
                ok = false;
                break;
            }
            facet_vhs[static_cast<size_t>(k)] = vh;
        }
        if (!ok) {
            continue;
        }

        count_cycle_edge(facet_vhs[0], facet_vhs[1]);
        count_cycle_edge(facet_vhs[1], facet_vhs[2]);
        count_cycle_edge(facet_vhs[2], facet_vhs[0]);
    }

    const int v0_index = v0->info().index;
    int candidates_ge4 = 0;
    for (const auto& kv : edge_multiplicity) {
        const EdgeKey& edge_key = kv.first;
        const CycleEdgeMultiplicityEntry& rec = kv.second;
        if (rec.count < 4 || !rec.has_endpoints) {
            continue;
        }

        ++candidates_ge4;

        const bool incident_on_v0 = (edge_key.v0 == v0_index || edge_key.v1 == v0_index);
        int other_vidx = -1;
        Vector3 dir_candidate(0.0, 0.0, 0.0);
        if (incident_on_v0) {
            if (edge_key.v0 == v0_index) {
                other_vidx = edge_key.v1;
                dir_candidate = Vector3(v0->point(), rec.p_high);
            } else {
                other_vidx = edge_key.v0;
                dir_candidate = Vector3(v0->point(), rec.p_low);
            }
        } else {
            // Non-incident case is rare; orient edge direction consistently away from v0.
            dir_candidate = Vector3(rec.p_low, rec.p_high);
            const Point midpoint(
                0.5 * (rec.p_low.x() + rec.p_high.x()),
                0.5 * (rec.p_low.y() + rec.p_high.y()),
                0.5 * (rec.p_low.z() + rec.p_high.z()));
            if (vec_dot(dir_candidate, Vector3(v0->point(), midpoint)) < 0.0) {
                dir_candidate = -dir_candidate;
            }
        }

        const double dlen = vec_norm(dir_candidate);
        if (dlen < 1e-12) {
            continue;
        }

        bool better = false;
        if (!out.found) {
            better = true;
        } else if (incident_on_v0 != out.incident_on_v0) {
            // Prefer an edge that actually includes v0.
            better = incident_on_v0;
        } else if (rec.count != out.count) {
            // Prefer stronger multiplicity evidence.
            better = (rec.count > out.count);
        } else if (incident_on_v0) {
            // Deterministic tie-break.
            better = (other_vidx < out.other_vidx);
        } else {
            // Deterministic tie-break for non-incident edges.
            better = (edge_key.v0 < out.v_low) ||
                     (edge_key.v0 == out.v_low && edge_key.v1 < out.v_high);
        }

        if (!better) {
            continue;
        }

        out.found = true;
        out.count = rec.count;
        out.v_low = edge_key.v0;
        out.v_high = edge_key.v1;
        out.incident_on_v0 = incident_on_v0;
        out.other_vidx = other_vidx;
        out.unit_dir = dir_candidate / dlen;
    }

    if (out_candidates_ge4) {
        *out_candidates_ge4 = candidates_ge4;
    }
    return out;
}

// Helper: extract the 3 vertex points of a Delaunay facet (opposite to facet_idx vertex).
// Returns false if any vertex is invalid/infinite.
static bool get_facet_triangle_points(
    const Delaunay& dt,
    Cell_handle cell,
    int facet_idx,
    Point out_pts[3]
) {
    for (int k = 0; k < 3; ++k) {
        const int vi = (facet_idx + 1 + k) % 4;
        Vertex_handle vh = cell->vertex(vi);
        if (vh == Vertex_handle() || dt.is_infinite(vh)) return false;
        out_pts[k] = vh->point();
    }
    return true;
}

// Helper: compute outward unit normal of a facet, oriented toward the positive cell.
static Vector3 compute_facet_outward_normal(
    const Delaunay& dt,
    Cell_handle cell,
    int facet_idx,
    const Point pts[3]
) {
    const Vector3 e1(pts[0], pts[1]);
    const Vector3 e2(pts[0], pts[2]);
    Vector3 n = vec_cross(e1, e2);
    const double len = vec_norm(n);
    if (len < 1e-20) return Vector3(0, 0, 0);
    n = n / len;
    // Orient toward the cell side (positive cell by convention).
    Vertex_handle opp = cell->vertex(facet_idx);
    if (opp != Vertex_handle() && !dt.is_infinite(opp)) {
        if (vec_dot(n, Vector3(pts[0], opp->point())) < 0) {
            n = -n;
        }
    }
    return n;
}

static std::vector<std::pair<int,int>> simplify_cycle_facets(
    Vertex_handle v0,
    const std::vector<std::pair<int,int>>& cycle_facets,
    const std::vector<Cell_handle>& cell_by_index,
    const Delaunay& dt,
    int cycle_index
) {
    if (v0 == Vertex_handle()) {
        return cycle_facets;
    }
    const auto& all_cycles = v0->info().facet_cycles;

    // Build a set of FacetKeys belonging to OTHER cycles of v0 for quick lookup.
    std::set<std::pair<int,int>> other_cycle_facets;
    for (int c = 0; c < static_cast<int>(all_cycles.size()); ++c) {
        if (c == cycle_index) continue;
        for (const auto& fk : all_cycles[static_cast<size_t>(c)]) {
            other_cycle_facets.insert(fk);
        }
    }

    const bool trace = vdc_debug::trace_selfi_enabled(v0->info().index);
    const std::string trace_dir = trace ? vdc_debug::trace_selfi_dump_dir() : std::string();
    const bool dump_file = trace && !trace_dir.empty();

    // Open trace file if tracing enabled.
    std::ofstream trace_out;
    if (dump_file) {
        std::error_code ec;
        std::filesystem::path out_dir(trace_dir);
        std::filesystem::create_directories(out_dir, ec);
        const int vidx = v0->info().index;
        const std::filesystem::path path =
            out_dir / ("sim_cyc_v" + std::to_string(vidx) + "_c" + std::to_string(cycle_index) + ".txt");
        trace_out.open(path);
        if (trace_out) {
            trace_out << std::setprecision(17);
            const Point center = v0->point();
            trace_out << "vertex_index " << vidx << "\n";
            trace_out << "cycle_index " << cycle_index << "\n";
            trace_out << "center " << center.x() << " " << center.y() << " " << center.z() << "\n";
            trace_out << "initial_facet_count " << cycle_facets.size() << "\n\n";
        }
    }

    // Lambda: dump the current state of the cycle (all facets with geometry + normals).
    auto dump_iteration_state = [&](int iteration, const std::vector<std::pair<int,int>>& facets,
                                    int sel_i, int sel_j, Cell_handle shared_tet,
                                    const std::vector<Cell_handle>& mc_cells,
                                    const std::vector<int>& mc_fis) {
        if (!trace_out) return;
        trace_out << "iteration " << iteration << "\n";
        trace_out << "facet_count " << facets.size() << "\n";

        // Dump all facets in this iteration.
        for (int i = 0; i < static_cast<int>(facets.size()); ++i) {
            const auto& [cell_idx, facet_idx] = facets[static_cast<size_t>(i)];
            Cell_handle cell = lookup_cell(cell_by_index, cell_idx);
            bool valid = (cell != Cell_handle() && !dt.is_infinite(cell)
                          && facet_idx >= 0 && facet_idx < 4);
            Point pts[3] = {Point(0,0,0), Point(0,0,0), Point(0,0,0)};
            Vector3 normal(0,0,0);
            if (valid) {
                valid = get_facet_triangle_points(dt, cell, facet_idx, pts);
                if (valid) {
                    normal = compute_facet_outward_normal(dt, cell, facet_idx, pts);
                }
            }

            // Mark if this facet is one of the selected pair.
            const char* role = "";
            if (i == sel_i) role = " SELECTED_F";
            else if (i == sel_j) role = " SELECTED_F'";

            // Mirror cell info (tet on non-cycle side).
            int mc_idx_val = -1;
            if (i < static_cast<int>(mc_cells.size()) && mc_cells[static_cast<size_t>(i)] != Cell_handle()) {
                mc_idx_val = mc_cells[static_cast<size_t>(i)]->info().index;
            }

            trace_out << "facet idx=" << i
                      << " cell=" << cell_idx << " fi=" << facet_idx
                      << " mirror_cell=" << mc_idx_val
                      << " p0=(" << pts[0].x() << "," << pts[0].y() << "," << pts[0].z() << ")"
                      << " p1=(" << pts[1].x() << "," << pts[1].y() << "," << pts[1].z() << ")"
                      << " p2=(" << pts[2].x() << "," << pts[2].y() << "," << pts[2].z() << ")"
                      << " n=(" << normal.x() << "," << normal.y() << "," << normal.z() << ")"
                      << role << "\n";
        }

        // Dump the shared tet geometry if a replacement is happening.
        if (shared_tet != Cell_handle() && !dt.is_infinite(shared_tet)) {
            trace_out << "shared_tet cell=" << shared_tet->info().index;
            for (int k = 0; k < 4; ++k) {
                Vertex_handle vh = shared_tet->vertex(k);
                if (vh != Vertex_handle() && !dt.is_infinite(vh)) {
                    const Point& p = vh->point();
                    trace_out << " v" << k << "=(" << p.x() << "," << p.y() << "," << p.z() << ")";
                } else {
                    trace_out << " v" << k << "=(inf)";
                }
            }
            trace_out << "\n";
        }
        trace_out << "\n";
    };

    auto compute_cycle_facet_normal = [&](int cell_idx, int facet_idx, Vector3* out_normal) -> bool {
        Cell_handle cell = lookup_cell(cell_by_index, cell_idx);
        if (cell == Cell_handle() || dt.is_infinite(cell)) return false;
        if (facet_idx < 0 || facet_idx >= 4) return false;
        Point pts[3];
        if (!get_facet_triangle_points(dt, cell, facet_idx, pts)) return false;
        const Vector3 normal = compute_facet_outward_normal(dt, cell, facet_idx, pts);
        if (vec_norm(normal) < 1e-20) return false;
        *out_normal = normal;
        return true;
    };

    std::vector<std::pair<int,int>> result = cycle_facets;
    bool changed = true;
    int iteration = 0;

    while (changed) {
        changed = false;
        const int n = static_cast<int>(result.size());
        if (n < 2) break;

        // For each facet, compute its mirror cell (the tet on the non-cycle side).
        // Group facets by mirror cell index.
        std::unordered_map<int, std::vector<int>> mirror_cell_to_facet_indices;
        std::vector<Cell_handle> mirror_cells(static_cast<size_t>(n), Cell_handle());
        std::vector<int> mirror_facet_indices(static_cast<size_t>(n), -1);

        for (int i = 0; i < n; ++i) {
            const auto& [cell_idx, facet_idx] = result[static_cast<size_t>(i)];
            Cell_handle cell = lookup_cell(cell_by_index, cell_idx);
            if (cell == Cell_handle() || dt.is_infinite(cell)) continue;
            if (facet_idx < 0 || facet_idx >= 4) continue;

            const Facet mirror = dt.mirror_facet(Facet(cell, facet_idx));
            Cell_handle mirror_cell = mirror.first;
            if (mirror_cell == Cell_handle() || dt.is_infinite(mirror_cell)) continue;

            mirror_cells[static_cast<size_t>(i)] = mirror_cell;
            mirror_facet_indices[static_cast<size_t>(i)] = mirror.second;
            const int mc_idx = mirror_cell->info().index;
            mirror_cell_to_facet_indices[mc_idx].push_back(i);
        }

        // Find a valid pair: co-tet, not in other cycle.
        // Choose the replacement facet orientation that best matches the normals
        // of the two facets being merged, while falling back to the old
        // positive-cell preference when the comparison is ambiguous.
        int replace_i = -1, replace_j = -1;
        Cell_handle shared_tet;
        int new_cell_idx = -1, new_facet_idx = -1;

        for (const auto& [mc_idx, indices] : mirror_cell_to_facet_indices) {
            if (indices.size() < 2) continue;
            for (size_t a = 0; a < indices.size() && replace_i < 0; ++a) {
                for (size_t b = a + 1; b < indices.size() && replace_i < 0; ++b) {
                    const int ii = indices[a];
                    const int jj = indices[b];
                    const auto& fi = result[static_cast<size_t>(ii)];
                    const auto& fj = result[static_cast<size_t>(jj)];
                    if (other_cycle_facets.count(fi) || other_cycle_facets.count(fj)) continue;

                    // Validate: find the third v0-incident facet and check it's an iso-facet.
                    Cell_handle cand_tet = mirror_cells[static_cast<size_t>(ii)];
                    int cand_mirror_fi = mirror_facet_indices[static_cast<size_t>(ii)];
                    int cand_mirror_fj = mirror_facet_indices[static_cast<size_t>(jj)];

                    const int v0_local = find_vertex_index_in_cell(cand_tet, v0);
                    if (v0_local < 0) continue;

                    int repl_fi = -1;
                    for (int k = 0; k < 4; ++k) {
                        if (k == v0_local) continue;
                        if (k == cand_mirror_fi) continue;
                        if (k == cand_mirror_fj) continue;
                        repl_fi = k;
                        break;
                    }
                    if (repl_fi < 0) continue;

                    struct ReplacementOption {
                        int cell_idx = -1;
                        int facet_idx = -1;
                        Vector3 normal = Vector3(0, 0, 0);
                    };

                    std::array<ReplacementOption, 2> options;
                    int option_count = 0;
                    auto append_option = [&](Cell_handle option_cell, int option_fi) {
                        if (option_cell == Cell_handle() || dt.is_infinite(option_cell)) return;
                        if (option_fi < 0 || option_fi >= 4) return;
                        const int option_cell_idx = option_cell->info().index;
                        for (int opt = 0; opt < option_count; ++opt) {
                            if (options[static_cast<size_t>(opt)].cell_idx == option_cell_idx &&
                                options[static_cast<size_t>(opt)].facet_idx == option_fi) {
                                return;
                            }
                        }
                        Vector3 option_normal;
                        if (!compute_cycle_facet_normal(option_cell_idx, option_fi, &option_normal)) return;
                        options[static_cast<size_t>(option_count)].cell_idx = option_cell_idx;
                        options[static_cast<size_t>(option_count)].facet_idx = option_fi;
                        options[static_cast<size_t>(option_count)].normal = option_normal;
                        ++option_count;
                    };

                    const Facet replacement_facet(cand_tet, repl_fi);
                    const Facet replacement_mirror = dt.mirror_facet(replacement_facet);
                    if (cand_tet->info().flag_positive) {
                        append_option(cand_tet, repl_fi);
                        append_option(replacement_mirror.first, replacement_mirror.second);
                    } else {
                        append_option(replacement_mirror.first, replacement_mirror.second);
                        append_option(cand_tet, repl_fi);
                    }
                    if (option_count == 0) continue;

                    Vector3 reference_normal(0, 0, 0);
                    bool has_reference_normal = false;
                    Vector3 merged_normal;
                    if (compute_cycle_facet_normal(fi.first, fi.second, &merged_normal)) {
                        reference_normal = reference_normal + merged_normal;
                        has_reference_normal = true;
                    }
                    if (compute_cycle_facet_normal(fj.first, fj.second, &merged_normal)) {
                        reference_normal = reference_normal + merged_normal;
                        has_reference_normal = true;
                    }

                    int best_option = 0;
                    if (has_reference_normal && vec_norm(reference_normal) >= 1e-12 && option_count > 1) {
                        double best_score =
                            vec_dot(options[0].normal, reference_normal);
                        for (int opt = 1; opt < option_count; ++opt) {
                            const double score =
                                vec_dot(options[static_cast<size_t>(opt)].normal, reference_normal);
                            if (score > best_score) {
                                best_score = score;
                                best_option = opt;
                            }
                        }
                    }

                    const int cand_new_cell_idx =
                        options[static_cast<size_t>(best_option)].cell_idx;
                    const int cand_new_facet_idx =
                        options[static_cast<size_t>(best_option)].facet_idx;

                    // All validation passed — accept this pair.
                    replace_i = ii;
                    replace_j = jj;
                    shared_tet = cand_tet;
                    new_cell_idx = cand_new_cell_idx;
                    new_facet_idx = cand_new_facet_idx;
                }
            }
            if (replace_i >= 0) break;
        }

        // Dump iteration state BEFORE the replacement (so we see the pair being combined).
        if (dump_file) {
            dump_iteration_state(iteration, result, replace_i, replace_j, shared_tet,
                                 mirror_cells, mirror_facet_indices);
        }

        if (replace_i < 0) break;

        if (trace) {
            const auto& fi = result[static_cast<size_t>(replace_i)];
            const auto& fj = result[static_cast<size_t>(replace_j)];
            DEBUG_PRINT("[DEL-SIMPLIFY-CYCLE] v=" << v0->info().index
                        << " cycle=" << cycle_index
                        << " replacing facets (" << fi.first << "," << fi.second << ")"
                        << " and (" << fj.first << "," << fj.second << ")"
                        << " via shared tet (cell " << shared_tet->info().index << ")"
                        << " with (" << new_cell_idx << "," << new_facet_idx << ")"
                        << " facet_count " << result.size() << " -> " << (result.size() - 1));
        }

        if (dump_file && trace_out) {
            const auto& fi = result[static_cast<size_t>(replace_i)];
            const auto& fj = result[static_cast<size_t>(replace_j)];
            trace_out << "replace f=(" << fi.first << "," << fi.second << ")"
                      << " f'=(" << fj.first << "," << fj.second << ")"
                      << " shared_tet=" << shared_tet->info().index
                      << " f''=(" << new_cell_idx << "," << new_facet_idx << ")\n\n";
        }

        // Perform the replacement.
        int hi = std::max(replace_i, replace_j);
        int lo = std::min(replace_i, replace_j);
        result.erase(result.begin() + hi);
        result.erase(result.begin() + lo);
        result.emplace_back(new_cell_idx, new_facet_idx);
        changed = true;
        ++iteration;
    }

    // Always dump the final state with incident tet info.
    if (dump_file && trace_out) {
        // Recompute mirror cells for the final state.
        const int fn = static_cast<int>(result.size());
        std::vector<Cell_handle> final_mirror_cells(static_cast<size_t>(fn), Cell_handle());
        std::vector<int> final_mirror_fi(static_cast<size_t>(fn), -1);
        for (int i = 0; i < fn; ++i) {
            const auto& [cell_idx, facet_idx] = result[static_cast<size_t>(i)];
            Cell_handle cell = lookup_cell(cell_by_index, cell_idx);
            if (cell == Cell_handle() || dt.is_infinite(cell)) continue;
            if (facet_idx < 0 || facet_idx >= 4) continue;
            const Facet mirror = dt.mirror_facet(Facet(cell, facet_idx));
            if (mirror.first != Cell_handle() && !dt.is_infinite(mirror.first)) {
                final_mirror_cells[static_cast<size_t>(i)] = mirror.first;
                final_mirror_fi[static_cast<size_t>(i)] = mirror.second;
            }
        }
        dump_iteration_state(iteration, result, -1, -1, Cell_handle(),
                             final_mirror_cells, final_mirror_fi);
        trace_out << "final_facet_count " << result.size() << "\n";
        trace_out << "total_iterations " << iteration << "\n";
        trace_out << "simplified " << (result.size() < cycle_facets.size() ? "yes" : "no") << "\n";
    }

    if (trace && result.size() < cycle_facets.size()) {
        DEBUG_PRINT("[DEL-SIMPLIFY-CYCLE] v=" << v0->info().index
                    << " cycle=" << cycle_index
                    << " simplified " << cycle_facets.size()
                    << " -> " << result.size() << " facets");
    }

    return result;
}

//! @brief Compute a direction that separates new cycle facets from old cycle facets.
//!
//! For each facet in the cycle, computes the unit outward normal (pointing toward the
//! positive cell, i.e., where circumcenter_scalar >= isovalue). Then uses CGAL's
//! Min_sphere_of_spheres_d to find the smallest enclosing sphere of these normals
//! (treated as points). The center of this sphere gives the separation direction.
//! Before min-sphere, normals may be globally flipped for this cycle if the non-symmetric
//! cycle-separation matrix indicates this cycle is not separated from any other cycle.
//!
//! @param v0 The vertex around which the cycle is centered.
//! @param cycle_facets The facets in the cycle as (cell_index, facet_index) pairs.
//! @param cell_by_index Lookup table for Cell_handle by index.
//! @param dt The Delaunay triangulation.
//! @param[out] separation_direction The computed separation direction (unit vector).
//! @return true if a valid separation direction was found, false otherwise.
bool compute_cycle_separating_direction(
    Vertex_handle v0,
    const std::vector<std::pair<int, int>>& cycle_facets,
    const std::vector<Cell_handle>& cell_by_index,
    const Delaunay& dt,
    int cycle_index,
    float separation_direction[3],
    const std::string& trace_suffix
) {
    if (v0 == Vertex_handle()) {
        return false;
    }
    if (cycle_facets.empty()) {
        return false;
    }

    const bool trace = vdc_debug::trace_selfi_enabled(v0->info().index);
    const std::string trace_dir = trace ? vdc_debug::trace_selfi_dump_dir() : std::string();
    const bool dump_normals = trace && !trace_dir.empty();
    const auto& all_cycles = v0->info().facet_cycles;
    const int num_cycles = static_cast<int>(all_cycles.size());

    // Non-symmetric cycle separation matrix (user-requested behavior):
    // For this cycle row i, if any j != i is not separated, flip all facet normals
    // before min-sphere so the resulting direction moves cycle i away from others.
    std::vector<std::vector<char>> is_separated;
    std::vector<int> separation_fail_targets;
    bool flip_normals_for_cycle = false;
    if (cycle_index >= 0 && cycle_index < num_cycles && num_cycles >= 2) {
        is_separated = determine_cycle_separation_matrix(v0, all_cycles, cell_by_index, dt);
        for (int j = 0; j < num_cycles; ++j) {
            if (j == cycle_index) {
                continue;
            }
            if (!is_separated.empty() &&
                !is_separated[static_cast<size_t>(cycle_index)][static_cast<size_t>(j)]) {
                flip_normals_for_cycle = true;
                separation_fail_targets.push_back(j);
            }
        }
    }

    // Collect unit facet normals for each facet in the cycle.
    // Each facet is stored with its positive cell (where scalar >= isovalue).
    // Orient each normal to point toward the positive cell side of the facet plane.
    std::vector<Point> normals_as_points;
    normals_as_points.reserve(cycle_facets.size());
    int foldover_edge_candidates_ge4 = 0;
    const CycleFoldoverEdgeInfo crossover_edge_dir = detect_cycle_crossover_edge_direction(
        v0, cycle_facets, cell_by_index, dt, &foldover_edge_candidates_ge4);


    using vdc_debug::FacetNormalDumpEntry;

    std::vector<FacetNormalDumpEntry> dump_entries;
    if (dump_normals) {
        dump_entries.reserve(cycle_facets.size());
    }

    int circumcenter_across_count = 0;
    int normals_point_to_pos_count = 0;
    int normals_point_to_neg_count = 0;

    for (const auto& [cell_idx, facet_idx] : cycle_facets) {
        Cell_handle cell = lookup_cell(cell_by_index, cell_idx);
        if (cell == Cell_handle()) {
            continue;
        }
        if (dt.is_infinite(cell)) {
            continue;
        }
        if (facet_idx < 0 || facet_idx >= 4) {
            continue;
        }

        // Get the three vertices of the facet (opposite to facet_idx vertex).
        std::array<Vertex_handle, 3> facet_vhs;
        std::array<Point, 3> facet_verts;
        bool ok = true;
        for (int k = 0; k < 3; ++k) {
            const int cell_vertex_idx = (facet_idx + 1 + k) % 4;
            Vertex_handle vh = cell->vertex(cell_vertex_idx);
            if (vh == Vertex_handle() || dt.is_infinite(vh)) {
                ok = false;
                break;
            }
            facet_vhs[static_cast<size_t>(k)] = vh;
            facet_verts[static_cast<size_t>(k)] = vh->point();
        }
        if (!ok) {
            continue;
        }

        // Compute the facet normal using cross product
        const Vector3 edge1(facet_verts[0], facet_verts[1]);
        const Vector3 edge2(facet_verts[0], facet_verts[2]);
        Vector3 normal = vec_cross(edge1, edge2);

        const double normal_len = vec_norm(normal);
        if (normal_len < 1e-12) {
            continue;  // Degenerate facet
        }
        normal = normal / normal_len;  // Normalize

        // Orient the normal so it points toward `cell` (by FacetKey convention this
        // is the positive cell). Use a point guaranteed to be on the cell side of
        // the facet plane: the vertex opposite this facet within the cell.
        // (Using the cell circumcenter is incorrect for obtuse/sliver cells where
        // the circumcenter may lie outside and across the facet plane.)
        Vertex_handle opposite_vh = cell->vertex(facet_idx);
        if (opposite_vh == Vertex_handle() || dt.is_infinite(opposite_vh)) {
            continue;
        }
        const Vector3 to_opposite(facet_verts[0], opposite_vh->point());
        double dot_to_opposite = vec_dot(normal, to_opposite);
        if (dot_to_opposite < 0) {
            normal = -normal;
            dot_to_opposite = -dot_to_opposite;
        }

        const Vector3 to_circumcenter(facet_verts[0], cell->info().circumcenter);
        const double dot_to_circumcenter = vec_dot(normal, to_circumcenter);

        if (trace) {
            // Diagnostics: count how often the circumcenter lies across the facet plane
            // from the (guaranteed) cell side. This was a root cause of misoriented
            // normals in the earlier circumcenter-based implementation.
            if (dot_to_circumcenter < -1e-12) {
                ++circumcenter_across_count;
            }
        }

        const int cell_sign = cell->info().flag_positive ? +1 : -1;
        int neighbor_cell_index = -1;
        int neighbor_facet_index = -1;
        int neighbor_sign = 0;
        const Facet mirror = dt.mirror_facet(Facet(cell, facet_idx));
        Cell_handle neighbor_cell = mirror.first;
        neighbor_facet_index = mirror.second;
        if (neighbor_cell != Cell_handle() && !dt.is_infinite(neighbor_cell)) {
            neighbor_cell_index = neighbor_cell->info().index;
            neighbor_sign = neighbor_cell->info().flag_positive ? +1 : -1;
        }

        bool normal_flip_applied = false;
        if (flip_normals_for_cycle) {
            normal = -normal;
            dot_to_opposite = -dot_to_opposite;
            normal_flip_applied = true;
        }

        int normal_points_to_cell_index = cell_idx;
        int normal_points_to_sign = cell_sign;
        int normal_points_away_cell_index = neighbor_cell_index;
        int normal_points_away_sign = neighbor_sign;
        if (normal_flip_applied) {
            normal_points_to_cell_index = neighbor_cell_index;
            normal_points_to_sign = neighbor_sign;
            normal_points_away_cell_index = cell_idx;
            normal_points_away_sign = cell_sign;
        }

        if (normal_points_to_sign > 0) {
            ++normals_point_to_pos_count;
        } else if (normal_points_to_sign < 0) {
            ++normals_point_to_neg_count;
        }

        // Store the unit normal as a point (vector from origin)
        normals_as_points.emplace_back(normal.x(), normal.y(), normal.z());

        if (dump_normals) {
            FacetNormalDumpEntry e;
            e.cell_index = cell_idx;
            e.facet_index = facet_idx;
            e.cell_sign = cell_sign;
            e.neighbor_cell_index = neighbor_cell_index;
            e.neighbor_facet_index = neighbor_facet_index;
            e.neighbor_sign = neighbor_sign;
            e.normal_points_to_cell_index = normal_points_to_cell_index;
            e.normal_points_to_sign = normal_points_to_sign;
            e.normal_points_away_cell_index = normal_points_away_cell_index;
            e.normal_points_away_sign = normal_points_away_sign;
            e.normal_flip_applied = normal_flip_applied;
            for (int k = 0; k < 3; ++k) {
                e.facet_p[k] = facet_verts[static_cast<size_t>(k)];
                e.facet_vidx[k] = facet_vhs[static_cast<size_t>(k)]->info().index;
            }
            e.opposite_vidx = opposite_vh->info().index;
            e.opposite_p = opposite_vh->point();
            e.nx = normal.x();
            e.ny = normal.y();
            e.nz = normal.z();
            e.dot_to_opposite = dot_to_opposite;
            e.dot_to_circumcenter = dot_to_circumcenter;
            dump_entries.push_back(e);
        }
    }

    if (normals_as_points.empty()) {
        return false;
    }

    // Use CGAL Min_sphere to compute the smallest enclosing sphere of the normals.
    // The center of this sphere gives the "average" direction.
    MinSphere min_sphere(normals_as_points.begin(), normals_as_points.end());

    if (min_sphere.is_empty()) {
        return false;
    }

    // Extract center coordinates
    auto center_it = min_sphere.center_cartesian_begin();
    const double cx = *center_it++;
    const double cy = *center_it++;
    const double cz = *center_it;

    const double center_len = std::sqrt(cx * cx + cy * cy + cz * cz);
    constexpr double eps = 1e-6;

    if (dump_normals) {
        vdc_debug::SepDirDumpContext dump_ctx;
        dump_ctx.v0 = v0;
        dump_ctx.cycle_index = cycle_index;
        dump_ctx.num_cycles = num_cycles;
        dump_ctx.cycle_facets_size = cycle_facets.size();
        dump_ctx.normals_used = normals_as_points.size();
        dump_ctx.circumcenter_across_count = circumcenter_across_count;
        dump_ctx.normals_point_to_pos_count = normals_point_to_pos_count;
        dump_ctx.normals_point_to_neg_count = normals_point_to_neg_count;
        dump_ctx.flip_normals_for_cycle = flip_normals_for_cycle;
        dump_ctx.is_separated = &is_separated;
        dump_ctx.separation_fail_targets = &separation_fail_targets;
        dump_ctx.cx = cx;
        dump_ctx.cy = cy;
        dump_ctx.cz = cz;
        dump_ctx.center_len = center_len;
        dump_ctx.foldover_edge_candidates_ge4 = foldover_edge_candidates_ge4;
        dump_ctx.crossover_found = crossover_edge_dir.found;
        dump_ctx.crossover_v_low = crossover_edge_dir.v_low;
        dump_ctx.crossover_v_high = crossover_edge_dir.v_high;
        dump_ctx.crossover_count = crossover_edge_dir.count;
        dump_ctx.crossover_incident_on_v0 = crossover_edge_dir.incident_on_v0;
        dump_ctx.crossover_other_vidx = crossover_edge_dir.other_vidx;
        dump_ctx.crossover_dir_x = crossover_edge_dir.unit_dir.x();
        dump_ctx.crossover_dir_y = crossover_edge_dir.unit_dir.y();
        dump_ctx.crossover_dir_z = crossover_edge_dir.unit_dir.z();
        dump_ctx.normals_as_points = &normals_as_points;
        vdc_debug::dump_sep_dir_normals_trace(trace_dir, dump_ctx, dump_entries, trace_suffix);
    }

    // If center is very close to origin, normals span all directions - no separation
    if (center_len < eps) {
        if (trace) {
            DEBUG_PRINT("[DEL-SEP-DIR] v=" << v0->info().index
                        << " facets=" << cycle_facets.size()
                        << " used=" << normals_as_points.size()
                        << " circumcenter_across=" << circumcenter_across_count
                        << "/" << normals_as_points.size()
                        << " flip_normals=" << (flip_normals_for_cycle ? "YES" : "NO")
                        << " center_len=" << center_len
                        << " -> NO_DIRECTION");
        }
        return false;
    }

    // Normalize and output
    separation_direction[0] = static_cast<float>(cx / center_len);
    separation_direction[1] = static_cast<float>(cy / center_len);
    separation_direction[2] = static_cast<float>(cz / center_len);

    if (trace) {
        DEBUG_PRINT("[DEL-SEP-DIR] v=" << v0->info().index
                    << " facets=" << cycle_facets.size()
                    << " used=" << normals_as_points.size()
                    << " circumcenter_across=" << circumcenter_across_count
                    << "/" << normals_as_points.size()
                    << " flip_normals=" << (flip_normals_for_cycle ? "YES" : "NO")
                    << " center_len=" << center_len
                    << " dir=("
                    << separation_direction[0] << ","
                    << separation_direction[1] << ","
                    << separation_direction[2] << ")");
    }

    return true;
}

static bool is_cycle_in_positive_half_space(
    const Point& p,
    const Vector3& normal,
    const std::vector<Vertex_handle>& cycle_vertices
) {
    constexpr double eps = 1e-12;
    for (Vertex_handle vh : cycle_vertices) {
        const Vector3 u(p, vh->point());
        if (vec_dot(u, normal) <= eps) {
            return false;
        }
    }
    return true;
}

static bool is_cycle_in_negative_half_space(
    const Point& p,
    const Vector3& normal,
    const std::vector<Vertex_handle>& cycle_vertices
) {
    return is_cycle_in_positive_half_space(p, -normal, cycle_vertices);
}

static bool bisecting_plane_separates_cycle_boundaries(
    const CycleBoundaryInfo& cycleA,
    const CycleBoundaryInfo& cycleB,
    const Point& vposA,
    const Point& vposB
) {
    const Vector3 normal_raw(vposB, vposA);
    const double len = vec_norm(normal_raw);
    if (len < 1e-12) {
        return false;
    }

    const Vector3 normalA = normal_raw / len;
    const Point posC(
        0.5 * (vposA.x() + vposB.x()),
        0.5 * (vposA.y() + vposB.y()),
        0.5 * (vposA.z() + vposB.z()));

    if (!is_cycle_in_positive_half_space(posC, normalA, cycleA.boundary_verts)) {
        return false;
    }
    if (!is_cycle_in_negative_half_space(posC, normalA, cycleB.boundary_verts)) {
        return false;
    }

    return true;
}

static Point reflect_through_center(const Point& center, const Point& p) {
    return Point(
        2.0 * center.x() - p.x(),
        2.0 * center.y() - p.y(),
        2.0 * center.z() - p.z());
}

static bool try_separate_cycle_pair_positions_A(
    Vertex_handle v0,
    const CycleBoundaryInfo& cycleA,
    const CycleBoundaryInfo& cycleB,
    int cycleA_index,
    int cycleB_index,
    const Point& vposA,
    const Point& vposB,
    Point& new_vposA,
    Point& new_vposB
) {
    new_vposA = vposA;
    new_vposB = vposB;

    const bool trace = vdc_debug::trace_selfi_enabled(v0->info().index);

    // First check: do current positions already satisfy bisecting plane separation?
    if (bisecting_plane_separates_cycle_boundaries(cycleA, cycleB, vposA, vposB)) {
        if (trace) {
            DEBUG_PRINT("[DEL-SELFI-A] v=" << v0->info().index
                        << " pair(" << cycleA_index << "," << cycleB_index << ") bisect-plane: PASS");
        }
        return true;
    }

    if (trace) {
        DEBUG_PRINT("[DEL-SELFI-A] v=" << v0->info().index
                    << " pair(" << cycleA_index << "," << cycleB_index << ") bisect-plane: FAIL");

        const Vector3 normal_raw(vposB, vposA);
        const double len = vec_norm(normal_raw);
        if (len < 1e-12) {
            DEBUG_PRINT("[DEL-SELFI-A]   reason: |vposA-vposB| ~ 0 (degenerate normal)");
        } else {
            const Vector3 normalA = normal_raw / len;
            const Point posC(
                0.5 * (vposA.x() + vposB.x()),
                0.5 * (vposA.y() + vposB.y()),
                0.5 * (vposA.z() + vposB.z()));

            constexpr double eps = 1e-12;
            double min_dot_A = std::numeric_limits<double>::infinity();
            int min_dot_A_vidx = -1;
            int bad_A = 0;
            for (Vertex_handle vh : cycleA.boundary_verts) {
                const Vector3 u(posC, vh->point());
                const double d = vec_dot(u, normalA);
                if (d <= eps) {
                    ++bad_A;
                }
                if (d < min_dot_A) {
                    min_dot_A = d;
                    min_dot_A_vidx = vh->info().index;
                }
            }

            double max_dot_B = -std::numeric_limits<double>::infinity();
            int max_dot_B_vidx = -1;
            int bad_B = 0;
            for (Vertex_handle vh : cycleB.boundary_verts) {
                const Vector3 u(posC, vh->point());
                const double d = vec_dot(u, normalA);
                if (d >= -eps) { // should be strictly negative
                    ++bad_B;
                }
                if (d > max_dot_B) {
                    max_dot_B = d;
                    max_dot_B_vidx = vh->info().index;
                }
            }

            DEBUG_PRINT("[DEL-SELFI-A]   cycleA boundary halfspace: bad=" << bad_A
                        << "/" << cycleA.boundary_verts.size()
                        << ", min_dot=" << min_dot_A
                        << " at delv=" << min_dot_A_vidx);
            DEBUG_PRINT("[DEL-SELFI-A]   cycleB boundary halfspace: bad=" << bad_B
                        << "/" << cycleB.boundary_verts.size()
                        << ", max_dot=" << max_dot_B
                        << " at delv=" << max_dot_B_vidx);
        }
    }

    const Point center = v0->point();

    // Try simple reflection through center (PositionMultiIsov.txt: find_vertex_positions_separating_cycles)
    const Point vposB2 = reflect_through_center(center, vposA);
    const Point vposA2 = reflect_through_center(center, vposB);
    const bool sepB = bisecting_plane_separates_cycle_boundaries(cycleA, cycleB, vposA, vposB2);
    const bool sepA = bisecting_plane_separates_cycle_boundaries(cycleA, cycleB, vposA2, vposB);

    if (trace) {
        DEBUG_PRINT("[DEL-SELFI-A]   reflect-B: bisect-plane-after=" << (sepB ? "PASS" : "FAIL"));
        DEBUG_PRINT("[DEL-SELFI-A]   reflect-A: bisect-plane-after=" << (sepA ? "PASS" : "FAIL"));
    }

    // Prefer a reflection that makes the bisect-plane separation test pass, but do not
    // enforce it as a necessary condition (the in-code self-intersection check is the
    // ground truth acceptance test).
    if (sepB || !sepA) {
        new_vposB = vposB2;
    } else {
        new_vposA = vposA2;
    }

    return true;
}

// ============================================================================
// Multi-cycle resolution
// ============================================================================

static ResolutionResult try_resolve_multicycle_by_cycle_separation_tests(
    Vertex_handle v,
    std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    std::vector<Point>* geometric_attempt_positions_out,
    bool use_sep_dir,
    const StrictTetPreservationContext& strict_context
) {
    const auto& cycles = v->info().facet_cycles;
    const int num_cycles = static_cast<int>(cycles.size());
    if (num_cycles < 2 || static_cast<int>(cycle_isovertices.size()) < 2) {
        return {ResolutionStatus::NOT_NEEDED, ResolutionStrategy::NONE};
    }

    std::vector<Point> baseline_positions = cycle_isovertices;
    if (apply_dihedral_orientation_gate_to_positions(
            dt, cell_by_index, v, baseline_positions, strict_context) > 0) {
        cycle_isovertices = baseline_positions;
    }

    std::vector<CycleBoundaryInfo> cycle_data;
    cycle_data.reserve(static_cast<size_t>(num_cycles));
    for (int c = 0; c < num_cycles; ++c) {
        cycle_data.push_back(gather_cycle_boundary_info(dt, cell_by_index, v, cycles[c]));
    }

    // PositionMultiIsov(A) is intentionally restrictive (diametric reflections only).
    // For 2-cycle vertices, try both reflection directions and pick the one that actually
    // resolves (or strictly reduces) the local self-intersection count (tie-break by minimal movement).
	    if (num_cycles == 2) {
	        const bool trace = vdc_debug::trace_selfi_enabled(v->info().index);
	        const Point center = v->point();

	        struct CandidateA2 {
	            std::vector<Point> positions;
	            const char* label = nullptr;
	            double move2 = 0.0;
	            int intersections = std::numeric_limits<int>::max();
	        };

		        std::vector<CandidateA2> candidates;
		        if (use_sep_dir) {
		            // One directional candidate per cycle.
		            candidates.reserve(2);
		        } else {
		            // Legacy reflection mode for 2-cycle vertices.
		            candidates.reserve(6);
		        }

        const Point vpos0 = baseline_positions[0];
        const Point vpos1 = baseline_positions[1];

        auto consider_positions = [&](const char* label, std::vector<Point> positions) {
            CandidateA2 cand;
            cand.positions = std::move(positions);
            (void)apply_dihedral_orientation_gate_to_positions(
                dt, cell_by_index, v, cand.positions, strict_context);
            cand.label = label;
            cand.move2 = 0.0;
            for (size_t i = 0; i < cand.positions.size() && i < baseline_positions.size(); ++i) {
                cand.move2 += squared_distance(cand.positions[i], baseline_positions[i]);
            }
            candidates.push_back(std::move(cand));
        };

        auto add_reflection_candidates = [&]() {
            // reflect-B: move cycle1 to the diametric position of cycle0
            {
                std::vector<Point> cand = baseline_positions;
                cand[1] = reflect_through_center(center, vpos0);
                consider_positions("reflect_B", std::move(cand));
            }
            // reflect-A: move cycle0 to the diametric position of cycle1
            {
                std::vector<Point> cand = baseline_positions;
                cand[0] = reflect_through_center(center, vpos1);
                consider_positions("reflect_A", std::move(cand));
            }

            // Also try self reflections (still diametric through the Delaunay site).
            {
                std::vector<Point> cand = baseline_positions;
                cand[0] = reflect_through_center(center, vpos0);
                consider_positions("reflect_self0", std::move(cand));
            }
            {
                std::vector<Point> cand = baseline_positions;
                cand[1] = reflect_through_center(center, vpos1);
                consider_positions("reflect_self1", std::move(cand));
            }
            {
                std::vector<Point> cand = baseline_positions;
                cand[0] = reflect_through_center(center, vpos0);
                cand[1] = reflect_through_center(center, vpos1);
                consider_positions("reflect_self_both", std::move(cand));
            }
            {
                std::vector<Point> cand = baseline_positions;
                cand[0] = reflect_through_center(center, vpos1);
                cand[1] = reflect_through_center(center, vpos0);
                consider_positions("reflect_cross_both", std::move(cand));
            }
        };

	        auto add_sep_dir_candidates = [&]() {
	            struct SepCycleState {
	                bool has_foldover_edge = false;
	                CycleFoldoverEdgeInfo foldover;
	                bool has_sep_dir = false;
	                float sep_dir[3] = {0, 0, 0};
	            };

	            std::array<SepCycleState, 2> st;
	            int fold_candidates0 = 0;
	            int fold_candidates1 = 0;
	            st[0].foldover = detect_cycle_crossover_edge_direction(
	                v, cycles[0], cell_by_index, dt, &fold_candidates0);
	            st[1].foldover = detect_cycle_crossover_edge_direction(
	                v, cycles[1], cell_by_index, dt, &fold_candidates1);
	            st[0].has_foldover_edge = st[0].foldover.found;
	            st[1].has_foldover_edge = st[1].foldover.found;

	            if (!st[0].has_foldover_edge) {
	                st[0].has_sep_dir = compute_cycle_separating_direction(
	                    v, cycles[0], cell_by_index, dt, 0, st[0].sep_dir);
	            }
	            if (!st[1].has_foldover_edge) {
	                st[1].has_sep_dir = compute_cycle_separating_direction(
	                    v, cycles[1], cell_by_index, dt, 1, st[1].sep_dir);
	            }

	            // Keep the same per-cycle radius convention as existing sep-dir candidates.
	            const double r0 = std::sqrt(squared_distance(center, vpos0));
	            const double r1 = std::sqrt(squared_distance(center, vpos1));

	            auto move_in_dir = [&center](const float dir[3], double r) -> Point {
	                return Point(center.x() + r * dir[0], center.y() + r * dir[1], center.z() + r * dir[2]);
	            };
	            auto move_in_vec = [&center](const Vector3& dir, double r) -> Point {
	                return Point(center.x() + r * dir.x(), center.y() + r * dir.y(), center.z() + r * dir.z());
	            };

	            const int fold_count =
	                static_cast<int>(st[0].has_foldover_edge) + static_cast<int>(st[1].has_foldover_edge);

	            if (trace) {
	                DEBUG_PRINT("[DEL-SELFI-A2] v=" << v->info().index
	                            << " sep_dir_precheck: "
	                            << "c0 foldover=" << (st[0].has_foldover_edge ? "yes" : "no")
	                            << " candidates_ge4=" << fold_candidates0
	                            << " has_sep_dir=" << (st[0].has_sep_dir ? "yes" : "no")
	                            << " | c1 foldover=" << (st[1].has_foldover_edge ? "yes" : "no")
	                            << " candidates_ge4=" << fold_candidates1
	                            << " has_sep_dir=" << (st[1].has_sep_dir ? "yes" : "no")
	                            << " | fold_count=" << fold_count);
	            }

	            if (fold_count == 2) {
	                // Both cycles have foldover edges: move both along their own shared-edge directions.
	                std::vector<Point> cand = baseline_positions;
	                cand[0] = move_in_vec(st[0].foldover.unit_dir, r0);
	                cand[1] = move_in_vec(st[1].foldover.unit_dir, r1);
	                consider_positions("foldover_both_edges", std::move(cand));
	                return;
	            }

	            if (fold_count == 1) {
	                // Exactly one foldover cycle: keep it at DelV center, move the other by sep-dir.
	                const int fold_idx = st[0].has_foldover_edge ? 0 : 1;
	                const int other_idx = 1 - fold_idx;
	                if (st[other_idx].has_sep_dir) {
	                    std::vector<Point> cand = baseline_positions;
	                    cand[static_cast<size_t>(fold_idx)] = center;
	                    const double r = (other_idx == 0) ? r0 : r1;
	                    cand[static_cast<size_t>(other_idx)] = move_in_dir(st[other_idx].sep_dir, r);
	                    consider_positions("foldover_one_keep_center_sep_other", std::move(cand));
	                }
	                return;
	            }

	            // No crossover edges detected: standard sep-dir candidates.
	            if (st[0].has_sep_dir) {
	                std::vector<Point> cand = baseline_positions;
	                cand[0] = move_in_dir(st[0].sep_dir, r0);
	                consider_positions("sep_dir0", std::move(cand));
	            }
	            if (st[1].has_sep_dir) {
	                std::vector<Point> cand = baseline_positions;
	                cand[1] = move_in_dir(st[1].sep_dir, r1);
	                consider_positions("sep_dir1", std::move(cand));
	            }
	        };

        // ====================================================================
        // Candidate generation
        // ====================================================================
        if (use_sep_dir) {
            add_sep_dir_candidates();
        } else {
            add_reflection_candidates();
        }

        const int baseline_intersections = count_self_intersection_pairs(
            v, baseline_positions, dt, cell_by_index);
        if (trace) {
            DEBUG_PRINT("[DEL-SELFI-A2] v=" << v->info().index
                        << " baseline intersections=" << baseline_intersections);
        }

        CandidateA2* best_resolved = nullptr;
        CandidateA2* best_improved = nullptr;
        CandidateA2* best_any = nullptr;

        auto evaluate_candidate = [&](CandidateA2& cand) {
            const bool bisect_ok = bisecting_plane_separates_cycle_boundaries(
                cycle_data[0], cycle_data[1], cand.positions[0], cand.positions[1]);
            cand.intersections = count_self_intersection_pairs(
                v, cand.positions, dt, cell_by_index);
            const bool resolved = (cand.intersections == 0);
            const bool improved =
                (!resolved) && (baseline_intersections > 0) && (cand.intersections < baseline_intersections);
            if (trace) {
                DEBUG_PRINT("[DEL-SELFI-A2] v=" << v->info().index
                            << " " << cand.label << ": "
                            << "bisect-plane-after=" << (bisect_ok ? "PASS" : "FAIL")
                            << " intersections=" << cand.intersections
                            << "/" << baseline_intersections
                            << " "
                            << (resolved ? "RESOLVED" :
                                    (improved ? "IMPROVED" : "NO_IMPROVEMENT"))
                            << " move2=" << cand.move2);
            }
            if (!best_any ||
                cand.intersections < best_any->intersections ||
                (cand.intersections == best_any->intersections &&
                 cand.move2 < best_any->move2)) {
                best_any = &cand;
            }
            if (resolved) {
                if (!best_resolved || cand.move2 < best_resolved->move2) {
                    best_resolved = &cand;
                }
                return;
            }
            if (improved) {
                if (!best_improved ||
                    cand.intersections < best_improved->intersections ||
                    (cand.intersections == best_improved->intersections &&
                     cand.move2 < best_improved->move2)) {
                    best_improved = &cand;
                }
            }
        };

        for (auto& cand : candidates) {
            evaluate_candidate(cand);
        }

        if (best_resolved) {
            cycle_isovertices = best_resolved->positions;
            if (geometric_attempt_positions_out) {
                *geometric_attempt_positions_out = best_resolved->positions;
            }
            return {ResolutionStatus::RESOLVED, ResolutionStrategy::GEOMETRIC_SEPARATION};
        }

        if (best_improved) {
            cycle_isovertices = best_improved->positions;
            if (geometric_attempt_positions_out) {
                *geometric_attempt_positions_out = best_improved->positions;
            }
            return {ResolutionStatus::UNRESOLVED, ResolutionStrategy::GEOMETRIC_SEPARATION};
        }

        if (geometric_attempt_positions_out) {
            if (best_any) {
                *geometric_attempt_positions_out = best_any->positions;
            } else {
                *geometric_attempt_positions_out = baseline_positions;
            }
        }

        cycle_isovertices = baseline_positions;
        return {ResolutionStatus::UNRESOLVED, ResolutionStrategy::GEOMETRIC_SEPARATION};
    }

    bool any_changed = false;
    const double change_eps2 = 1e-6;
    const int MAX_PAIRWISE_PASSES = 6;

    for (int pass = 0; pass < MAX_PAIRWISE_PASSES; ++pass) {
        bool pass_changed = false;

        for (int a = 0; a < num_cycles; ++a) {
            for (int b = a + 1; b < num_cycles; ++b) {
                Point newA = cycle_isovertices[static_cast<size_t>(a)];
                Point newB = cycle_isovertices[static_cast<size_t>(b)];

                if (!try_separate_cycle_pair_positions_A(
                        v,
                        cycle_data[static_cast<size_t>(a)],
                        cycle_data[static_cast<size_t>(b)],
                        a,
                        b,
                        cycle_isovertices[static_cast<size_t>(a)],
                        cycle_isovertices[static_cast<size_t>(b)],
                        newA,
                        newB)) {
                    continue;
                }

                (void)apply_dihedral_orientation_gate_to_cycle_position(
                    dt, cell_by_index, v, a, newA, strict_context);
                (void)apply_dihedral_orientation_gate_to_cycle_position(
                    dt, cell_by_index, v, b, newB, strict_context);

                if (squared_distance(newA, cycle_isovertices[static_cast<size_t>(a)]) > change_eps2) {
                    cycle_isovertices[static_cast<size_t>(a)] = newA;
                    pass_changed = true;
                }
                if (squared_distance(newB, cycle_isovertices[static_cast<size_t>(b)]) > change_eps2) {
                    cycle_isovertices[static_cast<size_t>(b)] = newB;
                    pass_changed = true;
                }
            }
        }

        if (pass_changed) {
            any_changed = true;
        } else {
            break;
        }
    }

    // Validation: accept if this resolves all local self-intersections.
    if (!check_self_intersection(v, cycle_isovertices, dt, cell_by_index)) {
        if (geometric_attempt_positions_out) {
            *geometric_attempt_positions_out = cycle_isovertices;
        }
        if (!any_changed) {
            return {ResolutionStatus::NOT_NEEDED, ResolutionStrategy::NONE};
        }
        return {ResolutionStatus::RESOLVED, ResolutionStrategy::GEOMETRIC_SEPARATION};
    }

	    // If pairwise reflections don't fully resolve, also try per-cycle reflection candidates
	    // (still diametric through the Delaunay site) and keep a candidate that resolves or strictly
	    // reduces the number of local self-intersection pairs.
	        const bool trace = vdc_debug::trace_selfi_enabled(v->info().index);
	        const Point center = v->point();

        struct CandidateAN {
            std::vector<Point> positions;
            std::string label;
            double move2 = 0.0;
            int intersections = std::numeric_limits<int>::max();
        };

        const int baseline_intersections = count_self_intersection_pairs(
            v, baseline_positions, dt, cell_by_index);
        if (trace) {
            DEBUG_PRINT("[DEL-SELFI-A] v=" << v->info().index
                        << " baseline intersections=" << baseline_intersections);
        }

        auto compute_move2 = [&](const std::vector<Point>& pos) -> double {
            double move2 = 0.0;
            const size_t n = std::min(pos.size(), baseline_positions.size());
            for (size_t i = 0; i < n; ++i) {
                move2 += squared_distance(pos[i], baseline_positions[i]);
            }
            return move2;
        };

	        const size_t n_cycles = static_cast<size_t>(num_cycles);
	        const size_t cross_count = (num_cycles > 1) ? (n_cycles * (n_cycles - 1)) : 0;

	        std::vector<CandidateAN> candidates;
	        candidates.reserve(1 + 3 * n_cycles + cross_count);

	        // Candidate 0: the current pairwise-pass result.
	        {
	            CandidateAN cand;
	            cand.positions = cycle_isovertices;
	            cand.label = "pairwise_attempt";
	            cand.move2 = compute_move2(cand.positions);
	            candidates.push_back(std::move(cand));
	        }

	        // Candidates: move one cycle to the Delaunay site or halfway toward it.
	        for (int c = 0; c < num_cycles; ++c) {
	            {
	                CandidateAN cand;
	                cand.positions = baseline_positions;
	                cand.positions[static_cast<size_t>(c)] = center;
	                std::ostringstream label;
	                label << "center_cycle" << c;
	                cand.label = label.str();
	                cand.move2 = compute_move2(cand.positions);
	                candidates.push_back(std::move(cand));
	            }
	            {
	                CandidateAN cand;
	                cand.positions = baseline_positions;
	                const Point p0 = baseline_positions[static_cast<size_t>(c)];
	                cand.positions[static_cast<size_t>(c)] = Point(
	                    0.5 * (center.x() + p0.x()),
	                    0.5 * (center.y() + p0.y()),
	                    0.5 * (center.z() + p0.z()));
	                std::ostringstream label;
	                label << "shrink_half_cycle" << c;
	                cand.label = label.str();
	                cand.move2 = compute_move2(cand.positions);
	                candidates.push_back(std::move(cand));
	            }
	        }

        // Candidates: self-reflect each cycle individually from the baseline.
        for (int c = 0; c < num_cycles; ++c) {
            CandidateAN cand;
            cand.positions = baseline_positions;
            cand.positions[static_cast<size_t>(c)] =
                reflect_through_center(center, baseline_positions[static_cast<size_t>(c)]);
            std::ostringstream label;
            label << "reflect_cycle" << c;
            cand.label = label.str();
            cand.move2 = compute_move2(cand.positions);
            candidates.push_back(std::move(cand));
        }

        // Candidates: move one cycle to the diametric position of another cycle (baseline).
        for (int a = 0; a < num_cycles; ++a) {
            for (int b = 0; b < num_cycles; ++b) {
                if (a == b) continue;

                CandidateAN cand;
                cand.positions = baseline_positions;
                cand.positions[static_cast<size_t>(a)] =
                    reflect_through_center(center, baseline_positions[static_cast<size_t>(b)]);
                std::ostringstream label;
                label << "reflect_cycle" << a << "_from" << b;
                cand.label = label.str();
                cand.move2 = compute_move2(cand.positions);
                candidates.push_back(std::move(cand));
            }
        }

        CandidateAN* best_resolved = nullptr;
        CandidateAN* best_improved = nullptr;
        CandidateAN* best_any = nullptr;
        for (auto& cand : candidates) {
            (void)apply_dihedral_orientation_gate_to_positions(
                dt, cell_by_index, v, cand.positions, strict_context);
            cand.move2 = compute_move2(cand.positions);
            cand.intersections = count_self_intersection_pairs(v, cand.positions, dt, cell_by_index);
            const bool resolved = (cand.intersections == 0);
            const bool improved =
                (!resolved) && (baseline_intersections > 0) && (cand.intersections < baseline_intersections);
            if (trace) {
                DEBUG_PRINT("[DEL-SELFI-A] v=" << v->info().index
                            << " " << cand.label
                            << ": intersections=" << cand.intersections
                            << "/" << baseline_intersections
                            << " "
                            << (resolved ? "RESOLVED" :
                                    (improved ? "IMPROVED" : "NO_IMPROVEMENT"))
                            << " move2=" << cand.move2);
            }
            if (!best_any ||
                cand.intersections < best_any->intersections ||
                (cand.intersections == best_any->intersections &&
                 cand.move2 < best_any->move2)) {
                best_any = &cand;
            }
            if (resolved) {
                if (!best_resolved || cand.move2 < best_resolved->move2) {
                    best_resolved = &cand;
                }
                continue;
            }
            if (improved) {
                if (!best_improved ||
                    cand.intersections < best_improved->intersections ||
                    (cand.intersections == best_improved->intersections &&
                     cand.move2 < best_improved->move2)) {
                    best_improved = &cand;
                }
            }
        }

        if (best_resolved) {
            cycle_isovertices = best_resolved->positions;
            if (geometric_attempt_positions_out) {
                *geometric_attempt_positions_out = best_resolved->positions;
            }
            return {ResolutionStatus::RESOLVED, ResolutionStrategy::GEOMETRIC_SEPARATION};
        }

        if (best_improved) {
            cycle_isovertices = best_improved->positions;
            if (geometric_attempt_positions_out) {
                *geometric_attempt_positions_out = best_improved->positions;
            }
            return {ResolutionStatus::UNRESOLVED, ResolutionStrategy::GEOMETRIC_SEPARATION};
        }

        if (geometric_attempt_positions_out) {
            if (best_any) {
                *geometric_attempt_positions_out = best_any->positions;
            } else {
                *geometric_attempt_positions_out = cycle_isovertices;
            }
        }

        cycle_isovertices = baseline_positions;
        return {ResolutionStatus::UNRESOLVED, ResolutionStrategy::GEOMETRIC_SEPARATION};

    if (geometric_attempt_positions_out) {
        *geometric_attempt_positions_out = cycle_isovertices;
    }

    cycle_isovertices = baseline_positions;
    return {ResolutionStatus::UNRESOLVED, ResolutionStrategy::GEOMETRIC_SEPARATION};
}

ResolutionResult try_resolve_multicycle_by_cycle_separation_tests(
    Vertex_handle v,
    std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    std::vector<Point>* geometric_attempt_positions_out,
    bool use_sep_dir
) {
    StrictTetPreservationContext relaxed_context;
    return try_resolve_multicycle_by_cycle_separation_tests(
        v, cycle_isovertices, dt, cell_by_index,
        geometric_attempt_positions_out, use_sep_dir, relaxed_context);
}

static bool multicycle_is_definitely_separated_by_bisecting_planes(
    Vertex_handle v,
    const std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index
) {
    const auto& cycles = v->info().facet_cycles;
    const int num_cycles = static_cast<int>(cycles.size());
    if (num_cycles < 2 || static_cast<int>(cycle_isovertices.size()) < 2) {
        return true;
    }

    // Within-cycle self-intersections are only possible if a cycle fan has >= 4 triangles.
    for (int c = 0; c < num_cycles; ++c) {
        if (cycles[static_cast<size_t>(c)].size() >= 4) {
            return false;
        }
    }

    std::vector<CycleBoundaryInfo> cycle_data;
    cycle_data.reserve(static_cast<size_t>(num_cycles));
    for (int c = 0; c < num_cycles; ++c) {
        cycle_data.push_back(gather_cycle_boundary_info(dt, cell_by_index, v, cycles[c]));
    }

    for (int a = 0; a < num_cycles; ++a) {
        for (int b = a + 1; b < num_cycles; ++b) {
            if (!bisecting_plane_separates_cycle_boundaries(
                    cycle_data[static_cast<size_t>(a)],
                    cycle_data[static_cast<size_t>(b)],
                    cycle_isovertices[static_cast<size_t>(a)],
                    cycle_isovertices[static_cast<size_t>(b)])) {
                return false;
            }
        }
    }

    return true;
}

static ResolutionResult try_resolve_cycle_fan_foldover_at_vertex_cycle(
    Vertex_handle v,
    int cycle_idx,
    std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    double foldover_sphere_radius,
    const StrictTetPreservationContext& strict_context
);

static ResolutionResult resolve_multicycle_self_intersection_at_vertex(
    Vertex_handle v,
    std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    const CycleIsovertexOptions& options,
    double foldover_sphere_radius,
    std::unordered_set<int>* local_unresolved_dumped_vertices,
    const StrictTetPreservationContext& strict_context
) {
    const int num_cycles = static_cast<int>(v->info().facet_cycles.size());
    if (num_cycles < 2) {
        return {ResolutionStatus::NOT_NEEDED, ResolutionStrategy::NONE};
    }
    const bool trace = vdc_debug::trace_selfi_enabled(v->info().index);

    // Cheap sufficient-condition filter: if all cycles are separated by bisecting planes and
    // each cycle fan is too small to self-intersect internally, skip expensive triangle tests.
    if (multicycle_is_definitely_separated_by_bisecting_planes(
            v, cycle_isovertices, dt, cell_by_index)) {
        if (trace) {
            DEBUG_PRINT("[DEL-SELFI-TRACE] Vertex " << v->info().index
                        << ": skipped (bisect-plane sufficient condition).");
        }
        return {ResolutionStatus::NOT_NEEDED, ResolutionStrategy::NONE};
    }

    if (!check_self_intersection(v, cycle_isovertices, dt, cell_by_index)) {
        if (trace) {
            DEBUG_PRINT("[DEL-SELFI-TRACE] Vertex " << v->info().index
                        << ": no self-intersection detected.");
        }
        return {ResolutionStatus::NOT_NEEDED, ResolutionStrategy::NONE};
    }

    const std::vector<Point> baseline_positions = cycle_isovertices;
    const auto witness = find_self_intersection_witness(v, baseline_positions, dt, cell_by_index);
    if (trace) {
        DEBUG_PRINT("[DEL-SELFI-TRACE] Vertex " << v->info().index
                    << ": detected self-intersection (cycles=" << num_cycles << ").");
        for (int c = 0; c < num_cycles; ++c) {
            const Point& p = baseline_positions[static_cast<size_t>(c)];
            DEBUG_PRINT("[DEL-SELFI-TRACE]   baseline[" << c << "]=("
                        << p.x() << ", " << p.y() << ", " << p.z() << ")");
        }
        if (witness) {
            DEBUG_PRINT("[DEL-SELFI-TRACE]   witness: cycle(" << witness->cycle0 << "," << witness->cycle1
                        << ") tri(" << witness->tri0 << "," << witness->tri1 << ")");
        }
    }
    vdc_debug::maybe_dump_selfi_cycle_metadata(dt, v, cell_by_index);
    vdc_debug::maybe_dump_selfi_stage(dt, v, baseline_positions, cell_by_index, "baseline");

    // If the first detected witness is within a single cycle, treat this as a fan foldover
    // and attempt a deterministic, direction-based correction *before* running the between-cycle
    // separation logic. This targets cases where the multi-cycle separation candidates are
    // irrelevant (the self-intersection is entirely within one cycle's triangle fan).
    //
    // Note: The foldover resolver may still return UNRESOLVED if the vertex also has between-cycle
    // intersections (or if the deterministic move cannot eliminate the within-cycle foldover).
    // In that case we continue into the standard separation resolver using the updated positions.
    if (options.foldover && witness && witness->cycle0 == witness->cycle1) {
        const int fold_cycle = witness->cycle0;
        if (trace) {
            DEBUG_PRINT("[DEL-SELFI-TRACE] Vertex " << v->info().index
                        << ": within-cycle witness in cycle " << fold_cycle
                        << " -> attempting deterministic fan foldover resolution.");
        }

        const ResolutionResult foldover_result = try_resolve_cycle_fan_foldover_at_vertex_cycle(
            v, fold_cycle, cycle_isovertices, dt, cell_by_index, foldover_sphere_radius, strict_context);
        if (foldover_result.status == ResolutionStatus::RESOLVED) {
            DEBUG_PRINT("[DEL-SELFI] Vertex " << v->info().index
                        << " resolved by deterministic fan foldover (cycles=" << num_cycles << ").");
            vdc_debug::maybe_dump_selfi_stage(dt, v, cycle_isovertices, cell_by_index, "fan");
            return foldover_result;
        }
        if (foldover_result.status == ResolutionStatus::UNRESOLVED) {
            vdc_debug::maybe_dump_selfi_stage(dt, v, cycle_isovertices, cell_by_index, "fan");
        }
    }

    std::vector<Point> geometric_attempt_positions;
    const ResolutionResult separation_result =
        try_resolve_multicycle_by_cycle_separation_tests(
            v,
            cycle_isovertices,
            dt,
            cell_by_index,
            &geometric_attempt_positions,
            options.use_sep_dir,
            strict_context);
    if (separation_result.status == ResolutionStatus::RESOLVED) {
        DEBUG_PRINT("[DEL-SELFI] Vertex " << v->info().index
                    << " resolved by geometric separation (cycles=" << num_cycles << ").");
        if (trace) {
            for (int c = 0; c < num_cycles; ++c) {
                const Point& p = cycle_isovertices[static_cast<size_t>(c)];
                DEBUG_PRINT("[DEL-SELFI-TRACE]   geometric[" << c << "]=("
                            << p.x() << ", " << p.y() << ", " << p.z() << ")");
            }
            if (const auto w = find_self_intersection_witness(v, cycle_isovertices, dt, cell_by_index)) {
                DEBUG_PRINT("[DEL-SELFI-TRACE]   geometric witness: cycle(" << w->cycle0 << "," << w->cycle1
                            << ") tri(" << w->tri0 << "," << w->tri1 << ")");
            }
        }
        vdc_debug::maybe_dump_selfi_stage(dt, v, cycle_isovertices, cell_by_index, "A");
        return separation_result;
    }
    if (trace) {
        DEBUG_PRINT("[DEL-SELFI-TRACE] Vertex " << v->info().index
                    << ": geometric separation failed.");
        if (const auto w = find_self_intersection_witness(v, cycle_isovertices, dt, cell_by_index)) {
            DEBUG_PRINT("[DEL-SELFI-TRACE]   geometric witness: cycle(" << w->cycle0 << "," << w->cycle1
                        << ") tri(" << w->tri0 << "," << w->tri1 << ")");
        }
    }
    vdc_debug::maybe_dump_selfi_stage(dt, v, cycle_isovertices, cell_by_index, "A");

    if (trace || vdc_debug::log_unresolved_A_vertices()) {
        DEBUG_PRINT("[DEL-SELFI] Vertex " << v->info().index
                    << " unresolved by A-separation (cycles=" << num_cycles << ").");
    }

    // -multi_isov_trace: dump the first time this vertex fails locally
    // (per-run, not per-pass) into <trace_dir>/local/.
    if (options.multi_isov_trace &&
        !options.multi_isov_trace_dir.empty() &&
        separation_result.status == ResolutionStatus::UNRESOLVED &&
        local_unresolved_dumped_vertices) {
        const int vidx = v->info().index;
        if (local_unresolved_dumped_vertices->insert(vidx).second) {
            const std::filesystem::path out_dir =
                std::filesystem::path(options.multi_isov_trace_dir) / "local";

            // Also dump the deterministic A-mode candidate set for the active strategy.
            vdc_debug::dump_multi_isov_trace_case(
                out_dir, v, baseline_positions, dt, cell_by_index, options.use_sep_dir);

            // Match the historical "simple_multi_failures" trace format: baseline + attempt.
            // (baseline is already emitted by vdc_debug::dump_multi_isov_trace_case.)
            if (!geometric_attempt_positions.empty()) {
                vdc_debug::dump_simple_multi_failure_stage(
                    out_dir, v, geometric_attempt_positions, dt, cell_by_index, "geometric_attempt");
            }
        }
    }

    return separation_result;
}

// ============================================================================
// Helper: Update statistics based on resolution result
// ============================================================================

/**
 * @brief Updates statistics counters based on a self-intersection resolution result.
 * 
 * This helper consolidates the repetitive switch logic for updating stats
 * after each resolution attempt.
 * 
 * @param result The resolution result from resolution attempts.
 * @param stats [in/out] Statistics to update.
 * @param any_changed [out] Set to true if resolution changed isovertices.
 * @param pass_detected [out] Incremented if self-intersection was detected.
 * @param pass_resolved [out] Incremented if successfully resolved.
 * @param pass_unresolved [out] Incremented if resolution failed.
 */
static void update_resolution_stats(
    const ResolutionResult& result,
    IsovertexComputationStats& stats,
    bool& any_changed,
    int& pass_detected,
    int& pass_resolved,
    int& pass_unresolved
) {
    const ResolutionStatus status = result.status;
    
    if (status != ResolutionStatus::NOT_NEEDED) {
        pass_detected++;
    }
    
    switch (status) {
        case ResolutionStatus::NOT_NEEDED:
            break;
        case ResolutionStatus::RESOLVED:
            any_changed = true;
            pass_resolved++;
            if (result.strategy == ResolutionStrategy::GEOMETRIC_SEPARATION) {
                stats.strat_geometric_separation++;
            }
            break;
        case ResolutionStatus::UNRESOLVED:
            stats.had_unresolved = true;
            pass_unresolved++;
            stats.self_intersection_unresolved++;
            break;
    }
}

// ============================================================================
// Sub-routine: Initialize isovertices per cycle
// ============================================================================

static std::vector<Vertex_handle> initialize_cycle_isovertices(
    Delaunay& dt,
    const UnifiedGrid& grid,
    float isovalue,
    const std::vector<Cell_handle>& cell_by_index,
    IsovertexComputationStats& stats,
    bool position_multi_isov_on_delv,
    bool move_cap,
    bool move_cap_strict,
    int sep_split,
    int supersample_r,
    const StrictTetPreservationContext& strict_context
) {
    std::vector<Vertex_handle> multi_cycle_vertices;
    
    // Iterate over all finite vertices in the Delaunay triangulation.
    // For each active vertex with assigned cycles, compute initial isovertex positions.
    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
        // Skip inactive vertices and boundary dummy vertices.
        if (!vit->info().active) continue;
        if (vit->info().is_dummy) continue;

        auto& cycles = vit->info().facet_cycles;
        auto& isovertices = vit->info().cycle_isovertices;

        // Skip vertices with no cycles (shouldn't happen for active vertices).
        if (cycles.empty()) continue;

        isovertices.resize(cycles.size());

        Point cube_center = vit->point();

        if (cycles.size() == 1) {
            // ----------------------------------------------------------------
            // Single cycle vertex: position isovertex at the isosurface sample
            // ----------------------------------------------------------------
            // Prefer the precomputed isosurface sample point when available.
            // Refinement-inserted vertices may not have a cube-derived sample.
            if (vit->info().has_isov_sample) {
                // The isosurface sample point is stored in vit->info().isov.
                // When -position_delv_on_isov is enabled, this equals cube_center.
                isovertices[0] = vit->info().isov;
            } else {
                // Fallback: compute centroid of Voronoi-edge/isovalue intersections for this cycle.
                isovertices[0] = compute_centroid_of_voronoi_edge_and_isosurface(
                    dt, grid, isovalue, vit, cycles[0], cell_by_index);
            }
            stats.single_cycle_count++;
        } else {
            // ----------------------------------------------------------------
            // Multi-cycle vertex: min-sphere separation direction + sphere projection
            // ----------------------------------------------------------------
            // For vertices where the isosurface passes through multiple times
            // (creating multiple cycles), we compute a separate isovertex per cycle.
            // Each isovertex is positioned by:
            //   1. Computing cycle separating direction using min-sphere of facet normals
            //   2. Projecting from the Delaunay site along that direction onto a sphere
            // This matches the directional model used by -use_sep_dir.
            if (position_multi_isov_on_delv) {
                for (size_t c = 0; c < cycles.size(); ++c) {
                    isovertices[c] = cube_center;
                }
            } else {
                const double sphere_radius = compute_sphere_radius(vit, dt, grid, sep_split, supersample_r);
                const bool trace_cap = vdc_debug::trace_selfi_enabled(vit->info().index);
                if (trace_cap) {
                    // For traced multi-cycle vertices, always dump cycle metadata even when
                    // this vertex has no local Stage-2 self-intersection.
                    vdc_debug::maybe_dump_selfi_cycle_metadata(dt, vit, cell_by_index);
                }
                struct InitialCycleDirPlan {
                    bool has_foldover_edge = false;
                    CycleFoldoverEdgeInfo foldover;
                    int foldover_candidates_ge4 = 0;
                    bool skip_sep_dir_for_foldover = false;
                    bool sep_dir_attempted = false;
                    bool has_sep_dir = false;
                    float sep_dir[3] = {0.0f, 0.0f, 0.0f};
                    double maxdist = std::numeric_limits<double>::infinity();
                    int acute_dihedral_constraints = 0;
                    int strict_preserve_hits = 0;
                    int strict_preserve_rejects = 0;
                    int direction_keep = 0;
                    int direction_skip = 0;
                    bool directional_cap_applied = false;
                    double cap_radius = 0.0;
                    Point projected = Point(0.0, 0.0, 0.0);
                };

                std::vector<InitialCycleDirPlan> plans(cycles.size());
                auto project_from_center_vec = [&cube_center](const Vector3& dir, double r) -> Point {
                    return Point(
                        cube_center.x() + r * dir.x(),
                        cube_center.y() + r * dir.y(),
                        cube_center.z() + r * dir.z());
                };

                for (size_t c = 0; c < cycles.size(); ++c) {
                    InitialCycleDirPlan& plan = plans[c];
                    plan.projected = cube_center;

                    int fold_candidates = 0;
                    plan.foldover = detect_cycle_crossover_edge_direction(
                        vit, cycles[c], cell_by_index, dt, &fold_candidates);
                    plan.foldover_candidates_ge4 = fold_candidates;
                    plan.has_foldover_edge = plan.foldover.found;

                    // Requested logic (2-cycle case): if this cycle has a 4-facet shared edge,
                    // skip sep-dir computation for that cycle and handle placement after both
                    // cycles are classified.
                    plan.skip_sep_dir_for_foldover =
                        (cycles.size() == 2 && plan.has_foldover_edge);
                    if (!plan.skip_sep_dir_for_foldover) {
                        plan.sep_dir_attempted = true;
                        plan.has_sep_dir = compute_cycle_separating_direction(
                            vit, cycles[c], cell_by_index, dt, static_cast<int>(c), plan.sep_dir);
                    }
                }

                // Only run simplify pass when more than one cycle fails sep-dir, and only
                // on those failed cycles.
                std::vector<size_t> sep_dir_failed_cycles;
                sep_dir_failed_cycles.reserve(cycles.size());
                for (size_t c = 0; c < cycles.size(); ++c) {
                    const InitialCycleDirPlan& plan = plans[c];
                    if (plan.sep_dir_attempted && !plan.has_sep_dir) {
                        sep_dir_failed_cycles.push_back(c);
                    }
                }
                if (sep_dir_failed_cycles.size() > 1) {
                    for (size_t c : sep_dir_failed_cycles) {
                        InitialCycleDirPlan& plan = plans[c];
                        const auto simplified = simplify_cycle_facets(
                            vit, cycles[c], cell_by_index, dt, static_cast<int>(c));
                        if (simplified.size() < cycles[c].size()) {
                            plan.has_sep_dir = compute_cycle_separating_direction(
                                vit, simplified, cell_by_index, dt, static_cast<int>(c), plan.sep_dir, "_final");
                        }
                    }
                }

                for (size_t c = 0; c < cycles.size(); ++c) {
                    InitialCycleDirPlan& plan = plans[c];
                    if (move_cap) {
                        const bool run_relax_cap = !move_cap_strict;
                        const bool run_strict_directional_cap = move_cap_strict && plan.has_sep_dir;
                        if (run_relax_cap || run_strict_directional_cap) {
                            const Vector3 move_dir(
                                static_cast<double>(plan.sep_dir[0]),
                                static_cast<double>(plan.sep_dir[1]),
                                static_cast<double>(plan.sep_dir[2]));
                            const Vector3* move_dir_ptr =
                                (run_strict_directional_cap ? &move_dir : nullptr);
                            const CycleMoveCapBound cap = compute_cycle_maxdist_to_opposite_face_planes(
                                dt,
                                cell_by_index,
                                vit,
                                static_cast<int>(c),
                                strict_context,
                                move_dir_ptr);
                            plan.maxdist = cap.maxdist;
                            plan.acute_dihedral_constraints = cap.acute_dihedral_constraints;
                            plan.strict_preserve_hits = cap.strict_preserve_hits;
                            plan.strict_preserve_rejects = cap.strict_preserve_rejects;
                            plan.direction_keep = cap.direction_keep;
                            plan.direction_skip = cap.direction_skip;
                            plan.directional_cap_applied = run_strict_directional_cap;
                        }
                    }
                    plan.cap_radius =
                        (move_cap && std::isfinite(plan.maxdist)
                            ? std::min(sphere_radius, 0.5 * plan.maxdist)
                            : sphere_radius);

                    if (trace_cap) {
                        DEBUG_PRINT("[DEL-MOVE-CAP] v=" << vit->info().index
                                    << " cycle=" << c
                                    << " r1=" << sphere_radius
                                    << " maxdist=" << plan.maxdist
                                    << " cap_r=" << plan.cap_radius
                                    << " acute_dihedral=" << plan.acute_dihedral_constraints
                                    << " strict_keep=" << plan.strict_preserve_hits
                                    << " strict_skip=" << plan.strict_preserve_rejects
                                    << " direction_keep=" << plan.direction_keep
                                    << " direction_skip=" << plan.direction_skip
                                    << " foldover_candidates_ge4=" << plan.foldover_candidates_ge4
                                    << " foldover_edge=" << (plan.has_foldover_edge ? "yes" : "no")
                                    << " has_sep_dir=" << (plan.has_sep_dir ? "yes" : "no")
                                    << " directional_cap=" << (plan.directional_cap_applied ? "yes" : "no")
                                    << (plan.skip_sep_dir_for_foldover ? " (skipped_by_foldover)" : "")
                                    << (move_cap ? (move_cap_strict ? " (strict)" : " (relax)") : " (disabled)"));
                    }
                }

                if (cycles.size() == 2) {
                    const bool f0 = plans[0].has_foldover_edge;
                    const bool f1 = plans[1].has_foldover_edge;
                    const int fold_count = static_cast<int>(f0) + static_cast<int>(f1);

                    if (fold_count == 2) {
                        // Both cycles have foldover edges: move each cycle along its own foldover edge.
                        plans[0].projected = project_from_center_vec(
                            plans[0].foldover.unit_dir, plans[0].cap_radius);
                        plans[1].projected = project_from_center_vec(
                            plans[1].foldover.unit_dir, plans[1].cap_radius);
                    } else if (fold_count == 1) {
                        // Exactly one foldover cycle: keep it on DelV center, move the other by sep-dir.
                        const int fold_idx = f0 ? 0 : 1;
                        const int other_idx = 1 - fold_idx;
                        plans[fold_idx].projected = cube_center;
                        if (plans[other_idx].has_sep_dir) {
                            plans[other_idx].projected = Point(
                                cube_center.x() + plans[other_idx].cap_radius *
                                                    static_cast<double>(plans[other_idx].sep_dir[0]),
                                cube_center.y() + plans[other_idx].cap_radius *
                                                    static_cast<double>(plans[other_idx].sep_dir[1]),
                                cube_center.z() + plans[other_idx].cap_radius *
                                                    static_cast<double>(plans[other_idx].sep_dir[2]));
                        }
                    } else {
                        // No foldover edges detected: standard sep-dir placement.
                        for (int c = 0; c < 2; ++c) {
                            if (!plans[static_cast<size_t>(c)].has_sep_dir) {
                                continue;
                            }
                            plans[static_cast<size_t>(c)].projected = Point(
                                cube_center.x() + plans[static_cast<size_t>(c)].cap_radius *
                                                    static_cast<double>(plans[static_cast<size_t>(c)].sep_dir[0]),
                                cube_center.y() + plans[static_cast<size_t>(c)].cap_radius *
                                                    static_cast<double>(plans[static_cast<size_t>(c)].sep_dir[1]),
                                cube_center.z() + plans[static_cast<size_t>(c)].cap_radius *
                                                    static_cast<double>(plans[static_cast<size_t>(c)].sep_dir[2]));
                        }
                    }

                    if (trace_cap) {
                        DEBUG_PRINT("[DEL-MOVE-CAP] v=" << vit->info().index
                                    << " two-cycle placement: fold_count=" << fold_count
                                    << " f0=" << (f0 ? "1" : "0")
                                    << " f1=" << (f1 ? "1" : "0"));
                    }
                } else {
                    // Existing behavior for >=3 cycles.
                    for (size_t c = 0; c < cycles.size(); ++c) {
                        if (!plans[c].has_sep_dir) {
                            continue;
                        }
                        plans[c].projected = Point(
                            cube_center.x() + plans[c].cap_radius * static_cast<double>(plans[c].sep_dir[0]),
                            cube_center.y() + plans[c].cap_radius * static_cast<double>(plans[c].sep_dir[1]),
                            cube_center.z() + plans[c].cap_radius * static_cast<double>(plans[c].sep_dir[2]));
                    }
                }

                for (size_t c = 0; c < cycles.size(); ++c) {
                    if (!plans[c].has_sep_dir &&
                        squared_distance(plans[c].projected, cube_center) <= 1e-24 &&
                        trace_cap) {
                        DEBUG_PRINT("[DEL-MOVE-CAP] v=" << vit->info().index
                                    << " cycle=" << c
                                    << " no separating direction; fallback to DelV center");
                    }
                    (void)apply_dihedral_orientation_gate_to_cycle_position(
                        dt, cell_by_index, vit, static_cast<int>(c), plans[c].projected, strict_context);
                    isovertices[c] = plans[c].projected;
                }
            }
            multi_cycle_vertices.push_back(vit);
            stats.multi_cycle_count++;
        }
    }
    
    return multi_cycle_vertices;
}

// ============================================================================
// Sub-routine: Resolve self-intersections on multi-cycle vertices
// ============================================================================

/**
 * @brief Collect boundary vertices around multi-cycle cycles.
 *
 * Motivation:
 * - Stage 2 primarily targets multi-cycle vertices. However, we have observed global mesh
 *   self-intersections that originate from *single-cycle* triangle fans adjacent to a
 *   multi-cycle vertex (a small multi-cycle cycle can "poke into" a neighboring fan
 *   without triggering a local multi-cycle self-intersection at the multi-cycle vertex).
 * - Scanning all single-cycle vertices for foldovers is expensive. Instead, we collect a
 *   targeted set of boundary vertices adjacent to multi-cycle cycles and only
 *   attempt foldover correction there when `options.foldover` is enabled.
 *
 * This keeps overhead low while covering the aneurysm foldover pattern:
 * the problematic single-cycle fan is consistently adjacent to a small multi-cycle cycle.
 */
static std::vector<Vertex_handle> collect_foldover_boundary_vertices_around_multicycles(
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    const std::vector<Vertex_handle>& multi_cycle_vertices
) {
    std::vector<Vertex_handle> out;
    out.reserve(1024);

    std::unordered_set<int> seen;
    seen.reserve(2048);

    for (Vertex_handle mh : multi_cycle_vertices) {
        if (mh == Vertex_handle()) {
            continue;
        }
        if (!mh->info().active || mh->info().is_dummy) {
            continue;
        }
        const auto& mcycles = mh->info().facet_cycles;
        if (mcycles.size() < 2) {
            continue;
        }

        for (int mc = 0; mc < static_cast<int>(mcycles.size()); ++mc) {
            const auto& facets = mcycles[static_cast<size_t>(mc)];
            if (facets.empty()) {
                continue;
            }

            const CycleBoundaryInfo boundary = gather_cycle_boundary_info(
                dt, cell_by_index, mh, facets);
            for (Vertex_handle bh : boundary.boundary_verts) {
                if (bh == Vertex_handle()) {
                    continue;
                }
                if (!bh->info().active || bh->info().is_dummy) {
                    continue;
                }
                const int idx = bh->info().index;
                if (seen.insert(idx).second) {
                    out.push_back(bh);
                }
            }
        }
    }

    return out;
}

static std::vector<Vertex_handle> resolve_multicycle_self_intersections(
    Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    const std::vector<Vertex_handle>& multi_cycle_vertices,
    IsovertexComputationStats& stats,
    const CycleIsovertexOptions& options,
    double foldover_sphere_radius,
    const StrictTetPreservationContext& strict_context
) {
    std::vector<Vertex_handle> modified_multi_cycle_vertices;
    std::unordered_set<int> modified_multi_cycle_indices;
    std::unordered_set<int> local_unresolved_dumped_vertices;
    std::unordered_set<int>* local_unresolved_dumped_vertices_ptr = nullptr;

    // Targeted set of boundary vertices where we also allow foldover resolution
    // (these are often single-cycle culprits for global selfI).
    const std::vector<Vertex_handle> foldover_boundary_vertices =
        options.foldover
            ? collect_foldover_boundary_vertices_around_multicycles(
                  dt, cell_by_index, multi_cycle_vertices)
            : std::vector<Vertex_handle>();

    if (options.multi_isov_trace &&
        !options.multi_isov_trace_dir.empty()) {
        local_unresolved_dumped_vertices_ptr = &local_unresolved_dumped_vertices;
    }

    // Run multiple stabilization passes to handle cascading conflicts.
    // Changing one multi-cycle vertex may introduce a new conflict at another.
    const int MAX_MULTI_CYCLE_PASSES = 3;
    
        for (int pass = 0; pass < MAX_MULTI_CYCLE_PASSES; ++pass) {
            bool any_changed = false;
            int pass_detected = 0;
            int pass_resolved = 0;
            int pass_unresolved = 0;

	        // Process each multi-cycle vertex for self-intersection.
	        for (Vertex_handle vh : multi_cycle_vertices) {
	            if (!vh->info().active) continue;
	            if (vh->info().is_dummy) continue;
	            if (vh->info().facet_cycles.size() < 2) continue;

	            auto& isovertices = vh->info().cycle_isovertices;
	            const std::vector<Point> isovertices_before = isovertices;

            const ResolutionResult result = resolve_multicycle_self_intersection_at_vertex(
                vh,
                isovertices,
                dt,
                cell_by_index,
                options,
                foldover_sphere_radius,
                local_unresolved_dumped_vertices_ptr,
                strict_context);

	            update_resolution_stats(result, stats, any_changed,
	                pass_detected, pass_resolved, pass_unresolved);

	            // Even if we couldn't fully resolve, keep stabilization passes going when
	            // the vertex moved (partial improvement in A-only mode).
	            {
	                bool moved = false;
	                if (isovertices_before.size() != isovertices.size()) {
	                    moved = true;
	                } else {
	                    constexpr double eps2 = 1e-18;
	                    for (size_t i = 0; i < isovertices.size(); ++i) {
	                        if (squared_distance(isovertices_before[i], isovertices[i]) > eps2) {
	                            moved = true;
	                            break;
	                        }
	                    }
	                }
	                if (moved) {
	                    any_changed = true;
	                }
	            }

	            // Track which vertices were modified for targeted cleanup passes.
	            const ResolutionStatus status = result.status;
	            if (status == ResolutionStatus::RESOLVED ||
	                status == ResolutionStatus::UNRESOLVED) {
                if (modified_multi_cycle_indices.insert(vh->info().index).second) {
                    modified_multi_cycle_vertices.push_back(vh);
                }
            }
        }

        // --------------------------------------------------------------------
        // Foldover resolution on a targeted set of boundary vertices.
        // --------------------------------------------------------------------
        //
        // This extends Stage 2 beyond multi-cycle vertices *only when enabled*.
        // We focus on boundary vertices adjacent to multi-cycle cycles (any size).
        // This targets global self-intersections caused by single-cycle fan foldovers
        // near a multi-cycle site without scanning all single-cycle vertices.
        if (options.foldover && !foldover_boundary_vertices.empty()) {
            for (Vertex_handle vh : foldover_boundary_vertices) {
                if (vh == Vertex_handle()) {
                    continue;
                }
                if (!vh->info().active || vh->info().is_dummy) {
                    continue;
                }

                // Only process single-cycle vertices here. Multi-cycle vertices are already handled
                // by the main Stage 2 loop above (including within-cycle witnesses).
                if (vh->info().facet_cycles.size() != 1) {
                    continue;
                }
                if (vh->info().facet_cycles[0].size() < 4) {
                    continue;
                }

                auto& isovertices = vh->info().cycle_isovertices;
                if (isovertices.empty()) {
                    continue;
                }
                const std::vector<Point> isovertices_before = isovertices;

                const ResolutionResult result = try_resolve_cycle_fan_foldover_at_vertex_cycle(
                    vh, /*cycle_idx=*/0, isovertices, dt, cell_by_index, foldover_sphere_radius, strict_context);

                update_resolution_stats(result, stats, any_changed,
                    pass_detected, pass_resolved, pass_unresolved);

                // Keep stabilization passes going if we moved any isovertex, even if the local
                // self-intersection was not fully resolved (UNRESOLVED). This mirrors the multi-cycle loop.
                {
                    bool moved = false;
                    if (isovertices_before.size() != isovertices.size()) {
                        moved = true;
                    } else {
                        constexpr double eps2 = 1e-18;
                        for (size_t i = 0; i < isovertices.size(); ++i) {
                            if (squared_distance(isovertices_before[i], isovertices[i]) > eps2) {
                                moved = true;
                                break;
                            }
                        }
                    }
                    if (moved) {
                        any_changed = true;
                    }
                }
            }
        }

        // Accumulate pass statistics into overall stats.
        stats.self_intersection_detected += pass_detected;
        stats.self_intersection_resolved += pass_resolved;

        if (pass_detected > 0) {
            DEBUG_PRINT("[DEL-SELFI] Stage2 pass " << (pass + 1) << "/" << MAX_MULTI_CYCLE_PASSES
                        << ": detected=" << pass_detected
                        << ", resolved=" << pass_resolved
                        << ", unresolved=" << pass_unresolved);
        }

        // Early exit if no changes were made this pass.
        if (!any_changed) {
            break;
        }
    }
    vdc_debug::finalize_multi_isov_trace(
        options,
        local_unresolved_dumped_vertices,
        multi_cycle_vertices,
        dt,
        cell_by_index);
    
    return modified_multi_cycle_vertices;
}

// ============================================================================
// Sub-routine: Resolve within-cycle foldovers on suspect cycle fans
// ============================================================================

struct VertexCycleRef {
    Vertex_handle v;
    int cycle = -1;
};

static std::vector<VertexCycleRef> collect_suspect_cycle_fans_by_normal_inconsistency(
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index
) {
    std::vector<VertexCycleRef> suspects;
    suspects.reserve(1024);

    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
        if (!vit->info().active) continue;
        if (vit->info().is_dummy) continue;

        const auto& cycles = vit->info().facet_cycles;
        if (cycles.empty()) continue;

        for (int c = 0; c < static_cast<int>(cycles.size()); ++c) {
            const auto& cycle_facets = cycles[static_cast<size_t>(c)];
            if (cycle_facets.size() < 4) {
                continue;
            }

            std::vector<Vector3> normals;
            normals.reserve(cycle_facets.size());

            Vector3 sum(0, 0, 0);
            for (const auto& [cell_idx, facet_idx] : cycle_facets) {
                Cell_handle cell = lookup_cell(cell_by_index, cell_idx);
                if (cell == Cell_handle()) {
                    continue;
                }
                if (facet_idx < 0 || facet_idx >= 4) {
                    continue;
                }

                // Reconstruct the oriented isosurface triangle for this facet
                // (same parity convention as generate_isosurface_triangles()).
                std::array<Point, 3> p;
                bool ok = true;
                for (int k = 0; k < 3; ++k) {
                    const int cell_vertex_idx = (facet_idx + 1 + k) % 4;
                    Vertex_handle vh = cell->vertex(cell_vertex_idx);
                    if (vh == Vertex_handle()) {
                        ok = false;
                        break;
                    }
                    if (dt.is_infinite(vh) || vh->info().is_dummy || !vh->info().active) {
                        ok = false;
                        break;
                    }

                    int cycle_idx = -1;
                    const int anchor = (facet_idx + 1) % 4;
                    const int slot = (cell_vertex_idx - anchor + 4) % 4;
                    if (slot >= 0 && slot < 3) {
                        cycle_idx = cell->info().facet_info[facet_idx].dualCellEdgeIndex[slot];
                    }
                    if (cycle_idx < 0) {
                        cycle_idx = find_cycle_containing_facet(vh, cell_idx, facet_idx);
                    }
                    if (cycle_idx < 0 ||
                        cycle_idx >= static_cast<int>(vh->info().cycle_isovertices.size())) {
                        ok = false;
                        break;
                    }
                    p[static_cast<size_t>(k)] = vh->info().cycle_isovertices[static_cast<size_t>(cycle_idx)];
                }

                if (!ok) {
                    continue;
                }

                if (facet_idx % 2 != 0) {
                    std::swap(p[1], p[2]);
                }

                const Vector3 n = vec_cross(Vector3(p[0], p[1]), Vector3(p[0], p[2]));
                if (n.squared_length() < 1e-24) {
                    continue;
                }

                normals.push_back(n);
                sum = sum + n;
            }

            // Only meaningful when there are enough non-degenerate triangles.
            if (normals.size() < 4) {
                continue;
            }
            if (sum.squared_length() < 1e-24) {
                continue;
            }

            // Heuristic: look for strongly inconsistent normals (likely foldover).
            // Use a cosine threshold to avoid over-triggering on near-perpendicular noise.
            const double sum_len = std::sqrt(sum.squared_length());
            bool inconsistent = false;
            int neg_count = 0;
            double min_cos = 1.0;
            for (const Vector3& n : normals) {
                const double n_len2 = n.squared_length();
                if (n_len2 < 1e-24) {
                    continue;
                }
                const double cos = vec_dot(n, sum) / (std::sqrt(n_len2) * sum_len);
                min_cos = std::min(min_cos, cos);
                if (cos < 0.0) {
                    ++neg_count;
                    if (cos < -0.25) {
                        inconsistent = true;
                        break;
                    }
                }
            }
            if (!inconsistent && neg_count >= 2 && min_cos < -0.1) {
                inconsistent = true;
            }

            if (inconsistent) {
                suspects.push_back(VertexCycleRef{vit, c});
            }
        }
    }

    std::sort(
        suspects.begin(),
        suspects.end(),
        [](const VertexCycleRef& a, const VertexCycleRef& b) {
            const int ai = (a.v == Vertex_handle()) ? -1 : a.v->info().index;
            const int bi = (b.v == Vertex_handle()) ? -1 : b.v->info().index;
            if (ai != bi) return ai < bi;
            return a.cycle < b.cycle;
        });

    return suspects;
}

static bool cycle_fan_has_within_cycle_intersection(
    const Delaunay& dt,
    Vertex_handle v,
    int cycle_idx,
    const Point& cycle_isovertex,
    const std::vector<Cell_handle>& cell_by_index
) {
    const auto tris = collect_cycle_triangles(dt, v, cycle_idx, cycle_isovertex, cell_by_index);
    if (tris.size() < 4) {
        return false;
    }
    for (size_t i = 0; i < tris.size(); ++i) {
        for (size_t j = i + 1; j < tris.size(); ++j) {
            if (triangles_have_nontrivial_intersection(tris[i], tris[j])) {
                return true;
            }
        }
    }
    return false;
}

static ResolutionResult try_resolve_cycle_fan_foldover_at_vertex_cycle(
    Vertex_handle v,
    int cycle_idx,
    std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    double foldover_sphere_radius,
    const StrictTetPreservationContext& strict_context
) {
    const auto& cycles = v->info().facet_cycles;
    if (cycle_idx < 0 || cycle_idx >= static_cast<int>(cycles.size())) {
        return {ResolutionStatus::NOT_NEEDED, ResolutionStrategy::NONE};
    }
    if (cycle_idx >= static_cast<int>(cycle_isovertices.size())) {
        return {ResolutionStatus::NOT_NEEDED, ResolutionStrategy::NONE};
    }
    if (cycles[static_cast<size_t>(cycle_idx)].size() < 4) {
        return {ResolutionStatus::NOT_NEEDED, ResolutionStrategy::NONE};
    }

    const Point baseline_p = cycle_isovertices[static_cast<size_t>(cycle_idx)];

    // Collect the local fan once (used for both detection and candidate construction).
    const auto tris = collect_cycle_triangles(dt, v, cycle_idx, baseline_p, cell_by_index);
    if (tris.size() < 4) {
        return {ResolutionStatus::NOT_NEEDED, ResolutionStrategy::NONE};
    }

    bool has_within = false;
    for (size_t i = 0; i < tris.size() && !has_within; ++i) {
        for (size_t j = i + 1; j < tris.size(); ++j) {
            if (triangles_have_nontrivial_intersection(tris[i], tris[j])) {
                has_within = true;
                break;
            }
        }
    }

    if (!has_within) {
        return {ResolutionStatus::NOT_NEEDED, ResolutionStrategy::NONE};
    }

    // ------------------------------------------------------------------------
    // Deterministic foldover resolution (no candidate-set search)
    // ------------------------------------------------------------------------
    //
    // The failures we want to target here are "within-cycle triangle fan foldovers":
    // triangles incident to this cycle isovertex overlap in a non-trivial way
    // (often detected as a shared-vertex + opposite-edge penetration).
    //
    // The previous implementation tried a small menu of fixed candidate positions.
    // While deterministic, it was still a *search* over hand-picked candidates.
    //
    // This implementation is closer in spirit to `-use_sep_dir`:
    //  - compute a direction directly from local geometry,
    //  - take a single step along that direction (optionally trying the opposite sign),
    //  - and use `check_self_intersection()` as the acceptance test.
    //
    // This is intentionally conservative:
    //  - we only apply it when there is already a detected within-cycle intersection,
    //  - we only move one cycle isovertex at a time,
    //  - and we revert to baseline if the move does not eliminate the local self-intersection.
    const bool trace = vdc_debug::trace_selfi_enabled(v->info().index);
    const Point center = v->point();
    const std::vector<Point> baseline_positions = cycle_isovertices;

    // Direction derived from the fan.
    //
    // We keep this intentionally minimal:
    //   - primary: normalized sum of fan triangle normals,
    //   - fallback (only if the sum cancels / is degenerate): direction to the fan boundary centroid.
    //
    Vector3 sum_n(0, 0, 0);
    double sx = 0.0, sy = 0.0, sz = 0.0;
    int boundary_cnt = 0;
    double min_boundary_r2 = std::numeric_limits<double>::max();
    for (const auto& tri : tris) {
        const Vector3 n = vec_cross(
            Vector3(tri.vertex(0), tri.vertex(1)),
            Vector3(tri.vertex(0), tri.vertex(2)));
        if (n.squared_length() >= 1e-24) {
            sum_n = sum_n + n;
        }

        // In collect_cycle_triangles(), vertex(0) is the cycle isovertex, and vertices 1/2
        // are the boundary vertices (isovertices of neighboring Delaunay vertices).
        for (int k = 1; k < 3; ++k) {
            const Point q = tri.vertex(k);
            sx += q.x();
            sy += q.y();
            sz += q.z();
            ++boundary_cnt;
            min_boundary_r2 = std::min(min_boundary_r2, squared_distance(center, q));
        }
    }

    const double sum_n_len = vec_norm(sum_n);
    Vector3 dir(0, 0, 0);
    const char* dir_label = "fan_normal";
    if (sum_n_len > 1e-12) {
        dir = sum_n / sum_n_len;
    } else if (boundary_cnt > 0) {
        const Point centroid(sx / boundary_cnt, sy / boundary_cnt, sz / boundary_cnt);
        dir = vec_normalize_or(Vector3(center, centroid), Vector3(1, 0, 0));
        dir_label = "boundary_centroid_fallback";
    } else {
        // Extremely degenerate fan: fall back to a fixed axis to avoid "no direction".
        dir = Vector3(1, 0, 0);
        dir_label = "fixed_axis_fallback";
    }

    // Determine the local scale for step size.
    //
    // We intentionally keep the step schedule tiny (2 radii, 2 signs):
    // - Most aneurysm foldovers resolve with a small fraction of the boundary scale.
    // - A single global radius (foldover_sphere_radius) is too large in some cases and can fail
    //   the within-fan acceptance test.
    //
    // Local boundary scale: r_min = min ||q-center|| over all fan boundary vertices q.
    // If it is unavailable (degenerate), fall back to foldover_sphere_radius.
    double r_min = 0.0;
    if (min_boundary_r2 < std::numeric_limits<double>::max() && min_boundary_r2 > 1e-24) {
        r_min = std::sqrt(min_boundary_r2);
    } else if (foldover_sphere_radius > 1e-12) {
        r_min = foldover_sphere_radius;
    }
    if (r_min <= 1e-12) {
        cycle_isovertices = baseline_positions;
        return {ResolutionStatus::UNRESOLVED, ResolutionStrategy::GEOMETRIC_SEPARATION};
    }

    // Two deterministic radii (smallest-first).
    // This is the minimal schedule that still covers our aneurysm regression cases.
    const double r_small = 0.20 * r_min;
    const double r_large = 0.50 * r_min;

    // Choose a deterministic sign preference that is stable under tiny perturbations:
    // try to keep directions consistent with the current isovertex radial direction (if any).
    const Vector3 baseline_dir = Vector3(center, baseline_p);
    if (baseline_dir.squared_length() > 1e-24 && vec_dot(dir, baseline_dir) < 0.0) {
        dir = -dir;
    }

    auto try_candidate = [&](double sign, double r) -> std::optional<ResolutionResult> {
        Point cand(
            center.x() + sign * r * dir.x(),
            center.y() + sign * r * dir.y(),
            center.z() + sign * r * dir.z());
        (void)apply_dihedral_orientation_gate_to_cycle_position(
            dt, cell_by_index, v, cycle_idx, cand, strict_context);
        const bool within_ok = !cycle_fan_has_within_cycle_intersection(
            dt, v, cycle_idx, cand, cell_by_index);
        if (trace) {
            DEBUG_PRINT("[DEL-SELFI-FAN] v=" << v->info().index
                        << " cycle=" << cycle_idx
                        << " " << dir_label
                        << " sign=" << (sign > 0.0 ? "+" : "-")
                        << " r=" << r
                        << ": within_ok=" << (within_ok ? "YES" : "NO")
                        << " dir=(" << dir.x() << "," << dir.y() << "," << dir.z() << ")"
                        << " baseline=(" << baseline_p.x() << "," << baseline_p.y() << "," << baseline_p.z() << ")"
                        << " cand=(" << cand.x() << "," << cand.y() << "," << cand.z() << ")");
        }
        if (!within_ok) {
            return std::nullopt;
        }

        cycle_isovertices = baseline_positions;
        cycle_isovertices[static_cast<size_t>(cycle_idx)] = cand;
        const bool all_ok = !check_self_intersection(v, cycle_isovertices, dt, cell_by_index);
        return ResolutionResult{
            all_ok ? ResolutionStatus::RESOLVED : ResolutionStatus::UNRESOLVED,
            ResolutionStrategy::GEOMETRIC_SEPARATION};
    };

    // Minimal deterministic candidate list (4 tries):
    //   (1) +r_small, (2) -r_small, (3) +r_large, (4) -r_large
    for (double r : {r_small, r_large}) {
        if (r <= 1e-12) {
            continue;
        }
        if (auto res = try_candidate(+1.0, r)) {
            return *res;
        }
        if (auto res = try_candidate(-1.0, r)) {
            return *res;
        }
    }

    // No resolution: revert to baseline.
    cycle_isovertices = baseline_positions;
    return {ResolutionStatus::UNRESOLVED, ResolutionStrategy::GEOMETRIC_SEPARATION};
}

static void resolve_suspect_cycle_fan_foldovers(
    Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    const std::vector<Vertex_handle>& multi_cycle_vertices,
    const std::vector<Vertex_handle>& modified_multi_cycle_vertices,
    IsovertexComputationStats& stats,
    const CycleIsovertexOptions& options,
    double foldover_sphere_radius,
    const StrictTetPreservationContext& strict_context
) {
    if (!options.foldover) {
        return;
    }

    const int MAX_FAN_PASSES = 2;
    for (int pass = 0; pass < MAX_FAN_PASSES; ++pass) {
        auto suspects = collect_suspect_cycle_fans_by_normal_inconsistency(dt, cell_by_index);

        // Ensure we don't miss within-cycle foldovers on multi-cycle vertices that were involved
        // in multi-cycle repositioning. This list is typically small and avoids relying solely on
        // the normal-inconsistency heuristic.
        if (!modified_multi_cycle_vertices.empty()) {
            for (Vertex_handle vh : modified_multi_cycle_vertices) {
                if (vh == Vertex_handle()) {
                    continue;
                }
                if (!vh->info().active || vh->info().is_dummy) {
                    continue;
                }
                const auto& cycles = vh->info().facet_cycles;
                if (cycles.size() < 2) {
                    continue;
                }

                bool has_large_cycle = false;
                for (const auto& cycle_facets : cycles) {
                    if (cycle_facets.size() >= 4) {
                        has_large_cycle = true;
                        break;
                    }
                }
                if (!has_large_cycle) {
                    continue;
                }

                const auto& isovertices = vh->info().cycle_isovertices;
                if (!check_self_intersection(vh, isovertices, dt, cell_by_index)) {
                    continue;
                }

                for (int c = 0; c < static_cast<int>(cycles.size()); ++c) {
                    if (cycles[static_cast<size_t>(c)].size() < 4) {
                        continue;
                    }
                    suspects.push_back(VertexCycleRef{vh, c});
                }
            }
        }

        // Also include boundary-vertex stars around small multi-cycle cycles (<=4 facets).
        // This targets rare global self-intersections where a multi-cycle isovertex pierces a
        // neighboring single-cycle triangle fan without causing local multi-cycle intersections.
        for (Vertex_handle mh : multi_cycle_vertices) {
            if (mh == Vertex_handle()) {
                continue;
            }
            if (!mh->info().active || mh->info().is_dummy) {
                continue;
            }
            const auto& mcycles = mh->info().facet_cycles;
            if (mcycles.size() < 2) {
                continue;
            }

            for (int mc = 0; mc < static_cast<int>(mcycles.size()); ++mc) {
                const auto& facets = mcycles[static_cast<size_t>(mc)];
                if (facets.empty() || facets.size() > 4) {
                    continue;
                }

                const CycleBoundaryInfo boundary = gather_cycle_boundary_info(
                    dt, cell_by_index, mh, facets);
                for (Vertex_handle bh : boundary.boundary_verts) {
                    if (bh == Vertex_handle()) {
                        continue;
                    }
                    if (!bh->info().active || bh->info().is_dummy) {
                        continue;
                    }

                    const auto& bcycles = bh->info().facet_cycles;
                    for (int bc = 0; bc < static_cast<int>(bcycles.size()); ++bc) {
                        if (bcycles[static_cast<size_t>(bc)].size() < 4) {
                            continue;
                        }
                        suspects.push_back(VertexCycleRef{bh, bc});
                    }
                }
            }
        }

        if (!suspects.empty()) {
            std::sort(
                suspects.begin(),
                suspects.end(),
                [](const VertexCycleRef& a, const VertexCycleRef& b) {
                    const int ai = (a.v == Vertex_handle()) ? -1 : a.v->info().index;
                    const int bi = (b.v == Vertex_handle()) ? -1 : b.v->info().index;
                    if (ai != bi) return ai < bi;
                    return a.cycle < b.cycle;
                });
            suspects.erase(
                std::unique(
                    suspects.begin(),
                    suspects.end(),
                    [](const VertexCycleRef& a, const VertexCycleRef& b) {
                        const int ai = (a.v == Vertex_handle()) ? -1 : a.v->info().index;
                        const int bi = (b.v == Vertex_handle()) ? -1 : b.v->info().index;
                        return ai == bi && a.cycle == b.cycle;
                    }),
                suspects.end());
        }
        if (suspects.empty()) {
            break;
        }

        bool any_changed = false;
        int pass_detected = 0;
        int pass_resolved = 0;
        int pass_unresolved = 0;

        for (const auto& ref : suspects) {
            if (ref.v == Vertex_handle()) {
                continue;
            }
            if (!ref.v->info().active || ref.v->info().is_dummy) {
                continue;
            }

            auto& isovertices = ref.v->info().cycle_isovertices;
            const std::vector<Point> before = isovertices;

            const ResolutionResult result = try_resolve_cycle_fan_foldover_at_vertex_cycle(
                ref.v, ref.cycle, isovertices, dt, cell_by_index, foldover_sphere_radius, strict_context);

            if (result.status != ResolutionStatus::NOT_NEEDED) {
                pass_detected++;
            }
            if (result.status == ResolutionStatus::RESOLVED) {
                pass_resolved++;
                stats.strat_geometric_separation++;
            } else if (result.status == ResolutionStatus::UNRESOLVED) {
                pass_unresolved++;
                stats.self_intersection_unresolved++;
                stats.had_unresolved = true;
            }

            // Keep passes going on any movement (even partial improvements).
            if (before.size() != isovertices.size()) {
                any_changed = true;
            } else {
                constexpr double eps2 = 1e-18;
                for (size_t i = 0; i < isovertices.size(); ++i) {
                    if (squared_distance(before[i], isovertices[i]) > eps2) {
                        any_changed = true;
                        break;
                    }
                }
            }
        }

        stats.self_intersection_detected += pass_detected;
        stats.self_intersection_resolved += pass_resolved;

        if (pass_detected > 0) {
            DEBUG_PRINT("[DEL-SELFI] fan pass " << (pass + 1) << "/" << MAX_FAN_PASSES
                        << ": suspects=" << suspects.size()
                        << ", detected=" << pass_detected
                        << ", resolved=" << pass_resolved
                        << ", unresolved=" << pass_unresolved);
        }

        if (!any_changed) {
            break;
        }
    }
}


// ============================================================================
// Sub-routine: Report isovertex computation statistics
// ============================================================================

static void report_isovertex_statistics(const IsovertexComputationStats& stats) {
    DEBUG_PRINT("[DEL-ISOV] Computed isovertices for "
                << stats.single_cycle_count << " single-cycle and "
                << stats.multi_cycle_count << " multi-cycle vertices");

    if (stats.self_intersection_detected > 0) {
        DEBUG_PRINT("[DEL-SELFI] Self-intersection stats: "
                    << "detected=" << stats.self_intersection_detected
                    << ", resolved=" << stats.self_intersection_resolved
                    << ", unresolved=" << stats.self_intersection_unresolved);

        const int64_t resolved_total = stats.strat_geometric_separation;
        if (resolved_total > 0) {
            auto pct = [&](int64_t count) -> int {
                return static_cast<int>((count * 100 + resolved_total / 2) / resolved_total);
            };

            DEBUG_PRINT("[DEL-SELFI] Resolution strategy breakdown:");
            DEBUG_PRINT("[DEL-SELFI]   geometric_separation: " << stats.strat_geometric_separation
                        << " (" << pct(stats.strat_geometric_separation) << "%)");
        }
    }
}

// ============================================================================
// Main orchestrator: compute_cycle_isovertices
// ============================================================================

/**
 * @brief Compute isovertex positions for all cycles around active Delaunay vertices.
 *
 * This is the main entry point for Step 9 of the Delaunay-based VDC algorithm.
 * The computation proceeds in two stages:
 *
 *   Stage 1 (initialize_cycle_isovertices):
 *     Compute initial positions for all active vertices.
 *     - Single-cycle vertices: use the isosurface sample point
 *     - Multi-cycle vertices: compute min-sphere separation directions and project to sphere
 *
 *   Stage 2 (resolve_multicycle_self_intersections):
 *     Iteratively resolve local self-intersections on multi-cycle vertices, and (when enabled)
 *     also resolve a targeted set of within-cycle foldovers on boundary single-cycle vertices.
 *     - Between-cycle conflicts: deterministic hyperplane separation tests.
 *     - Within-cycle fan foldovers (`-foldover`): deterministic direction-based fan correction.
 *
 * @param dt The Delaunay triangulation with cycles computed.
 * @param grid The scalar field grid.
 * @param isovalue The isovalue threshold.
 * @param position_on_isov If true, Delaunay sites are positioned on isosurface samples.
 */
void compute_cycle_isovertices(
    Delaunay& dt,
    const UnifiedGrid& grid,
    float isovalue,
    bool position_on_isov,
    const CycleIsovertexOptions& options
) {
    DEBUG_PRINT("[DEL-ISOV] Computing isovertex positions for cycles...");
    
    // Suppress unused parameter warning - the choice is encoded in vertex info.
    (void)position_on_isov;

    // Build a direct lookup table for efficient cell access by index.
    const std::vector<Cell_handle> cell_by_index = build_cell_index_vector(dt);

    IsosurfaceEdgeSet strict_isosurface_edges;
    StrictTetPreservationContext strict_context;
    strict_context.strict_mode = options.move_cap_strict;
    if (strict_context.strict_mode) {
        strict_isosurface_edges = build_isosurface_edge_membership_set(dt);
        strict_context.isosurface_edges = &strict_isosurface_edges;
        DEBUG_PRINT("[DEL-MOVE-CAP] strict-preserve-edge-set size="
                    << strict_isosurface_edges.size());
    }

    // Initialize statistics tracking.
    IsovertexComputationStats stats;

    // A single global "projection radius" is used throughout:
    // - Stage 1 uses it to project multi-cycle separation directions onto an inscribed sphere.
    // - Stage 2 uses it (when -foldover is enabled) as a deterministic step size
    //   for within-cycle fan foldover correction.
    //
    // This depends only on grid spacing + sep_split + supersample
    const double sphere_radius = 0.5 *
        (std::min({grid.spacing[0], grid.spacing[1], grid.spacing[2]}) /
         (static_cast<double>(std::max(1, options.sep_split + 1)) *
          static_cast<double>(std::max(1, options.supersample_r))));

    // ========================================================================
    // Stage 1: Compute initial isovertex positions for all active vertices.
    // ========================================================================
    std::vector<Vertex_handle> multi_cycle_vertices = initialize_cycle_isovertices(
        dt, grid, isovalue, cell_by_index, stats, options.position_multi_isov_on_delv,
        options.move_cap, options.move_cap_strict, options.sep_split, options.supersample_r,
        strict_context);

    if (options.position_multi_isov_on_delv) {
        DEBUG_PRINT("[DEL-ISOV] Skipping multi-cycle repositioning (-position_multi_isov_on_delv).");
        report_isovertex_statistics(stats);
        return;
    }

    // ========================================================================
    // Stage 2: Resolve self-intersections
    // ========================================================================
    (void)resolve_multicycle_self_intersections(
        dt, cell_by_index, multi_cycle_vertices, stats, options, sphere_radius, strict_context);

    // ========================================================================
    // Report final statistics.
    // ========================================================================
    report_isovertex_statistics(stats);
}
// ============================================================================
// Helper: Find cycle containing a facet
// ============================================================================

int find_cycle_containing_facet(
    Vertex_handle v_handle,
    int cell_index,
    int facet_index
) {
    const auto& cycles = v_handle->info().facet_cycles;

    // Straight forward case: single-cycle vertices only have one cycle, so every incident
    // isosurface facet must belong to cycle 0.
    if (cycles.size() == 1) {
        return 0;
    }

    for (size_t c = 0; c < cycles.size(); ++c) {
        for (const auto& [ci, fi] : cycles[c]) {
            if (ci == cell_index && fi == facet_index) {
                return static_cast<int>(c);
            }
        }
    }

    return -1; // Not found
}

// ============================================================================
// Modify Cycles Pass: Fix non-manifolds by flipping problematic matchings
// ============================================================================
// Note: EdgeMatchingState is now declared in vdc_del_cycles.h

//! @brief Flip the matching method for an edge
static void flip_edge_matching(EdgeMatchingState& state) {
    state.flip_count++;
    state.flipped = true;

    // After being flipped twice, switch to UNCONSTRAINED mode
    // This prevents oscillation between SEP_NEG and SEP_POS
    if (state.flip_count >= 2 && state.method != BIPOLAR_MATCH_METHOD::UNCONSTRAINED_MATCH) {
        state.method = BIPOLAR_MATCH_METHOD::UNCONSTRAINED_MATCH;
        state.unconstrained_offset = 0;
        return;
    }

    switch (state.method) {
        case BIPOLAR_MATCH_METHOD::SEP_POS:
            state.method = BIPOLAR_MATCH_METHOD::SEP_NEG;
            break;
        case BIPOLAR_MATCH_METHOD::SEP_NEG:
            state.method = BIPOLAR_MATCH_METHOD::SEP_POS;
            break;
        case BIPOLAR_MATCH_METHOD::UNCONSTRAINED_MATCH:
            ++state.unconstrained_offset;
            break;
        default:
            state.method = BIPOLAR_MATCH_METHOD::SEP_NEG;
            break;
    }
}

//! @brief Recompute cycles for a single vertex using current edge matching states
static void recompute_cycles_for_vertex(
    Delaunay& dt,
    Vertex_handle v_handle,
    const std::vector<Cell_handle>& cell_by_index,
    const std::unordered_map<EdgeKey, EdgeMatchingState, EdgeKeyHash>& edge_matching_states
) {
    if (!v_handle->info().active || v_handle->info().is_dummy) {
        return;
    }

    // Clear existing cycles
    v_handle->info().facet_cycles.clear();

    // Collect all isosurface facets incident on this vertex (same as compute_facet_cycles)
    std::vector<std::pair<int, int>> isosurface_facets;
    std::vector<Cell_handle> incident_cells;
    dt.incident_cells(v_handle, std::back_inserter(incident_cells));

    for (Cell_handle ch : incident_cells) {
        if (dt.is_infinite(ch)) continue;

        for (int i = 0; i < 4; ++i) {
            if (!ch->info().facet_is_isosurface[i]) continue;

            const int v_idx = find_vertex_index_in_cell(ch, v_handle);
            if (v_idx < 0 || v_idx == i) continue;

            // Ignore facets involving dummy sites
            bool has_dummy = false;
            for (int j = 0; j < 4; ++j) {
                if (j == i) continue;
                if (ch->vertex(j)->info().is_dummy) {
                    has_dummy = true;
                    break;
                }
            }
            if (has_dummy) continue;

            isosurface_facets.push_back({ch->info().index, i});
        }
    }

    if (isosurface_facets.empty()) {
        return;
    }

    // Build facet -> id map
    const int n = static_cast<int>(isosurface_facets.size());
    std::unordered_map<FacetKey, int, FacetKeyHash> facetkey_to_id;
    facetkey_to_id.reserve(isosurface_facets.size());
    for (int i = 0; i < n; ++i) {
        facetkey_to_id.emplace(FacetKey{isosurface_facets[i].first, isosurface_facets[i].second}, i);
    }

    // Collect neighbor -> facets mapping
    const int v_global_index = v_handle->info().index;
    std::unordered_map<int, std::vector<int>> neighbor_to_facets;

    for (int fid = 0; fid < n; ++fid) {
        const int cell_index = isosurface_facets[fid].first;
        const int facet_index = isosurface_facets[fid].second;
        Cell_handle cell = lookup_cell(cell_by_index, cell_index);
        if (cell == Cell_handle()) {
            continue;
        }
        Vertex_handle other0 = nullptr, other1 = nullptr;
        for (int j = 0; j < 4; ++j) {
            if (j == facet_index) continue;
            Vertex_handle vh = cell->vertex(j);
            if (vh == v_handle) continue;
            if (other0 == nullptr) other0 = vh;
            else { other1 = vh; break; }
        }

        if (other0 != Vertex_handle() && other1 != Vertex_handle()) {
            neighbor_to_facets[other0->info().index].push_back(fid);
            neighbor_to_facets[other1->info().index].push_back(fid);
        }
    }

    // Build adjacency using stored edge matching states
    std::vector<std::vector<int>> adj(n);
    std::vector<Edge> incident_edges;
    dt.incident_edges(v_handle, std::back_inserter(incident_edges));
    std::unordered_map<int, Edge> neighbor_edge;
    for (const Edge& e : incident_edges) {
        Vertex_handle a = e.first->vertex(e.second);
        Vertex_handle b = e.first->vertex(e.third);
        if (a != v_handle && b != v_handle) continue;
        Vertex_handle other = (a == v_handle) ? b : a;
        neighbor_edge[other->info().index] = e;
    }

    for (const auto& [neighbor_index, facet_ids] : neighbor_to_facets) {
        if (facet_ids.size() < 2) continue;

        if (facet_ids.size() == 2) {
            adj[facet_ids[0]].push_back(facet_ids[1]);
            adj[facet_ids[1]].push_back(facet_ids[0]);
            continue;
        }

        // Ambiguous edge: use stored matching state
        EdgeKey ekey = EdgeKey::FromVertexIndices(v_global_index, neighbor_index);
        bool applied_matching = false;

        BIPOLAR_MATCH_METHOD method = BIPOLAR_MATCH_METHOD::SEP_NEG;
        auto state_it = edge_matching_states.find(ekey);
        if (state_it != edge_matching_states.end()) {
            method = state_it->second.method;
        }

        auto edge_it = neighbor_edge.find(neighbor_index);
        if (edge_it != neighbor_edge.end()) {
            auto pairs = compute_edge_bipolar_matching(dt, edge_it->second, method, cell_by_index);
            if (!pairs.empty()) {
                std::unordered_set<int> group_set(facet_ids.begin(), facet_ids.end());
                std::unordered_set<int> covered;
                int connections = 0;

                for (const auto& [k0, k1] : pairs) {
                    auto it0 = facetkey_to_id.find(k0);
                    auto it1 = facetkey_to_id.find(k1);
                    if (it0 == facetkey_to_id.end() || it1 == facetkey_to_id.end()) continue;

                    int id0 = it0->second, id1 = it1->second;
                    if (group_set.count(id0) && group_set.count(id1)) {
                        adj[id0].push_back(id1);
                        adj[id1].push_back(id0);
                        covered.insert(id0);
                        covered.insert(id1);
                        connections++;
                    }
                }

                if (covered.size() == facet_ids.size() && connections == static_cast<int>(facet_ids.size() / 2)) {
                    applied_matching = true;
                }
            }
        }

    }

    // Extract connected components as cycles
    std::vector<bool> visited(n, false);
    std::vector<std::vector<std::pair<int, int>>> cycles;

    for (int start = 0; start < n; ++start) {
        if (visited[start]) continue;

        std::vector<std::pair<int, int>> cycle;
        std::stack<int> stack;
        stack.push(start);

        while (!stack.empty()) {
            int curr = stack.top();
            stack.pop();
            if (visited[curr]) continue;
            visited[curr] = true;
            cycle.push_back(isosurface_facets[curr]);
            for (int neighbor : adj[curr]) {
                if (!visited[neighbor]) stack.push(neighbor);
            }
        }
        cycles.push_back(cycle);
    }

    v_handle->info().facet_cycles = cycles;
    write_facet_cycle_indices_for_vertex(v_handle, cell_by_index);
}

ModifyCyclesResult modify_cycles_pass(Delaunay& dt) {
    ModifyCyclesResult result;
    const int MAX_ITERATIONS = 10;

    const std::vector<Cell_handle> cell_by_index = build_cell_index_vector(dt);

    // Store only edges whose matching has been modified; all other edges implicitly use SEP_NEG.
    std::unordered_map<EdgeKey, EdgeMatchingState, EdgeKeyHash> edge_matching_states;

    // Precompute ambiguous edges (>=4 incident valid isosurface facets) since only these can change under flips.
    std::vector<Edge> ambiguous_edges;
    ambiguous_edges.reserve(1024);
    for (auto eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); ++eit) {
        const Edge& e = *eit;
        Vertex_handle a = e.first->vertex(e.second);
        Vertex_handle b = e.first->vertex(e.third);
        if (dt.is_infinite(a) || dt.is_infinite(b)) continue;
        if (a->info().is_dummy || b->info().is_dummy) continue;
        if (!a->info().active || !b->info().active) continue;
        if (!edge_has_at_least_n_valid_isosurface_facets(dt, e, cell_by_index, 4)) continue;
        ambiguous_edges.push_back(e);
    }

    if (ambiguous_edges.empty()) {
        return result;
    }

    std::vector<Vertex_handle> dirty_vertices;
    dirty_vertices.reserve(1024);

    // Iterate until no conflicts or max iterations.
    for (int iter = 0; iter < MAX_ITERATIONS; ++iter) {
        result.iterations = iter + 1;
        bool found_conflict = false;
        int conflicts_this_iteration = 0;
        dirty_vertices.clear();

        for (const Edge& e : ambiguous_edges) {
            Vertex_handle a = e.first->vertex(e.second);
            Vertex_handle b = e.first->vertex(e.third);
            if (dt.is_infinite(a) || dt.is_infinite(b)) continue;
            if (a->info().is_dummy || b->info().is_dummy) continue;
            if (!a->info().active || !b->info().active) continue;

            EdgeKey ekey = EdgeKey::FromVertexIndices(a->info().index, b->info().index);
            auto state_it = edge_matching_states.find(ekey);
            if (state_it != edge_matching_states.end() && state_it->second.flipped) {
                continue;
            }

            BIPOLAR_MATCH_METHOD method = BIPOLAR_MATCH_METHOD::SEP_NEG;
            if (state_it != edge_matching_states.end()) {
                method = state_it->second.method;
            }

            auto pairs = compute_edge_bipolar_matching(dt, e, method, cell_by_index);
            if (pairs.size() < 2) continue;

            // Check for duplicate cycle pairs (same endpoints in terms of (vertex-cycle, vertex-cycle)).
            Vertex_handle v0 = a;
            Vertex_handle v1 = b;
            if (v0->info().index > v1->info().index) {
                std::swap(v0, v1);
            }

            bool has_duplicate = false;
            std::vector<std::pair<int, int>> seen_cycle_pairs;
            seen_cycle_pairs.reserve(pairs.size());

            for (const auto& [k0, k1] : pairs) {
                auto cycle_index_for_vertex_on_facet = [&](Vertex_handle vh, const FacetKey& fk) -> int {
                    Cell_handle cell = lookup_cell(cell_by_index, fk.cell_index);
                    if (cell == Cell_handle()) {
                        return -1;
                    }
                    if (fk.facet_index < 0 || fk.facet_index >= 4) {
                        return -1;
                    }
                    const int v_local = find_vertex_index_in_cell(cell, vh);
                    if (v_local < 0 || v_local >= 4) {
                        return -1;
                    }
                    const int anchor = (fk.facet_index + 1) % 4;
                    const int slot = (v_local - anchor + 4) % 4;
                    if (slot < 0 || slot >= 3) {
                        return -1;
                    }
                    return cell->info().facet_info[fk.facet_index].dualCellEdgeIndex[slot];
                };

                int c0 = cycle_index_for_vertex_on_facet(v0, k0);
                int c1 = cycle_index_for_vertex_on_facet(v1, k0);
                if (c0 < 0 || c1 < 0) {
                    c0 = cycle_index_for_vertex_on_facet(v0, k1);
                    c1 = cycle_index_for_vertex_on_facet(v1, k1);
                }

                if (c0 < 0 || c1 < 0) {
                    continue;
                }

                const std::pair<int, int> cycle_pair = std::make_pair(c0, c1);
                for (const auto& prev : seen_cycle_pairs) {
                    if (prev == cycle_pair) {
                        has_duplicate = true;
                        break;
                    }
                }
                if (has_duplicate) {
                    break;
                }
                seen_cycle_pairs.push_back(cycle_pair);
            }

            if (!has_duplicate) {
                continue;
            }

            EdgeMatchingState& state = edge_matching_states[ekey];
            flip_edge_matching(state);
            found_conflict = true;
            ++conflicts_this_iteration;
            ++result.total_flips;
            dirty_vertices.push_back(a);
            dirty_vertices.push_back(b);
        }

        result.problematic_edges = conflicts_this_iteration;

        if (!found_conflict) {
            break;
        }

        std::sort(dirty_vertices.begin(), dirty_vertices.end(),
                  [](Vertex_handle lhs, Vertex_handle rhs) {
                      return lhs->info().index < rhs->info().index;
                  });
        dirty_vertices.erase(
            std::unique(dirty_vertices.begin(), dirty_vertices.end(),
                        [](Vertex_handle lhs, Vertex_handle rhs) {
                            return lhs->info().index == rhs->info().index;
                        }),
            dirty_vertices.end());

        for (Vertex_handle vh : dirty_vertices) {
            recompute_cycles_for_vertex(dt, vh, cell_by_index, edge_matching_states);
        }

        // Reset flipped flags for next iteration.
        for (auto& [k, state] : edge_matching_states) {
            state.flipped = false;
        }
    }

    if (indicator && result.total_flips > 0) {
        std::cout << "  ModCyc: " << result.total_flips << " edge flips ("
                  << result.iterations << " iterations)" << std::endl;
    }

    DEBUG_PRINT("[CYC-MOD] mod_cyc: flips=" << result.total_flips
                << ", iterations=" << result.iterations);

    return result;
}
