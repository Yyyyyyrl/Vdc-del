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
#include <cmath>
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
    constexpr double NEAR_TOL2 = 1e-16; // ~1e-8 distance squared
    TriIntersectionTrace trace;
    auto segment_interior_intersects_triangle = [](
        const Point& shared,
        const Point& other,
        const Triangle_3& tri
    ) -> bool {
        // Treat intersections at the shared endpoint as trivial by testing an open segment
        // starting slightly away from `shared`.
        constexpr double t = 1e-6;
        const double dx = other.x() - shared.x();
        const double dy = other.y() - shared.y();
        const double dz = other.z() - shared.z();
        const double len2 = dx * dx + dy * dy + dz * dz;
        if (len2 < 1e-24) {
            return false;
        }

        const Point nudged(
            shared.x() + t * dx,
            shared.y() + t * dy,
            shared.z() + t * dz);
        if (nudged == other) {
            return false;
        }

        return CGAL::do_intersect(Segment3(nudged, other), tri);
    };

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

    // Shared single vertex: triangles always intersect at that vertex, so avoid triangle/triangle
    // intersection tests. A non-trivial intersection must involve the opposite edge of one triangle
    // intersecting the other triangle away from the shared vertex.
    const Point& s = shared_point;

    std::array<Point, 2> t1_other;
    std::array<Point, 2> t2_other;
    int t1_count = 0;
    int t2_count = 0;
    for (const Point& p : t1v) {
        if (p != s && t1_count < 2) {
            t1_other[t1_count++] = p;
        }
    }
    for (const Point& p : t2v) {
        if (p != s && t2_count < 2) {
            t2_other[t2_count++] = p;
        }
    }

    if (t1_count != 2 || t2_count != 2) {
        return false;
    }

    const Segment3 t1_opposite(t1_other[0], t1_other[1]);
    const Segment3 t2_opposite(t2_other[0], t2_other[1]);

    if (CGAL::do_intersect(t1_opposite, t2)) {
        trace.kind = TriIntersectionKind::SHARED_VERTEX_T1_OPPOSITE_SEGMENT;
        if (trace_out) {
            *trace_out = trace;
        }
        return true;
    }
    if (CGAL::do_intersect(t2_opposite, t1)) {
        trace.kind = TriIntersectionKind::SHARED_VERTEX_T2_OPPOSITE_SEGMENT;
        if (trace_out) {
            *trace_out = trace;
        }
        return true;
    }

    // Handle cases where the intersection segment touches the shared vertex and extends into
    // the interior of the opposite triangle (a common "foldover" failure mode in a fan).
    if (segment_interior_intersects_triangle(s, t1_other[0], t2)) {
        trace.kind = TriIntersectionKind::SHARED_VERTEX_INTERIOR_SEGMENT_T1_0;
        if (trace_out) {
            *trace_out = trace;
        }
        return true;
    }
    if (segment_interior_intersects_triangle(s, t1_other[1], t2)) {
        trace.kind = TriIntersectionKind::SHARED_VERTEX_INTERIOR_SEGMENT_T1_1;
        if (trace_out) {
            *trace_out = trace;
        }
        return true;
    }
    if (segment_interior_intersects_triangle(s, t2_other[0], t1)) {
        trace.kind = TriIntersectionKind::SHARED_VERTEX_INTERIOR_SEGMENT_T2_0;
        if (trace_out) {
            *trace_out = trace;
        }
        return true;
    }
    if (segment_interior_intersects_triangle(s, t2_other[1], t1)) {
        trace.kind = TriIntersectionKind::SHARED_VERTEX_INTERIOR_SEGMENT_T2_1;
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

    auto kind_name = [](TriIntersectionKind k) -> const char* {
        switch (k) {
            case TriIntersectionKind::NONE: return "NONE";
            case TriIntersectionKind::DISJOINT_TRIANGLE_TRIANGLE: return "DISJOINT_TRIANGLE_TRIANGLE";
            case TriIntersectionKind::SHARED_VERTEX_T1_OPPOSITE_SEGMENT: return "SHARED_VERTEX_T1_OPPOSITE_SEGMENT";
            case TriIntersectionKind::SHARED_VERTEX_T2_OPPOSITE_SEGMENT: return "SHARED_VERTEX_T2_OPPOSITE_SEGMENT";
            case TriIntersectionKind::SHARED_VERTEX_INTERIOR_SEGMENT_T1_0: return "SHARED_VERTEX_INTERIOR_SEGMENT_T1_0";
            case TriIntersectionKind::SHARED_VERTEX_INTERIOR_SEGMENT_T1_1: return "SHARED_VERTEX_INTERIOR_SEGMENT_T1_1";
            case TriIntersectionKind::SHARED_VERTEX_INTERIOR_SEGMENT_T2_0: return "SHARED_VERTEX_INTERIOR_SEGMENT_T2_0";
            case TriIntersectionKind::SHARED_VERTEX_INTERIOR_SEGMENT_T2_1: return "SHARED_VERTEX_INTERIOR_SEGMENT_T2_1";
            default: return "UNKNOWN";
        }
    };

    auto dump_tri_pair = [&](const Triangle_3& a, const Triangle_3& b) {
        DEBUG_PRINT("[DEL-SELFI-INT]   t1=("
                    << a.vertex(0) << ", " << a.vertex(1) << ", " << a.vertex(2) << ")");
        DEBUG_PRINT("[DEL-SELFI-INT]   t2=("
                    << b.vertex(0) << ", " << b.vertex(1) << ", " << b.vertex(2) << ")");
    };

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
                    if (log_within) {
                        static std::unordered_set<int> logged_vertices;
                        if (logged_vertices.insert(v_idx).second) {
                            DEBUG_PRINT("[DEL-SELFI-WITHIN] Vertex " << v_idx
                                        << " cycle " << c
                                        << " tri(" << i << "," << j << ")"
                                        << " kind=" << kind_name(info.kind)
                                        << " shared_count=" << info.shared_count);
                        }
                    }
                    if (trace_details) {
                        DEBUG_PRINT("[DEL-SELFI-INT] Vertex " << v_idx
                                    << " within-cycle cycle " << c
                                    << " tri(" << i << "," << j << ")"
                                    << " kind=" << kind_name(info.kind)
                                    << " shared_count=" << info.shared_count);
                        dump_tri_pair(tris[i], tris[j]);
                    }
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
    const UnifiedGrid& grid
) {
    // Sphere radius for multi-cycle isovertex projection: proportional to the
    // minimum incident Delaunay edge length at the site.
    std::vector<Cell_handle> incident_cells;
    incident_cells.reserve(128);
    dt.incident_cells(v, std::back_inserter(incident_cells));

    // Compute minimum incident edge length (excluding dummy/infinite neighbors).
    double min_edge_length = std::numeric_limits<double>::infinity();
    {
        std::set<Vertex_handle> neighbor_vertices;
        for (Cell_handle ch : incident_cells) {
            if (ch == Cell_handle() || dt.is_infinite(ch)) continue;
            for (int i = 0; i < 4; ++i) {
                Vertex_handle vh = ch->vertex(i);
                if (vh != v && !vh->info().is_dummy && !dt.is_infinite(vh)) {
                    neighbor_vertices.insert(vh);
                }
            }
        }

        const Point p = v->point();
        for (Vertex_handle vh : neighbor_vertices) {
            const Point q = vh->point();
            const double dx = p.x() - q.x();
            const double dy = p.y() - q.y();
            const double dz = p.z() - q.z();
            const double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
            min_edge_length = std::min(min_edge_length, dist);
        }
    }

    // Also consider grid spacing as a fallback (isolated vertex).
    const double grid_spacing = std::min({grid.spacing[0], grid.spacing[1], grid.spacing[2]});
    if (!std::isfinite(min_edge_length)) {
        min_edge_length = grid_spacing;
    }

    return 0.1 * min_edge_length;
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

//! @brief Compute a direction that separates new cycle facets from old cycle facets.
//!
//! For each facet in the cycle, computes the unit outward normal (pointing toward the
//! positive cell, i.e., where circumcenter_scalar >= isovalue). Then uses CGAL's
//! Min_sphere_of_spheres_d to find the smallest enclosing sphere of these normals
//! (treated as points). The center of this sphere gives the separation direction.
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
    float separation_direction[3]
) {
    if (cycle_facets.empty()) {
        return false;
    }

    // Collect outward unit normals for each facet in the cycle.
    // Each facet is stored with its positive cell (where scalar >= isovalue).
    // The outward normal points from the facet toward the positive cell's circumcenter.
    std::vector<Point> normals_as_points;
    normals_as_points.reserve(cycle_facets.size());

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

        // Get the three vertices of the facet (opposite to facet_idx vertex)
        std::array<Point, 3> facet_verts;
        bool ok = true;
        for (int k = 0; k < 3; ++k) {
            const int cell_vertex_idx = (facet_idx + 1 + k) % 4;
            Vertex_handle vh = cell->vertex(cell_vertex_idx);
            if (vh == Vertex_handle() || dt.is_infinite(vh)) {
                ok = false;
                break;
            }
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

        // Determine if normal points toward the positive cell or away.
        // The positive cell is `cell` (by FacetKey convention).
        // The facet lies between `cell` (positive) and its neighbor (negative).
        // We want the normal pointing toward `cell`, i.e., toward circumcenter of `cell`.
        const Point& circumcenter = cell->info().circumcenter;
        const Vector3 to_circumcenter(facet_verts[0], circumcenter);
        if (vec_dot(normal, to_circumcenter) < 0) {
            normal = -normal;  // Flip to point toward positive cell
        }

        // Store the unit normal as a point (vector from origin)
        normals_as_points.emplace_back(normal.x(), normal.y(), normal.z());
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

    // If center is very close to origin, normals span all directions - no separation
    constexpr double eps = 1e-6;
    if (center_len < eps) {
        return false;
    }

    // Normalize and output
    separation_direction[0] = static_cast<float>(cx / center_len);
    separation_direction[1] = static_cast<float>(cy / center_len);
    separation_direction[2] = static_cast<float>(cz / center_len);

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

ResolutionResult try_resolve_multicycle_by_cycle_separation_tests(
    Vertex_handle v,
    std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    std::vector<Point>* geometric_attempt_positions_out,
    bool use_sep_dir
) {
    const auto& cycles = v->info().facet_cycles;
    const int num_cycles = static_cast<int>(cycles.size());
    if (num_cycles < 2 || static_cast<int>(cycle_isovertices.size()) < 2) {
        return {ResolutionStatus::NOT_NEEDED, ResolutionStrategy::NONE};
    }

    const std::vector<Point> baseline_positions = cycle_isovertices;

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
	        candidates.reserve(6);

        const Point vpos0 = baseline_positions[0];
        const Point vpos1 = baseline_positions[1];

        auto consider_positions = [&](const char* label, std::vector<Point> positions) {
            CandidateA2 cand;
            cand.positions = std::move(positions);
            cand.label = label;
            cand.move2 = 0.0;
            for (size_t i = 0; i < cand.positions.size() && i < baseline_positions.size(); ++i) {
                cand.move2 += squared_distance(cand.positions[i], baseline_positions[i]);
            }
            candidates.push_back(std::move(cand));
        };

        // ====================================================================
        // Old approach: reflection-based candidates (only when !use_sep_dir)
        // ====================================================================
        if (!use_sep_dir) {
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
        }

        // ====================================================================
        // New approach: Separation-direction-based candidates (only when use_sep_dir)
        // ====================================================================
        // Compute separation direction for each cycle and place isovertices along
        // those directions at various radii from the Delaunay vertex.
        if (use_sep_dir) {
            float sep_dir0[3] = {0, 0, 0};
            float sep_dir1[3] = {0, 0, 0};
            const bool has_sep0 = compute_cycle_separating_direction(
                v, cycles[0], cell_by_index, dt, sep_dir0);
            const bool has_sep1 = compute_cycle_separating_direction(
                v, cycles[1], cell_by_index, dt, sep_dir1);

            if (has_sep0 || has_sep1) {
                // Compute a reasonable radius: distance from center to baseline positions
                const double r0 = std::sqrt(squared_distance(center, vpos0));
                const double r1 = std::sqrt(squared_distance(center, vpos1));
                const double avg_r = (r0 + r1) * 0.5;

                // Helper: move from center in direction dir by distance r
                auto move_in_dir = [&center](const float dir[3], double r) -> Point {
                    return Point(
                        center.x() + r * dir[0],
                        center.y() + r * dir[1],
                        center.z() + r * dir[2]
                    );
                };

                // Candidate: move cycle0 along its separation direction
                if (has_sep0) {
                    std::vector<Point> cand = baseline_positions;
                    cand[0] = move_in_dir(sep_dir0, avg_r);
                    consider_positions("sep_dir0", std::move(cand));

                    // Also try opposite direction
                    float neg_sep_dir0[3] = {-sep_dir0[0], -sep_dir0[1], -sep_dir0[2]};
                    std::vector<Point> cand_neg = baseline_positions;
                    cand_neg[0] = move_in_dir(neg_sep_dir0, avg_r);
                    consider_positions("sep_dir0_neg", std::move(cand_neg));
                }

                // Candidate: move cycle1 along its separation direction
                if (has_sep1) {
                    std::vector<Point> cand = baseline_positions;
                    cand[1] = move_in_dir(sep_dir1, avg_r);
                    consider_positions("sep_dir1", std::move(cand));

                    // Also try opposite direction
                    float neg_sep_dir1[3] = {-sep_dir1[0], -sep_dir1[1], -sep_dir1[2]};
                    std::vector<Point> cand_neg = baseline_positions;
                    cand_neg[1] = move_in_dir(neg_sep_dir1, avg_r);
                    consider_positions("sep_dir1_neg", std::move(cand_neg));
                }

                // Candidate: move both cycles along their separation directions
                if (has_sep0 && has_sep1) {
                    std::vector<Point> cand = baseline_positions;
                    cand[0] = move_in_dir(sep_dir0, avg_r);
                    cand[1] = move_in_dir(sep_dir1, avg_r);
                    consider_positions("sep_dir_both", std::move(cand));

                    // Try opposite directions
                    float neg_sep_dir0[3] = {-sep_dir0[0], -sep_dir0[1], -sep_dir0[2]};
                    float neg_sep_dir1[3] = {-sep_dir1[0], -sep_dir1[1], -sep_dir1[2]};
                    std::vector<Point> cand_neg = baseline_positions;
                    cand_neg[0] = move_in_dir(neg_sep_dir0, avg_r);
                    cand_neg[1] = move_in_dir(neg_sep_dir1, avg_r);
                    consider_positions("sep_dir_both_neg", std::move(cand_neg));
                }
            }
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

static ResolutionResult resolve_multicycle_self_intersection_at_vertex(
    Vertex_handle v,
    std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    const CycleIsovertexOptions& options,
    std::unordered_set<int>* local_unresolved_dumped_vertices
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
    if (trace) {
        DEBUG_PRINT("[DEL-SELFI-TRACE] Vertex " << v->info().index
                    << ": detected self-intersection (cycles=" << num_cycles << ").");
        for (int c = 0; c < num_cycles; ++c) {
            const Point& p = baseline_positions[static_cast<size_t>(c)];
            DEBUG_PRINT("[DEL-SELFI-TRACE]   baseline[" << c << "]=("
                        << p.x() << ", " << p.y() << ", " << p.z() << ")");
        }
        if (const auto w = find_self_intersection_witness(v, baseline_positions, dt, cell_by_index)) {
            DEBUG_PRINT("[DEL-SELFI-TRACE]   witness: cycle(" << w->cycle0 << "," << w->cycle1
                        << ") tri(" << w->tri0 << "," << w->tri1 << ")");
        }
    }
    vdc_debug::maybe_dump_selfi_cycle_metadata(dt, v, cell_by_index);
    vdc_debug::maybe_dump_selfi_stage(dt, v, baseline_positions, cell_by_index, "baseline");

    std::vector<Point> geometric_attempt_positions;
	    const ResolutionResult separation_result =
	        try_resolve_multicycle_by_cycle_separation_tests(
	            v,
	            cycle_isovertices,
	            dt,
	            cell_by_index,
	            &geometric_attempt_positions,
	            options.use_sep_dir);
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

            // Also dump the deterministic A-only candidate set (both reflections for 2-cycle,
            // per-cycle reflections for 3+ cycles, plus the pairwise-pass attempt).
            vdc_debug::dump_multi_isov_trace_case(
                out_dir, v, baseline_positions, dt, cell_by_index);

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
    bool position_multi_isov_on_delv
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
            // Multi-cycle vertex: compute centroids and project to sphere
            // ----------------------------------------------------------------
            // For vertices where the isosurface passes through multiple times
            // (creating multiple cycles), we compute a separate isovertex per cycle.
            // Each isovertex is positioned by:
            //   1. Computing the centroid of Voronoi edge / isosurface intersections
            //   2. Projecting this centroid onto a sphere around the Delaunay site
            // This provides directional separation between cycles.
            if (position_multi_isov_on_delv) {
                for (size_t c = 0; c < cycles.size(); ++c) {
                    isovertices[c] = cube_center;
                }
            } else {
                const double sphere_radius = compute_sphere_radius(vit, dt, grid);
                for (size_t c = 0; c < cycles.size(); ++c) {
                    Point centroid = compute_centroid_of_voronoi_edge_and_isosurface(
                        dt, grid, isovalue, vit, cycles[c], cell_by_index);

                    isovertices[c] = project_to_sphere(centroid, cube_center, sphere_radius);
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

static std::vector<Vertex_handle> resolve_multicycle_self_intersections(
    Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    const std::vector<Vertex_handle>& multi_cycle_vertices,
    IsovertexComputationStats& stats,
    const CycleIsovertexOptions& options
) {
    std::vector<Vertex_handle> modified_multi_cycle_vertices;
    std::unordered_set<int> modified_multi_cycle_indices;
    std::unordered_set<int> local_unresolved_dumped_vertices;
    std::unordered_set<int>* local_unresolved_dumped_vertices_ptr = nullptr;

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
                local_unresolved_dumped_vertices_ptr);

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

        // Accumulate pass statistics into overall stats.
        stats.self_intersection_detected += pass_detected;
        stats.self_intersection_resolved += pass_resolved;

        if (pass_detected > 0) {
            DEBUG_PRINT("[DEL-SELFI] multi-cycle pass " << (pass + 1) << "/" << MAX_MULTI_CYCLE_PASSES
                        << ": detected=" << pass_detected
                        << ", resolved=" << pass_resolved
                        << ", unresolved=" << pass_unresolved);
        }

        // Early exit if no changes were made this pass.
        if (!any_changed) {
            break;
        }
    }

    // ========================================================================
    // Joint resolution pass for coupled multi-cycle sites
    // ========================================================================
    // After individual resolution attempts, identify still-unresolved multi-cycle
    // vertices that share triangle geometry (coupled sites). For these, attempt
    // joint resolution by evaluating combinations of reflection candidates.
    {
        // Collect still-unresolved vertices
        std::vector<Vertex_handle> unresolved;
        std::unordered_set<int> unresolved_indices;
        for (Vertex_handle vh : multi_cycle_vertices) {
            if (!vh->info().active) continue;
            if (vh->info().is_dummy) continue;
            if (vh->info().facet_cycles.size() < 2) continue;
            const auto& isovertices = vh->info().cycle_isovertices;
            if (check_self_intersection(vh, isovertices, dt, cell_by_index)) {
                unresolved.push_back(vh);
                unresolved_indices.insert(vh->info().index);
            }
        }

        if (unresolved.size() >= 2) {
            // For each pair of unresolved vertices, check if they share triangle geometry
            for (size_t i = 0; i < unresolved.size(); ++i) {
                for (size_t j = i + 1; j < unresolved.size(); ++j) {
                    Vertex_handle v1 = unresolved[i];
                    Vertex_handle v2 = unresolved[j];

                    // Check if v1's triangles use v2's isovertices or vice versa
                    auto shares_geometry = [&](Vertex_handle va, Vertex_handle vb) -> bool {
                        const auto& cycles = va->info().facet_cycles;
                        for (size_t c = 0; c < cycles.size(); ++c) {
                            for (const auto& [cell_idx, facet_idx] : cycles[c]) {
                                Cell_handle ch = lookup_cell(cell_by_index, cell_idx);
                                if (ch == Cell_handle()) continue;
                                for (int t = 0; t < 3; ++t) {
                                    int cv = (facet_idx + 1 + t) % 4;
                                    if (ch->vertex(cv) == vb) return true;
                                }
                            }
                        }
                        return false;
                    };

                    if (!shares_geometry(v1, v2) && !shares_geometry(v2, v1)) {
                        continue; // Not coupled
                    }

                    // These are coupled vertices - attempt joint resolution
                    const bool trace = vdc_debug::trace_selfi_enabled(v1->info().index) ||
                                       vdc_debug::trace_selfi_enabled(v2->info().index);

                    const std::vector<Point> baseline1 = v1->info().cycle_isovertices;
                    const std::vector<Point> baseline2 = v2->info().cycle_isovertices;
                    const Point center1 = v1->point();
                    const Point center2 = v2->point();

                    if (baseline1.size() < 2 || baseline2.size() < 2) continue;

                    // Generate reflection candidates for each vertex
                    auto generate_candidates = [](const std::vector<Point>& baseline, const Point& center)
                        -> std::vector<std::vector<Point>> {
                        std::vector<std::vector<Point>> cands;
                        if (baseline.size() != 2) return cands;
                        const Point& p0 = baseline[0];
                        const Point& p1 = baseline[1];
                        auto reflect = [&](const Point& p) {
                            return Point(2.0*center.x() - p.x(), 2.0*center.y() - p.y(), 2.0*center.z() - p.z());
                        };
                        // baseline, reflect_A, reflect_B, reflect_self0, reflect_self1, reflect_self_both, reflect_cross_both
                        cands.push_back({p0, p1});
                        cands.push_back({reflect(p1), p1});
                        cands.push_back({p0, reflect(p0)});
                        cands.push_back({reflect(p0), p1});
                        cands.push_back({p0, reflect(p1)});
                        cands.push_back({reflect(p0), reflect(p1)});
                        cands.push_back({reflect(p1), reflect(p0)});
                        return cands;
                    };

                    auto cands1 = generate_candidates(baseline1, center1);
                    auto cands2 = generate_candidates(baseline2, center2);

                    if (cands1.empty() || cands2.empty()) continue;

                    // Evaluate all joint candidates
                    int best_total = std::numeric_limits<int>::max();
                    int best_i = -1, best_j = -1;

                    for (size_t ci = 0; ci < cands1.size(); ++ci) {
                        for (size_t cj = 0; cj < cands2.size(); ++cj) {
                            // Build position override map
                            PositionOverrideMap overrides;
                            overrides[v1->info().index] = cands1[ci];
                            overrides[v2->info().index] = cands2[cj];

                            // Count total intersections for both vertices with these positions
                            int total = 0;
                            total += count_self_intersection_pairs_with_overrides(
                                v1, cands1[ci], dt, cell_by_index, overrides, 10);
                            total += count_self_intersection_pairs_with_overrides(
                                v2, cands2[cj], dt, cell_by_index, overrides, 10);

                            if (total < best_total) {
                                best_total = total;
                                best_i = static_cast<int>(ci);
                                best_j = static_cast<int>(cj);
                            }

                            if (total == 0) break; // Found perfect solution
                        }
                        if (best_total == 0) break;
                    }

                    if (trace) {
                        DEBUG_PRINT("[DEL-SELFI-JOINT] Coupled v" << v1->info().index
                                    << " + v" << v2->info().index
                                    << ": best_total=" << best_total
                                    << " (cand " << best_i << "," << best_j << ")");
                    }

                    // Accept if the joint candidate resolves (or at least improves)
                    if (best_total == 0 && best_i >= 0 && best_j >= 0) {
                        v1->info().cycle_isovertices = cands1[static_cast<size_t>(best_i)];
                        v2->info().cycle_isovertices = cands2[static_cast<size_t>(best_j)];
                        stats.self_intersection_resolved += 2;
                        if (trace) {
                            DEBUG_PRINT("[DEL-SELFI-JOINT] RESOLVED coupled v" << v1->info().index
                                        << " + v" << v2->info().index);
                        }
                    }
                }
            }
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
    const std::vector<Cell_handle>& cell_by_index
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

    const bool trace = vdc_debug::trace_selfi_enabled(v->info().index);
    const Point center = v->point();
    const std::vector<Point> baseline_positions = cycle_isovertices;
    const int baseline_intersections = count_self_intersection_pairs(
        v, baseline_positions, dt, cell_by_index);

    struct CandidateFan {
        Point p;
        const char* label = nullptr;
        double move2 = 0.0;
        int intersections = std::numeric_limits<int>::max();
    };

    std::vector<CandidateFan> candidates;
    candidates.reserve(8);

    auto add_candidate = [&](const Point& p, const char* label) {
        // Avoid duplicates (exact coordinates are stable for deterministic candidates).
        for (const auto& cand : candidates) {
            if (cand.p == p) {
                return;
            }
        }
        CandidateFan cand;
        cand.p = p;
        cand.label = label;
        cand.move2 = squared_distance(p, baseline_p);
        candidates.push_back(std::move(cand));
    };

    // Deterministic finite candidates: shrink toward the site and optional reflections.
    add_candidate(Point(
        0.75 * baseline_p.x() + 0.25 * center.x(),
        0.75 * baseline_p.y() + 0.25 * center.y(),
        0.75 * baseline_p.z() + 0.25 * center.z()),
        "shrink_75");
    add_candidate(Point(
        0.5 * (baseline_p.x() + center.x()),
        0.5 * (baseline_p.y() + center.y()),
        0.5 * (baseline_p.z() + center.z())),
        "shrink_50");
    add_candidate(Point(
        0.25 * baseline_p.x() + 0.75 * center.x(),
        0.25 * baseline_p.y() + 0.75 * center.y(),
        0.25 * baseline_p.z() + 0.75 * center.z()),
        "shrink_25");
    add_candidate(center, "center");

    const Point reflect = reflect_through_center(center, baseline_p);
    add_candidate(reflect, "reflect");
    add_candidate(Point(
        0.5 * (reflect.x() + center.x()),
        0.5 * (reflect.y() + center.y()),
        0.5 * (reflect.z() + center.z())),
        "reflect_shrink_50");

    // Candidate: steer toward the boundary centroid of the fan (changes direction, not only radius).
    {
        double sx = 0.0, sy = 0.0, sz = 0.0;
        int cnt = 0;
        double min_boundary_r2 = std::numeric_limits<double>::max();
        for (const auto& tri : tris) {
            // In collect_cycle_triangles(), vertex(0) is the cycle isovertex.
            for (int k = 1; k < 3; ++k) {
                const Point q = tri.vertex(k);
                sx += q.x();
                sy += q.y();
                sz += q.z();
                ++cnt;
                min_boundary_r2 = std::min(min_boundary_r2, squared_distance(center, q));
            }
        }
        if (cnt > 0) {
            const Point centroid(sx / cnt, sy / cnt, sz / cnt);
            add_candidate(centroid, "boundary_centroid");

            const Vector3 dir(center, centroid);
            const double len = vec_norm(dir);
            if (len > 1e-12) {
                // Rotate deterministically by projecting onto a sphere around the Delaunay site.
                // If the baseline is at the site (radius ~0), use a small radius derived from
                // the fan boundary to allow a directional move.
                const double baseline_r2 = squared_distance(center, baseline_p);
                const bool baseline_is_center = (baseline_r2 <= 1e-24);

                auto add_sphere_candidate = [&](double r, const char* label_pos, const char* label_neg) {
                    if (r <= 1e-12) {
                        return;
                    }
                    const Point proj_pos(
                        center.x() + r * dir.x() / len,
                        center.y() + r * dir.y() / len,
                        center.z() + r * dir.z() / len);
                    const Point proj_neg(
                        center.x() - r * dir.x() / len,
                        center.y() - r * dir.y() / len,
                        center.z() - r * dir.z() / len);
                    add_candidate(proj_pos, label_pos);
                    add_candidate(proj_neg, label_neg);
                };

                if (!baseline_is_center) {
                    const double r = std::sqrt(baseline_r2);
                    add_sphere_candidate(r, "boundary_centroid_sphere", "boundary_centroid_sphere_neg");
                } else if (min_boundary_r2 < std::numeric_limits<double>::max() && min_boundary_r2 > 1e-24) {
                    const double boundary_r = std::sqrt(min_boundary_r2);
                    add_sphere_candidate(0.10 * boundary_r, "boundary_centroid_sphere_r10", "boundary_centroid_sphere_r10_neg");
                    add_sphere_candidate(0.20 * boundary_r, "boundary_centroid_sphere_r20", "boundary_centroid_sphere_r20_neg");
                    add_sphere_candidate(0.35 * boundary_r, "boundary_centroid_sphere_r35", "boundary_centroid_sphere_r35_neg");
                    add_sphere_candidate(0.50 * boundary_r, "boundary_centroid_sphere_r50", "boundary_centroid_sphere_r50_neg");
                }
            }
        }
    }

    // Candidate: steer along the aggregate fan normal (often resolves foldovers when the isovertex
    // is at the Delaunay site).
    {
        Vector3 sum_n(0, 0, 0);
        for (const auto& tri : tris) {
            const Vector3 n = vec_cross(Vector3(tri.vertex(0), tri.vertex(1)), Vector3(tri.vertex(0), tri.vertex(2)));
            if (n.squared_length() < 1e-24) {
                continue;
            }
            sum_n = sum_n + n;
        }

        const double len = vec_norm(sum_n);
        if (len > 1e-12) {
            const double baseline_r2 = squared_distance(center, baseline_p);
            const bool baseline_is_center = (baseline_r2 <= 1e-24);

            auto add_dir_candidate = [&](double r, const char* label_pos, const char* label_neg) {
                if (r <= 1e-12) {
                    return;
                }
                const Point pos(
                    center.x() + r * sum_n.x() / len,
                    center.y() + r * sum_n.y() / len,
                    center.z() + r * sum_n.z() / len);
                const Point neg(
                    center.x() - r * sum_n.x() / len,
                    center.y() - r * sum_n.y() / len,
                    center.z() - r * sum_n.z() / len);
                add_candidate(pos, label_pos);
                add_candidate(neg, label_neg);
            };

            if (!baseline_is_center) {
                const double r = std::sqrt(baseline_r2);
                add_dir_candidate(r, "fan_normal_sphere", "fan_normal_sphere_neg");
            } else {
                double min_boundary_r2 = std::numeric_limits<double>::max();
                for (const auto& tri : tris) {
                    min_boundary_r2 = std::min(min_boundary_r2, squared_distance(center, tri.vertex(1)));
                    min_boundary_r2 = std::min(min_boundary_r2, squared_distance(center, tri.vertex(2)));
                }
                if (min_boundary_r2 < std::numeric_limits<double>::max() && min_boundary_r2 > 1e-24) {
                    const double boundary_r = std::sqrt(min_boundary_r2);
                    add_dir_candidate(0.10 * boundary_r, "fan_normal_sphere_r10", "fan_normal_sphere_r10_neg");
                    add_dir_candidate(0.20 * boundary_r, "fan_normal_sphere_r20", "fan_normal_sphere_r20_neg");
                    add_dir_candidate(0.35 * boundary_r, "fan_normal_sphere_r35", "fan_normal_sphere_r35_neg");
                    add_dir_candidate(0.50 * boundary_r, "fan_normal_sphere_r50", "fan_normal_sphere_r50_neg");
                }
            }
        }
    }

    CandidateFan* best_resolved = nullptr;
    CandidateFan* best_improved = nullptr;
    CandidateFan* best_any = nullptr;
    for (auto& cand : candidates) {
        std::vector<Point> positions = baseline_positions;
        positions[static_cast<size_t>(cycle_idx)] = cand.p;
        cand.intersections = count_self_intersection_pairs(v, positions, dt, cell_by_index);

        const bool resolved = (cand.intersections == 0);
        const bool improved =
            (!resolved) && (baseline_intersections > 0) && (cand.intersections < baseline_intersections);

        if (trace) {
            DEBUG_PRINT("[DEL-SELFI-FAN] v=" << v->info().index
                        << " cycle=" << cycle_idx
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
            (cand.intersections == best_any->intersections && cand.move2 < best_any->move2)) {
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
                (cand.intersections == best_improved->intersections && cand.move2 < best_improved->move2)) {
                best_improved = &cand;
            }
        }
    }

    if (best_resolved) {
        cycle_isovertices[static_cast<size_t>(cycle_idx)] = best_resolved->p;
        return {ResolutionStatus::RESOLVED, ResolutionStrategy::GEOMETRIC_SEPARATION};
    }

    if (best_improved) {
        cycle_isovertices[static_cast<size_t>(cycle_idx)] = best_improved->p;
        return {ResolutionStatus::UNRESOLVED, ResolutionStrategy::GEOMETRIC_SEPARATION};
    }

    // No improvement: keep baseline.
    cycle_isovertices = baseline_positions;
    return {ResolutionStatus::UNRESOLVED, ResolutionStrategy::GEOMETRIC_SEPARATION};
}

static void resolve_suspect_cycle_fan_foldovers(
    Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    const std::vector<Vertex_handle>& multi_cycle_vertices,
    const std::vector<Vertex_handle>& modified_multi_cycle_vertices,
    IsovertexComputationStats& stats,
    const CycleIsovertexOptions& options
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
                ref.v, ref.cycle, isovertices, dt, cell_by_index);

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
 *     - Multi-cycle vertices: compute Voronoi/isosurface centroids and project to sphere
 *
 *   Stage 2 (resolve_multicycle_self_intersections):
 *     Iteratively resolve self-intersections on multi-cycle vertices.
 *     Uses only deterministic hyperplane separation tests + diametric reflections through the Delaunay site.
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

    // Initialize statistics tracking.
    IsovertexComputationStats stats;

    // ========================================================================
    // Stage 1: Compute initial isovertex positions for all active vertices.
    // ========================================================================
    std::vector<Vertex_handle> multi_cycle_vertices = initialize_cycle_isovertices(
        dt, grid, isovalue, cell_by_index, stats, options.position_multi_isov_on_delv);

    if (options.position_multi_isov_on_delv) {
        DEBUG_PRINT("[DEL-ISOV] Skipping multi-cycle repositioning (-position_multi_isov_on_delv).");
        report_isovertex_statistics(stats);
        return;
    }

    // ========================================================================
    // Stage 2: Resolve self-intersections on multi-cycle vertices.
    // ========================================================================
    const std::vector<Vertex_handle> modified_multi_cycle_vertices = resolve_multicycle_self_intersections(
        dt, cell_by_index, multi_cycle_vertices, stats, options);

    // ========================================================================
    // Stage 3: Resolve within-cycle foldovers (targeted, deterministic).
    // ========================================================================
    if (options.foldover) {
        resolve_suspect_cycle_fan_foldovers(
            dt, cell_by_index, multi_cycle_vertices, modified_multi_cycle_vertices, stats, options);
    }

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
