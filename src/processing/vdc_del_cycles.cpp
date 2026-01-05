//! @file vdc_del_cycles.cpp
//! @brief Implementation of cycle detection for Delaunay-based isosurface extraction.
//!
//! This file implements Steps 8-9 of the Delaunay-based VDC algorithm:
//! - Step 8: Compute cycles of isosurface facets around active vertices
//! - Step 9: Compute isosurface vertex positions for each cycle

#include "processing/vdc_del_isosurface.h"
#include "processing/vdc_del_cycles.h"
#include "core/vdc_debug.h"
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <optional>
#include <array>
#include <limits>
#include <cstdint>
#include <stdexcept>
#include <CGAL/Triangle_3.h>
#include <CGAL/intersections.h>

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

        // Delaunay facet is only marked isosurface on the positive-cell side, so
        // (sA - sB) should be strictly positive.
        const float denom = sA - sB;
        if (denom <= 0.0f) {
            throw std::runtime_error("Invalid circumcenter scalar ordering (denom <= 0).");
        }

        // Linear interpolation to find isosurface crossing point
        // p = ((isovalue - sB) * wA + (sA - isovalue) * wB) / (sA - sB)
        //
        // Using p = wA + t * (wB - wA) gives t = (sA - isovalue) / (sA - sB).
        float t = (sA - isovalue) / denom;

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

//! @brief Check if two triangles intersect (excluding shared edges/vertices)
/*!
 * Uses CGAL's do_intersect for Triangle_3 objects.
 * Returns true if the triangles have a non-trivial intersection
 * (i.e., more than just touching at shared vertices).
 */
static bool triangles_have_nontrivial_intersection(
    const Triangle_3& t1,
    const Triangle_3& t2
) {
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

    // If the triangles are disjoint in vertex set, any intersection is non-trivial.
    if (shared_count == 0) {
        return CGAL::do_intersect(t1, t2);
    }

    // Shared edge/triangle: adjacent triangles are expected to intersect along their shared simplex.
    // Overlap beyond that would indicate degeneracy, which we ignore here for performance.
    if (shared_count >= 2) {
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
        return true;
    }
    if (CGAL::do_intersect(t2_opposite, t1)) {
        return true;
    }

    return false;
}

bool check_self_intersection(
    Vertex_handle v,
    const std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index
) {
    const auto& cycles = v->info().facet_cycles;
    size_t num_cycles = cycles.size();

    if (num_cycles < 1) {
        return false;
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
                if (triangles_have_nontrivial_intersection(tris[i], tris[j])) {
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

double compute_sphere_radius(
    Vertex_handle v,
    const Delaunay& dt,
    const UnifiedGrid& grid
) {
    // Compute minimum incident edge length
    double min_edge_length = std::numeric_limits<double>::max();

    // Get all edges incident to v by examining incident cells
    std::set<Vertex_handle> neighbor_vertices;
    std::vector<Cell_handle> incident_cells;
    dt.incident_cells(v, std::back_inserter(incident_cells));

    for (Cell_handle ch : incident_cells) {
        if (dt.is_infinite(ch)) continue;

        for (int i = 0; i < 4; ++i) {
            Vertex_handle vh = ch->vertex(i);
            if (vh != v && !vh->info().is_dummy && !dt.is_infinite(vh)) {
                neighbor_vertices.insert(vh);
            }
        }
    }

    Point p = v->point();
    for (Vertex_handle vh : neighbor_vertices) {
        Point q = vh->point();
        double dx = p.x() - q.x();
        double dy = p.y() - q.y();
        double dz = p.z() - q.z();
        double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
        min_edge_length = std::min(min_edge_length, dist);
    }

    // Also consider grid spacing as a fallback
    double grid_spacing = std::min({grid.spacing[0], grid.spacing[1], grid.spacing[2]});

    if (min_edge_length == std::numeric_limits<double>::max()) {
        min_edge_length = grid_spacing;
    }

    // Use a fraction of the minimum edge length
    // 0.1 provides enough separation while keeping isovertices close to the Delaunay vertex
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

static CGAL::Orientation opposite_orientation(const CGAL::Orientation o) {
    switch (o) {
        case CGAL::POSITIVE:
            return CGAL::NEGATIVE;
        case CGAL::NEGATIVE:
            return CGAL::POSITIVE;
        default:
            return o;
    }
}

static CGAL::Orientation cell_orientation(const Cell_handle cell) {
    return CGAL::orientation(
        cell->vertex(0)->point(),
        cell->vertex(1)->point(),
        cell->vertex(2)->point(),
        cell->vertex(3)->point());
}

static CGAL::Orientation orientation_with_replaced_delv(
    const Cell_handle cell,
    const Vertex_handle v,
    const Point& p
) {
    Point verts[4];
    for (int i = 0; i < 4; ++i) {
        Vertex_handle vi = cell->vertex(i);
        verts[i] = (vi == v) ? p : vi->point();
    }
    return CGAL::orientation(verts[0], verts[1], verts[2], verts[3]);
}


static StarCellSetPositionMultiIsov build_star_cell_set_position_multi_isov(
    const Delaunay& dt,
    Vertex_handle v0
) {
    StarCellSetPositionMultiIsov out;

    std::vector<Cell_handle> incident_cells;
    dt.incident_cells(v0, std::back_inserter(incident_cells));
    out.cells.reserve(incident_cells.size());
    out.local_by_cell_idx.reserve(incident_cells.size() * 2);

    for (Cell_handle ch : incident_cells) {
        if (dt.is_infinite(ch)) {
            continue;
        }
        const int idx = ch->info().index;
        if (idx < 0) {
            continue;
        }
        const int local = static_cast<int>(out.cells.size());
        if (out.local_by_cell_idx.emplace(idx, local).second) {
            out.cells.push_back(ch);
        }
    }

    return out;
}

static CycleDataPositionMultiIsov build_cycle_data_position_multi_isov(
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    Vertex_handle v0,
    const std::vector<std::pair<int, int>>& cycle_facets
) {
    CycleDataPositionMultiIsov out;
    out.facets.reserve(cycle_facets.size());

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

        Facet f{cell, facet_idx};
        out.facets.push_back(f);

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

static bool are_cycle_cells_separated_from_cellA(
    const Delaunay& dt,
    const StarCellSetPositionMultiIsov& star,
    Vertex_handle v0,
    const std::vector<Facet>& cycle_facets,
    Cell_handle cellA
) {
    if (cycle_facets.empty()) {
        return true;
    }
    if (dt.is_infinite(cellA)) {
        return true;
    }

    const auto itA = star.local_by_cell_idx.find(cellA->info().index);
    if (itA == star.local_by_cell_idx.end()) {
        return true;
    }

    std::vector<char> visited(star.cells.size(), 0);

    std::vector<Cell_handle> stackB;
    stackB.reserve(cycle_facets.size() * 2);

    // Initialize barrier: mark the mirror-side cells as visited, seed with cycle-side cells.
    for (const Facet& f0 : cycle_facets) {
        if (dt.is_infinite(f0.first)) {
            continue;
        }
        const auto it0 = star.local_by_cell_idx.find(f0.first->info().index);
        if (it0 != star.local_by_cell_idx.end()) {
            stackB.push_back(f0.first);
        }

        const Facet f1 = dt.mirror_facet(f0);
        if (!dt.is_infinite(f1.first)) {
            const auto it1 = star.local_by_cell_idx.find(f1.first->info().index);
            if (it1 != star.local_by_cell_idx.end()) {
                visited[static_cast<size_t>(it1->second)] = 1;
            }
        }
    }

    while (!stackB.empty()) {
        Cell_handle cellB = stackB.back();
        stackB.pop_back();

        if (cellB == cellA) {
            return false;
        }

        const auto itB = star.local_by_cell_idx.find(cellB->info().index);
        if (itB == star.local_by_cell_idx.end()) {
            continue;
        }
        const int localB = itB->second;
        if (visited[static_cast<size_t>(localB)]) {
            continue;
        }
        visited[static_cast<size_t>(localB)] = 1;

        // Traverse only within the star of v0: cross facets that contain v0.
        for (int fi = 0; fi < 4; ++fi) {
            if (cellB->vertex(fi) == v0) {
                continue; // facet fi does NOT contain v0
            }
            Cell_handle nb = cellB->neighbor(fi);
            if (dt.is_infinite(nb)) {
                continue;
            }
            const auto itNB = star.local_by_cell_idx.find(nb->info().index);
            if (itNB == star.local_by_cell_idx.end()) {
                continue;
            }
            if (nb == cellA) {
                return false;
            }
            const int localNB = itNB->second;
            if (!visited[static_cast<size_t>(localNB)]) {
                stackB.push_back(nb);
            }
        }
    }

    return true;
}

static bool does_starA_separate_cycles(
    const Delaunay& dt,
    const StarCellSetPositionMultiIsov& star,
    Vertex_handle v0,
    const std::vector<Facet>& cycleA_facets,
    const std::vector<Facet>& cycleB_facets,
    const Point& vposB
) {
    if (cycleA_facets.empty() || cycleB_facets.empty()) {
        return true;
    }

    Cell_handle cellA = cycleA_facets.front().first;
    if (dt.is_infinite(cellA)) {
        return true;
    }

    const bool flag_separateA = are_cycle_cells_separated_from_cellA(dt, star, v0, cycleB_facets, cellA);

    for (Facet fB : cycleB_facets) {
        if (!flag_separateA) {
            fB = dt.mirror_facet(fB);
        }
        Cell_handle cellB = fB.first;
        if (dt.is_infinite(cellB)) {
            continue;
        }

        const int facet_index = fB.second;
        if (facet_index < 0 || facet_index >= 4) {
            continue;
        }

        Vertex_handle vB = cellB->vertex(facet_index);

        const CGAL::Orientation ori_cell = cell_orientation(cellB);
        if (ori_cell == CGAL::COPLANAR) {
            return false;
        }

        const CGAL::Orientation target = opposite_orientation(ori_cell);
        const CGAL::Orientation ori_new = orientation_with_replaced_delv(cellB, vB, vposB);

        if (ori_new != target) {
            return false;
        }
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

static bool does_bisecting_plane_separate_cycles(
    const CycleDataPositionMultiIsov& cycleA,
    const CycleDataPositionMultiIsov& cycleB,
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

static bool find_vertex_positions_separating_cycles(
    const Delaunay& dt,
    const StarCellSetPositionMultiIsov& star,
    Vertex_handle v0,
    const CycleDataPositionMultiIsov& cycleA,
    const CycleDataPositionMultiIsov& cycleB,
    const Point& vposA,
    const Point& vposB,
    Point& new_vposA,
    Point& new_vposB
) {
    new_vposA = vposA;
    new_vposB = vposB;

    if (does_bisecting_plane_separate_cycles(cycleA, cycleB, vposA, vposB)) {
        return true;
    }

    const Point center = v0->point();
    const Point vposB2 = reflect_through_center(center, vposA);
    if (does_starA_separate_cycles(dt, star, v0, cycleA.facets, cycleB.facets, vposB2)) {
        new_vposB = vposB2;
        return true;
    }

    const Point vposA2 = reflect_through_center(center, vposB);
    if (does_starA_separate_cycles(dt, star, v0, cycleB.facets, cycleA.facets, vposA2)) {
        new_vposA = vposA2;
        return true;
    }

    return false;
}

static uint32_t xorshift32(uint32_t& state) {
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
}

static double rand01(uint32_t& state) {
    return (xorshift32(state) & 0xFFFFFFu) / static_cast<double>(0x1000000u);
}

static Vector3 rotate_axis_angle(const Vector3& vec, const Vector3& axis_unit, double angle_rad) {
    const double c = std::cos(angle_rad);
    const double s = std::sin(angle_rad);

    const Vector3 term0 = vec * c;
    const Vector3 term1 = vec_cross(axis_unit, vec) * s;
    const Vector3 term2 = axis_unit * (vec_dot(axis_unit, vec) * (1.0 - c));
    return term0 + term1 + term2;
}

struct IsovertexCandidateEvaluator {
    Vertex_handle v;
    const Delaunay& dt;
    const std::vector<Cell_handle>& cell_by_index;
    const std::vector<Point>& original_positions;
    double min_sep2 = 0.0;

    bool found_candidate = false;
    double best_cost = std::numeric_limits<double>::infinity();
    std::vector<Point> best_candidate;
    ResolutionStatus best_status = ResolutionStatus::UNRESOLVED;
    ResolutionStrategy best_strategy = ResolutionStrategy::NONE;

    bool accept_if_better(
        const std::vector<Point>& candidate,
        ResolutionStatus status,
        ResolutionStrategy strategy
    ) {
        if (!have_min_separation(candidate, min_sep2)) {
            return false;
        }
        if (check_self_intersection(v, candidate, dt, cell_by_index)) {
            return false;
        }

        const int n = static_cast<int>(candidate.size());
        double cost = 0.0;
        for (int i = 0; i < n; ++i) {
            cost += squared_distance(candidate[i], original_positions[i]);
        }

        if (!found_candidate || cost < best_cost) {
            found_candidate = true;
            best_cost = cost;
            best_candidate = candidate;
            best_status = status;
            best_strategy = strategy;
            return true;
        }

        return false;
    }
};

static void add_two_cycle_diametric_opposite_candidates(
    const Point& center,
    double base_radius,
    const std::vector<Point>& centroids,
    IsovertexCandidateEvaluator& evaluator
) {
    if (centroids.size() != 2) {
        return;
    }

    auto diametric_opposite = [&](const Point& p) -> Point {
        return Point(
            2.0 * center.x() - p.x(),
            2.0 * center.y() - p.y(),
            2.0 * center.z() - p.z());
    };

    // First try reflecting the current positions 
    if (evaluator.original_positions.size() == 2) {
        std::vector<Point> cand = evaluator.original_positions;
        cand[1] = diametric_opposite(evaluator.original_positions[0]);
        evaluator.accept_if_better(
            cand,
            ResolutionStatus::RESOLVED_FALLBACK,
            ResolutionStrategy::TWO_CYCLE_DIAMETRIC);

        cand = evaluator.original_positions;
        cand[0] = diametric_opposite(evaluator.original_positions[1]);
        evaluator.accept_if_better(
            cand,
            ResolutionStatus::RESOLVED_FALLBACK,
            ResolutionStrategy::TWO_CYCLE_DIAMETRIC);
    }

    const double radius_scales[] = {0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0};
    for (double scale : radius_scales) {
        const double r = base_radius * scale;

        const Point vposA = project_to_sphere(centroids[0], center, r);
        std::vector<Point> cand(2);
        cand[0] = vposA;
        cand[1] = diametric_opposite(vposA);
        evaluator.accept_if_better(
            cand,
            ResolutionStatus::RESOLVED_FALLBACK,
            ResolutionStrategy::TWO_CYCLE_DIAMETRIC);

        const Point vposB = project_to_sphere(centroids[1], center, r);
        cand[1] = vposB;
        cand[0] = diametric_opposite(vposB);
        evaluator.accept_if_better(
            cand,
            ResolutionStatus::RESOLVED_FALLBACK,
            ResolutionStrategy::TWO_CYCLE_DIAMETRIC);
    }
}

static void add_fibonacci_fallback_candidates(
    const Point& center,
    double base_radius,
    const std::vector<Vector3>& centroid_dirs,
    IsovertexCandidateEvaluator& evaluator
) {
    const int n = static_cast<int>(centroid_dirs.size());
    if (n < 1) {
        return;
    }

    const double golden_angle = M_PI * (3.0 - std::sqrt(5.0));
    const std::vector<double> fallback_scales = (n == 1)
        ? std::vector<double>{0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0}
        : std::vector<double>{3.0, 4.0, 6.0, 8.0, 10.0};
    const int fallback_attempts = (n == 1) ? 32 : 16;

    for (double s : fallback_scales) {
        const double r = base_radius * s;

        for (int attempt = 0; attempt < fallback_attempts; ++attempt) {
            uint32_t prng = static_cast<uint32_t>(evaluator.v->info().index) * 747796405u +
                            static_cast<uint32_t>(attempt) * 2891336453u + 277803737u;

            const Vector3 axis = vec_normalize_or(
                Vector3(2.0 * rand01(prng) - 1.0, 2.0 * rand01(prng) - 1.0, 2.0 * rand01(prng) - 1.0),
                Vector3(1, 0, 0));
            const double angle = 2.0 * M_PI * rand01(prng);

            // Generate evenly-distributed unit directions.
            std::vector<Vector3> dirs(n);
            for (int i = 0; i < n; ++i) {
                const double t = (static_cast<double>(i) + 0.5) / static_cast<double>(n);
                const double z = 1.0 - 2.0 * t;
                const double xy = std::sqrt(std::max(0.0, 1.0 - z * z));
                const double theta = golden_angle * static_cast<double>(i);
                const Vector3 base_dir(xy * std::cos(theta), xy * std::sin(theta), z);
                dirs[i] = vec_normalize_or(rotate_axis_angle(base_dir, axis, angle), base_dir);
            }

            // Greedy assignment of generated directions to cycles to stay near centroid directions.
            std::vector<int> assign(n, -1);
            std::vector<bool> used(n, false);
            for (int c = 0; c < n; ++c) {
                double best_score = -std::numeric_limits<double>::infinity();
                int best_i = -1;
                for (int i = 0; i < n; ++i) {
                    if (used[i]) continue;
                    const double score = vec_dot(dirs[i], centroid_dirs[c]);
                    if (score > best_score) {
                        best_score = score;
                        best_i = i;
                    }
                }
                if (best_i < 0) {
                    best_i = c;
                }
                assign[c] = best_i;
                used[best_i] = true;
            }

            std::vector<Point> candidate(n);
            for (int c = 0; c < n; ++c) {
                candidate[c] = center + dirs[assign[c]] * r;
            }

            evaluator.accept_if_better(
                candidate,
                ResolutionStatus::RESOLVED_FALLBACK,
                ResolutionStrategy::FIBONACCI_FALLBACK);

            if (n == 2) {
                std::swap(candidate[0], candidate[1]);
                evaluator.accept_if_better(
                    candidate,
                    ResolutionStatus::RESOLVED_FALLBACK,
                    ResolutionStrategy::FIBONACCI_FALLBACK);
            }
        }
    }
}

//! @brief Attempt to resolve self-intersection for a multi-cycle vertex
/*!
 * Tries multiple strategies to eliminate self-intersections between cycles:
 * 1. Increase sphere radius (tries 1x, 2x, 3x of base radius)
 * 2. uniform spherical distribution with different rotations
 *
 * @param v The vertex with potential self-intersection
 * @param cycle_isovertices The isovertex positions to modify
 * @param dt The Delaunay triangulation
 * @param grid The volumetric grid
 * @param isovalue The isovalue for centroid computation
 * @param cell_by_index Cell handles indexed by `cell->info().index`
 * @return Resolution status indicating what happened
 */
ResolutionResult resolve_self_intersection(
    Vertex_handle v,
    std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    const UnifiedGrid& grid,
    float isovalue,
    const std::vector<Cell_handle>& cell_by_index
) {
    if (!check_self_intersection(v, cycle_isovertices, dt, cell_by_index)) {
        return {ResolutionStatus::NOT_NEEDED, ResolutionStrategy::NONE};
    }

    const auto& cycles = v->info().facet_cycles;
    const Point center = v->point();
    const int n = static_cast<int>(cycle_isovertices.size());

    if (n < 1) {
        return {ResolutionStatus::NOT_NEEDED, ResolutionStrategy::NONE};
    }

    const std::vector<Point> original_positions = cycle_isovertices;
    const double base_radius = compute_sphere_radius(v, dt, grid);

    // Precompute centroids per-cycle once.
    std::vector<Point> centroids(n);
    for (int c = 0; c < n; ++c) {
        centroids[c] = compute_centroid_of_voronoi_edge_and_isosurface(
            dt, grid, isovalue, v, cycles[c], cell_by_index);
    }

    std::vector<Vector3> centroid_dirs(n);
    for (int i = 0; i < n; ++i) {
        centroid_dirs[i] = unit_dir_to_or(center, centroids[i], Vector3(1, 0, 0));
    }

    const double min_sep = std::max(1e-8, base_radius * 1e-3);
    const double min_sep2 = min_sep * min_sep;

    IsovertexCandidateEvaluator evaluator{v, dt, cell_by_index, original_positions, min_sep2};

    if (n == 2) {
        add_two_cycle_diametric_opposite_candidates(center, base_radius, centroids, evaluator);
    }

    if (!evaluator.found_candidate) {
        add_fibonacci_fallback_candidates(center, base_radius, centroid_dirs, evaluator);
    }

    if (evaluator.found_candidate) {
        cycle_isovertices = evaluator.best_candidate;
        return {evaluator.best_status, evaluator.best_strategy};
    }

    cycle_isovertices = original_positions;

    DEBUG_PRINT("[DEBUG-UNRESOLVED] Vertex " << v->info().index
                << " at (" << center.x() << ", " << center.y() << ", " << center.z() << ")"
                << " has " << n << " cycles, could not resolve self-intersection");
    for (int i = 0; i < n; ++i) {
        DEBUG_PRINT("[DEBUG-UNRESOLVED]   Isovertex " << i << ": ("
                    << cycle_isovertices[i].x() << ", "
                    << cycle_isovertices[i].y() << ", "
                    << cycle_isovertices[i].z() << ")");
        DEBUG_PRINT("[DEBUG-UNRESOLVED]     Cycle " << i << " contains "
                    << cycles[i].size() << " facets:");
        for (const auto& fkey : cycles[i]) {
            DEBUG_PRINT("[DEBUG-UNRESOLVED]       Facet (cell=" << fkey.first
                        << ", local_idx=" << fkey.second << ")");
        }
    }

    return {ResolutionStatus::UNRESOLVED, ResolutionStrategy::NONE};
}

// ============================================================================
// Multi-cycle resolution
// ============================================================================

static ResolutionResult resolve_multi_cycle_position_multi_isov(
    Vertex_handle v,
    std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    const UnifiedGrid& grid,
    float isovalue,
    const std::vector<Cell_handle>& cell_by_index
) {
    const auto& cycles = v->info().facet_cycles;
    const int num_cycles = static_cast<int>(cycles.size());
    if (num_cycles < 2 || static_cast<int>(cycle_isovertices.size()) < 2) {
        return {ResolutionStatus::NOT_NEEDED, ResolutionStrategy::NONE};
    }

    // attempt only when a local self-intersection is present.
    const std::vector<Point> original_positions = cycle_isovertices;
    if (!check_self_intersection(v, original_positions, dt, cell_by_index)) {
        return {ResolutionStatus::NOT_NEEDED, ResolutionStrategy::NONE};
    }

    std::vector<CycleDataPositionMultiIsov> cycle_data;
    cycle_data.reserve(static_cast<size_t>(num_cycles));
    for (int c = 0; c < num_cycles; ++c) {
        cycle_data.push_back(build_cycle_data_position_multi_isov(dt, cell_by_index, v, cycles[c]));
    }

    const StarCellSetPositionMultiIsov star = build_star_cell_set_position_multi_isov(dt, v);

    bool any_changed = false;
    const double change_eps2 = 1e-6;
    const int MAX_PAIRWISE_PASSES = 6;

    for (int pass = 0; pass < MAX_PAIRWISE_PASSES; ++pass) {
        bool pass_changed = false;

        for (int a = 0; a < num_cycles; ++a) {
            for (int b = a + 1; b < num_cycles; ++b) {
                Point newA = cycle_isovertices[static_cast<size_t>(a)];
                Point newB = cycle_isovertices[static_cast<size_t>(b)];

                if (!find_vertex_positions_separating_cycles(
                        dt,
                        star,
                        v,
                        cycle_data[static_cast<size_t>(a)],
                        cycle_data[static_cast<size_t>(b)],
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
        if (!any_changed) {
            return {ResolutionStatus::NOT_NEEDED, ResolutionStrategy::NONE};
        }
        return {ResolutionStatus::RESOLVED_FALLBACK, ResolutionStrategy::TWO_CYCLE_DIAMETRIC};
    }

    // Fallback: retain the previous heuristic resolver when PositionMultiIsov
    // separation tests are insufficient (e.g., >2 cycles or degenerate cases).
    cycle_isovertices = original_positions;
    return resolve_self_intersection(v, cycle_isovertices, dt, grid, isovalue, cell_by_index);
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
 * @param result The resolution result from resolve_self_intersection.
 * @param stats [in/out] Statistics to update.
 * @param any_changed [out] Set to true if resolution changed isovertices.
 * @param pass_detected [out] Incremented if self-intersection was detected.
 * @param pass_resolved [out] Incremented if resolved via sphere projection.
 * @param pass_fallback [out] Incremented if resolved via fallback strategy.
 * @param pass_unresolved [out] Incremented if resolution failed.
 */
static void update_resolution_stats(
    const ResolutionResult& result,
    IsovertexComputationStats& stats,
    bool& any_changed,
    int& pass_detected,
    int& pass_resolved,
    int& pass_fallback,
    int& pass_unresolved
) {
    const ResolutionStatus status = result.status;
    
    if (status != ResolutionStatus::NOT_NEEDED) {
        pass_detected++;
    }
    
    switch (status) {
        case ResolutionStatus::NOT_NEEDED:
            break;
        case ResolutionStatus::RESOLVED_SPHERE:
            any_changed = true;
            pass_resolved++;
            switch (result.strategy) {
                case ResolutionStrategy::TWO_CYCLE_DIAMETRIC:
                    stats.strat_two_cycle_diametric++;
                    break;
                case ResolutionStrategy::FIBONACCI_FALLBACK:
                    stats.strat_fibonacci_fallback++;
                    break;
                default:
                    break;
            }
            break;
        case ResolutionStatus::RESOLVED_FALLBACK:
            any_changed = true;
            pass_fallback++;
            switch (result.strategy) {
                case ResolutionStrategy::TWO_CYCLE_DIAMETRIC:
                    stats.strat_two_cycle_diametric++;
                    break;
                case ResolutionStrategy::FIBONACCI_FALLBACK:
                    stats.strat_fibonacci_fallback++;
                    break;
                default:
                    break;
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
// Sub-routine: Compute initial cycle isovertices (Pass 1)
// ============================================================================

std::vector<Vertex_handle> compute_initial_cycle_isovertices(
    Delaunay& dt,
    const UnifiedGrid& grid,
    float isovalue,
    const std::vector<Cell_handle>& cell_by_index,
    IsovertexComputationStats& stats
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
            // The isosurface sample point is stored in vit->info().isov.
            // When -position_delv_on_isov is enabled, this equals cube_center.
            isovertices[0] = vit->info().isov;
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
            const double sphere_radius = compute_sphere_radius(vit, dt, grid);
            for (size_t c = 0; c < cycles.size(); ++c) {
                Point centroid = compute_centroid_of_voronoi_edge_and_isosurface(
                    dt, grid, isovalue, vit, cycles[c], cell_by_index);

                isovertices[c] = project_to_sphere(centroid, cube_center, sphere_radius);
            }
            multi_cycle_vertices.push_back(vit);
            stats.multi_cycle_count++;
        }
    }
    
    return multi_cycle_vertices;
}

// ============================================================================
// Sub-routine: Resolve self-intersections on multi-cycle vertices (Pass 2)
// ============================================================================

std::vector<Vertex_handle> resolve_multi_cycle_self_intersections(
    Delaunay& dt,
    const UnifiedGrid& grid,
    float isovalue,
    const std::vector<Cell_handle>& cell_by_index,
    const std::vector<Vertex_handle>& multi_cycle_vertices,
    IsovertexComputationStats& stats
) {
    std::vector<Vertex_handle> modified_multi_cycle_vertices;
    std::unordered_set<int> modified_multi_cycle_indices;

    // Run multiple stabilization passes to handle cascading conflicts.
    // Changing one multi-cycle vertex may introduce a new conflict at another.
    const int MAX_MULTI_CYCLE_PASSES = 3;
    
    for (int pass = 0; pass < MAX_MULTI_CYCLE_PASSES; ++pass) {
        bool any_changed = false;
        int pass_detected = 0;
        int pass_resolved = 0;
        int pass_fallback = 0;
        int pass_unresolved = 0;

        // Process each multi-cycle vertex for self-intersection.
        for (Vertex_handle vh : multi_cycle_vertices) {
            if (!vh->info().active) continue;
            if (vh->info().is_dummy) continue;
            if (vh->info().facet_cycles.size() < 2) continue;

            auto& isovertices = vh->info().cycle_isovertices;

            // First try the PositionMultiIsov orientation/half-space separation logic.
            // Fall back to the legacy intersection-based resolver only if needed.
            const ResolutionResult result = resolve_multi_cycle_position_multi_isov(
                vh, isovertices, dt, grid, isovalue, cell_by_index);

            update_resolution_stats(result, stats, any_changed,
                pass_detected, pass_resolved, pass_fallback, pass_unresolved);

            // Track which vertices were modified for targeted cleanup passes.
            const ResolutionStatus status = result.status;
            if (status == ResolutionStatus::RESOLVED_SPHERE ||
                status == ResolutionStatus::RESOLVED_FALLBACK ||
                status == ResolutionStatus::UNRESOLVED) {
                if (modified_multi_cycle_indices.insert(vh->info().index).second) {
                    modified_multi_cycle_vertices.push_back(vh);
                }
            }
        }

        // Accumulate pass statistics into overall stats.
        stats.self_intersection_detected += pass_detected;
        stats.self_intersection_resolved += pass_resolved;
        stats.self_intersection_fallback += pass_fallback;

        if (pass_detected > 0) {
            DEBUG_PRINT("[DEL-SELFI] multi-cycle pass " << (pass + 1) << "/" << MAX_MULTI_CYCLE_PASSES
                        << ": detected=" << pass_detected
                        << ", resolved=" << pass_resolved
                        << ", fallback=" << pass_fallback
                        << ", unresolved=" << pass_unresolved);
        }

        // Early exit if no changes were made this pass.
        if (!any_changed) {
            break;
        }
    }
    
    return modified_multi_cycle_vertices;
}

// ============================================================================
// Sub-routine: Resolve self-intersections on nearby single-cycle vertices (Pass 3)
// ============================================================================

void resolve_single_cycle_self_intersections(
    Delaunay& dt,
    const UnifiedGrid& grid,
    float isovalue,
    const std::vector<Cell_handle>& cell_by_index,
    const std::vector<Vertex_handle>& seed_vertices,
    IsovertexComputationStats& stats
) {
    const int max_hops = 2;  // Hard-coded hop radius for single-cycle cleanup
    if (seed_vertices.empty()) {
        return;
    }

    std::vector<Vertex_handle> single_cycle_candidates;
    std::vector<uint8_t> is_candidate(static_cast<size_t>(dt.number_of_vertices() + 1), 0);
    std::vector<uint8_t> visited(static_cast<size_t>(dt.number_of_vertices() + 1), 0);
    std::vector<Edge> incident_edges;
    incident_edges.reserve(32);

    auto ensure_index_capacity = [&](int idx) {
        if (idx < 0) {
            return;
        }
        const size_t needed = static_cast<size_t>(idx) + 1;
        if (needed > visited.size()) {
            visited.resize(needed, 0);
            is_candidate.resize(needed, 0);
        }
    };

    // Lambda to gather single-cycle vertices within a certain hop distance.
    auto add_single_cycle_candidates_within_hops = [&](
        const std::vector<Vertex_handle>& seeds,
        int max_hops
    ) {
        std::fill(visited.begin(), visited.end(), 0);

        std::vector<Vertex_handle> frontier;
        frontier.reserve(seeds.size());
        for (Vertex_handle vh : seeds) {
            if (vh == Vertex_handle()) {
                continue;
            }
            if (!vh->info().active || vh->info().is_dummy) {
                continue;
            }
            const int seed_idx = vh->info().index;
            if (seed_idx < 0) {
                continue;
            }
            ensure_index_capacity(seed_idx);
            if (visited[static_cast<size_t>(seed_idx)] == 0) {
                visited[static_cast<size_t>(seed_idx)] = 1;
                frontier.push_back(vh);
            }
        }

        for (int depth = 0; depth < max_hops; ++depth) {
            std::vector<Vertex_handle> next_frontier;
            next_frontier.reserve(frontier.size() * 4);

            for (Vertex_handle src : frontier) {
                incident_edges.clear();
                dt.incident_edges(src, std::back_inserter(incident_edges));

                for (const Edge& e : incident_edges) {
                    Vertex_handle a = e.first->vertex(e.second);
                    Vertex_handle b = e.first->vertex(e.third);
                    Vertex_handle other = (a == src) ? b : a;

                    if (dt.is_infinite(other)) continue;
                    if (other->info().is_dummy) continue;
                    if (!other->info().active) continue;

                    const int other_idx = other->info().index;
                    if (other_idx < 0) {
                        continue;
                    }
                    ensure_index_capacity(other_idx);
                    if (visited[static_cast<size_t>(other_idx)] == 0) {
                        visited[static_cast<size_t>(other_idx)] = 1;
                        next_frontier.push_back(other);
                    }

                    // Collect single-cycle candidates.
                    if (other->info().facet_cycles.size() == 1 &&
                        other->info().facet_cycles[0].size() >= 4) {
                        if (static_cast<size_t>(other_idx) < is_candidate.size() &&
                            is_candidate[static_cast<size_t>(other_idx)] == 0) {
                            is_candidate[static_cast<size_t>(other_idx)] = 1;
                            single_cycle_candidates.push_back(other);
                        }
                    }
                }
            }

            frontier.swap(next_frontier);
            if (frontier.empty()) {
                break;
            }
        }
    };

    // Targeted cleanup: check single-cycle vertices within a small hop radius of seed vertices.
    add_single_cycle_candidates_within_hops(seed_vertices, max_hops);

    const int MAX_SINGLE_CYCLE_PASSES = 2;
    const int MAX_EXPANSION_ROUNDS = 2;
    
    for (int round = 0; round < MAX_EXPANSION_ROUNDS; ++round) {
        std::vector<Vertex_handle> unresolved_vertices;
        unresolved_vertices.reserve(32);

        for (int pass = 0; pass < MAX_SINGLE_CYCLE_PASSES; ++pass) {
            bool any_changed = false;
            int pass_detected = 0;
            int pass_resolved = 0;
            int pass_fallback = 0;
            int pass_unresolved = 0;

            unresolved_vertices.clear();

            for (Vertex_handle vh : single_cycle_candidates) {
                if (!vh->info().active) continue;
                if (vh->info().is_dummy) continue;
                if (vh->info().facet_cycles.size() != 1) continue;

                auto& isovertices = vh->info().cycle_isovertices;
                const ResolutionResult result = resolve_self_intersection(
                    vh, isovertices, dt, grid, isovalue, cell_by_index);

                update_resolution_stats(result, stats, any_changed,
                    pass_detected, pass_resolved, pass_fallback, pass_unresolved);

                if (result.status == ResolutionStatus::UNRESOLVED) {
                    unresolved_vertices.push_back(vh);
                }
            }

            stats.self_intersection_detected += pass_detected;
            stats.self_intersection_resolved += pass_resolved;
            stats.self_intersection_fallback += pass_fallback;

            if (pass_detected > 0) {
                DEBUG_PRINT("[DEL-SELFI] single-cycle pass " << (pass + 1) << "/" << MAX_SINGLE_CYCLE_PASSES
                            << ": detected=" << pass_detected
                            << ", resolved=" << pass_resolved
                            << ", fallback=" << pass_fallback
                            << ", unresolved=" << pass_unresolved);
            }

            if (!any_changed) {
                break;
            }
        }

        if (unresolved_vertices.empty()) {
            break;
        }

        // Expand the candidate region around unresolved single-cycle vertices.
        const size_t before = single_cycle_candidates.size();
        add_single_cycle_candidates_within_hops(unresolved_vertices, max_hops);
        if (single_cycle_candidates.size() == before) {
            break;
        }
    }
}


// ============================================================================
// Sub-routine: Report isovertex computation statistics
// ============================================================================

void report_isovertex_statistics(const IsovertexComputationStats& stats) {
    DEBUG_PRINT("[DEL-ISOV] Computed isovertices for "
                << stats.single_cycle_count << " single-cycle and "
                << stats.multi_cycle_count << " multi-cycle vertices");

    if (stats.self_intersection_detected > 0) {
        DEBUG_PRINT("[DEL-SELFI] Self-intersection stats: "
                    << "detected=" << stats.self_intersection_detected
                    << ", resolved_sphere=" << stats.self_intersection_resolved
                    << ", resolved_fallback=" << stats.self_intersection_fallback
                    << ", unresolved=" << stats.self_intersection_unresolved);

        const int64_t resolved_total = stats.strat_two_cycle_diametric +
                                       stats.strat_fibonacci_fallback;
        if (resolved_total > 0) {
            auto pct = [&](int64_t count) -> int {
                return static_cast<int>((count * 100 + resolved_total / 2) / resolved_total);
            };

            DEBUG_PRINT("[DEL-SELFI] Resolution strategy breakdown:");
            DEBUG_PRINT("[DEL-SELFI]   two_cycle_diametric: " << stats.strat_two_cycle_diametric
                        << " (" << pct(stats.strat_two_cycle_diametric) << "%)");
            DEBUG_PRINT("[DEL-SELFI]   fibonacci_fallback: " << stats.strat_fibonacci_fallback
                        << " (" << pct(stats.strat_fibonacci_fallback) << "%)");
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
 * The computation proceeds in three passes:
 *
 *   Pass 1 (compute_initial_cycle_isovertices):
 *     Compute initial positions for all active vertices.
 *     - Single-cycle vertices: use the isosurface sample point
 *     - Multi-cycle vertices: compute Voronoi/isosurface centroids and project to sphere
 *
 *   Pass 2 (resolve_multi_cycle_self_intersections):
 *     Iteratively resolve self-intersections on multi-cycle vertices.
 *     Uses multiple strategies: centroid projection, diametric placement, Fibonacci fallback.
 *
 *   Pass 3 (resolve_single_cycle_self_intersections):
 *     Targeted cleanup of single-cycle vertices near modified vertices.
 *     Checks vertices within 2 hops of modified multi-cycle vertices.
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
    bool position_on_isov
) {
    DEBUG_PRINT("[DEL-ISOV] Computing isovertex positions for cycles...");
    
    // Suppress unused parameter warning - the choice is encoded in vertex info.
    (void)position_on_isov;

    // Build a direct lookup table for efficient cell access by index.
    const std::vector<Cell_handle> cell_by_index = build_cell_index_vector(dt);

    // Initialize statistics tracking.
    IsovertexComputationStats stats;

    // ========================================================================
    // Pass 1: Compute initial isovertex positions for all active vertices.
    // ========================================================================
    std::vector<Vertex_handle> multi_cycle_vertices = compute_initial_cycle_isovertices(
        dt, grid, isovalue, cell_by_index, stats);

    // ========================================================================
    // Pass 2: Resolve self-intersections on multi-cycle vertices.
    // ========================================================================
    std::vector<Vertex_handle> modified_vertices = resolve_multi_cycle_self_intersections(
        dt, grid, isovalue, cell_by_index, multi_cycle_vertices, stats);

    // ========================================================================
    // Pass 3: Targeted cleanup of single-cycle vertices near modified multi-cycle vertices that might be affected.
    // ========================================================================
    resolve_single_cycle_self_intersections(
        dt, grid, isovalue, cell_by_index, modified_vertices, stats);

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

    if (result.total_flips > 0) {
        std::cout << "  ModCyc: " << result.total_flips << " edge flips ("
                  << result.iterations << " iterations)" << std::endl;
    }

    DEBUG_PRINT("[CYC-MOD] mod_cyc: flips=" << result.total_flips
                << ", iterations=" << result.iterations);

    return result;
}
