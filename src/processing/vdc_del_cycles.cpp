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
#include <CGAL/Triangle_3.h>
#include <CGAL/intersections.h>

// ============================================================================
// Helper: Build cell index lookup
// ============================================================================

std::unordered_map<int, Cell_handle> build_cell_index_map(const Delaunay& dt) {
    std::unordered_map<int, Cell_handle> map;
    for (auto cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit) {
        map[cit->info().index] = cit;
    }
    return map;
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
    const std::unordered_map<int, Cell_handle>& cell_map
) {
    auto it = cell_map.find(key.cell_index);
    if (it == cell_map.end()) {
        return false;
    }

    Cell_handle cell = it->second;
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
    const std::unordered_map<int, Cell_handle>& cell_map
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
        if (!facet_is_valid_for_cycle_matching(dt, key, cell_map)) {
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


    // Build a map for quick cell lookup by index
    std::unordered_map<int, Cell_handle> cell_map = build_cell_index_map(dt);

    // =========================================================================
    // Step 1: Build global matching data for all facets
    // =========================================================================
    
    std::unordered_map<FacetKey, FacetMatchingData, FacetKeyHash> facet_matching;
    
    // Track which edges we've already processed (keyed by min,max vertex indices)
    std::unordered_set<EdgeKey, EdgeKeyHash> processed_edges;
    
    // For each active Delaunay vertex v0, match facets on its incident edges
    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
        if (!vit->info().active) continue;
        if (vit->info().is_dummy) continue;
        
        const Vertex_handle v0 = vit;
        const int v0_idx = v0->info().index;
        
        // Get all incident edges
        std::vector<Edge> incident_edges;
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
            auto pairs = compute_edge_bipolar_matching(dt, e, BIPOLAR_MATCH_METHOD::SEP_NEG, cell_map);
            if (pairs.empty()) continue;
            
            // For each matched pair (f0, f1), store the matching per vertex slot
            // Per algorithm: f0.cycleMatchingFacet[j0] = f1, f1.cycleMatchingFacet[j1] = f0
            for (const auto& [fkey0, fkey1] : pairs) {
                auto cell0_it = cell_map.find(fkey0.cell_index);
                auto cell1_it = cell_map.find(fkey1.cell_index);
                if (cell0_it == cell_map.end() || cell1_it == cell_map.end()) continue;
                
                Cell_handle cell0 = cell0_it->second;
                Cell_handle cell1 = cell1_it->second;
                
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
    
    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
        if (!vit->info().active) continue;
        if (vit->info().is_dummy) continue;
        const Vertex_handle v_handle = vit;
        const int v_idx = v_handle->info().index;

        // Clear any existing cycles
        vit->info().facet_cycles.clear();

        // Collect all isosurface facets incident on this vertex
        std::vector<FacetKey> isosurface_facets;
        std::unordered_set<FacetKey, FacetKeyHash> facet_set;
        
        std::vector<Cell_handle> incident_cells;
        dt.incident_cells(v_handle, std::back_inserter(incident_cells));

        for (Cell_handle ch : incident_cells) {
            if (dt.is_infinite(ch)) continue;

            for (int i = 0; i < 4; ++i) {
                if (!ch->info().facet_is_isosurface[i]) continue;

                const int v_local = find_vertex_index_in_cell(ch, v_handle);
                if (v_local < 0 || v_local == i) continue;
                
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

                FacetKey fkey{ch->info().index, i};
                if (facet_set.insert(fkey).second) {
                    isosurface_facets.push_back(fkey);
                }
            }
        }


        // Map facet key to vertex index within that facet (for this vertex v)
        std::unordered_map<FacetKey, int, FacetKeyHash> facet_to_v_index;
        for (const FacetKey& fkey : isosurface_facets) {
            auto it = cell_map.find(fkey.cell_index);
            if (it == cell_map.end()) continue;
            Cell_handle cell = it->second;
            int j = get_vertex_index_in_facet(cell, fkey.facet_index, v_handle);
            if (j >= 0) {
                facet_to_v_index[fkey] = j;
            }
        }

        // Traverse matchings to form cycles
        std::unordered_set<FacetKey, FacetKeyHash> visited;
        std::vector<std::vector<std::pair<int, int>>> cycles;

        for (const FacetKey& start_fkey : isosurface_facets) {
            if (visited.count(start_fkey)) continue;

            std::vector<std::pair<int, int>> cycle;
            FacetKey current_fkey = start_fkey;

            while (true) {
                if (visited.count(current_fkey)) {
                    break;  // Completed cycle or hit already-visited facet
                }
                visited.insert(current_fkey);
                cycle.push_back({current_fkey.cell_index, current_fkey.facet_index});

                // Find the next facet in this cycle using cycleMatchingFacet[j]
                auto v_idx_it = facet_to_v_index.find(current_fkey);
                if (v_idx_it == facet_to_v_index.end()) {
                    break;  // No vertex index found
                }
                int j = v_idx_it->second;

                auto match_it = facet_matching.find(current_fkey);
                if (match_it == facet_matching.end() || !match_it->second.hasMatch[j]) {
                    break;  // No matching for this vertex slot
                }

                FacetKey next_fkey = match_it->second.cycleMatchingFacet[j];
                
                // Verify next facet is incident on v (should be in our facet list)
                if (facet_set.count(next_fkey) == 0) {
                    break;  // Next facet is not incident on this vertex
                }

                current_fkey = next_fkey;
            }

            if (!cycle.empty()) {
                cycles.push_back(cycle);
            }
        }

        // Store cycles in vertex info
        vit->info().facet_cycles = cycles;


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
    const std::unordered_map<int, Cell_handle>& cell_map
) {

    double cx = 0, cy = 0, cz = 0;
    int numI = 0;

    for (const auto& [cell_idx, facet_idx] : cycle_facets) {
        // Get the cell and its neighbor
        auto it = cell_map.find(cell_idx);
        if (it == cell_map.end()) continue;

        Cell_handle cellA = it->second;
        Cell_handle cellB = cellA->neighbor(facet_idx);

        if (dt.is_infinite(cellB)) continue;

        // Get circumcenters (dual Voronoi vertices)
        Point wA = cellA->info().circumcenter;
        Point wB = cellB->info().circumcenter;

        // Get scalar values at circumcenters
        float sA = cellA->info().circumcenter_scalar;
        float sB = cellB->info().circumcenter_scalar;

        // Avoid division by zero
        float denom = sA - sB;
        if (std::abs(denom) < 1e-10f) continue;

        // Linear interpolation to find isosurface crossing point
        // p = ((isovalue - sB) * wA + (sA - isovalue) * wB) / (sA - sB)
        //
        // Using p = wA + t * (wB - wA) gives t = (sA - isovalue) / (sA - sB).
        float t = (sA - isovalue) / denom;

        // Clamp t to [0, 1] for safety
        t = std::max(0.0f, std::min(1.0f, t));

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
        // Fallback to vertex position if no valid intersections
        return v_handle->point();
    }

    return Point(cx / numI, cy / numI, cz / numI);
}

Point compute_centroid_of_voronoi_edge_and_isosurface(
    const Delaunay& dt,
    const UnifiedGrid& grid,
    float isovalue,
    Vertex_handle v_handle,
    const std::vector<std::pair<int, int>>& cycle_facets
) {
    const std::unordered_map<int, Cell_handle> cell_map = build_cell_index_map(dt);
    return compute_centroid_of_voronoi_edge_and_isosurface(
        dt, grid, isovalue, v_handle, cycle_facets, cell_map);
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
 * @param cell_map Map from cell indices to cell handles
 * @return Vector of Triangle_3 objects for this cycle
 */
static std::vector<Triangle_3> collect_cycle_triangles(
    const Delaunay& dt,
    Vertex_handle v,
    int cycle_idx,
    const Point& cycle_isovertex,
    const std::unordered_map<int, Cell_handle>& cell_map
) {
    std::vector<Triangle_3> triangles;
    const auto& cycles = v->info().facet_cycles;

    if (cycle_idx < 0 || cycle_idx >= static_cast<int>(cycles.size())) {
        return triangles;
    }

    const auto& cycle_facets = cycles[cycle_idx];

    for (const auto& [cell_idx, facet_idx] : cycle_facets) {
        auto it = cell_map.find(cell_idx);
        if (it == cell_map.end()) continue;

        Cell_handle ch = it->second;

        // Get the 3 vertices of this facet (excluding the one opposite to facet_idx)
        std::vector<Vertex_handle> facet_verts;
        for (int j = 0; j < 4; ++j) {
            if (j != facet_idx) {
                facet_verts.push_back(ch->vertex(j));
            }
        }

        // Find the two vertices that are not v
        std::vector<Vertex_handle> other_verts;
        for (auto vh : facet_verts) {
            if (vh != v) {
                other_verts.push_back(vh);
            }
        }

        if (other_verts.size() != 2) continue;

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
        int c1 = find_cycle_containing_facet(v1, cell_idx, facet_idx);
        if (c1 >= 0 && c1 < static_cast<int>(v1->info().cycle_isovertices.size())) {
            p2 = v1->info().cycle_isovertices[c1];
            found_p2 = true;
        }

        // Find cycle for v2 containing this facet
        int c2 = find_cycle_containing_facet(v2, cell_idx, facet_idx);
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
    std::vector<Point> shared;
    shared.reserve(3);
    for (const Point& a : t1v) {
        for (const Point& b : t2v) {
            if (a == b) {
                shared.push_back(a);
                break;
            }
        }
    }

    // If the triangles are disjoint in vertex set, any intersection is non-trivial.
    if (shared.empty()) {
        return CGAL::do_intersect(t1, t2);
    }

    // Shared edge/triangle: adjacent triangles are expected to intersect along their shared simplex.
    // Overlap beyond that would indicate degeneracy, which we ignore here for performance.
    if (shared.size() >= 2) {
        return false;
    }

    // Shared single vertex: triangles always intersect at that vertex, so avoid triangle/triangle
    // intersection tests. A non-trivial intersection must involve the opposite edge of one triangle
    // intersecting the other triangle away from the shared vertex.
    const Point& s = shared[0];

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
    const std::unordered_map<int, Cell_handle>& cell_map
) {
    const auto& cycles = v->info().facet_cycles;
    size_t num_cycles = cycles.size();

    if (num_cycles < 1) {
        return false;
    }

    std::vector<Point> original_isovertices = v->info().cycle_isovertices;

    auto& mutable_isovertices = const_cast<std::vector<Point>&>(v->info().cycle_isovertices);
    mutable_isovertices = cycle_isovertices;

    // Collect triangles for each cycle
    std::vector<std::vector<Triangle_3>> cycle_triangles(num_cycles);
    for (size_t c = 0; c < num_cycles; ++c) {
        cycle_triangles[c] = collect_cycle_triangles(
            dt, v, static_cast<int>(c), cycle_isovertices[c], cell_map);
    }

    // Restore original isovertices
    mutable_isovertices = original_isovertices;

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
// Note: ResolutionStatus, ResolutionStrategy, ResolutionResult are now declared in vdc_del_cycles.h

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
    const std::unordered_map<int, Cell_handle>& cell_map;
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
        if (check_self_intersection(v, candidate, dt, cell_map)) {
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

static void add_centroid_projection_candidates(
    const Point& center,
    double base_radius,
    const std::vector<Point>& centroids,
    IsovertexCandidateEvaluator& evaluator
) {
    const double radius_scales[] = {0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0};
    const int n = static_cast<int>(centroids.size());
    for (double scale : radius_scales) {
        const double r = base_radius * scale;
        std::vector<Point> candidate(n);
        for (int c = 0; c < n; ++c) {
            candidate[c] = project_to_sphere(centroids[c], center, r);
        }
        evaluator.accept_if_better(
            candidate,
            ResolutionStatus::RESOLVED_SPHERE,
            ResolutionStrategy::CENTROID_PROJECTION);
    }
}

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
 * @param cell_map Map from cell indices to cell handles
 * @return Resolution status indicating what happened
 */
ResolutionResult resolve_self_intersection(
    Vertex_handle v,
    std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    const UnifiedGrid& grid,
    float isovalue,
    const std::unordered_map<int, Cell_handle>& cell_map
) {
    if (!check_self_intersection(v, cycle_isovertices, dt, cell_map)) {
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
            dt, grid, isovalue, v, cycles[c], cell_map);
    }

    std::vector<Vector3> centroid_dirs(n);
    for (int i = 0; i < n; ++i) {
        centroid_dirs[i] = unit_dir_to_or(center, centroids[i], Vector3(1, 0, 0));
    }

    const double min_sep = std::max(1e-8, base_radius * 1e-3);
    const double min_sep2 = min_sep * min_sep;

    // Two-cycle special case 
    // keep one isovertex on the current sphere and place the other diametrically opposite.
    if (n == 2) {
        auto diametric_opposite = [&](const Point& p) -> Point {
            return Point(
                2.0 * center.x() - p.x(),
                2.0 * center.y() - p.y(),
                2.0 * center.z() - p.z());
        };

        // vposB2 = 2*(v0 - vposA) + vposA = 2*v0 - vposA
        {
            std::vector<Point> candidate = original_positions;
            candidate[1] = diametric_opposite(original_positions[0]);
            if (have_min_separation(candidate, min_sep2) &&
                !check_self_intersection(v, candidate, dt, cell_map)) {
                cycle_isovertices = candidate;
                return {ResolutionStatus::RESOLVED_FALLBACK, ResolutionStrategy::TWO_CYCLE_DIAMETRIC};
            }
        }

        // vposA2 = 2*(v0 - vposB) + vposB = 2*v0 - vposB
        {
            std::vector<Point> candidate = original_positions;
            candidate[0] = diametric_opposite(original_positions[1]);
            if (have_min_separation(candidate, min_sep2) &&
                !check_self_intersection(v, candidate, dt, cell_map)) {
                cycle_isovertices = candidate;
                return {ResolutionStatus::RESOLVED_FALLBACK, ResolutionStrategy::TWO_CYCLE_DIAMETRIC};
            }
        }
    }

    IsovertexCandidateEvaluator evaluator{v, dt, cell_map, original_positions, min_sep2};

    add_centroid_projection_candidates(center, base_radius, centroids, evaluator);

    // If simple centroid-based placements resolve the issue, avoid expensive fallbacks.
    if (evaluator.found_candidate) {
        cycle_isovertices = evaluator.best_candidate;
        return {evaluator.best_status, evaluator.best_strategy};
    }

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
                case ResolutionStrategy::CENTROID_PROJECTION:
                    stats.strat_centroid_projection++;
                    break;
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
                case ResolutionStrategy::CENTROID_PROJECTION:
                    stats.strat_centroid_projection++;
                    break;
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
    const std::unordered_map<int, Cell_handle>& cell_map,
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
                    dt, grid, isovalue, vit, cycles[c], cell_map);

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
    const std::unordered_map<int, Cell_handle>& cell_map,
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

            // Attempt to resolve self-intersection using multiple strategies.
            const ResolutionResult result = resolve_self_intersection(
                vh, isovertices, dt, grid, isovalue, cell_map);

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
    const std::unordered_map<int, Cell_handle>& cell_map,
    const std::vector<Vertex_handle>& modified_multi_cycle_vertices,
    IsovertexComputationStats& stats
) {
    // Skip if no multi-cycle vertices were modified.
    if (modified_multi_cycle_vertices.empty()) {
        return;
    }

    std::vector<Vertex_handle> single_cycle_candidates;
    std::unordered_set<int> single_cycle_candidate_indices;

    // Lambda to gather single-cycle vertices within a certain hop distance.
    auto add_single_cycle_candidates_within_hops = [&](
        const std::vector<Vertex_handle>& seeds,
        int max_hops
    ) {
        std::unordered_set<int> visited_indices;
        visited_indices.reserve(seeds.size() * 8);

        std::vector<Vertex_handle> frontier = seeds;
        for (Vertex_handle vh : frontier) {
            visited_indices.insert(vh->info().index);
        }

        for (int depth = 0; depth < max_hops; ++depth) {
            std::vector<Vertex_handle> next_frontier;
            next_frontier.reserve(frontier.size() * 4);

            for (Vertex_handle src : frontier) {
                std::vector<Edge> incident_edges;
                dt.incident_edges(src, std::back_inserter(incident_edges));

                for (const Edge& e : incident_edges) {
                    Vertex_handle a = e.first->vertex(e.second);
                    Vertex_handle b = e.first->vertex(e.third);
                    Vertex_handle other = (a == src) ? b : a;

                    if (dt.is_infinite(other)) continue;
                    if (other->info().is_dummy) continue;
                    if (!other->info().active) continue;

                    const int other_idx = other->info().index;
                    if (visited_indices.insert(other_idx).second) {
                        next_frontier.push_back(other);
                    }

                    // Collect single-cycle candidates.
                    if (other->info().facet_cycles.size() == 1) {
                        if (single_cycle_candidate_indices.insert(other_idx).second) {
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

    // Targeted cleanup: check single-cycle vertices within 2 hops of modified multi-cycle vertices.
    add_single_cycle_candidates_within_hops(modified_multi_cycle_vertices, /*max_hops=*/2);

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
                    vh, isovertices, dt, grid, isovalue, cell_map);

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
        add_single_cycle_candidates_within_hops(unresolved_vertices, /*max_hops=*/2);
        if (single_cycle_candidates.size() == before) {
            break;
        }
    }
}

// ============================================================================
// Sub-routine: Global fallback for single-cycle vertices (Pass 4)
// ============================================================================

void resolve_global_single_cycle_fallback(
    Delaunay& dt,
    const UnifiedGrid& grid,
    float isovalue,
    const std::unordered_map<int, Cell_handle>& cell_map,
    IsovertexComputationStats& stats
) {
    // Small-mesh fallback: if unresolved cases remain and mesh is small enough,
    // run a global pass over all single-cycle vertices.
    const int MAX_GLOBAL_SINGLE_CYCLE_VERTICES = 20000;
    
    if (!stats.had_unresolved || stats.single_cycle_count > MAX_GLOBAL_SINGLE_CYCLE_VERTICES) {
        return;
    }

    const int MAX_GLOBAL_SINGLE_CYCLE_PASSES = 2;
    
    for (int pass = 0; pass < MAX_GLOBAL_SINGLE_CYCLE_PASSES; ++pass) {
        bool any_changed = false;
        int pass_detected = 0;
        int pass_resolved = 0;
        int pass_fallback = 0;
        int pass_unresolved = 0;

        for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
            if (!vit->info().active) continue;
            if (vit->info().is_dummy) continue;
            if (vit->info().facet_cycles.size() != 1) continue;

            auto& isovertices = vit->info().cycle_isovertices;
            const ResolutionResult result = resolve_self_intersection(
                vit, isovertices, dt, grid, isovalue, cell_map);

            update_resolution_stats(result, stats, any_changed,
                pass_detected, pass_resolved, pass_fallback, pass_unresolved);
        }

        stats.self_intersection_detected += pass_detected;
        stats.self_intersection_resolved += pass_resolved;
        stats.self_intersection_fallback += pass_fallback;

        if (pass_detected > 0) {
            DEBUG_PRINT("[DEL-SELFI] global single-cycle pass " << (pass + 1) << "/" << MAX_GLOBAL_SINGLE_CYCLE_PASSES
                        << ": detected=" << pass_detected
                        << ", resolved=" << pass_resolved
                        << ", fallback=" << pass_fallback
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
                                       stats.strat_centroid_projection +
                                       stats.strat_fibonacci_fallback;
        if (resolved_total > 0) {
            auto pct = [&](int64_t count) -> int {
                return static_cast<int>((count * 100 + resolved_total / 2) / resolved_total);
            };

            DEBUG_PRINT("[DEL-SELFI] Resolution strategy breakdown:");
            DEBUG_PRINT("[DEL-SELFI]   centroid_projection: " << stats.strat_centroid_projection
                        << " (" << pct(stats.strat_centroid_projection) << "%)");
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
 * The computation proceeds in four passes:
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
 *     Targeted cleanup of single-cycle vertices near modified multi-cycle vertices.
 *     Only checks vertices within 2 hops to limit overhead.
 *
 *   Pass 4 (resolve_global_single_cycle_fallback):
 *     For small meshes with unresolved cases, run a global check on all single-cycle vertices.
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

    // Build cell lookup map once for efficient cell access by index.
    std::unordered_map<int, Cell_handle> cell_map = build_cell_index_map(dt);

    // Initialize statistics tracking.
    IsovertexComputationStats stats;

    // ========================================================================
    // Pass 1: Compute initial isovertex positions for all active vertices.
    // ========================================================================
    std::vector<Vertex_handle> multi_cycle_vertices = compute_initial_cycle_isovertices(
        dt, grid, isovalue, cell_map, stats);

    // ========================================================================
    // Pass 2: Resolve self-intersections on multi-cycle vertices.
    // ========================================================================
    std::vector<Vertex_handle> modified_vertices = resolve_multi_cycle_self_intersections(
        dt, grid, isovalue, cell_map, multi_cycle_vertices, stats);

    // ========================================================================
    // Pass 3: Targeted cleanup of single-cycle vertices near modified vertices.
    // ========================================================================
    resolve_single_cycle_self_intersections(
        dt, grid, isovalue, cell_map, modified_vertices, stats);

    // ========================================================================
    // Pass 4: Global fallback for small meshes with unresolved cases.
    // ========================================================================
    resolve_global_single_cycle_fallback(dt, grid, isovalue, cell_map, stats);

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
    const std::unordered_map<int, Cell_handle>& cell_map,
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
        auto it = cell_map.find(cell_index);
        if (it == cell_map.end()) continue;

        Cell_handle cell = it->second;
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
        auto state_it = edge_matching_states.find(ekey);
        bool applied_matching = false;

        if (state_it != edge_matching_states.end()) {
            auto edge_it = neighbor_edge.find(neighbor_index);
            if (edge_it != neighbor_edge.end()) {
                auto pairs = compute_edge_bipolar_matching(dt, edge_it->second, state_it->second.method, cell_map);
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

        if (!applied_matching) {
            // Clique fallback
            for (size_t i = 0; i < facet_ids.size(); ++i) {
                for (size_t j = i + 1; j < facet_ids.size(); ++j) {
                    adj[facet_ids[i]].push_back(facet_ids[j]);
                    adj[facet_ids[j]].push_back(facet_ids[i]);
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
}

ModifyCyclesResult modify_cycles_pass(Delaunay& dt) {
    ModifyCyclesResult result;
    const int MAX_ITERATIONS = 10;

    std::unordered_map<int, Cell_handle> cell_map = build_cell_index_map(dt);

    // Build initial edge matching states (all start with SEP_NEG as in compute_facet_cycles)
    std::unordered_map<EdgeKey, EdgeMatchingState, EdgeKeyHash> edge_matching_states;

    // Collect all ambiguous edges (edges with >2 incident isosurface facets)
    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
        if (!vit->info().active || vit->info().is_dummy) continue;

        std::vector<Edge> incident_edges;
        dt.incident_edges(vit, std::back_inserter(incident_edges));

        for (const Edge& e : incident_edges) {
            Vertex_handle a = e.first->vertex(e.second);
            Vertex_handle b = e.first->vertex(e.third);
            if (&*a != &*vit && &*b != &*vit) continue;
            Vertex_handle other = (&*a == &*vit) ? b : a;

            EdgeKey ekey = EdgeKey::FromVertexIndices(vit->info().index, other->info().index);
            if (edge_matching_states.find(ekey) == edge_matching_states.end()) {
                edge_matching_states[ekey] = EdgeMatchingState{};
            }
        }
    }

    // Build vertex index -> handle map for recomputation
    std::unordered_map<int, Vertex_handle> vertex_index_to_handle;
    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
        if (!vit->info().active || vit->info().is_dummy) continue;
        vertex_index_to_handle[vit->info().index] = vit;
    }

    // Iterate until no conflicts or max iterations
    for (int iter = 0; iter < MAX_ITERATIONS; ++iter) {
        result.iterations = iter + 1;
        bool found_conflict = false;
        std::set<int> dirty_vertex_indices;
        int conflicts_this_iteration = 0;

        // Check each vertex for problematic iso-segment assignments ( not just multi-cycle ones, including single-cycles)
        for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
            if (!vit->info().active || vit->info().is_dummy) continue;

            // Process all edges, not just those around multi-cycle vertices
            std::vector<Edge> incident_edges;
            dt.incident_edges(vit, std::back_inserter(incident_edges));

            for (const Edge& e : incident_edges) {
                Vertex_handle a = e.first->vertex(e.second);
                Vertex_handle b = e.first->vertex(e.third);
                if (&*a != &*vit && &*b != &*vit) continue;
                Vertex_handle other = (&*a == &*vit) ? b : a;

                EdgeKey ekey = EdgeKey::FromVertexIndices(vit->info().index, other->info().index);
                auto state_it = edge_matching_states.find(ekey);
                if (state_it == edge_matching_states.end()) continue;
                if (state_it->second.flipped) continue; // Already flipped this iteration

                auto pairs = compute_edge_bipolar_matching(dt, e, state_it->second.method, cell_map);
                if (pairs.size() < 2) continue;

                // Check for duplicate cycle pairs
                std::set<std::pair<int, int>> seen_cycle_pairs;
                bool has_duplicate = false;

                for (const auto& [k0, k1] : pairs) {
                    // k0 and k1 are matched facets around this edge
                    // Find which cycle each facet belongs to on its respective endpoint vertex
                    // c_vit = cycle on vit, c_other = cycle on other
                    int c_vit = find_cycle_containing_facet(vit, k0.cell_index, k0.facet_index);
                    if (c_vit < 0) {
                        c_vit = find_cycle_containing_facet(vit, k1.cell_index, k1.facet_index);
                    }
                    
                    int c_other = -1;
                    if (!other->info().is_dummy && other->info().active) {
                        c_other = find_cycle_containing_facet(other, k0.cell_index, k0.facet_index);
                        if (c_other < 0) {
                            c_other = find_cycle_containing_facet(other, k1.cell_index, k1.facet_index);
                        }
                    }

                    if (c_vit >= 0 && c_other >= 0) {
                        // The cycle pair is (vit.cycle, other.cycle)
                        auto cycle_pair = std::make_pair(c_vit, c_other);
                        if (seen_cycle_pairs.count(cycle_pair)) {
                            has_duplicate = true;
                            break;
                        }
                        seen_cycle_pairs.insert(cycle_pair);
                    }
                }

                if (has_duplicate) {
                    flip_edge_matching(state_it->second);
                    found_conflict = true;
                    conflicts_this_iteration++;
                    result.total_flips++;
                    dirty_vertex_indices.insert(vit->info().index);
                    if (!other->info().is_dummy && other->info().active) {
                        dirty_vertex_indices.insert(other->info().index);
                    }
                }
            }
        }

        result.problematic_edges = conflicts_this_iteration;

        if (!found_conflict) {
            break;
        }

        // Recompute cycles for affected vertices
        for (int v_idx : dirty_vertex_indices) {
            auto it = vertex_index_to_handle.find(v_idx);
            if (it != vertex_index_to_handle.end()) {
                recompute_cycles_for_vertex(dt, it->second, cell_map, edge_matching_states);
            }
        }

        // Reset flipped flags for next iteration
        for (auto& [k, state] : edge_matching_states) {
            state.flipped = false;
        }
    }

    // ========================================================================
    // Phase 2: Cross-edge conflict detection
    // Detect when two different edges map to the same global (vertex-cycle, vertex-cycle) pair
    // ========================================================================
    
    // Helper to make a unique component key from (vertex_index, cycle_index)
    auto make_component_key = [](int vertex_idx, int cycle_idx) -> int64_t {
        return (static_cast<int64_t>(vertex_idx) << 32) | static_cast<uint32_t>(cycle_idx);
    };

    const int MAX_CROSS_ITERATIONS = 10;
    for (int cross_iter = 0; cross_iter < MAX_CROSS_ITERATIONS; ++cross_iter) {
        bool found_cross_conflict = false;
        
        // Map from (compKeyA, compKeyB) -> EdgeKey that produces this pairing
        std::map<std::pair<int64_t, int64_t>, EdgeKey> globalCompPairToEdge;
        std::set<int> dirty_vertex_indices_cross;

        // Scan all edges and their current matchings
        for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
            if (!vit->info().active || vit->info().is_dummy) continue;

            const auto& cycles = vit->info().facet_cycles;
            if (cycles.size() <= 1) continue;

            std::vector<Edge> incident_edges;
            dt.incident_edges(vit, std::back_inserter(incident_edges));

            for (const Edge& e : incident_edges) {
                Vertex_handle a = e.first->vertex(e.second);
                Vertex_handle b = e.first->vertex(e.third);
                if (&*a != &*vit && &*b != &*vit) continue;
                Vertex_handle other = (&*a == &*vit) ? b : a;

                if (other->info().is_dummy) continue;

                EdgeKey ekey = EdgeKey::FromVertexIndices(vit->info().index, other->info().index);
                auto state_it = edge_matching_states.find(ekey);
                if (state_it == edge_matching_states.end()) continue;
                if (state_it->second.flipped) continue;

                auto pairs = compute_edge_bipolar_matching(dt, e, state_it->second.method, cell_map);
                
                for (const auto& [k0, k1] : pairs) {
                    // Get cycle indices from both endpoints
                    int c_vit = find_cycle_containing_facet(vit, k0.cell_index, k0.facet_index);
                    int c_other = -1;
                    if (!other->info().is_dummy && other->info().active) {
                        c_other = find_cycle_containing_facet(other, k1.cell_index, k1.facet_index);
                    }
                    
                    if (c_vit < 0 || c_other < 0) continue;

                    int64_t keyA = make_component_key(vit->info().index, c_vit);
                    int64_t keyB = make_component_key(other->info().index, c_other);
                    if (keyA > keyB) std::swap(keyA, keyB);
                    auto globalKey = std::make_pair(keyA, keyB);

                    auto git = globalCompPairToEdge.find(globalKey);
                    if (git == globalCompPairToEdge.end()) {
                        globalCompPairToEdge.emplace(globalKey, ekey);
                        continue;
                    }

                    EdgeKey otherEdge = git->second;
                    if (otherEdge == ekey) continue;

                    // Conflict: two different edges map to same component pair
                    // Flip the current edge if not already flipped
                    auto other_state_it = edge_matching_states.find(otherEdge);
                    EdgeKey target_edge = ekey;
                    
                    if (state_it->second.flipped && other_state_it != edge_matching_states.end() 
                        && !other_state_it->second.flipped) {
                        target_edge = otherEdge;
                    }

                    auto& target_state = edge_matching_states[target_edge];
                    flip_edge_matching(target_state);
                    result.total_flips++;
                    found_cross_conflict = true;

                    // Mark affected vertices dirty
                    dirty_vertex_indices_cross.insert(target_edge.v0);
                    dirty_vertex_indices_cross.insert(target_edge.v1);
                    break;
                }
                if (found_cross_conflict) break;
            }
            if (found_cross_conflict) break;
        }

        if (!found_cross_conflict) break;

        // Recompute cycles for affected vertices
        for (int v_idx : dirty_vertex_indices_cross) {
            auto it = vertex_index_to_handle.find(v_idx);
            if (it != vertex_index_to_handle.end()) {
                recompute_cycles_for_vertex(dt, it->second, cell_map, edge_matching_states);
            }
        }

        // Reset flipped flags
        for (auto& [k, state] : edge_matching_states) {
            state.flipped = false;
        }

        result.iterations++;
    }

    if (result.total_flips > 0) {
        std::cout << "  ModCyc: " << result.total_flips << " edge flips ("
                  << result.iterations << " iterations)" << std::endl;
    }

    DEBUG_PRINT("[CYC-MOD] mod_cyc: flips=" << result.total_flips
                << ", iterations=" << result.iterations);

    return result;
}
