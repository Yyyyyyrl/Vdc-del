//! @file vdc_del_cycles.cpp
//! @brief Implementation of cycle detection for Delaunay-based isosurface extraction.
//!
//! This file implements Steps 8-9 of the Delaunay-based VDC algorithm:
//! - Step 8: Compute cycles of isosurface facets around active vertices
//! - Step 9: Compute isosurface vertex positions for each cycle

#include "processing/vdc_del_isosurface.h"
#include "core/vdc_debug.h"
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <optional>
#include <CGAL/Triangle_3.h>
#include <CGAL/intersections.h>

// ============================================================================
// Helper: Build cell index lookup
// ============================================================================

//! @brief Build a mapping from cell indices to cell handles
static std::unordered_map<int, Cell_handle> build_cell_index_map(const Delaunay& dt) {
    std::unordered_map<int, Cell_handle> map;
    for (auto cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit) {
        map[cit->info().index] = cit;
    }
    return map;
}

// ============================================================================
// Helpers: Matching + adjacency on a Delaunay edge
// ============================================================================

struct FacetKey {
    int cell_index = -1;  // positive cell index
    int facet_index = -1; // facet index in that positive cell

    bool operator==(const FacetKey& other) const {
        return cell_index == other.cell_index && facet_index == other.facet_index;
    }
};

struct FacetKeyHash {
    size_t operator()(const FacetKey& key) const noexcept {
        // 32-bit pack is enough for typical indices; still fine on 64-bit.
        return (static_cast<size_t>(static_cast<uint32_t>(key.cell_index)) << 32) ^
               static_cast<size_t>(static_cast<uint32_t>(key.facet_index));
    }
};

struct EdgeKey {
    int v0 = -1;
    int v1 = -1;

    static EdgeKey FromVertexIndices(int a, int b) {
        if (a < b) {
            return EdgeKey{a, b};
        }
        return EdgeKey{b, a};
    }

    bool operator==(const EdgeKey& other) const {
        return v0 == other.v0 && v1 == other.v1;
    }
};

struct EdgeKeyHash {
    size_t operator()(const EdgeKey& key) const noexcept {
        return (static_cast<size_t>(static_cast<uint32_t>(key.v0)) << 32) ^
               static_cast<size_t>(static_cast<uint32_t>(key.v1));
    }
};

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
            desired_between_positive = false;
            break;
        case BIPOLAR_MATCH_METHOD::SEP_POS:
            desired_between_positive = true;
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
// Step 8: Compute Facet Cycles
// ============================================================================

void compute_facet_cycles(Delaunay& dt) {
    DEBUG_PRINT("[DEL-CYCLE] Computing facet cycles around active vertices...");

    int total_cycles = 0;
    int multi_cycle_vertices = 0;
    int ambiguous_edges_total = 0;
    int ambiguous_edges_matched = 0;

    // Build a map for quick cell lookup by index
    std::unordered_map<int, Cell_handle> cell_map = build_cell_index_map(dt);

    // Cache Delaunay-edge bipolar matchings (Voronoi-facet matchings).
    // Keyed by (min(v0,v1), max(v0,v1)).
    std::unordered_map<EdgeKey, std::vector<std::pair<FacetKey, FacetKey>>, EdgeKeyHash> edge_matching_cache;

    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
        if (!vit->info().active) continue;
        if (vit->info().is_dummy) continue;
        const Vertex_handle v_handle = vit;

        // Clear any existing cycles
        vit->info().facet_cycles.clear();

        // Collect all isosurface facets incident on this vertex
        // Store as (cell_index, facet_index) pairs
        std::vector<std::pair<int, int>> isosurface_facets;

        // Iterate through all cells incident on this vertex
        std::vector<Cell_handle> incident_cells;
        dt.incident_cells(v_handle, std::back_inserter(incident_cells));

        for (Cell_handle ch : incident_cells) {
            if (dt.is_infinite(ch)) continue;

            // Check each facet of this cell
            for (int i = 0; i < 4; ++i) {
                if (!ch->info().facet_is_isosurface[i]) continue;

                // Check if this facet is incident on vertex vit
                // Facet i is opposite vertex i, so vit must not be at position i
                const int v_idx = find_vertex_index_in_cell(ch, v_handle);
                if (v_idx < 0) {
                    continue;
                }
                if (v_idx != i) {
                    // Ignore facets involving any dummy Delaunay sites.
                    bool has_dummy = false;
                    for (int j = 0; j < 4; ++j) {
                        if (j == i) continue;
                        if (ch->vertex(j)->info().is_dummy) {
                            has_dummy = true;
                            break;
                        }
                    }
                    if (has_dummy) {
                        continue;
                    }

                    // This isosurface facet is incident on vit
                    isosurface_facets.push_back({ch->info().index, i});
                }
            }
        }

        if (isosurface_facets.empty()) {
            // Should not happen if vertex is marked active
            DEBUG_PRINT("[DEL-CYCLE] Warning: Active vertex "
                        << vit->info().index << " has no isosurface facets");
            continue;
        }

        // Build facet->id map for adjacency and matching lookup.
        const int n = static_cast<int>(isosurface_facets.size());
        std::unordered_map<FacetKey, int, FacetKeyHash> facetkey_to_id;
        facetkey_to_id.reserve(isosurface_facets.size());
        for (int i = 0; i < n; ++i) {
            facetkey_to_id.emplace(FacetKey{isosurface_facets[i].first, isosurface_facets[i].second}, i);
        }

        // For each Delaunay edge (v,u), collect the incident isosurface facets in this star.
        const int v_global_index = vit->info().index;
        std::unordered_map<int, std::vector<int>> neighbor_to_facets;
        neighbor_to_facets.reserve(isosurface_facets.size());

        for (int fid = 0; fid < n; ++fid) {
            const int cell_index = isosurface_facets[fid].first;
            const int facet_index = isosurface_facets[fid].second;
            auto it = cell_map.find(cell_index);
            if (it == cell_map.end()) {
                continue;
            }

            Cell_handle cell = it->second;
            Vertex_handle other0 = nullptr;
            Vertex_handle other1 = nullptr;
            for (int j = 0; j < 4; ++j) {
                if (j == facet_index) continue;
                Vertex_handle vh = cell->vertex(j);
                if (vh == v_handle) continue;
                if (other0 == nullptr) {
                    other0 = vh;
                } else {
                    other1 = vh;
                    break;
                }
            }

            if (other0 == nullptr || other1 == nullptr) {
                continue;
            }
            neighbor_to_facets[other0->info().index].push_back(fid);
            neighbor_to_facets[other1->info().index].push_back(fid);
        }

        // Map neighbor vertex index -> an Edge object representing (v, neighbor).
        std::unordered_map<int, Edge> neighbor_edge;
        {
            std::vector<Edge> incident_edges;
            dt.incident_edges(v_handle, std::back_inserter(incident_edges));
            neighbor_edge.reserve(incident_edges.size());

            for (const Edge& e : incident_edges) {
                Vertex_handle a = e.first->vertex(e.second);
                Vertex_handle b = e.first->vertex(e.third);
                if (a != v_handle && b != v_handle) {
                    continue;
                }
                Vertex_handle other = (a == v_handle) ? b : a;
                neighbor_edge[other->info().index] = e;
            }
        }

        // Build adjacency graph for facets:
        // - Baseline: clique connections among facets incident to the same Delaunay edge (v,u).
        // - Refinement: for ambiguous edges with >2 incident facets, replace clique with
        //   Voronoi-facet bipolar matching derived from the Delaunay-edge ring.
        std::vector<std::vector<int>> adj(n);
        for (const auto& [neighbor_index, facet_ids] : neighbor_to_facets) {
            if (facet_ids.size() < 2) {
                continue;
            }

            if (facet_ids.size() == 2) {
                const int a = facet_ids[0];
                const int b = facet_ids[1];
                adj[a].push_back(b);
                adj[b].push_back(a);
                continue;
            }

            ambiguous_edges_total++;

            // Try to apply matching on this ambiguous edge.
            bool applied_matching = false;
            const auto edge_it = neighbor_edge.find(neighbor_index);
            if (edge_it != neighbor_edge.end()) {
                const Edge& e = edge_it->second;
                const EdgeKey ekey = EdgeKey::FromVertexIndices(v_global_index, neighbor_index);

                auto cache_it = edge_matching_cache.find(ekey);
                if (cache_it == edge_matching_cache.end()) {
                    edge_matching_cache[ekey] = compute_edge_bipolar_matching(
                        dt, e, BIPOLAR_MATCH_METHOD::SEP_NEG, cell_map);
                    cache_it = edge_matching_cache.find(ekey);
                }

                const auto& pairs = cache_it->second;
                if (!pairs.empty()) {
                    std::unordered_set<int> group_set(facet_ids.begin(), facet_ids.end());
                    std::unordered_set<int> covered;
                    covered.reserve(facet_ids.size());

                    int edge_connections_added = 0;
                    for (const auto& [k0, k1] : pairs) {
                        auto it0 = facetkey_to_id.find(k0);
                        auto it1 = facetkey_to_id.find(k1);
                        if (it0 == facetkey_to_id.end() || it1 == facetkey_to_id.end()) {
                            continue;
                        }

                        const int id0 = it0->second;
                        const int id1 = it1->second;
                        if (group_set.count(id0) == 0 || group_set.count(id1) == 0) {
                            continue;
                        }

                        adj[id0].push_back(id1);
                        adj[id1].push_back(id0);
                        covered.insert(id0);
                        covered.insert(id1);
                        edge_connections_added++;
                    }

                    if (covered.size() == facet_ids.size() &&
                        edge_connections_added == static_cast<int>(facet_ids.size() / 2)) {
                        applied_matching = true;
                        ambiguous_edges_matched++;
                    } else {
                        // Roll back partial matching edges by clearing and falling back to clique, should never happen though.
                        applied_matching = false;
                        DEBUG_PRINT("Rolling back partial matching edges\n");
                    }
                }
            }

            if (!applied_matching) {
                // Clique fallback to avoid seams/cracks if matching is incomplete.
                for (size_t i = 0; i < facet_ids.size(); ++i) {
                    for (size_t j = i + 1; j < facet_ids.size(); ++j) {
                        const int a = facet_ids[i];
                        const int b = facet_ids[j];
                        adj[a].push_back(b);
                        adj[b].push_back(a);
                    }
                }
            }
        }

        // Find connected components (cycles) using DFS
        std::vector<bool> visited(static_cast<size_t>(n), false);
        std::vector<std::vector<std::pair<int, int>>> cycles;

        for (int start = 0; start < n; ++start) {
            if (visited[static_cast<size_t>(start)]) continue;

            // Start a new cycle (connected component)
            std::vector<std::pair<int, int>> cycle;
            std::stack<int> stack;
            stack.push(start);

            while (!stack.empty()) {
                int curr = stack.top();
                stack.pop();

                if (visited[static_cast<size_t>(curr)]) continue;
                visited[static_cast<size_t>(curr)] = true;

                cycle.push_back(isosurface_facets[curr]);

                for (int neighbor : adj[curr]) {
                    if (!visited[static_cast<size_t>(neighbor)]) {
                        stack.push(neighbor);
                    }
                }
            }

            cycles.push_back(cycle);
        }

        // Store cycles in vertex info
        vit->info().facet_cycles = cycles;
        total_cycles += cycles.size();

        if (cycles.size() > 1) {
            multi_cycle_vertices++;
            DEBUG_PRINT("[DEL-CYCLE] Vertex " << vit->info().index
                        << " has " << cycles.size() << " cycles");
        }
    }

    DEBUG_PRINT("[DEL-CYCLE] Found " << total_cycles << " cycles total, "
                << multi_cycle_vertices << " vertices with multiple cycles");
    DEBUG_PRINT("[DEL-CYCLE] Ambiguous Delaunay edges: " << ambiguous_edges_total
                << " (matched " << ambiguous_edges_matched << ")");
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
    // Check if triangles share any vertices
    std::set<Point, std::less<Point>> t1_verts = {t1.vertex(0), t1.vertex(1), t1.vertex(2)};
    int shared_count = 0;
    for (int i = 0; i < 3; ++i) {
        if (t1_verts.count(t2.vertex(i)) > 0) {
            shared_count++;
        }
    }

    // If triangles share 1+ vertices, they touch along a mesh vertex/edge.
    if (shared_count >= 1) {
        return false;
    }

    // Use CGAL's intersection test
    return CGAL::do_intersect(t1, t2);
}

bool check_self_intersection(
    Vertex_handle v,
    const std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    const std::unordered_map<int, Cell_handle>& cell_map
) {
    const auto& cycles = v->info().facet_cycles;
    size_t num_cycles = cycles.size();

    if (num_cycles < 2) {
        return false; // No self-intersection possible with single cycle
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

    // Check all pairs of cycles for intersection
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

//! @brief Compute minimum angular separation between points on a sphere
static double compute_min_angular_separation(
    const std::vector<Point>& points,
    const Point& center
) {
    if (points.size() < 2) return 2.0 * M_PI; // Maximum separation for single point

    double min_angle = 2.0 * M_PI;

    for (size_t i = 0; i < points.size(); ++i) {
        for (size_t j = i + 1; j < points.size(); ++j) {
            // Compute vectors from center to each point
            Vector3 v1(center, points[i]);
            Vector3 v2(center, points[j]);

            double len1 = std::sqrt(v1.squared_length());
            double len2 = std::sqrt(v2.squared_length());

            if (len1 < 1e-10 || len2 < 1e-10) continue;

            // Compute angle between vectors
            double dot = (v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z()) / (len1 * len2);
            dot = std::max(-1.0, std::min(1.0, dot)); // Clamp for numerical stability
            double angle = std::acos(dot);

            min_angle = std::min(min_angle, angle);
        }
    }

    return min_angle;
}

//! @brief Redistribute points evenly on a sphere
/*!
 * Places points at evenly distributed positions on the sphere.
 * Uses a simple approach: distribute along a great circle for 2 points,
 * or use vertices of a regular polyhedron for more points.
 */
static void redistribute_on_sphere(
    std::vector<Point>& points,
    const Point& center,
    double radius
) {
    int n = static_cast<int>(points.size());
    if (n < 2) return;

    // For simplicity, distribute points evenly around a great circle
    // More sophisticated methods could use spherical Fibonacci or other distributions
    for (int i = 0; i < n; ++i) {
        double theta = 2.0 * M_PI * i / n;
        double x = center.x() + radius * std::cos(theta);
        double y = center.y() + radius * std::sin(theta);
        double z = center.z();
        points[i] = Point(x, y, z);
    }
}

void resolve_self_intersection(
    Vertex_handle v,
    std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    double sphere_radius,
    const std::unordered_map<int, Cell_handle>& cell_map
) {
    Point center = v->point();
    int n = static_cast<int>(cycle_isovertices.size());

    if (n < 2) return;

    // Step 1: Project all centroids onto sphere
    for (int i = 0; i < n; ++i) {
        cycle_isovertices[i] = project_to_sphere(cycle_isovertices[i], center, sphere_radius);
    }

    // Step 2: Check minimum angular separation
    double min_angle = compute_min_angular_separation(cycle_isovertices, center);
    double required_angle = 2.0 * M_PI / (3.0 * n); // Heuristic threshold

    // Step 3: If too close, redistribute on sphere
    if (min_angle < required_angle) {
        redistribute_on_sphere(cycle_isovertices, center, sphere_radius);
    }

    // Caller decides what to do if this is still intersecting.
}

void compute_cycle_isovertices(
    Delaunay& dt,
    const UnifiedGrid& grid,
    float isovalue,
    bool position_on_isov
) {
    DEBUG_PRINT("[DEL-ISOV] Computing isovertex positions for cycles...");

    int single_cycle_count = 0;
    int multi_cycle_count = 0;
    int self_intersection_detected = 0;
    int self_intersection_resolved = 0;
    int self_intersection_fallback = 0;

    // Build cell lookup map ONCE for the entire function (performance fix)
    std::unordered_map<int, Cell_handle> cell_map = build_cell_index_map(dt);

    // Pass 1: compute initial isovertices for all active vertices.
    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
        if (!vit->info().active) continue;
        if (vit->info().is_dummy) continue;

        auto& cycles = vit->info().facet_cycles;
        auto& isovertices = vit->info().cycle_isovertices;

        if (cycles.empty()) continue;

        isovertices.resize(cycles.size());

        Point cube_center = vit->point();

        if (cycles.size() == 1) {
            // Single cycle: use the isosurface sample point associated with this site.
            // When `-position_delv_on_isov` is enabled, this equals `cube_center`.
            (void)position_on_isov; // the choice is encoded in `vit->info().isov`
            isovertices[0] = vit->info().isov;
            single_cycle_count++;
        } else {
            // Multiple cycles: compute centroids and project them to a sphere around the site.
            // This gives directional separation between cycles.
            const double sphere_radius = compute_sphere_radius(vit, dt, grid);
            for (size_t c = 0; c < cycles.size(); ++c) {
                Point centroid = compute_centroid_of_voronoi_edge_and_isosurface(
                    dt, grid, isovalue, vit, cycles[c], cell_map);

                isovertices[c] = project_to_sphere(centroid, cube_center, sphere_radius);
            }
            multi_cycle_count++;
        }
    }

    // Pass 2: resolve self-intersections after all isovertices are assigned.
    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
        if (!vit->info().active) continue;
        if (vit->info().is_dummy) continue;

        auto& cycles = vit->info().facet_cycles;
        auto& isovertices = vit->info().cycle_isovertices;

        if (cycles.size() < 2) continue;

        if (check_self_intersection(vit, isovertices, dt, cell_map)) {
            self_intersection_detected++;

            const double sphere_radius = compute_sphere_radius(vit, dt, grid);
            resolve_self_intersection(vit, isovertices, dt, sphere_radius, cell_map);

        }
    }

    DEBUG_PRINT("[DEL-ISOV] Computed isovertices for "
                << single_cycle_count << " single-cycle and "
                << multi_cycle_count << " multi-cycle vertices");

    if (self_intersection_detected > 0) {
        DEBUG_PRINT("[DEL-SELFI] Self-intersection stats: "
                    << self_intersection_detected << " detected, "
                    << self_intersection_resolved << " resolved by sphere projection, "
                    << self_intersection_fallback << " using fallback (center position)");
    }
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

//! @brief Stores matching state per edge for potential flipping
struct EdgeMatchingState {
    BIPOLAR_MATCH_METHOD method = BIPOLAR_MATCH_METHOD::SEP_NEG;
    int unconstrained_offset = 0;
    bool flipped = false;
    int flip_count = 0;  // Permanent flip count to prevent oscillation
};

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

        // Check each vertex for problematic iso-segment assignments
        // Note: We need to check ALL vertices (not just multi-cycle ones) because
        // non-manifolds can occur when multiple facets on the same edge map to
        // the same (vertex0.cycle0, vertex1.cycle0) pair even with single-cycle endpoints
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

