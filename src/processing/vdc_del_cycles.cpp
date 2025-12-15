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
// Step 8: Compute Facet Cycles
// ============================================================================

//! @brief Check if two facets share an edge (excluding the central vertex v)
/*!
 * Two facets around vertex v share an edge if they have two vertices in common
 * (v itself plus one other vertex).
 */
static bool facets_share_edge_around_vertex(
    const Delaunay& dt,
    Cell_handle cell1, int facet_idx1,
    Cell_handle cell2, int facet_idx2,
    Vertex_handle v
) {
    // Get the vertices of each facet (excluding the opposite vertex)
    std::set<Vertex_handle> verts1, verts2;

    for (int k = 0; k < 4; ++k) {
        if (k != facet_idx1) {
            verts1.insert(cell1->vertex(k));
        }
        if (k != facet_idx2) {
            verts2.insert(cell2->vertex(k));
        }
    }

    // Count common vertices
    int common_count = 0;
    for (auto vh : verts1) {
        if (verts2.count(vh) > 0) {
            common_count++;
        }
    }

    // Two facets sharing vertex v share an edge if they have exactly 2 common vertices
    // (v plus one edge vertex)
    return common_count >= 2;
}

void compute_facet_cycles(Delaunay& dt) {
    DEBUG_PRINT("[DEL-CYCLE] Computing facet cycles around active vertices...");

    int total_cycles = 0;
    int multi_cycle_vertices = 0;

    // Build a map for quick cell lookup by index
    std::unordered_map<int, Cell_handle> cell_map = build_cell_index_map(dt);

    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
        if (!vit->info().active) continue;

        // Clear any existing cycles
        vit->info().facet_cycles.clear();

        // Collect all isosurface facets incident on this vertex
        // Store as (cell_index, facet_index) pairs
        std::vector<std::pair<int, int>> isosurface_facets;
        std::vector<std::pair<Cell_handle, int>> facet_handles;

        // Iterate through all cells incident on this vertex
        std::vector<Cell_handle> incident_cells;
        dt.incident_cells(vit, std::back_inserter(incident_cells));

        for (Cell_handle ch : incident_cells) {
            if (dt.is_infinite(ch)) continue;

            // Check each facet of this cell
            for (int i = 0; i < 4; ++i) {
                if (!ch->info().facet_is_isosurface[i]) continue;

                // Check if this facet is incident on vertex vit
                // Facet i is opposite vertex i, so vit must not be at position i
                int v_idx = ch->index(vit);
                if (v_idx != i) {
                    // This isosurface facet is incident on vit
                    isosurface_facets.push_back({ch->info().index, i});
                    facet_handles.push_back({ch, i});
                }
            }
        }

        if (isosurface_facets.empty()) {
            // Should not happen if vertex is marked active
            DEBUG_PRINT("[DEL-CYCLE] Warning: Active vertex "
                        << vit->info().index << " has no isosurface facets");
            continue;
        }

        // Build adjacency graph for facets
        // Two facets are adjacent if they share an edge around this vertex
        int n = isosurface_facets.size();
        std::vector<std::vector<int>> adj(n);

        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (facets_share_edge_around_vertex(
                        dt,
                        facet_handles[i].first, facet_handles[i].second,
                        facet_handles[j].first, facet_handles[j].second,
                        vit)) {
                    adj[i].push_back(j);
                    adj[j].push_back(i);
                }
            }
        }

        // Find connected components (cycles) using DFS
        std::vector<bool> visited(n, false);
        std::vector<std::vector<std::pair<int, int>>> cycles;

        for (int start = 0; start < n; ++start) {
            if (visited[start]) continue;

            // Start a new cycle (connected component)
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
                    if (!visited[neighbor]) {
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
}

// ============================================================================
// Step 9: Compute Cycle Isovertices
// ============================================================================

Point compute_centroid_of_voronoi_edge_and_isosurface(
    const Delaunay& dt,
    const UnifiedGrid& grid,
    float isovalue,
    Vertex_handle v_handle,
    const std::vector<std::pair<int, int>>& cycle_facets
) {
    // Build cell lookup map
    std::unordered_map<int, Cell_handle> cell_map;
    for (auto cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit) {
        cell_map[cit->info().index] = cit;
    }

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
        float t = (isovalue - sB) / denom;

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

void compute_cycle_isovertices(
    Delaunay& dt,
    const UnifiedGrid& grid,
    float isovalue,
    bool position_on_isov
) {
    DEBUG_PRINT("[DEL-ISOV] Computing isovertex positions for cycles...");

    int single_cycle_count = 0;
    int multi_cycle_count = 0;

    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
        if (!vit->info().active) continue;

        auto& cycles = vit->info().facet_cycles;
        auto& isovertices = vit->info().cycle_isovertices;

        if (cycles.empty()) continue;

        isovertices.resize(cycles.size());

        if (cycles.size() == 1) {
            // Single cycle: isovertex at vertex location (or isov if position_on_isov)
            if (position_on_isov) {
                // Use the pre-computed accurate iso-crossing point
                isovertices[0] = vit->info().isov;
            } else {
                // Use the Delaunay vertex position
                isovertices[0] = vit->point();
            }
            single_cycle_count++;
        } else {
            // Multiple cycles: compute centroids for each cycle
            for (size_t c = 0; c < cycles.size(); ++c) {
                Point centroid = compute_centroid_of_voronoi_edge_and_isosurface(
                    dt, grid, isovalue, vit, cycles[c]);

                // For multi-cycle vertices, we could optionally project the centroid
                // onto a sphere around the vertex. For now, just use the centroid directly.
                // In the future, self-intersection check (Step 9.b-c) would be added here.
                isovertices[c] = centroid;
            }
            multi_cycle_count++;
        }
    }

    DEBUG_PRINT("[DEL-ISOV] Computed isovertices for "
                << single_cycle_count << " single-cycle and "
                << multi_cycle_count << " multi-cycle vertices");
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
