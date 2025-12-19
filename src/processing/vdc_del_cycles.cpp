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
    // (This is a conservative approximation; true circumscribed radius is sqrt(3)/2 * side)
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
        // These vertices may also have cycles - we need to find which cycle
        // contains this facet for each vertex
        Point p1 = cycle_isovertex; // Our cycle's isovertex

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

        // Only create triangle if we found all three vertices
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

    // If triangles share 2 or more vertices, they share an edge - this is expected
    // and not a self-intersection in the problematic sense
    if (shared_count >= 2) {
        return false;
    }

    // Use CGAL's intersection test
    return CGAL::do_intersect(t1, t2);
}

bool check_self_intersection(
    Vertex_handle v,
    const std::vector<Point>& cycle_isovertices,
    const Delaunay& dt
) {
    const auto& cycles = v->info().facet_cycles;
    size_t num_cycles = cycles.size();

    if (num_cycles < 2) {
        return false; // No self-intersection possible with single cycle
    }

    // Build cell map for quick lookup
    std::unordered_map<int, Cell_handle> cell_map;
    for (auto cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit) {
        cell_map[cit->info().index] = cit;
    }

    // Temporarily set the cycle isovertices for triangle collection
    // (Save original values to restore later if needed)
    std::vector<Point> original_isovertices = v->info().cycle_isovertices;

    // Note: We need to use a const_cast here because we're checking with proposed positions
    // In the actual integration, the values will already be set
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
    double sphere_radius
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

    // Step 4: Verify no self-intersection; if still intersecting, fall back to vertex position
    if (check_self_intersection(v, cycle_isovertices, dt)) {
        // Fallback: place all isovertices at vertex center
        DEBUG_PRINT("[DEL-SELFI] Vertex " << v->info().index
                    << ": fallback to center position after sphere projection failed");
        for (int i = 0; i < n; ++i) {
            cycle_isovertices[i] = center;
        }
    }
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

        if (check_self_intersection(vit, isovertices, dt)) {
            self_intersection_detected++;

            const double sphere_radius = compute_sphere_radius(vit, dt, grid);
            resolve_self_intersection(vit, isovertices, dt, sphere_radius);

            const Point cube_center = vit->point();
            bool all_at_center = true;
            for (const auto& pos : isovertices) {
                double dx = pos.x() - cube_center.x();
                double dy = pos.y() - cube_center.y();
                double dz = pos.z() - cube_center.z();
                if (dx*dx + dy*dy + dz*dz > 1e-10) {
                    all_at_center = false;
                    break;
                }
            }

            if (all_at_center) {
                self_intersection_fallback++;
            } else {
                self_intersection_resolved++;
            }
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
