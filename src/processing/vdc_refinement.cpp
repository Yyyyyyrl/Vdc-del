//! @file vdc_refinement.cpp
//! @brief Implementation of facet-centric Delaunay refinement.

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>
#include <set>
#include <CGAL/squared_distance_3.h>

#include "processing/vdc_refinement.h"
#include "core/vdc_utilities.h"

constexpr double kPi = 3.14159265358979323846;

// Helper: Compute minimum and maximum angles of a triangle in degrees
static std::pair<double, double> triangle_angle_extents_deg(const Point &a, const Point &b, const Point &c)
{
    const Vector3 ab = b - a;
    const Vector3 ac = c - a;
    const Vector3 bc = c - b;

    const double ab_sq = ab.squared_length();
    const double ac_sq = ac.squared_length();
    const double bc_sq = bc.squared_length();
    const double eps = std::numeric_limits<double>::epsilon();
    if (ab_sq <= eps || ac_sq <= eps || bc_sq <= eps)
    {
        return {0.0, 0.0};
    }

    const double ab_len = std::sqrt(ab_sq);
    const double ac_len = std::sqrt(ac_sq);
    const double bc_len = std::sqrt(bc_sq);

    const double dot_ab_ac = ab * ac;
    const double dot_ab_bc = ab * bc;
    const double dot_ac_bc = ac * bc;

    auto angle_deg = [](double dot, double len1, double len2) -> double
    {
        double cos_theta = dot / (len1 * len2);
        cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
        return std::acos(cos_theta) * 180.0 / kPi;
    };

    const double angle_a = angle_deg(dot_ab_ac, ab_len, ac_len);
    const double angle_b = angle_deg(-dot_ab_bc, ab_len, bc_len);
    const double angle_c = angle_deg(dot_ac_bc, ac_len, bc_len);

    const double min_angle = std::min({angle_a, angle_b, angle_c});
    const double max_angle = std::max({angle_a, angle_b, angle_c});
    return {min_angle, max_angle};
}

// Helper: Check if a facet has a dummy vertex
static bool facet_has_dummy_vertex(const Facet &facet)
{
    const Cell_handle &cell = facet.first;
    const int facet_index = facet.second;
    for (int k = 0; k < 3; ++k)
    {
        const int cell_vertex_index = CellInfo::FacetVertexIndex(facet_index, k);
        if (cell->vertex(cell_vertex_index)->info().is_dummy)
        {
            return true;
        }
    }
    return false;
}

// Helper: Get dual segment from facet
static bool dual_as_segment(const Delaunay &dt, const Facet &facet, Segment3 &out_segment)
{
    Object dual_obj = dt.dual(facet);
    return CGAL::assign(out_segment, dual_obj);
}

// Helper: Check if segment is bipolar
static bool is_bipolar_segment(const Segment3 &seg, const UnifiedGrid &grid, float isovalue)
{
    const Point &p0 = seg.source();
    const Point &p1 = seg.target();
    // Use trilinear interpolation to get values at endpoints
    float v0 = trilinear_interpolate(adjust_outside_bound_points(p0, grid, p0, p1), grid);
    float v1 = trilinear_interpolate(adjust_outside_bound_points(p1, grid, p0, p1), grid);
    return is_bipolar(v0, v1, isovalue);
}

// Helper: Find nearest active cube
static const Cube *find_nearest_active_cube(const std::vector<Cube> &cubes, const Point &p)
{
    const Cube *best = nullptr;
    double best_dist_sq = std::numeric_limits<double>::max();
    for (const Cube &cube : cubes)
    {
        const double dist_sq = CGAL::squared_distance(cube.cubeCenter, p);
        if (dist_sq < best_dist_sq)
        {
            best_dist_sq = dist_sq;
            best = &cube;
        }
    }
    return best;
}

// Helper: Snap to subcell center
static Point snap_to_subcell_center(const Cube &cube, const UnifiedGrid &grid, const Point &target, int resolution)
{
    const int res = std::max(1, std::min(3, resolution));
    const double dx = grid.spacing[0];
    const double dy = grid.spacing[1];
    const double dz = grid.spacing[2];
    const double base_x = cube.indices[0] * dx + grid.min_coord[0];
    const double base_y = cube.indices[1] * dy + grid.min_coord[1];
    const double base_z = cube.indices[2] * dz + grid.min_coord[2];

    const auto coord_center = [&](double base, double d, double value) -> double
    {
        const double local = (value - base) / d;
        const double scaled = local * res;
        int idx = static_cast<int>(std::floor(scaled));
        idx = std::max(0, std::min(res - 1, idx));
        const double cell_size = d / static_cast<double>(res);
        return base + (static_cast<double>(idx) + 0.5) * cell_size;
    };

    return Point(coord_center(base_x, dx, target.x()),
                 coord_center(base_y, dy, target.y()),
                 coord_center(base_z, dz, target.z()));
}

// Helper: Get next vertex index
static int next_vertex_index(const Delaunay &dt)
{
    int next_index = 0;
    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit)
    {
        next_index = std::max(next_index, vit->info().index + 1);
    }
    return next_index;
}

// Helper: Reindex cells
static void reindex_cells(Delaunay &dt)
{
    int cell_index = 0;
    for (auto cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit)
    {
        cit->info().index = cell_index++;
        cit->info().dualVoronoiVertexIndex = -1;
    }
}

RefinementStats refine_delaunay(Delaunay &dt,
                                const UnifiedGrid &grid,
                                const std::vector<Cube> &active_cubes,
                                const VdcParam &params)
{
    RefinementStats stats;
    if (!params.refine_small_angles || active_cubes.empty())
    {
        return stats;
    }

    bool use_min_angle = params.refine_min_angle_enabled;
    bool use_max_angle = params.refine_max_angle_enabled;
    double min_angle_threshold = params.refine_min_surface_angle_deg;
    double max_angle_threshold = params.refine_max_surface_angle_deg;

    if (use_min_angle && min_angle_threshold <= 0.0)
    {
        min_angle_threshold = 20.0;
    }
    if (use_max_angle && max_angle_threshold <= 0.0)
    {
        max_angle_threshold = 120.0;
    }
    if (!use_min_angle && !use_max_angle)
    {
        use_max_angle = true;
        if (max_angle_threshold <= 0.0)
        {
            max_angle_threshold = 120.0;
        }
    }
    const int insert_resolution = params.refine_insert_resolution;
    const float isovalue = params.isovalue;

    const int max_iterations = 10;

    int vertex_index = next_vertex_index(dt);

    for (int iter = 0; iter < max_iterations; ++iter)
    {
        std::vector<Point> candidates;
        std::set<Point> unique_candidates; // To avoid duplicates in one pass

        // 1. Scan all finite tetrahedra (via facets)
        for (auto fit = dt.finite_facets_begin(); fit != dt.finite_facets_end(); ++fit)
        {
            const Facet &facet = *fit;
            
            // Skip facets with dummy vertices
            if (facet_has_dummy_vertex(facet))
            {
                continue;
            }

            // Get facet vertices
            const Cell_handle &cell = facet.first;
            const int facet_index = facet.second;
            const Point p0 = cell->vertex(CellInfo::FacetVertexIndex(facet_index, 0))->point();
            const Point p1 = cell->vertex(CellInfo::FacetVertexIndex(facet_index, 1))->point();
            const Point p2 = cell->vertex(CellInfo::FacetVertexIndex(facet_index, 2))->point();

            const auto [min_angle, max_angle] = triangle_angle_extents_deg(p0, p1, p2);
            bool trigger = false;
            if (use_min_angle && use_max_angle)
            {
                trigger = (min_angle < min_angle_threshold) || (max_angle >= max_angle_threshold);
            }
            else if (use_min_angle)
            {
                trigger = (min_angle < min_angle_threshold);
            }
            else if (use_max_angle)
            {
                trigger = (max_angle >= max_angle_threshold);
            }
            if (!trigger)
            {
                continue;
            }

            ++stats.candidate_facets;

            // Compute dual Voronoi edge
            Segment3 dual_segment;
            if (!dual_as_segment(dt, facet, dual_segment))
            {
                
                continue;
            }

            // Check if dual edge is bipolar
            if (!is_bipolar_segment(dual_segment, grid, isovalue))
            {
                continue;
            }
            
            ++stats.bipolar_facets;

            // Find insertion point
            Point centroid = CGAL::midpoint(dual_segment.source(), dual_segment.target());
            const Cube *nearest_cube = find_nearest_active_cube(active_cubes, centroid);

            if (nearest_cube)
            {
                Point p_ref = snap_to_subcell_center(*nearest_cube, grid, centroid, insert_resolution);
                
                // Avoid inserting duplicates in the same batch
                if (unique_candidates.find(p_ref) == unique_candidates.end())
                {
                    unique_candidates.insert(p_ref);
                    candidates.push_back(p_ref);
                }
            }
        }

        if (candidates.empty())
        {
            break;
        }

        // 2. Insert points
        std::size_t inserted_this_iter = 0;
        for (const auto &p : candidates)
        {
            Vertex_handle vh = dt.insert(p);

            vh->info().index = vertex_index++;
            vh->info().is_dummy = false;
            vh->info().voronoiCellIndex = -1;
            
            ++inserted_this_iter;
            ++stats.inserted_points;
        }
        
        stats.iterations_run++;
        
        if (inserted_this_iter == 0)
        {
            break;
        }
    }

    reindex_cells(dt);
    return stats;
}
