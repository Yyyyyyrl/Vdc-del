//! @file vdc_refinement.cpp
//! @brief Implementation of facet-centric Delaunay refinement.

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <utility>
#include <vector>

#include "processing/vdc_refinement.h"
#include "core/vdc_utilities.h"

constexpr double kPi = 3.14159265358979323846;

// Helper: Compute minimum and maximum corner angles of a triangle in degrees.
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

static bool cell_has_dummy_vertex(const Cell_handle &cell)
{
    for (int i = 0; i < 4; ++i)
    {
        if (cell->vertex(i)->info().is_dummy)
        {
            return true;
        }
    }
    return false;
}

static bool is_bipolar_segment(const Segment3 &seg, const UnifiedGrid &grid, float isovalue)
{
    const Point &p0 = seg.source();
    const Point &p1 = seg.target();
    const float v0 = trilinear_interpolate(adjust_outside_bound_points(p0, grid, p0, p1), grid);
    const float v1 = trilinear_interpolate(adjust_outside_bound_points(p1, grid, p0, p1), grid);
    return is_bipolar(v0, v1, isovalue);
}

static int next_vertex_index(const Delaunay &dt)
{
    int next_index = 0;
    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit)
    {
        next_index = std::max(next_index, vit->info().index + 1);
    }
    return next_index;
}

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

    if (!use_min_angle && !use_max_angle)
    {
        use_min_angle = true;
    }
    if (use_min_angle && min_angle_threshold <= 0.0)
    {
        min_angle_threshold = 20.0;
    }
    if (use_max_angle && max_angle_threshold <= 0.0)
    {
        max_angle_threshold = 120.0;
    }

    const float isovalue = params.isovalue;
    const int max_iterations = 10;
    const std::size_t max_insert_per_iteration = 500;
    int vertex_index = next_vertex_index(dt);

    for (int iter = 0; iter < max_iterations; ++iter)
    {
        std::vector<std::pair<double, Point>> candidates;

        for (auto fit = dt.finite_facets_begin(); fit != dt.finite_facets_end(); ++fit)
        {
            const Facet &facet = *fit;
            if (facet_has_dummy_vertex(facet))
            {
                continue;
            }

            const Cell_handle &cell_a = facet.first;
            const int facet_index = facet.second;
            const Cell_handle cell_b = cell_a->neighbor(facet_index);
            if (dt.is_infinite(cell_b))
            {
                continue;
            }
            if (cell_has_dummy_vertex(cell_a) || cell_has_dummy_vertex(cell_b))
            {
                continue;
            }

            const Vertex_handle v0 = cell_a->vertex(CellInfo::FacetVertexIndex(facet_index, 0));
            const Vertex_handle v1 = cell_a->vertex(CellInfo::FacetVertexIndex(facet_index, 1));
            const Vertex_handle v2 = cell_a->vertex(CellInfo::FacetVertexIndex(facet_index, 2));

            // Angle proxy: use per-site isosurface samples when available.
            // This targets small angles in the eventual isosurface triangles while still
            // applying refinement as pure Delaunay point insertions (tetrahedron circumcenters).
            const Point p0 = v0->info().has_isov_sample ? v0->info().isov : v0->point();
            const Point p1 = v1->info().has_isov_sample ? v1->info().isov : v1->point();
            const Point p2 = v2->info().has_isov_sample ? v2->info().isov : v2->point();

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

            const Point circumcenter_a = cell_a->circumcenter();
            const Point circumcenter_b = cell_b->circumcenter();
            const Segment3 dual_segment(circumcenter_a, circumcenter_b);
            if (!is_bipolar_segment(dual_segment, grid, isovalue))
            {
                continue;
            }
            ++stats.bipolar_facets;

            double severity = min_angle;
            if (!use_min_angle && use_max_angle)
            {
                severity = 180.0 - max_angle;
            }
            else if (use_min_angle && use_max_angle)
            {
                severity = std::min(min_angle, 180.0 - max_angle);
            }

            if (is_point_inside_grid(circumcenter_a, grid))
            {
                candidates.emplace_back(severity, circumcenter_a);
            }
            if (is_point_inside_grid(circumcenter_b, grid))
            {
                candidates.emplace_back(severity, circumcenter_b);
            }
        }

        if (candidates.empty())
        {
            break;
        }

        std::sort(candidates.begin(), candidates.end(),
                  [](const auto &a, const auto &b) { return a.first < b.first; });

        std::set<Point> unique_candidates;
        std::size_t inserted_this_iter = 0;

        for (const auto &[severity, p] : candidates)
        {
            (void)severity;
            if (unique_candidates.size() >= max_insert_per_iteration)
            {
                break;
            }
            if (!unique_candidates.insert(p).second)
            {
                continue;
            }

            Vertex_handle vh = dt.insert(p);
            if (vh == Vertex_handle())
            {
                continue;
            }
            if (vh->info().index >= 0)
            {
                continue;
            }

            vh->info().index = vertex_index++;
            vh->info().is_dummy = false;
            vh->info().voronoiCellIndex = -1;
            vh->info().has_isov_sample = false;

            ++inserted_this_iter;
            ++stats.inserted_points;
        }

        if (inserted_this_iter == 0)
        {
            break;
        }
        ++stats.iterations_run;
    }

    reindex_cells(dt);
    return stats;
}
