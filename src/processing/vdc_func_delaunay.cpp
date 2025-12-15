#include "processing/vdc_func.h"
#include "core/vdc_timing.h"

//! @brief Adds dummy points from a facet for Voronoi diagram bounding.
std::vector<Point> add_dummy_from_facet(const GridFacets &facet,
                                        const UnifiedGrid &data_grid,
                                        double supersample_multiplier)
{
    std::vector<Point> points;

    // 2D slice dimension
    int dim0 = facet.axis_size[0];
    int dim1 = facet.axis_size[1];

    // For convenience
    int d = facet.orth_dir;
    int d1 = facet.axis_dir[0];
    int d2 = facet.axis_dir[1];

    // and localSize[] = (maxIndex[i] - minIndex[i] + 1)
    // The grid spacing in each dimension
    double dx[3] = {data_grid.spacing[0], data_grid.spacing[1], data_grid.spacing[2]};

    // Loop over the 2D slice
    for (int coord1 = 0; coord1 < dim1; coord1++)
    {
        for (int coord0 = 0; coord0 < dim0; coord0++)
        {
            if (!facet.CubeFlag(coord0, coord1))
                continue;

            // localX[d1] = coord0, localX[d2] = coord1
            int localX[3] = {0, 0, 0};
            localX[d1] = coord0;
            localX[d2] = coord1;

            // side=0 => localX[d] = 0, side=1 => localX[d] = localSize[d]-1
            localX[d] = (facet.side == 0) ? 0 : (facet.localSize[d] - 1);

            // Convert localX -> global indices
            int g[3];
            for (int i = 0; i < 3; i++)
            {
                g[i] = localX[i] + facet.minIndex[i];
            }

            // Compute center in real-world coordinates
            double cx = (g[0] + 0.5) * dx[0];
            double cy = (g[1] + 0.5) * dx[1];
            double cz = (g[2] + 0.5) * dx[2];

            // Offset by +/- dx[d]
            // The offset multiplier is for avoid bipolar edges touching dummy voronoi cells
            const double offsetMultiplier = supersample_multiplier > 0.0 ? supersample_multiplier : 1.0;
            double offset = ((facet.side == 0) ? -4 * dx[d] : 4 * dx[d]) * offsetMultiplier;
            if (d == 0)
                cx += offset;
            else if (d == 1)
                cy += offset;
            else
                cz += offset;

            points.emplace_back(cx, cy, cz);
        }
    }

    return points;
}

//! @brief Collects points for the Delaunay triangulation.
/*!
 * Gathers original points and dummy points from grid facets for multi-isovertex mode,
 * or uses only active cube centers for single-isovertex mode.
 *
 * @param grid The grid containing data.
 * @param grid_facets The grid facets for dummy point generation.
 * @param activeCubeCenters The list of cube centers of active cubes.
 * @param vdc_param The VdcParam instance containing user input options.
 * @param delaunay_points Output vector for all points (original + dummy).
 */
static int collect_delaunay_points(UnifiedGrid &grid,
                                   const std::vector<std::vector<GridFacets>> &grid_facets,
                                   const std::vector<Point> &activeCubeCenters,
                                   VdcParam &vdc_param,
                                   std::vector<Point> &delaunay_points)
{
    delaunay_points = activeCubeCenters;
    int first_dummy_index = delaunay_points.size(); // Dummies start here

    const double supersampleMultiplier = vdc_param.supersample
                                             ? static_cast<double>(vdc_param.supersample_r)
                                             : 1.0;

    if (vdc_param.multi_isov)
    {
        for (int d = 0; d < 3; ++d)
        { // Assuming 3 dimensions
            for (const auto &f : grid_facets[d])
            {
                auto pointsf = add_dummy_from_facet(f, grid, supersampleMultiplier);
                delaunay_points.insert(delaunay_points.end(), pointsf.begin(), pointsf.end());
            }
        }
    }
    return first_dummy_index;
}

//! @brief Inserts points into the Delaunay triangulation.
/*!
 * Inserts the collected points into the triangulation and writes debug output if enabled.
 *
 * @param dt The Delaunay triangulation to insert points into.
 * @param delaunay_points The points to insert.
 * @param activeCubeCenters The list of center points of active cubes.
 * @param vdc_param The VdcParam instance containing user input options.
 */
static Vertex_handle insert_point_into_delaunay_triangulation(Delaunay &dt,
                                                              const Point &p,
                                                              int index,
                                                              bool is_dummy)
{
    // Insert point and retrieve handle
    Vertex_handle vh = dt.insert(p);
    // Immediately assign index and dummy status
    vh->info().index = index;
    vh->info().is_dummy = is_dummy;
    return vh;
}

//! @brief Constructs a Delaunay triangulation from a grid and grid facets.
/*!
 * Constructs a 3D Delaunay triangulation using the grid's scalar values
 * and the facets of the active cubes.
 *
 * @param dt The Delaunay triangulation instance.
 * @param grid The grid containing scalar values.
 * @param grid_facets The grid facets to use in constructing the triangulation.
 * @param vdc_param The VdcParam instance holding user input options.
 * @param activeCubeCenters The list of cube centers of active cubes.
 */
void construct_delaunay_triangulation(Delaunay &dt, UnifiedGrid &grid, const std::vector<std::vector<GridFacets>> &grid_facets, VdcParam &vdc_param, std::vector<Point> &activeCubeCenters)
{
    TimingStats& timer = TimingStats::getInstance();

    timer.startTimer("Collect Delaunay points", "3. Delaunay Triangulation Construction");
    std::vector<Point> delaunay_points;
    size_t first_dummy_index = collect_delaunay_points(grid, grid_facets, activeCubeCenters, vdc_param, delaunay_points);
    timer.stopTimer("Collect Delaunay points", "3. Delaunay Triangulation Construction");

    dt.clear();

    // Batch insert all points
    timer.startTimer("Insert vertices", "3. Delaunay Triangulation Construction");
    dt.insert(delaunay_points.begin(), delaunay_points.end());
    timer.stopTimer("Insert vertices", "3. Delaunay Triangulation Construction");

    // Create map from point to original index
    timer.startTimer("Assign vertex info", "3. Delaunay Triangulation Construction");
    std::map<Point, size_t> point_to_index;
    for (size_t i = 0; i < delaunay_points.size(); ++i)
    {
        point_to_index[delaunay_points[i]] = i;
    }

    // Assign info to vertices
    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit)
    {
        const Point &p = vit->point();
        auto it = point_to_index.find(p);
        if (it == point_to_index.end())
        {
            std::cerr << "[ERROR] Vertex point not found in original points!" << std::endl;
            continue;
        }
        size_t original_index = it->second;
        vit->info().index = original_index;
        vit->info().is_dummy = (original_index >= first_dummy_index);
        vit->info().voronoiCellIndex = -1; // Initialize if needed
    }
    timer.stopTimer("Assign vertex info", "3. Delaunay Triangulation Construction");

    // Assign indices to cells for debugging and tracking
    timer.startTimer("Assign cell indices", "3. Delaunay Triangulation Construction");
    int cellIndex = 0;
    for (auto cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit)
    {
        cit->info().index = cellIndex++;
    }
    timer.stopTimer("Assign cell indices", "3. Delaunay Triangulation Construction");
}
