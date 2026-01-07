#include "processing/vdc_grid.h"
#include "core/vdc_timing.h"
#include "core/vdc_debug.h"
#include <algorithm>
#include <cmath>
#include <limits>


//! Constructor for UnifiedGrid
UnifiedGrid::UnifiedGrid(int nx, int ny, int nz, float dx, float dy, float dz, float min_x, float min_y, float min_z)
    : num_cells{nx, ny, nz},
      spacing{dx, dy, dz},
      physical_spacing{dx, dy, dz},
      min_coord{min_x, min_y, min_z},
      max_coord{0.0f, 0.0f, 0.0f}
{
    update_bounds();
    flat_data.resize(num_cells[0] * num_cells[1] * num_cells[2], 0.0f);
    data.resize(num_cells[0], std::vector<std::vector<float>>(num_cells[1], std::vector<float>(num_cells[2], 0.0f)));
}

// Retrieve a scalar value from the grid
float UnifiedGrid::get_value(int x, int y, int z) const
{
    if (x < 0 || x >= num_cells[0] || y < 0 || y >= num_cells[1] || z < 0 || z >= num_cells[2])
        return 0.0f;
    return data[x][y][z];
}

// Set a scalar value in the grid
void UnifiedGrid::set_value(int x, int y, int z, float value)
{
    if (x >= 0 && x < num_cells[0] && y >= 0 && y < num_cells[1] && z >= 0 && z < num_cells[2])
    {
        data[x][y][z] = value;
        flat_data[z * num_cells[0] * num_cells[1] + y * num_cells[0] + x] = value;
    }
}


// Convert point to grid index
std::tuple<int, int, int> UnifiedGrid::point_to_grid_index(const Point &point) const
{
    int x = static_cast<int>((point.x() - min_coord[0]) / spacing[0]);
    int y = static_cast<int>((point.y() - min_coord[1]) / spacing[1]);
    int z = static_cast<int>((point.z() - min_coord[2]) / spacing[2]);
    return {x, y, z};
}

// Get scalar value at a point using trilinear interpolation
float UnifiedGrid::get_scalar_value_at_point(const Point &point) const
{
    return trilinear_interpolate(point, *this);
}


//! @brief Constructor for the GridFacets structure.
GridFacets::GridFacets(int d, int s, const int minIdx[DIM3], const int maxIdx[DIM3])
    : orth_dir(d), side(s)
{
    //! Store `minIndex[]` and `maxIndex[]`, and compute `localSize[]`.
    for (int i = 0; i < DIM3; i++) {
        minIndex[i] = minIdx[i];
        maxIndex[i] = maxIdx[i];
        localSize[i] = (maxIndex[i] - minIndex[i] + 1);
    }

    //! Determine the 2D slice axes.
    axis_dir[0] = (orth_dir + 1) % DIM3;
    axis_dir[1] = (orth_dir + 2) % DIM3;

    //! Compute the size of the facet along the two axes.
    axis_size[0] = localSize[axis_dir[0]];
    axis_size[1] = localSize[axis_dir[1]];

    //! Allocate flags for the facet, initialized to `false`.
    cube_flag.resize(axis_size[0] * axis_size[1], false);
}

//! @brief Set the flag for a particular `(coord0, coord1)` in the facet.
void GridFacets::SetFlag(int coord0, int coord1, bool flag)
{
    cube_flag[index(coord0, coord1)] = flag;
}

//! @brief Get the flag for a particular `(coord0, coord1)` in the facet.
bool GridFacets::CubeFlag(int coord0, int coord1) const
{
    return cube_flag[index(coord0, coord1)];
}

//! @brief Convert a 2D coordinate `(coord0, coord1)` to a linear index.
int GridFacets::index(int coord0, int coord1) const
{
    return coord1 * axis_size[0] + coord0;
}

// Print grid metadata and data
void UnifiedGrid::print_grid() const
{
    const int nx = num_cells[0];
    const int ny = num_cells[1];
    const int nz = num_cells[2];
    const float dx = spacing[0];
    const float dy = spacing[1];
    const float dz = spacing[2];
    std::cout << "Unified Grid Information:\n";
    std::cout << "Dimensions: " << nx << "x" << ny << "x" << nz << "\n";
    std::cout << "Internal spacing (grid units): dx=" << dx << ", dy=" << dy << ", dz=" << dz << "\n";
    std::cout << "Physical spacing: dx=" << physical_spacing[0] << ", dy=" << physical_spacing[1] << ", dz=" << physical_spacing[2] << "\n";
    std::cout << "Bounds (grid units): [" << min_coord[0] << ", " << max_coord[0] << "] x [" << min_coord[1] << ", " << max_coord[1] << "] x [" << min_coord[2] << ", " << max_coord[2] << "]\n";
    const float phys_max_x = min_coord[0] + (nx - 1) * physical_spacing[0];
    const float phys_max_y = min_coord[1] + (ny - 1) * physical_spacing[1];
    const float phys_max_z = min_coord[2] + (nz - 1) * physical_spacing[2];
    std::cout << "Bounds (physical): [" << min_coord[0] << ", " << phys_max_x << "] x [" << min_coord[1] << ", " << phys_max_y << "] x [" << min_coord[2] << ", " << phys_max_z << "]\n\n";
    std::cout << "Data:\n";
    for (int z = 0; z < nz; ++z)
    {
        std::cout << "Slice z = " << z << ":\n";
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
                std::cout << std::setw(8) << data[x][y][z] << " ";
            std::cout << "\n";
        }
        std::cout << "\n";
    }
}

void UnifiedGrid::update_bounds()
{
    for (int d = 0; d < DIM3; ++d)
    {
        max_coord[d] = min_coord[d] + (num_cells[d] - 1) * spacing[d];
    }
}

void UnifiedGrid::force_unit_spacing()
{
    spacing[0] = spacing[1] = spacing[2] = 1.0f;
    update_bounds();
}

void UnifiedGrid::zero_boundary_shell()
{
    const int nx = num_cells[0];
    const int ny = num_cells[1];
    const int nz = num_cells[2];
    if (nx <= 0 || ny <= 0 || nz <= 0)
        return;

    const float fill_value = flat_data.empty()
                                 ? 0.0f
                                 : *std::min_element(flat_data.begin(), flat_data.end());

    auto assign = [&](int x, int y, int z) {
        if (x < 0 || x >= nx || y < 0 || y >= ny || z < 0 || z >= nz)
            return;
        data[x][y][z] = fill_value;
        flat_data[z * nx * ny + y * nx + x] = fill_value;
    };

    for (int x = 0; x < nx; ++x)
    {
        for (int y = 0; y < ny; ++y)
        {
            assign(x, y, 0);
            assign(x, y, nz - 1);
        }
    }

    for (int x = 0; x < nx; ++x)
    {
        for (int z = 0; z < nz; ++z)
        {
            assign(x, 0, z);
            assign(x, ny - 1, z);
        }
    }

    for (int y = 0; y < ny; ++y)
    {
        for (int z = 0; z < nz; ++z)
        {
            assign(0, y, z);
            assign(nx - 1, y, z);
        }
    }
}

bool UnifiedGrid::boundary_crosses_isovalue(float isovalue) const
{
    const int nx = num_cells[0];
    const int ny = num_cells[1];
    const int nz = num_cells[2];
    if (nx <= 0 || ny <= 0 || nz <= 0)
        return false;

    float boundary_min = std::numeric_limits<float>::infinity();
    float boundary_max = -std::numeric_limits<float>::infinity();

    auto consider = [&](int x, int y, int z)
    {
        if (x < 0 || x >= nx || y < 0 || y >= ny || z < 0 || z >= nz)
            return;
        float val = data[x][y][z];
        boundary_min = std::min(boundary_min, val);
        boundary_max = std::max(boundary_max, val);
    };

    for (int x = 0; x < nx; ++x)
    {
        for (int y = 0; y < ny; ++y)
        {
            consider(x, y, 0);
            consider(x, y, nz - 1);
        }
    }

    for (int x = 0; x < nx; ++x)
    {
        for (int z = 1; z < nz - 1; ++z)
        {
            consider(x, 0, z);
            consider(x, ny - 1, z);
        }
    }

    for (int y = 1; y < ny - 1; ++y)
    {
        for (int z = 1; z < nz - 1; ++z)
        {
            consider(0, y, z);
            consider(nx - 1, y, z);
        }
    }

    if (!std::isfinite(boundary_min) || !std::isfinite(boundary_max))
        return false;

    return (boundary_min < isovalue) && (isovalue < boundary_max);
}


//! @brief Converts raw data of type `T` into a float vector.
/*!
 * @tparam T The type of the input data.
 * @param data_ptr Pointer to the input data.
 * @param total_size Total number of elements in the input data.
 * @return A `std::vector<float>` containing the converted data.
 */
// Convert raw data to float vector
template <typename T>
std::vector<float> convert_to_float_vector(T *data_ptr, size_t total_size)
{
    std::vector<float> data(total_size);
    for (size_t i = 0; i < total_size; ++i)
        data[i] = static_cast<float>(data_ptr[i]);
    return data;
}

// Load NRRD data
UnifiedGrid load_nrrd_data(const std::string &file_path)
{
    TimingStats& timer = TimingStats::getInstance();
    timer.startTimer("Load NRRD file", "1. Load Data and Grid Formation");

    Nrrd *nrrd = nrrdNew();
    if (nrrdLoad(nrrd, file_path.c_str(), NULL))
    {
        char *err = biffGetDone(NRRD);
        std::cerr << "Error reading NRRD file: " << err << std::endl;
        free(err);
        nrrdNuke(nrrd);
        exit(1);
    }

    size_t total_size = nrrdElementNumber(nrrd);
    int nx = nrrd->axis[0].size;
    int ny = nrrd->axis[1].size;
    int nz = nrrd->axis[2].size;
    auto sanitize_spacing = [](double spacing) -> float {
        if (!std::isfinite(spacing) || spacing <= 0.0)
            return 1.0f;
        return static_cast<float>(spacing);
    };

    float dx = sanitize_spacing(nrrd->axis[0].spacing);
    float dy = sanitize_spacing(nrrd->axis[1].spacing);
    float dz = sanitize_spacing(nrrd->axis[2].spacing);
    float min_x = 0.0f, min_y = 0.0f, min_z = 0.0f;
    timer.stopTimer("Load NRRD file", "1. Load Data and Grid Formation");

    timer.startTimer("Grid initialization", "1. Load Data and Grid Formation");
    UnifiedGrid grid(nx, ny, nz, dx, dy, dz, min_x, min_y, min_z);

    if (nrrd->type == nrrdTypeFloat)
    {
        float *data_ptr = static_cast<float *>(nrrd->data);
        grid.flat_data = std::vector<float>(data_ptr, data_ptr + total_size);
    }
    else if (nrrd->type == nrrdTypeUChar)
    {
        unsigned char *data_ptr = static_cast<unsigned char *>(nrrd->data);
        grid.flat_data = convert_to_float_vector(data_ptr, total_size);
    }
    else
    {
        std::cerr << "Unsupported NRRD data type." << std::endl;
        nrrdNuke(nrrd);
        exit(1);
    }

    // Populate 3D data array
    for (int z = 0; z < nz; ++z)
        for (int y = 0; y < ny; ++y)
            for (int x = 0; x < nx; ++x)
                grid.data[x][y][z] = grid.flat_data[z * nx * ny + y * nx + x];

    nrrdNuke(nrrd);

    grid.force_unit_spacing();
    timer.stopTimer("Grid initialization", "1. Load Data and Grid Formation");

    const int grid_nx = grid.num_cells[0];
    const int grid_ny = grid.num_cells[1];
    const int grid_nz = grid.num_cells[2];
    const float phys_max_x = grid.min_coord[0] + (grid_nx - 1) * grid.physical_spacing[0];
    const float phys_max_y = grid.min_coord[1] + (grid_ny - 1) * grid.physical_spacing[1];
    const float phys_max_z = grid.min_coord[2] + (grid_nz - 1) * grid.physical_spacing[2];

    if (indicator) {
        std::cout << "Grid dimensions: " << grid_nx << "x" << grid_ny << "x" << grid_nz << "\n";
        std::cout << "Physical spacing: dx=" << grid.physical_spacing[0]
                  << ", dy=" << grid.physical_spacing[1]
                  << ", dz=" << grid.physical_spacing[2] << "\n";
        std::cout << "Internal spacing (grid units): dx=" << grid.spacing[0]
                  << ", dy=" << grid.spacing[1]
                  << ", dz=" << grid.spacing[2] << "\n";
        std::cout << "Bounds (grid units): [" << grid.min_coord[0] << ", " << grid.max_coord[0] << "] x ["
                  << grid.min_coord[1] << ", " << grid.max_coord[1] << "] x ["
                  << grid.min_coord[2] << ", " << grid.max_coord[2] << "]\n";
        std::cout << "Bounds (physical): [" << grid.min_coord[0] << ", " << phys_max_x << "] x ["
                  << grid.min_coord[1] << ", " << phys_max_y << "] x ["
                  << grid.min_coord[2] << ", " << phys_max_z << "]\n";
    }

    return grid;
}

// Supersample grid
UnifiedGrid supersample_grid(const UnifiedGrid &grid, int n)
{
    TimingStats& timer = TimingStats::getInstance();
    timer.startTimer("Supersample", "1. Load Data and Grid Formation");

    const int grid_nx = grid.num_cells[0];
    const int grid_ny = grid.num_cells[1];
    const int grid_nz = grid.num_cells[2];
    int nx2 = grid_nx * n - (n - 1);
    int ny2 = grid_ny * n - (n - 1);
    int nz2 = grid_nz * n - (n - 1);
    float dx2 = grid.physical_spacing[0] / n;
    float dy2 = grid.physical_spacing[1] / n;
    float dz2 = grid.physical_spacing[2] / n;

    UnifiedGrid new_grid(nx2, ny2, nz2, dx2, dy2, dz2, grid.min_coord[0], grid.min_coord[1], grid.min_coord[2]);

    for (int z = 0; z < nz2; ++z)
    {
        for (int y = 0; y < ny2; ++y)
        {
            for (int x = 0; x < nx2; ++x)
            {
                float px = grid.min_coord[0] + (static_cast<float>(x) / n) * grid.spacing[0];
                float py = grid.min_coord[1] + (static_cast<float>(y) / n) * grid.spacing[1];
                float pz = grid.min_coord[2] + (static_cast<float>(z) / n) * grid.spacing[2];
                float value = trilinear_interpolate(Point(px, py, pz), grid);
                new_grid.set_value(x, y, z, value);
            }
        }
    }

    new_grid.force_unit_spacing();
    timer.stopTimer("Supersample", "1. Load Data and Grid Formation");

    return new_grid;
}


// Check if cube is active
bool is_cube_active(const UnifiedGrid &grid, int x, int y, int z, float isovalue)
{
    static const std::vector<std::tuple<int,int,int>> vertex_offsets = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
        {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}
    };

    auto get_value = [&](int dx, int dy, int dz) -> float {
        return grid.get_value(x + dx, y + dy, z + dz);
    };

    float val0 = get_value(0, 0, 0);
    bool is_val0_negative = (val0 < isovalue);

    for (int i = 1; i < 8; i++) {
        auto [dx, dy, dz] = vertex_offsets[i];
        float vali = get_value(dx, dy, dz);
        bool is_vali_negative = (vali < isovalue);
        if (is_vali_negative != is_val0_negative)
            return true;
    }

    return false;
}

// Adjust points outside grid bounds
Point adjust_outside_bound_points(const Point &p, const UnifiedGrid &grid, const Point &v1, const Point &v2)
{
    const float min_x = grid.min_coord[0];
    const float min_y = grid.min_coord[1];
    const float min_z = grid.min_coord[2];
    const float dx = grid.spacing[0];
    const float dy = grid.spacing[1];
    const float dz = grid.spacing[2];
    const int nx = grid.num_cells[0];
    const int ny = grid.num_cells[1];
    const int nz = grid.num_cells[2];

    float gx = (p.x() - min_x) / dx;
    float gy = (p.y() - min_y) / dy;
    float gz = (p.z() - min_z) / dz;

    if (gx < 0 || gx >= nx || gy < 0 || gy >= ny || gz < 0 || gz >= nz)
    {
        float v1_gx = (v1.x() - min_x) / dx;
        float v1_gy = (v1.y() - min_y) / dy;
        float v1_gz = (v1.z() - min_z) / dz;

        float v2_gx = (v2.x() - min_x) / dx;
        float v2_gy = (v2.y() - min_y) / dy;
        float v2_gz = (v2.z() - min_z) / dz;

        float t = ((gx - v1_gx) * (v2_gx - v1_gx) + (gy - v1_gy) * (v2_gy - v1_gy) + (gz - v1_gz) * (v2_gz - v1_gz)) /
                  ((v2_gx - v1_gx) * (v2_gx - v1_gx) + (v2_gy - v1_gy) * (v2_gy - v1_gy) + (v2_gz - v1_gz) * (v2_gz - v1_gz));

        t = std::max(0.0f, std::min(t, 1.0f));

        float px = v1_gx + t * (v2_gx - v1_gx);
        float py = v1_gy + t * (v2_gy - v1_gy);
        float pz = v1_gz + t * (v2_gz - v1_gz);

        return Point(px * dx + min_x, py * dy + min_y, pz * dz + min_z);
    }

    return p;
}


// Helper function to compute iso-crossing point in an active cube
// Returns the centroid of edge/isovalue intersections (fallback: cube center).
// Note: The Delaunay site can still be chosen independently (e.g., cube center
// vs iso-crossing) via higher-level options such as `-position_delv_on_isov`.
Point compute_iso_crossing_point(const UnifiedGrid &grid, int i, int j, int k, float isovalue)
{
    static const int cubeVertices[8][3] = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
        {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}
    };
    static const int cubeEdges[12][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0},
        {4, 5}, {5, 6}, {6, 7}, {7, 4},
        {0, 4}, {1, 5}, {2, 6}, {3, 7}
    };

    const float dx = grid.spacing[0];
    const float dy = grid.spacing[1];
    const float dz = grid.spacing[2];
    const float min_x = grid.min_coord[0];
    const float min_y = grid.min_coord[1];
    const float min_z = grid.min_coord[2];

    const float base_x = i * dx + min_x;
    const float base_y = j * dy + min_y;
    const float base_z = k * dz + min_z;

    std::array<Point, 8> vertices;
    std::array<float, 8> scalarValues;
    for (int v = 0; v < 8; ++v)
    {
        vertices[v] = Point(
            base_x + cubeVertices[v][0] * dx,
            base_y + cubeVertices[v][1] * dy,
            base_z + cubeVertices[v][2] * dz);
        scalarValues[v] = grid.get_value(
            i + cubeVertices[v][0],
            j + cubeVertices[v][1],
            k + cubeVertices[v][2]);
    }

    std::vector<Point> intersectionPoints;
    intersectionPoints.reserve(12);
    for (const auto &edge : cubeEdges)
    {
        const int idx1 = edge[0];
        const int idx2 = edge[1];
        const float val1 = scalarValues[idx1];
        const float val2 = scalarValues[idx2];

        if ((val1 < isovalue && val2 >= isovalue) || (val1 >= isovalue && val2 < isovalue))
        {
            intersectionPoints.push_back(
                interpolate(vertices[idx1], vertices[idx2], val1, val2, isovalue, grid));
        }
    }

    if (!intersectionPoints.empty())
    {
        double cx = 0.0, cy = 0.0, cz = 0.0;
        for (const auto &pt : intersectionPoints)
        {
            cx += pt.x();
            cy += pt.y();
            cz += pt.z();
        }
        const double inv = 1.0 / static_cast<double>(intersectionPoints.size());
        return Point(cx * inv, cy * inv, cz * inv);
    }

    return Point(base_x + 0.5f * dx,
                 base_y + 0.5f * dy,
                 base_z + 0.5f * dz);
}

// Find active cubes
void find_active_cubes(const UnifiedGrid &grid, float isovalue, std::vector<Cube> &cubes)
{
    const int nx = grid.num_cells[0];
    const int ny = grid.num_cells[1];
    const int nz = grid.num_cells[2];
    const float dx = grid.spacing[0];
    const float dy = grid.spacing[1];
    const float dz = grid.spacing[2];
    const float min_x = grid.min_coord[0];
    const float min_y = grid.min_coord[1];
    const float min_z = grid.min_coord[2];
    cubes.clear();
    for (int i = 0; i < nx - 1; ++i)
    {
        for (int j = 0; j < ny - 1; ++j)
        {
            for (int k = 0; k < nz - 1; ++k)
            {
                if (is_cube_active(grid, i, j, k, isovalue))
                {
                    Point repVertex(i * dx + min_x, j * dy + min_y, k * dz + min_z);
                    Point cubeCenter((i + 0.5f) * dx + min_x, (j + 0.5f) * dy + min_y, (k + 0.5f) * dz + min_z);
                    Cube cube(repVertex, cubeCenter, i, j, k);
                    cube.accurateIsoCrossing = compute_iso_crossing_point(grid, i, j, k, isovalue);
                    cubes.push_back(cube);
                }
            }
        }
    }
}


// Load grid points
std::vector<Point> load_grid_points(const UnifiedGrid &grid)
{
    const int nx = grid.num_cells[0];
    const int ny = grid.num_cells[1];
    const int nz = grid.num_cells[2];
    const float dx = grid.spacing[0];
    const float dy = grid.spacing[1];
    const float dz = grid.spacing[2];
    const float min_x = grid.min_coord[0];
    const float min_y = grid.min_coord[1];
    const float min_z = grid.min_coord[2];
    std::vector<Point> points;
    points.reserve(static_cast<size_t>(nx) * ny * nz);
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j)
            for (int k = 0; k < nz; ++k)
                points.push_back(Point(i * dx + min_x, j * dy + min_y, k * dz + min_z));
    return points;
}

// Check if point is inside grid
bool is_point_inside_grid(const Point &p, const UnifiedGrid &grid)
{
    return (p.x() >= grid.min_coord[0] && p.x() <= grid.max_coord[0] &&
            p.y() >= grid.min_coord[1] && p.y() <= grid.max_coord[1] &&
            p.z() >= grid.min_coord[2] && p.z() <= grid.max_coord[2]);
}

// Interpolate along an edge
Point interpolate(const Point &p1, const Point &p2, float val1, float val2, float isovalue, const UnifiedGrid &grid)
{
    if (std::abs(val1 - val2) < 1e-6)
        return p1;
    float t = (isovalue - val1) / (val2 - val1);
    return Point(p1.x() + t * (p2.x() - p1.x()),
                p1.y() + t * (p2.y() - p1.y()),
                p1.z() + t * (p2.z() - p1.z()));
}

// Trilinear interpolation
float trilinear_interpolate(const Point &p, const UnifiedGrid &grid)
{
    const float min_x = grid.min_coord[0];
    const float min_y = grid.min_coord[1];
    const float min_z = grid.min_coord[2];
    const float dx = grid.spacing[0];
    const float dy = grid.spacing[1];
    const float dz = grid.spacing[2];
    const int nx = grid.num_cells[0];
    const int ny = grid.num_cells[1];
    const int nz = grid.num_cells[2];

    float gx = (p.x() - min_x) / dx;
    float gy = (p.y() - min_y) / dy;
    float gz = (p.z() - min_z) / dz;
    gx = std::max(0.0f, std::min(gx, static_cast<float>(nx - 1)));
    gy = std::max(0.0f, std::min(gy, static_cast<float>(ny - 1)));
    gz = std::max(0.0f, std::min(gz, static_cast<float>(nz - 1)));

    int x0 = static_cast<int>(std::floor(gx));
    int x1 = std::min(x0 + 1, nx - 1);
    int y0 = static_cast<int>(std::floor(gy));
    int y1 = std::min(y0 + 1, ny - 1);
    int z0 = static_cast<int>(std::floor(gz));
    int z1 = std::min(z0 + 1, nz - 1);

    float xd = gx - x0;
    float yd = gy - y0;
    float zd = gz - z0;

    float c000 = grid.get_value(x0, y0, z0);
    float c001 = grid.get_value(x0, y0, z1);
    float c010 = grid.get_value(x0, y1, z0);
    float c011 = grid.get_value(x0, y1, z1);
    float c100 = grid.get_value(x1, y0, z0);
    float c101 = grid.get_value(x1, y0, z1);
    float c110 = grid.get_value(x1, y1, z0);
    float c111 = grid.get_value(x1, y1, z1);

    float c00 = c000 * (1 - zd) + c001 * zd;
    float c01 = c010 * (1 - zd) + c011 * zd;
    float c10 = c100 * (1 - zd) + c101 * zd;
    float c11 = c110 * (1 - zd) + c111 * zd;

    float c0 = c00 * (1 - yd) + c01 * yd;
    float c1 = c10 * (1 - yd) + c11 * yd;

    return c0 * (1 - xd) + c1 * xd;
}


//! @brief Creates grid facets for active cubes.
std::vector<std::vector<GridFacets>> create_grid_facets(const std::vector<Cube> &activeCubes) {

    int minIdx[3];
    int maxIdx[3];

    minIdx[0] = minIdx[1] = minIdx[2] = INT_MAX;
    maxIdx[0] = maxIdx[1] = maxIdx[2] = INT_MIN;

    for (auto &cube : activeCubes)
    {
        if (cube.indices[0] < minIdx[0])
            minIdx[0] = cube.indices[0];
        if (cube.indices[0] > maxIdx[0])
            maxIdx[0] = cube.indices[0];

        if (cube.indices[1] < minIdx[1])
            minIdx[1] = cube.indices[1];
        if (cube.indices[1] > maxIdx[1])
            maxIdx[1] = cube.indices[1];

        if (cube.indices[2] < minIdx[2])
            minIdx[2] = cube.indices[2];
        if (cube.indices[2] > maxIdx[2])
            maxIdx[2] = cube.indices[2];
    }
    std::vector<std::vector<GridFacets>> grid_facets(3, std::vector<GridFacets>(2,
                                                                                  GridFacets(0, 0, minIdx, maxIdx)));

    // re-construct them properly with the correct (d, side):
    for (int d = 0; d < 3; d++)
    {
        for (int side = 0; side < 2; side++)
        {
            grid_facets[d][side] = GridFacets(d, side, minIdx, maxIdx);
        }
    }

    // Populate them
    for (auto &cube : activeCubes)
    {
        // Global index
        int g[3] = {cube.indices[0], cube.indices[1], cube.indices[2]};

        for (int d = 0; d < 3; d++)
        {
            int d1 = (d + 1) % 3;
            int d2 = (d + 2) % 3;

            for (int side = 0; side < 2; side++)
            {
                GridFacets &f = grid_facets[d][side];

                // Convert to local indices
                // localCoord = g - minIndex
                int coord0 = g[d1] - f.minIndex[d1];
                int coord1 = g[d2] - f.minIndex[d2];

                // Mark it
                f.SetFlag(coord0, coord1, true);
            }
        }
    }

    return grid_facets;
}


// Check if two cubes are adjacent in grid space
bool is_adjacent(const Cube &cubeA, const Cube &cubeB, const UnifiedGrid &grid)
{
    int di = std::abs(cubeA.indices[0] - cubeB.indices[0]);
    int dj = std::abs(cubeA.indices[1] - cubeB.indices[1]);
    int dk = std::abs(cubeA.indices[2] - cubeB.indices[2]);
    return (di <= 1 && dj <= 1 && dk <= 1) && !(di == 0 && dj == 0 && dk == 0);
}

// Calculate unique cube index
int get_cube_index(const Point &repVertex, const UnifiedGrid &grid)
{
    const float min_x = grid.min_coord[0];
    const float min_y = grid.min_coord[1];
    const float min_z = grid.min_coord[2];
    const float dx = grid.spacing[0];
    const float dy = grid.spacing[1];
    const float dz = grid.spacing[2];
    const int nx = grid.num_cells[0];
    const int ny = grid.num_cells[1];
    int i = static_cast<int>((repVertex.x() - min_x) / dx);
    int j = static_cast<int>((repVertex.y() - min_y) / dy);
    int k = static_cast<int>((repVertex.z() - min_z) / dz);
    return k * (nx - 1) * (ny - 1) + j * (nx - 1) + i;
}

// Find neighbor indices
std::vector<int> find_neighbor_indices(const Point &repVertex, const UnifiedGrid &grid)
{
    std::vector<int> neighbors;
    const float min_x = grid.min_coord[0];
    const float min_y = grid.min_coord[1];
    const float min_z = grid.min_coord[2];
    const float dx = grid.spacing[0];
    const float dy = grid.spacing[1];
    const float dz = grid.spacing[2];
    const int nx = grid.num_cells[0];
    const int ny = grid.num_cells[1];
    const int nz = grid.num_cells[2];
    int i = static_cast<int>((repVertex.x() - min_x) / dx);
    int j = static_cast<int>((repVertex.y() - min_y) / dy);
    int k = static_cast<int>((repVertex.z() - min_z) / dz);
    for (int di = -1; di <= 1; ++di)
        for (int dj = -1; dj <= 1; ++dj)
            for (int dk = -1; dk <= 1; ++dk)
                if (di != 0 || dj != 0 || dk != 0)
                {
                    int ni = i + di, nj = j + dj, nk = k + dk;
                    if (ni >= 0 && ni < nx - 1 && nj >= 0 && nj < ny - 1 && nk >= 0 && nk < nz - 1)
                        neighbors.push_back(nk * (nx - 1) * (ny - 1) + nj * (nx - 1) + ni);
                }
    return neighbors;
}

//! @brief Retrieves the centers of a list of cubes.
std::vector<Point> get_cube_centers(const std::vector<Cube> &cubes)
{
    std::vector<Point> cubeCenters;
    for (auto &cube : cubes)
    {
        cubeCenters.push_back(cube.cubeCenter);
    }
    return cubeCenters;
}

//! @brief Retrieves the accurate iso-crossing points of a list of cubes.
std::vector<Point> get_cube_accurate_iso_crossing_points(const std::vector<Cube> &cubes)
{
    std::vector<Point> accurateIsoCrossingPoints;
    for (auto &cube : cubes)
    {
        accurateIsoCrossingPoints.push_back(cube.accurateIsoCrossing);
    }
    return accurateIsoCrossingPoints;
}

// Separation routines have been moved to vdc_sep_isov.cpp
