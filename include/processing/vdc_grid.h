//! @file vdc_grid.h
//! @brief Header file for the grid structure and methods related
#ifndef VDC_GRID_H
#define VDC_GRID_H

#include "core/vdc_type.h"

// DIM = 3 for a 3D grid
static const int DIM3 = 3;

//! @brief Represents a cube in 3D space.
struct Cube {
    Point repVertex;        //!< Representative vertex (world coordinates).
    Point cubeCenter;       //!< Cube center used for Delaunay (provides robustness).
    
    Point accurateIsoCrossing; //!< Accurate iso-crossing point (centroid of edge intersections).
    int indices[DIM3];      //!< Grid indices of the cube: {i, j, k}.
    int isov_subgrid_index; //!< Subgrid index for refined separation grid.

    Cube()
        : repVertex(0, 0, 0),
          cubeCenter(0, 0, 0),
          accurateIsoCrossing(0, 0, 0),
          indices{0, 0, 0},
          isov_subgrid_index(0) {}

    Cube(Point v, Point center, int ix, int iy, int iz)
        : repVertex(v),
          cubeCenter(center),
          accurateIsoCrossing(center),
          indices{ix, iy, iz},
          isov_subgrid_index(0) {}

    //! @brief Print the cube information for debugging
    template <typename OSTREAM_TYPE>
    void Print(OSTREAM_TYPE & out) const {
        out << "Cube:\n";
        out << "  Representative vertex: " << repVertex << "\n";
        out << "  Cube center: " << cubeCenter << "\n";
        out << "  Accurate iso-crossing point: " << accurateIsoCrossing << "\n";
        out << "  Grid indices (i, j, k): (" << indices[0] << ", " << indices[1] << ", " << indices[2] << ")\n";
        out << "  Isov subgrid index: " << isov_subgrid_index << "\n";
    }
};

//! @brief A structure representing a 3D grid for storing scalar data.
/*!
 * Combines the functionality of the previous Grid and ScalarGrid structures.
 * Stores scalar values in both a flattened vector (for NRRD compatibility) and
 * a 3D array (for efficient access), along with grid dimensions, spacings, and bounds.
 */
struct UnifiedGrid
{
    //! @brief 1D vector storing the scalar values for NRRD compatibility.
    std::vector<float> flat_data;

    //! @brief 3D vector storing scalar values for efficient access.
    std::vector<std::vector<std::vector<float>>> data;

    //! @brief Number of grid cells along the x, y, and z axes.
    int num_cells[DIM3];

    //! @brief Spacing of the grid cells along the x, y, and z axes.
    float spacing[DIM3];

    //! @brief Physical spacing of the grid cells (retained when internal spacing is rescaled).
    float physical_spacing[DIM3];

    //! @brief Minimum and maximum coordinates of the grid.
    float min_coord[DIM3], max_coord[DIM3];

    //! @brief Constructor to initialize an empty grid.
    UnifiedGrid() : num_cells{0, 0, 0}, spacing{1.0f, 1.0f, 1.0f},
                    physical_spacing{1.0f, 1.0f, 1.0f},
                    min_coord{0.0f, 0.0f, 0.0f},
                    max_coord{0.0f, 0.0f, 0.0f} {}

    //! @brief Constructor to initialize a grid with specified dimensions and spacings.
    UnifiedGrid(int nx, int ny, int nz, float dx, float dy, float dz, float min_x, float min_y, float min_z);

    //! @brief Get a scalar value at the specified grid index.
    float get_value(int x, int y, int z) const;

    //! @brief Set a scalar value at the specified grid index.
    void set_value(int x, int y, int z, float value);

    //! @brief Get the scalar value at a given point in space using trilinear interpolation.
    float get_scalar_value_at_point(const Point &point) const;

    //! @brief Convert a 3D point to its corresponding grid index.
    std::tuple<int, int, int> point_to_grid_index(const Point &point) const;

    //! @brief Print the grid's metadata and data.
    void print_grid() const;

    //! @brief Recompute max bounds from current spacing.
    void update_bounds();

    //! @brief Force internal spacing to unit while keeping physical spacing metadata.
    void force_unit_spacing();

    //! @brief Clamp all boundary voxels to the global minimum scalar value.
    void zero_boundary_shell();

    //! @brief Detect whether boundary voxel range straddles the given isovalue.
    bool boundary_crosses_isovalue(float isovalue) const;

    //! @brief Print unified grid information for debugging
    template <typename OSTREAM_TYPE>
    void Print(OSTREAM_TYPE & out) const {
        out << "UnifiedGrid:\n";
        out << "  Dimensions (nx, ny, nz): (" << num_cells[0] << ", " << num_cells[1] << ", " << num_cells[2] << ")\n";
        out << "  Spacing (dx, dy, dz): (" << spacing[0] << ", " << spacing[1] << ", " << spacing[2] << ")\n";
        out << "  Physical spacing: (" << physical_spacing[0] << ", " << physical_spacing[1] << ", " << physical_spacing[2] << ")\n";
        out << "  Min coords: (" << min_coord[0] << ", " << min_coord[1] << ", " << min_coord[2] << ")\n";
        out << "  Max coords: (" << max_coord[0] << ", " << max_coord[1] << ", " << max_coord[2] << ")\n";
        out << "  Flat data size: " << flat_data.size() << " element(s)\n";
    }
};



//! @brief A structure representing a 2D "facet" of a 3D grid, orthogonal to one axis.
/*!
 * For example, if `orth_dir = 0` (the x-axis), then this facet stores data 
 * along the y- and z-axes. The `axis_dir[0]` corresponds to the first axis, 
 * and `axis_dir[1]` corresponds to the second axis.
 */
struct GridFacets
{
    //! @brief Orthogonal direction: 0=x, 1=y, 2=z.
    int orth_dir;

    //! @brief Side of the grid facet: 0=lower (min), 1=upper (max).
    int side;

    //! @brief The two axes defining this facet (e.g., y and z if orth_dir = 0).
    int axis_dir[2];

    //! @brief Size of the facet along the two axes.
    /*!
     * - `axis_size[0]`: Size along `axis_dir[0]`.
     * - `axis_size[1]`: Size along `axis_dir[1]`.
     */
    int axis_size[2];

    //! @brief Minimum global indices for the bounding box in each dimension.
    int minIndex[DIM3];

    //! @brief Maximum global indices for the bounding box in each dimension.
    int maxIndex[DIM3];

    //! @brief Local size of the facet in each dimension.
    int localSize[DIM3];

    //! @brief Boolean flags for each `(coord0, coord1)` in the facet.
    /*!
     * The flags are stored in a linearized format for efficient memory usage.
     */
    std::vector<bool> cube_flag;

    //! @brief Constructor to initialize a grid facet.
    /*!
     * @param d Orthogonal direction (0=x, 1=y, 2=z).
     * @param s Side of the grid (0=lower, 1=upper).
     * @param minIdx Global minimum indices for the bounding box.
     * @param maxIdx Global maximum indices for the bounding box.
     */
    GridFacets(int d, int s, const int minIdx[DIM3], const int maxIdx[DIM3]);

    //! @brief Default destructor (uses `std::vector<bool>` for automatic cleanup).
    ~GridFacets() = default;

    //! @brief Set the flag for a particular `(coord0, coord1)` in this facet.
    /*!
     * @param coord0 Coordinate along the first axis (`axis_dir[0]`).
     * @param coord1 Coordinate along the second axis (`axis_dir[1]`).
     * @param flag Boolean value to set (true or false).
     */
    void SetFlag(int coord0, int coord1, bool flag);

    //! @brief Get the flag for a particular `(coord0, coord1)` in this facet.
    /*!
     * @param coord0 Coordinate along the first axis (`axis_dir[0]`).
     * @param coord1 Coordinate along the second axis (`axis_dir[1]`).
     * @return The boolean flag at `(coord0, coord1)`.
     */
    bool CubeFlag(int coord0, int coord1) const;

    //! @brief Print grid facets information for debugging
    template <typename OSTREAM_TYPE>
    void Print(OSTREAM_TYPE & out) const {
        out << "GridFacets:\n";
        out << "  Orthogonal direction: " << orth_dir << " (0=x, 1=y, 2=z)\n";
        out << "  Side: " << side << " (0=lower, 1=upper)\n";
        out << "  Axis directions: [" << axis_dir[0] << ", " << axis_dir[1] << "]\n";
        out << "  Axis sizes: [" << axis_size[0] << ", " << axis_size[1] << "]\n";
        out << "  Min index: (" << minIndex[0] << ", " << minIndex[1] << ", " << minIndex[2] << ")\n";
        out << "  Max index: (" << maxIndex[0] << ", " << maxIndex[1] << ", " << maxIndex[2] << ")\n";
        out << "  Local size: (" << localSize[0] << ", " << localSize[1] << ", " << localSize[2] << ")\n";
        out << "  Cube flags: " << cube_flag.size() << " flag(s)\n";
    }

private:
    //! @brief Convert a 2D coordinate `(coord0, coord1)` to a linear index.
    /*!
     * @param coord0 Coordinate along the first axis (`axis_dir[0]`).
     * @param coord1 Coordinate along the second axis (`axis_dir[1]`).
     * @return The linearized index for accessing `cube_flag`.
     */
    int index(int coord0, int coord1) const;
};


// Functions for loading nrrd data

//! @brief Converts raw data of type `T` into a float vector.
/*!
 * @tparam T The type of the input data.
 * @param data_ptr Pointer to the input data.
 * @param total_size Total number of elements in the input data.
 * @return A `std::vector<float>` containing the converted data.
 */
template <typename T>
std::vector<float> convert_to_float_vector(T *data_ptr, size_t total_size);

//! @brief Adjusts a point that lies outside the bounds of the scalar grid.
/*!
 * @param p The input point to adjust.
 * @param grid The scalar grid defining the bounds.
 * @param v1 The start point of the segment.
 * @param v2 The end point of the segment.
 * @return A point adjusted to lie within the bounds of the grid.
 */
Point adjust_outside_bound_points(const Point &p, const UnifiedGrid &grid, const Point &v1, const Point &v2);

/*Preprocessing*/

//! @brief Loads NRRD data into a `Grid` structure.
/*!
 * @param file_path The path to the NRRD file.
 * @return A `Grid` object containing the loaded data.
 */
UnifiedGrid load_nrrd_data(const std::string &file_path);

//! @brief Supersamples a `Grid` by a factor of `n`.
/*!
 * @param grid The input grid to supersample.
 * @param n The supersampling factor.
 * @return A supersampled `Grid`.
 */
UnifiedGrid supersample_grid(const UnifiedGrid &grid, int n);

//! @brief Checks if a cube is active based on its scalar values.
/*!
 * @param grid The grid containing the scalar values.
 * @param x The x-coordinate of the cube.
 * @param y The y-coordinate of the cube.
 * @param z The z-coordinate of the cube.
 * @param isovalue The isovalue for activity determination.
 * @return `true` if the cube is active; `false` otherwise.
 */
bool is_cube_active(const UnifiedGrid &grid, int x, int y, int z, float isovalue);

//! @brief Finds all active cubes in the grid based on an isovalue.
/*!
 * @param grid The input grid.
 * @param isovalue The isovalue for activity determination.
 * @param cubes A vector to store the active cubes.
 */
void find_active_cubes(const UnifiedGrid &grid, float isovalue, std::vector<Cube> &cubes);

//! @brief Loads the grid points from a `Grid`.
/*!
 * @param grid The input grid.
 * @return A vector of points representing the grid points.
 */
std::vector<Point> load_grid_points(const UnifiedGrid &grid);

//! @brief Checks if a point is inside the bounds of the grid.
/*!
 * @param p The point to check.
 * @param grid The scalar grid defining the bounds.
 * @return `true` if the point is inside the grid; `false` otherwise.
 */
bool is_point_inside_grid(const Point &p, const UnifiedGrid &grid);

//! @brief Interpolates a point on a line segment based on scalar values.
/*!
 * @param p1 The first point of the segment.
 * @param p2 The second point of the segment.
 * @param val1 The scalar value at `p1`.
 * @param val2 The scalar value at `p2`.
 * @param isovalue The target scalar value for interpolation.
 * @param grid The grid containing the scalar data.
 * @return The interpolated point.
 */
Point interpolate(const Point &p1, const Point &p2, float val1, float val2, float isovalue, const UnifiedGrid &grid);


//! @brief Performs trilinear interpolation at a given point in a scalar grid.
/*!
 * @param p The point for interpolation.
 * @param grid The scalar grid containing the data.
 * @return The interpolated scalar value.
 */
float trilinear_interpolate(const Point &p, const UnifiedGrid &grid);


//! @brief Creates grid facets for a given set of active cubes.
/*!
 * @param activeCubes A vector of active cubes.
 * @return A 3D vector of `GridFacets`, grouped by direction and side.
 */
std::vector<std::vector<GridFacets>> create_grid_facets(const std::vector<Cube> &activeCubes);




//! @brief Tests 6-neighborhood adjacency between two cubes in grid space.
/*! 
 * Two cubes are adjacent if they share a face in the i/j/k lattice
 * (Manhattan distance of 1 between indices).
 *
 * @param cubeA First cube (with i/j/k indices)
 * @param cubeB Second cube (with i/j/k indices)
 * @param grid Grid providing dimensions for bounds checks
 * @return true if the cubes are face-adjacent; false otherwise
 */
bool is_adjacent(const Cube &cubeA, const Cube &cubeB, const UnifiedGrid &grid);

//! @brief Computes a unique linear index for a cube.
/*!
 * Uses the integer grid coordinates derived from the representative vertex and
 * maps them to a linear index in the (nx-1)×(ny-1)×(nz-1) active‑cell lattice.
 *
 * @param repVertex Representative vertex position of the cube (world coords)
 * @param grid Grid describing spacing and origin
 * @return Linear index in the active‑cell lattice
 */
int get_cube_index(const Point &repVertex, const UnifiedGrid &grid);

//! @brief Finds neighboring cube indices around a given cube.
/*!
 * Returns the linear indices of all valid 26‑neighbors around the cube that
 * contains the representative vertex.
 *
 * @param repVertex Representative vertex position of the reference cube
 * @param grid Grid describing spacing and size
 * @return Vector of neighbor indices in the active‑cell lattice
 */
std::vector<int> find_neighbor_indices(const Point &repVertex, const UnifiedGrid &grid);

//! @brief Extracts iso-crossing points for a list of cubes.
/*!
 * @param cubes Input cubes with precomputed iso-crossing points
 * @return Points where the isovalue is crossed within each cube
 */
std::vector<Point> get_cube_centers(const std::vector<Cube> &cubes);

//! @brief Extracts accurate iso-crossing points for a list of cubes.
/*!
 * @param cubes Input cubes with precomputed accurate iso-crossing points
 * @return Accurate iso-crossing points (centroid of edge intersections)
 */
std::vector<Point> get_cube_accurate_iso_crossing_points(const std::vector<Cube> &cubes);

//! @brief Filters active cubes to a set without vertex adjacency (greedy).
/*!
 * Iterates active cubes and keeps a cube only if none of its 26 neighbors
 * (sharing a face, edge, or vertex) was kept before. Useful for producing a
 * well‑separated subset of active regions.
 *
 * @param activeCubes Active cubes to filter
 * @param grid Grid for bounds/neighbor checks
 * @return A subset of input cubes with no vertex/edge/face adjacency among them
 */
std::vector<Cube> separate_active_cubes_greedy(std::vector<Cube> &activeCubes, const UnifiedGrid &grid);

#endif
