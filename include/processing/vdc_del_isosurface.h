//! @file vdc_del_isosurface.h
//! @brief Data structures and function declarations for Delaunay-based isosurface extraction.

#ifndef VDC_DEL_ISOSURFACE_H
#define VDC_DEL_ISOSURFACE_H

#include "core/vdc_type.h"
#include "processing/vdc_grid.h"
#include <map>
#include <array>

//! @brief Structure to hold the extracted isosurface from Delaunay-based algorithm.
/*!
 * Contains isosurface vertices (one or more per active Delaunay vertex),
 * triangles, and mappings between Delaunay structures and isosurface elements.
 */
struct DelaunayIsosurface
{
    //! @brief Isosurface vertex positions (one or more per active Delaunay vertex)
    std::vector<Point> isovertices;

    //! @brief Index of the Delaunay vertex that owns each isosurface vertex
    std::vector<int> isovertex_delaunay_vertex;

    //! @brief Index of the cycle within that Delaunay vertex for each isosurface vertex
    std::vector<int> isovertex_cycle_index;

    //! @brief Output triangles (indices into isovertices array)
    std::vector<std::array<int, 3>> triangles;

    //! @brief Mapping from (delaunay_vertex_index, cycle_index) -> isovertex index
    std::map<std::pair<int, int>, int> vertex_cycle_to_isovertex;

    //! @brief Clear all data
    void clear() {
        isovertices.clear();
        isovertex_delaunay_vertex.clear();
        isovertex_cycle_index.clear();
        triangles.clear();
        vertex_cycle_to_isovertex.clear();
    }

    //! @brief Get the number of isosurface vertices
    size_t num_vertices() const { return isovertices.size(); }

    //! @brief Get the number of triangles
    size_t num_triangles() const { return triangles.size(); }
};

// ============================================================================
// Function declarations for Delaunay-based isosurface extraction
// ============================================================================

//! @brief Compute circumcenters and scalar values for all finite Delaunay cells.
/*!
 * For each finite Delaunay cell:
 * - Computes the circumcenter (center of circumscribed sphere)
 * - Interpolates the scalar value at the circumcenter from the grid
 * - Sets flag_positive based on whether the scalar >= isovalue
 *
 * @param dt The Delaunay triangulation
 * @param grid The scalar field grid
 * @param isovalue The isovalue threshold
 */
void compute_cell_circumcenters_and_scalars(
    Delaunay& dt,
    const UnifiedGrid& grid,
    float isovalue
);

//! @brief Mark Delaunay facets that cross the isosurface.
/*!
 * A facet is marked as isosurface if it separates a positive cell from a negative cell.
 * Only marks from the positive side to avoid double processing.
 *
 * @param dt The Delaunay triangulation (must have circumcenters computed first)
 */
void mark_isosurface_facets(Delaunay& dt);

//! @brief Mark active Delaunay vertices.
/*!
 * A vertex is active if it is incident on any isosurface facet.
 *
 * @param dt The Delaunay triangulation (must have isosurface facets marked first)
 */
void mark_active_vertices(Delaunay& dt);

//! @brief Compute cycles of isosurface facets around each active vertex.
/*!
 * For each active Delaunay vertex, finds all isosurface facets incident on it
 * and groups them into cycles (connected components).
 *
 * @param dt The Delaunay triangulation (must have active vertices marked first)
 */
void compute_facet_cycles(Delaunay& dt);

//! @brief Compute the centroid of Voronoi edge / isosurface intersections for a cycle.
/*!
 * For each facet in the cycle:
 * - Gets the two incident cells and their circumcenters
 * - Computes the isosurface crossing point via linear interpolation
 * - Returns the centroid of all crossing points
 *
 * @param dt The Delaunay triangulation
 * @param grid The scalar field grid
 * @param isovalue The isovalue threshold
 * @param v_handle Handle to the Delaunay vertex
 * @param cycle_facets List of (cell_index, facet_index) pairs in the cycle
 * @return The centroid of all isosurface crossing points
 */
Point compute_centroid_of_voronoi_edge_and_isosurface(
    const Delaunay& dt,
    const UnifiedGrid& grid,
    float isovalue,
    Vertex_handle v_handle,
    const std::vector<std::pair<int, int>>& cycle_facets
);

//! @brief Compute isosurface vertex positions for each cycle around each active vertex.
/*!
 * For single-cycle vertices: position at vertex location (or isov if -position_delv_on_isov)
 * For multi-cycle vertices: compute centroid of Voronoi edge intersections for each cycle
 *
 * @param dt The Delaunay triangulation (must have cycles computed first)
 * @param grid The scalar field grid
 * @param isovalue The isovalue threshold
 * @param position_on_isov If true, use accurate iso-crossing point for single-cycle vertices
 */
void compute_cycle_isovertices(
    Delaunay& dt,
    const UnifiedGrid& grid,
    float isovalue,
    bool position_on_isov
);

//! @brief Generate isosurface triangles from marked facets.
/*!
 * For each isosurface facet, creates a triangle by looking up the isosurface
 * vertex for each of the facet's three Delaunay vertices.
 *
 * @param dt The Delaunay triangulation (must have cycle isovertices computed)
 * @param iso_surface Output structure to receive vertices and triangles
 */
void generate_isosurface_triangles(
    const Delaunay& dt,
    DelaunayIsosurface& iso_surface
);

//! @brief Find which cycle (if any) contains a given facet for a vertex.
/*!
 * @param v_handle Handle to the Delaunay vertex
 * @param cell_index Index of the cell containing the facet
 * @param facet_index Index of the facet within the cell (0-3)
 * @return Cycle index (0-based), or -1 if not found
 */
int find_cycle_containing_facet(
    Vertex_handle v_handle,
    int cell_index,
    int facet_index
);

//! @brief Project a point onto a sphere centered at a given point.
/*!
 * Used for projecting cycle centroids onto a sphere around the Delaunay vertex.
 *
 * @param point The point to project
 * @param center The center of the sphere
 * @param radius The radius of the sphere
 * @return The projected point
 */
Point project_to_sphere(const Point& point, const Point& center, double radius);

#endif // VDC_DEL_ISOSURFACE_H
