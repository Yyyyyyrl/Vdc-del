//! @file vdc_type.h
//! @brief Type definitions and includes for Voronoi Diagram and Delaunay Triangulation operations.

#ifndef VDC_TYPE_H
#define VDC_TYPE_H
// Note: moved to core/ subfolder; include via "core/vdc_type.h"

// Standard library headers for various utilities.
#include <iostream> // Input and output stream.
#include <fstream>  // File stream for reading and writing files.
#include <vector>   // STL vector container.
#include <array>    // STL array container.
#include <cmath>    // Math functions like abs(), pow(), etc.
#include <limits>   // To work with numeric limits of data types.
#include <set>
#include <map>
#include <stack>
#include <queue>
#include <unordered_map>
#include <tuple>
#include <vector>
#include <functional>
#include <iterator>    // for std::back_inserter
#include <algorithm>   // Common algorithms like sort(), min(), max().
#include <cstddef>     // Definitions for size_t, ptrdiff_t, etc.
#include <cstring>     // C-style string operations.
#include <teem/nrrd.h> // Teem library for handling NRRD data files.
#include <iomanip>     // For formatted output (e.g., precision control).
#include <variant>
#include <ctime>       // For time-related functions and data types. ( clock )

// CGAL headers for computational geometry operations.
#include <CGAL/bounding_box.h>                                  // Compute bounding boxes for geometric objects.
#include <CGAL/intersections.h>                                 // Perform geometric intersection tests.
#include <CGAL/Vector_3.h>                                      // Represents a vector in 3D space.
#include <CGAL/Polyhedron_3.h>                                  // Represents 3D polyhedral surfaces.
#include <CGAL/convex_hull_3.h>                                 // Compute convex hulls in 3D.
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>   // Compute Halfspace Intersections in 3D
#include <CGAL/centroid.h>                                      // Compute centroids of geometric objects.
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h> // Kernel for exact predicates and inexact constructions.
#include <CGAL/Delaunay_triangulation_3.h>                      // Delaunay triangulation in 3D.
#include <CGAL/Triangulation_vertex_base_with_info_3.h>         // Vertex base with additional user data.
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_with_circumcenter_3.h> // Cell base with circumcenter support.

//! @brief CGAL Kernel.
/*!
 * Provides geometric primitives like points, vectors, lines, and planes with
 * exact predicates and inexact constructions.
 */
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

//! @brief 3D point type from the kernel (forward declaration for use in VertexInfo)
typedef K::Point_3 Point3;


//! @brief Information associated with a vertex in the Delaunay triangulation.
/*!
 * Stores metadata for vertices in the Delaunay triangulation.
 * Extended for Delaunay-based isosurface extraction algorithm.
 */
struct VertexInfo
{
    bool is_dummy;        //!< Flag indicating if this is a dummy vertex added for bounding
    int voronoiCellIndex; //!< Index of the Voronoi cell dual to this vertex
    int index;            //!< Unique index of this vertex in the triangulation

    // ─── Fields for Delaunay-based isosurface extraction ───────────────────────
    bool active;          //!< true if vertex is incident on any isosurface facet
    Point3 isov;          //!< Isosurface vertex location (if -position_delv_on_isov)
    bool has_isov_sample; //!< true if `isov` stores a valid isosurface sample point

    //! Cycles of isosurface facets around this vertex.
    //! Each cycle is a list of (cell_index, facet_index) pairs representing facets in the cycle.
    std::vector<std::vector<std::pair<int, int>>> facet_cycles;

    //! Isosurface vertex positions for each cycle (one per cycle)
    std::vector<Point3> cycle_isovertices;

    //! @brief Default constructor
    VertexInfo()
        : is_dummy(false),
          voronoiCellIndex(-1),
          index(-1),
          active(false),
          isov(Point3(0, 0, 0)),
          has_isov_sample(false)
    {}

    //! @brief Print vertex info for debugging
    template <typename OSTREAM_TYPE>
    void Print(OSTREAM_TYPE & out) const {
        out << "VertexInfo:\n";
        out << "  Index: " << index << "\n";
        out << "  Is dummy: " << (is_dummy ? "true" : "false") << "\n";
        out << "  Voronoi cell index: " << voronoiCellIndex << "\n";
        out << "  Active: " << (active ? "true" : "false") << "\n";
        out << "  Has isov sample: " << (has_isov_sample ? "true" : "false") << "\n";
        out << "  Number of cycles: " << facet_cycles.size() << "\n";
    }
};

//! @brief Information associated with a facet in a Delaunay cell.
/*!
 * Each Delaunay facet is dual to a Voronoi edge. This structure stores
 * references to the dual Voronoi edge and the corresponding Voronoi cell edges
 * for quick lookup during Voronoi facet construction.
 */
struct DelaunayFacetInfo
{
    int dualEdgeIndex;           //!< Index of the Voronoi edge dual to this Delaunay facet
    int dualCellEdgeIndex[3];    //!< Indices of Voronoi cell edges for the 3 vertices of this facet
                                 //!< dualCellEdgeIndex[k] is the index of the cell edge dual to this facet
                                 //!< and in the Voronoi cell around the k'th vertex of the facet

    //! @brief Constructor to initialize DelaunayFacetInfo with default values
    DelaunayFacetInfo() : dualEdgeIndex(-1)
    {
        dualCellEdgeIndex[0] = -1;
        dualCellEdgeIndex[1] = -1;
        dualCellEdgeIndex[2] = -1;
    }

    //! @brief Check if a dual cell edge index is undefined
    /*!
     * @param k Index into dualCellEdgeIndex array (0, 1, or 2)
     * @return true if the index is undefined (-1), false otherwise
     */
    bool IsDualCellEdgeIndexUndefined(const int k) const
    {
        return (dualCellEdgeIndex[k] == -1);
    }

    //! @brief Print Delaunay facet info for debugging
    template <typename OSTREAM_TYPE>
    void Print(OSTREAM_TYPE & out) const {
        out << "DelaunayFacetInfo:\n";
        out << "  Dual edge index: " << dualEdgeIndex << "\n";
        out << "  Dual cell edge indices: ["
            << dualCellEdgeIndex[0] << ", "
            << dualCellEdgeIndex[1] << ", "
            << dualCellEdgeIndex[2] << "]\n";
    }
};

//! @brief Information associated with a cell in the Delaunay triangulation.
/*!
 * Stores metadata for cells (tetrahedra) in the Delaunay triangulation.
 * Each Delaunay cell (tetrahedron) has 4 facets indexed 0-3, where facet i
 * is opposite to vertex i. Each facet is dual to a Voronoi edge.
 * Extended for Delaunay-based isosurface extraction algorithm.
 */
struct CellInfo
{
    int dualVoronoiVertexIndex;     //!< Index of the Voronoi vertex dual to this cell
    int index;                      //!< Unique index of this cell in the triangulation
    DelaunayFacetInfo facet_info[4]; //!< Information for each of the 4 facets of this cell
    int edge_index[4][4];           //!< Global Delaunay edge indices: edge_index[i][j] is the
                                    //!< index of the edge (vertex i, vertex j). edge_index[i][i] is undefined.

    // ─── Fields for Delaunay-based isosurface extraction ───────────────────────
    Point3 circumcenter;            //!< Circumcenter of this cell (dual Voronoi vertex location)
    float circumcenter_scalar;      //!< Scalar value at the circumcenter
    bool flag_positive;             //!< true if circumcenter_scalar >= isovalue
    bool facet_is_isosurface[4];    //!< true if facet i crosses the isosurface

    //! @brief Constructor to initialize CellInfo with default values
    CellInfo() : dualVoronoiVertexIndex(-1), index(-1),
                 circumcenter_scalar(0.0f), flag_positive(false)
    {
        // Initialize edge_index to -1 (undefined)
        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < 4; ++j)
            {
                edge_index[i][j] = -1;
            }
            facet_is_isosurface[i] = false;
        }
    }

    //! @brief Return the index (0, 1, 2, or 3) of the k'th vertex of facet facet_index
    /*!
     * For a Delaunay cell with facet indexed as facet_index (opposite to vertex facet_index),
     * this computes the cell vertex index of the k'th vertex (k = 0, 1, 2) on that facet.
     *
     * @param facet_index The facet index (0-3)
     * @param k The vertex position on the facet (0, 1, or 2)
     * @return The cell vertex index (0-3)
     */
    static int FacetVertexIndex(const int facet_index, const int k)
    {
        return (facet_index + k + 1) % 4;
    }

    //! @brief Print cell info for debugging
    template <typename OSTREAM_TYPE>
    void Print(OSTREAM_TYPE & out) const {
        out << "CellInfo:\n";
        out << "  Index: " << index << "\n";
        out << "  Dual Voronoi vertex index: " << dualVoronoiVertexIndex << "\n";
        out << "  Facet info count: 4\n";
        out << "  Circumcenter scalar: " << circumcenter_scalar << "\n";
        out << "  Flag positive: " << (flag_positive ? "true" : "false") << "\n";
        out << "  Isosurface facets: ["
            << facet_is_isosurface[0] << ", "
            << facet_is_isosurface[1] << ", "
            << facet_is_isosurface[2] << ", "
            << facet_is_isosurface[3] << "]\n";
    }
};

//! @brief Enum indicating the method for matching bipolar edges on a facet.
enum class BIPOLAR_MATCH_METHOD {
    SEP_POS,               // Separate positive regions: pair neg->pos to next pos->neg.
    SEP_NEG,               // Separate negative regions: pair pos->neg to next neg->pos.
    UNCONSTRAINED_MATCH,   // Arbitrary pairing (e.g., consecutive pairs without sign consideration).
    UNDEFINED_MATCH_TYPE   // Undefined or error state.
};

//! @brief Vertex base for triangulations with additional information.
typedef CGAL::Triangulation_vertex_base_with_info_3<VertexInfo, K> Vb;

//! @brief Cell base for Delaunay triangulations.
typedef CGAL::Delaunay_triangulation_cell_base_with_circumcenter_3<K> Cb2;

//! @brief Cell base for Delaunay triangulations with additional information.
typedef CGAL::Triangulation_cell_base_with_info_3<CellInfo, K, Cb2> Cb;

//! @brief Data structure for triangulations.
/*!
 * Combines the vertex and cell bases into a complete triangulation data structure.
 */
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;

//! @brief 3D Delaunay triangulation.
/*!
 * A Delaunay triangulation is a geometric structure that divides 3D space
 * into tetrahedra such that no vertex lies inside the circumsphere of any tetrahedron.
 */
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;

//! @brief 3D polyhedral surface representation.
/*!
 * A polyhedron is used to represent Voronoi cells and other geometric structures.
 */
typedef CGAL::Polyhedron_3<K> Polyhedron_3;

//! @brief 3D point type used in Delaunay triangulations.
/*!
 * Points are the fundamental objects used for constructing the triangulation.
 */
typedef Delaunay::Point Point;

//! @brief Handle to a cell in the Delaunay triangulation.
/*!
 * A cell handle allows access to the tetrahedra in the triangulation.
 */
typedef Delaunay::Cell_handle Cell_handle;

//! @brief Handle to a vertex in the Delaunay triangulation.
/*!
 * A vertex handle allows access to a specific vertex in the triangulation.
 */
typedef Delaunay::Vertex_handle Vertex_handle;

//! @brief Circulator for iterating through cells around a vertex.
/*!
 * This circulator iterates through all cells sharing a specific vertex.
 */
typedef Delaunay::Cell_circulator Cell_circulator;

//! @brief CGAL object wrapper for geometric primitives.
/*!
 * Used to store generic geometric objects like points, segments, or lines.
 */
typedef CGAL::Object Object;

//! @brief Represents an edge in the Delaunay triangulation.
typedef Delaunay::Edge Edge;

//! @brief Represents a facet (triangle face) in the Delaunay triangulation.
typedef Delaunay::Facet Facet;

//! @brief Represents a triangle in 3D space.
typedef Delaunay::Triangle Triangle;

//! @brief Represents a line segment in 3D space.
typedef K::Segment_3 Segment3;

//! @brief Represents a ray (semi-infinite line) in 3D space.
typedef K::Ray_3 Ray3;

//! @brief Represents a line in 3D space.
typedef K::Line_3 Line3;

//! @brief Represents a point in 3D space.
typedef K::Point_3 Point3;

//! @brief Represents a vector in 3D space.
/*!
 * Vectors are used for translations and other vector operations.
 */
typedef K::Vector_3 Vector3;

//! @brief Represents a plane in 3D space.
/*!
 * A plane is defined by a point and a normal vector or by three non-collinear points.
 */
typedef K::Plane_3 Plane_3;

//! @brief Direction in 3D space.
/*!
 * A direction is a vector with norm 1, representing a line through the origin.
 */
typedef K::Direction_3 Direction3;

//
// Custom Output << Operators
//

//! @brief Output operator for Cell_handle
/*!
 * This operator allows printing of Cell_handle objects to output streams.
 * It displays the cell's vertices, circumcenter, and additional info.
 *
 * @param os The output stream to write to
 * @param ch The Cell_handle to output
 * @return Reference to the output stream
 */
template <typename OSTREAM_TYPE>
OSTREAM_TYPE &operator<<(OSTREAM_TYPE &os, const Cell_handle &ch)
{
    if (ch == nullptr)
    {
        os << "Cell_handle(nullptr)";
        return os;
    }

    // Print cell info
    os << "Cell_handle[" << &(*ch) << "]:\n";
    os << "  Info: { index: " << ch->info().index
       << ", dualVoronoiVertexIndex: " << ch->info().dualVoronoiVertexIndex << " }\n";

    // Print vertices
    os << "  Vertices:\n";
    for (int i = 0; i < 4; ++i)
    {
        Vertex_handle vh = ch->vertex(i);
        os << "    [" << i << "]: ("
           << vh->point().x() << ", "
           << vh->point().y() << ", "
           << vh->point().z() << ")";

        if (vh->info().is_dummy)
        {
            os << " (dummy)";
        }
        os << " index: " << vh->info().index;
        os << "\n";
    }

    // Print circumcenter
    try
    {
        Point cc = ch->circumcenter();
        os << "  Circumcenter: ("
           << cc.x() << ", "
           << cc.y() << ", "
           << cc.z() << ")\n";
    }
    catch (...)
    {
        os << "  Circumcenter: (calculation failed)\n";
    }

    return os;
}

//! @brief Output operator for Cell_circulator
/*!
 * This operator allows printing of Cell_circulator objects to output streams.
 * It displays the current cell the circulator points to.
 *
 * @param os The output stream to write to
 * @param cc The Cell_circulator to output
 * @return Reference to the output stream
 */
template <typename OSTREAM_TYPE>
OSTREAM_TYPE &operator<<(OSTREAM_TYPE &os, const Cell_circulator &cc)
{
    if (cc == Cell_circulator())
    {
        os << "Cell_circulator(default)";
        return os;
    }

    // Get the current cell and output it
    Cell_handle ch = cc;
    os << "Cell_circulator -> " << ch;

    return os;
}

//! @brief Output operator for Delaunay triangulation.
/*!
 * Prints comprehensive information about a Delaunay triangulation including:
 * - Basic statistics (dimension, validity)
 * - Vertex statistics (total, regular, dummy)
 * - Cell statistics (finite, infinite)
 * - Sample cell details (limited to 5 cells)
 * - Vertex degree distribution
 * - Sample vertex details (limited to 10 vertices)
 *
 * @tparam OSTREAM_TYPE Output stream type (e.g. std::ostream, std::ofstream)
 * @param os Output stream to write to
 * @param dt Delaunay triangulation to print
 * @return Reference to the output stream
 */
template <typename OSTREAM_TYPE>
OSTREAM_TYPE &operator<<(OSTREAM_TYPE &os, const Delaunay &dt)
{
    // Basic statistics
    os << "Delaunay Triangulation Statistics:\n"
       << "  Dimension: " << dt.dimension() << "\n"
       << "  Is valid: " << (dt.is_valid() ? "yes" : "no") << "\n\n";

    // Vertex statistics
    int total_vertices = dt.number_of_vertices();
    int dummy_vertices = 0;
    int regular_vertices = 0;

    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit)
    {
        if (vit->info().is_dummy)
        {
            dummy_vertices++;
        }
        else
        {
            regular_vertices++;
        }
    }

    os << "Vertex Statistics:\n"
       << "  Total vertices: " << total_vertices << "\n"
       << "  Regular vertices: " << regular_vertices << "\n"
       << "  Dummy vertices: " << dummy_vertices << "\n";

    // Cell statistics
    int finite_cells = dt.number_of_finite_cells();
    int infinite_cells = dt.number_of_cells() - finite_cells;

    os << "\nCell Statistics:\n"
       << "  Total cells: " << dt.number_of_cells() << "\n"
       << "  Finite cells: " << finite_cells << "\n"
       << "  Infinite cells: " << infinite_cells << "\n";

    // Output detailed information for finite cells
    os << "\nCell Details (limited to first 5 finite cells):\n";
    int cell_count = 0;
    for (auto cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit)
    {
        Cell_handle ch = cit;
        os << ch;

        // Limit the number of cells to display to avoid excessive output
        if (++cell_count >= 5)
        {
            os << "... and " << (finite_cells - 5) << " more finite cells\n";
            break;
        }

        // Add a separator line
        os << "  -----------------------------------------\n";
    }

    // Vertex degree distribution (number of adjacent cells per vertex)
    std::map<int, int> degree_distribution;
    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit)
    {
        int degree = dt.degree(vit);
        degree_distribution[degree]++;
    }

    os << "\nVertex Degree Distribution:\n"
       << "  (degree: number of vertices)\n";
    for (const auto &pair : degree_distribution)
    {
        os << "  " << pair.first << " cells: " << pair.second << " vertices\n";
    }

    // Vertex detailed information
    os << "\nVertex Details:\n";
    int vertex_count = 0;
    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit)
    {
        os << "  [" << vit->info().index << "] ("
           << std::fixed << std::setprecision(3)
           << vit->point().x() << ", "
           << vit->point().y() << ", "
           << vit->point().z() << ")"
           << (vit->info().is_dummy ? " (dummy)" : "")
           << " voronoi_cell: " << vit->info().voronoiCellIndex
           << " degree: " << dt.degree(vit)
           << "\n";

        // Limit the number of vertices to display to avoid excessive output
        if (++vertex_count >= 30)
        {
            os << "  ... and " << (total_vertices - 30) << " more vertices\n";
            break;
        }
    }

    return os;
}

//! @brief Helper function to print all cells in a Cell_circulator.
/*!
 * Iterates through all cells in a circulator and prints detailed information
 * about each cell. Includes safety checks to prevent infinite loops.
 * 
 * @tparam OSTREAM_TYPE Type of output stream (e.g. std::ostream, std::ofstream)
 * @param os Output stream to write to
 * @param start Starting cell circulator to begin iteration
 * @param vertex_index Index of the vertex being circulated around (for context)
 */
template <typename OSTREAM_TYPE>
void print_cell_circuit(OSTREAM_TYPE &os, Cell_circulator start, int vertex_index)
{
    if (start == Cell_circulator())
    {
        os << "Empty Cell_circulator\n";
        return;
    }

    os << "Cells around vertex " << vertex_index << ":\n";

    Cell_circulator current = start;
    int count = 0;
    do
    {
        os << "Cell " << count << ":\n"
           << current << "\n";
        ++current;
        ++count;
    } while (current != start && count < 100); // Safety limit to prevent infinite loops

    if (count >= 100)
    {
        os << "Warning: Circuit may be infinite, stopped after 100 cells\n";
    }

    os << "Total: " << count << " cells\n";
}

//! @brief Structure for comparing points for approximate equality.
/*!
 * This structure defines a custom comparator for points, allowing comparison
 * with a small epsilon tolerance to handle floating-point precision issues.
 */
struct PointApproxEqual
{
    //! @brief Compares two points for approximate equality.
    /*!
     * Uses squared distance comparison with a small epsilon (1e-6) threshold
     * to determine if points are effectively the same location.
     * 
     * @param p1 First point to compare
     * @param p2 Second point to compare
     * @return true if points are within epsilon distance, false otherwise
     */
    bool operator()(const Point &p1, const Point &p2) const
    {
        const double eps = 5e-5;
        return (std::abs(p1.x() - p2.x()) < eps) &&
           (std::abs(p1.y() - p2.y()) < eps) &&
           (std::abs(p1.z() - p2.z()) < eps);
    }
};

//! @brief Strict weak ordering for `Vector3` values.
/*!
 * Compares vectors lexicographically by x, then y, then z components.
 * Useful for ordered containers like `std::map` and `std::set`.
 */
struct Vector3Comparator
{
    bool operator()(const Vector3 &a, const Vector3 &b) const
    {
        if (a.x() < b.x())
            return true;
        if (a.x() > b.x())
            return false;
        if (a.y() < b.y())
            return true;
        if (a.y() > b.y())
            return false;
        return a.z() < b.z();
    }
};

//! @brief Ordering for `Direction3` using underlying vector components.
/*!
 * Converts directions to vectors and applies `Vector3Comparator`.
 * Note that opposite directions will be ordered based on raw component signs.
 */
struct Direction3Comparator
{
    bool operator()(const Direction3 &a, const Direction3 &b) const
    {
        // Convert directions to vectors for comparison
        Vector3 va = a.to_vector();
        Vector3 vb = b.to_vector();
        return Vector3Comparator()(va, vb);
    }
};

//! @brief Comparator for CGAL::Object types.
/*!
 * Provides a strict weak ordering for CGAL::Object types containing geometric primitives.
 * Used for sorting and storing objects in ordered containers.
 */
struct ObjectComparator
{
    //! @brief Compares two CGAL::Object instances for ordering.
    /*!
     * Handles comparison of Segment3, Ray3 and Line3 objects with consistent ordering:
     * 1. Segments (ordered by normalized endpoints)
     * 2. Rays (ordered by source point then direction)
     * 3. Lines (ordered by point then normalized direction)
     *
     * @param a First object to compare
     * @param b Second object to compare
     * @return true if a should be ordered before b, false otherwise
     */
    bool operator()(const CGAL::Object &a, const CGAL::Object &b) const
    {
        // Extract and compare geometric properties based on type
        if (const Segment3 *segA = CGAL::object_cast<Segment3>(&a))
        {
            if (const Segment3 *segB = CGAL::object_cast<Segment3>(&b))
            {
                // Compare segments by their endpoints (order-independent)
                Point sA1 = segA->source(), sA2 = segA->target();
                Point sB1 = segB->source(), sB2 = segB->target();
                std::pair<Point, Point> keyA = (sA1 < sA2) ? std::make_pair(sA1, sA2) : std::make_pair(sA2, sA1);
                std::pair<Point, Point> keyB = (sB1 < sB2) ? std::make_pair(sB1, sB2) : std::make_pair(sB2, sB1);
                return keyA < keyB;
            }
            return true; // Segments come before rays/lines
        }
        if (const Ray3 *rayA = CGAL::object_cast<Ray3>(&a))
        {
            if (const Ray3 *rayB = CGAL::object_cast<Ray3>(&b))
            {
                // Compare rays by source and direction
                auto pairA = std::make_pair(rayA->source(), rayA->direction());
                auto pairB = std::make_pair(rayB->source(), rayB->direction());
                if (pairA.first < pairB.first)
                    return true;
                if (pairB.first < pairA.first)
                    return false;
                return Direction3Comparator()(pairA.second, pairB.second);
            }
            if (CGAL::object_cast<Segment3>(&b))
                return false; // Rays after segments
            return true;      // Rays before lines
        }
        if (const Line3 *lineA = CGAL::object_cast<Line3>(&a))
        {
            if (const Line3 *lineB = CGAL::object_cast<Line3>(&b))
            {
                // Compare lines by a point and direction (normalize direction)
                Vector3 dA = lineA->direction().vector();
                Vector3 dB = lineB->direction().vector();
                if (dA * dB < 0)
                    dB = -dB; // Normalize direction
                auto pairA = std::make_pair(lineA->point(0), dA);
                auto pairB = std::make_pair(lineB->point(0), dB);
                if (pairA.first < pairB.first)
                    return true;
                if (pairB.first < pairA.first)
                    return false;
                return Vector3Comparator()(pairA.second, pairB.second);
            }
            return false; // Lines after segments and rays
        }
        return false; // Default case
    }
};

//! @brief Hash functor for `Point` to use in `unordered_*` containers.
/*!
 * Combines hashed Cartesian coordinates. Intended for exact‑coordinate keys; if
 * approximate equivalence is required consider `PointApproxEqual` and custom
 * bucketing.
 */
struct PointHash
{
    std::size_t operator()(const Point &p) const
    {
        std::hash<double> h;
        return h(p.x()) ^ (h(p.y()) << 1) ^ (h(p.z()) << 2);
    }
};

#endif // VDC_TYPE_H
