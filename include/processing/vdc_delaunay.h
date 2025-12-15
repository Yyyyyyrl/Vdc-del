//! @file vdc_delaunay.h
//! @brief Header file defining structures and functions related to Delaunay triangulations and isosurface triangles.

#ifndef VDC_DELAUNAY_H
#define VDC_DELAUNAY_H

#include "core/vdc_type.h"

//! @brief Represents a triangle in the Delaunay triangulation.
/*!
 * This structure is used to store the geometric representation of a triangle,
 * defined by three vertices.
 */
struct DelaunayTriangle {
    Vertex_handle vertex1; //!< First vertex of the triangle.
    Vertex_handle vertex2; //!< Second vertex of the triangle.
    Vertex_handle vertex3; //!< Third vertex of the triangle.

    //! @brief Constructor to initialize a DelaunayTriangle.
    /*!
     * @param v1 Handle to the first vertex.
     * @param v2 Handle to the second vertex.
     * @param v3 Handle to the third vertex.
     */
    DelaunayTriangle(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3)
        : vertex1(v1), vertex2(v2), vertex3(v3) {}

    //! @brief Print Delaunay triangle information for debugging
    template <typename OSTREAM_TYPE>
    void Print(OSTREAM_TYPE & out) const {
        out << "DelaunayTriangle:\n";
        out << "  Vertex 1: " << vertex1->point() << "\n";
        out << "  Vertex 2: " << vertex2->point() << "\n";
        out << "  Vertex 3: " << vertex3->point() << "\n";
    }
};

//! @brief Represents a triangle on an isosurface.
/*!
 * An isosurface triangle is defined by the indices of its three vertices,
 * which refer to positions in the isosurfaceVertices vector.
 */
struct IsoTriangle {
    int vertex_indices[3]; //!< Indices into the isosurfaceVertices vector.

    //! @brief Constructor to initialize an IsoTriangle.
    /*!
     * @param idx1 Index of the first vertex in the isosurfaceVertices vector.
     * @param idx2 Index of the second vertex in the isosurfaceVertices vector.
     * @param idx3 Index of the third vertex in the isosurfaceVertices vector.
     */
    IsoTriangle(int idx1, int idx2, int idx3) {
        vertex_indices[0] = idx1;
        vertex_indices[1] = idx2;
        vertex_indices[2] = idx3;
    }

    //! @brief Print iso-triangle information for debugging
    template <typename OSTREAM_TYPE>
    void Print(OSTREAM_TYPE & out) const {
        out << "IsoTriangle:\n";
        out << "  Vertex indices: [" << vertex_indices[0] << ", "
            << vertex_indices[1] << ", " << vertex_indices[2] << "]\n";
    }
};

//! @brief Overloaded output operator for IsoTriangle.
/*!
 * Provides a formatted string representation of an IsoTriangle for debugging
 * and logging purposes.
 *
 * @tparam OSTREAM_TYPE The type of the output stream (e.g., `std::ostream`).
 * @param os The output stream to write to.
 * @param it The IsoTriangle to output.
 * @return A reference to the output stream.
 */
template <typename OSTREAM_TYPE>
OSTREAM_TYPE& operator<<(OSTREAM_TYPE& os, const IsoTriangle& it)
{
    os << "IsoTriangle:\n";
    os << "  Vertex indices: " << it.vertex_indices[0] << ", "
       << it.vertex_indices[1] << ", " << it.vertex_indices[2] << "\n";

    return os;
}

#endif // VDC_DELAUNAY_H
