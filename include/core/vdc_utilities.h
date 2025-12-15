//! @file vdc_utilities.h
//! @brief Utility functions for Voronoi and Delaunay computations, including centroid calculations and geometric checks.

#ifndef VDC_UTILITIES_H
#define VDC_UTILITIES_H

#include "core/vdc_debug.h"
#include "core/vdc_type.h"
#include "processing/vdc_grid.h"
#include "processing/vdc_delaunay.h"
#include "core/vdc_commandline.h"


//! @brief Checks if two scalar values are bipolar with respect to an isovalue.
/*!
 * Determines if one value is above and the other is below the specified isovalue.
 * 
 * @param val1 First scalar value.
 * @param val2 Second scalar value.
 * @param isovalue The isovalue used for comparison (default is 0).
 * @return `true` if the values are bipolar, otherwise `false`.
 */
bool is_bipolar(float val1, float val2, float isovalue = 0);

//! @brief Checks if a CGAL::Object is degenerate.
/*!
 * A segment is degenerate if its start and end points are the same.
 *
 * @param obj The CGAL::Object to check.
 * @return `true` if the object is degenerate, otherwise `false`.
 */
bool is_degenerate(const Object &obj);

//! @brief Checks if a Delaunay cell is degenerate.
/*!
 * A Delaunay cell is degenerate if its volume is approximately zero.
 * 
 * @param cell The Delaunay cell handle.
 * @return `true` if the cell is degenerate, otherwise `false`.
 */
bool is_degenerate(Delaunay::Cell_handle cell);

//! @brief Checks if a value is positive relative to an isovalue.
/*!
 * @param value The scalar value to check.
 * @param isovalue The isovalue threshold for comparison.
 * @return `true` if the value is greater than or equal to the isovalue, otherwise `false`.
 */
bool is_positive(double value, float isovalue);

/* General Helper Functions */

//! @brief Computes the centroid of a set of points.
/*!
 * This function calculates the geometric centroid of a set of points.
 *
 * @param points Vector of points to compute the centroid for.
 * @param supersample Indicates if supersampling is applied (currently unused, reserved for future use).
 * @param ratio Supersampling ratio (currently unused, reserved for future use).
 * @return The computed centroid as a Point.
 */
Point compute_centroid(const std::vector<Point> &points, bool supersample, int ratio);

//! @brief Computes the centroid of a set of points using CGAL.
/*!
 * @param points Vector of points to compute the centroid for.
 * @return The computed centroid as a Point.
 */
Point compute_centroid(const std::vector<Point> &points);

//! @brief Gets the corners of a cube.
/*!
 * Computes the eight corner points of a cube given its center and side length.
 * 
 * @param center The center of the cube.
 * @param side_length The side length of the cube.
 * @return An array containing the eight corner points.
 */
std::array<Point, 8> get_cube_corners(const Point &center, float side_length);

//! @brief Determines the orientation of a facet.
/*!
 * Based on the facet index, vertex positions, and scalar values, this function determines
 * the orientation of a facet (positive or negative).
 * 
 * @param iFacet The facet index.
 * @param v1 First vertex of the facet.
 * @param v2 Second vertex of the facet.
 * @param f1 Scalar value at v1.
 * @param f2 Scalar value at v2.
 * @return 1 for positive orientation, -1 for negative.
 */
int get_orientation(const int iFacet, const Point v1, const Point v2, const float f1, const float f2);


//! @brief Loads 3D points from a text file.
/*!
 * Accepts simple whitespace/comma separated lines containing three numeric
 * coordinates per line. Lines that cannot be parsed are ignored.
 *
 * @param filename Path to the input file
 * @param points Output container receiving parsed points
 * @return true on success, false if the file could not be opened
 */
bool read_points_from_file(const std::string &filename, std::vector<Point> &points);

//! @brief Writes a Delaunay triangulation and its input points for debugging.
/*!
 * Produces a text dump next to the given input file name to aid inspection of
 * triangulation quality and topology.
 *
 * @param dt The triangulation to serialize
 * @param points The set of points used to build the triangulation
 * @param input_filename The original input filename used to derive an output stem
 */
void write_triangulation(Delaunay dt, std::vector<Point> &points, std::string &input_filename);

#endif // VDC_UTILITIES_H
