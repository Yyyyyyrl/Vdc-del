//! @file vdc_debug.h
//! @brief Header file for debugging functions and utilities.

#ifndef VDC_DEBUG_H
#define VDC_DEBUG_H

#include "core/vdc_type.h"
#include "processing/vdc_grid.h"

//! @brief Global variable to enable or disable debug output.
/*!
 * When `debug` is set to `true`, debug-related print commands are executed.
 */
extern bool debug;

//! @brief Global variable to enable or disable progress indicators.
/*!
 * When `indicator` is set to `true`, progress-related print commands are executed.
 */
extern bool indicator;

//! @brief Debug print macro that respects the global debug flag.
#define DEBUG_PRINT(msg) do { if (debug) { std::cerr << msg << std::endl; } } while(0)

//! @brief Prints information about a Delaunay facet.
/*!
 * This function outputs the vertices of a given Delaunay facet to the console.
 * 
 * @param f The facet to be printed. A facet is represented as a pair of a cell handle
 *          and an index of the facet within the cell.
 */
void print_facet(Facet f);

//! @brief Prints information about a Delaunay cell.
/*!
 * This function outputs the vertices of a given Delaunay cell to the standard error stream.
 * 
 * @param c The cell to be printed.
 */
void print_cell(Delaunay::Cell c);

//! @brief Writes dummy points to a CSV file for debugging purposes.
/*!
 * The output CSV file includes the grid dimensions, spacing, and the coordinates of the dummy points.
 * 
 * @param grid The grid containing the dimensions and spacing information.
 * @param dummy_points A vector of points to be written to the CSV file.
 */
void write_dummy_points(UnifiedGrid &grid, std::vector<Point> dummy_points);

#endif // VDC_DEBUG_H
