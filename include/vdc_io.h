//! @file vdc_io.h
//! @brief Header file for input/output operations, including mesh export.

#ifndef VDC_IO_H
#define VDC_IO_H

#include "core/vdc_utilities.h"
#include "processing/vdc_del_isosurface.h"

//! @brief Writes a DelaunayIsosurface mesh in OFF format.
/*!
 * @param filename The output file path.
 * @param iso_surface The isosurface container providing vertices and triangles.
 * @return `true` on success.
 */
bool write_off_delaunay(const std::string &filename, const DelaunayIsosurface &iso_surface);

//! @brief Writes a DelaunayIsosurface mesh in PLY format.
/*!
 * @param filename The output file path.
 * @param iso_surface The isosurface container providing vertices and triangles.
 * @return `true` on success.
 */
bool write_ply_delaunay(const std::string &filename, const DelaunayIsosurface &iso_surface);

#endif // VDC_IO_H
