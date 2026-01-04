//! @file vdc_delaunay.h
//! @brief Header file for Delaunay triangulation support types.
//!
//! This header provides type definitions used across Delaunay-based
//! isosurface extraction. The main data structures are defined in
//! vdc_type.h (VertexInfo, CellInfo) and vdc_del_isosurface.h (DelaunayIsosurface).

#ifndef VDC_DELAUNAY_H
#define VDC_DELAUNAY_H

#include "core/vdc_type.h"

// Note: This header is kept for compatibility. The main isosurface data
// structure (DelaunayIsosurface) is defined in vdc_del_isosurface.h.
// Legacy triangle structs (DelaunayTriangle, IsoTriangle) have been removed
// as they were not used in the current codebase.

#endif // VDC_DELAUNAY_H
