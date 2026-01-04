//! @file vdc_utilities.h
//! @brief Utility functions for Voronoi and Delaunay computations.

#ifndef VDC_UTILITIES_H
#define VDC_UTILITIES_H

#include "core/vdc_type.h"

// ============================================================================
// Bipolar Check
// ============================================================================

/**
 * @brief Checks if two scalar values are bipolar with respect to an isovalue.
 *
 * Determines if one value is above and the other is below the specified isovalue.
 *
 * @param val1 First scalar value.
 * @param val2 Second scalar value.
 * @param isovalue The isovalue used for comparison (default is 0).
 * @return `true` if the values are bipolar, otherwise `false`.
 */
bool is_bipolar(float val1, float val2, float isovalue = 0);

#endif // VDC_UTILITIES_H
