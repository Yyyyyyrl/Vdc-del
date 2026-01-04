//! @file vdc_utilities.cpp
//! @brief Implementation of utility functions for Voronoi and Delaunay computations.

#include "core/vdc_utilities.h"

// ============================================================================
// Bipolar Check
// ============================================================================

bool is_bipolar(float val1, float val2, float isovalue) {
    return ((val1 < isovalue) && (val2 >= isovalue)) || ((val1 >= isovalue) && (val2 < isovalue));
}
