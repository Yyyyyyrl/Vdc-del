//! @file vdc_del_cells.cpp
//! @brief Implementation of cell-level processing for Delaunay-based isosurface extraction.
//!
//! This file implements Steps 5-7 of the Delaunay-based VDC algorithm:
//! - Step 5: Compute circumcenters and scalar values for Delaunay cells
//! - Step 6: Mark isosurface facets based on cell signs
//! - Step 7: Mark active vertices incident on isosurface facets

#include "processing/vdc_del_isosurface.h"
#include "core/vdc_debug.h"

// ============================================================================
// Step 5: Compute Cell Circumcenters and Scalars
// ============================================================================

void compute_cell_circumcenters_and_scalars(
    Delaunay& dt,
    const UnifiedGrid& grid,
    float isovalue
) {
    DEBUG_PRINT("[DEL-CELL] Computing circumcenters and scalar values for "
                << dt.number_of_finite_cells() << " finite cells");

    int positive_count = 0;
    int negative_count = 0;

    for (auto cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit) {
        // Step 5a: Compute circumcenter using CGAL's built-in method
        // The Delaunay triangulation stores circumcenters if using the appropriate cell base
        Point circumcenter = cit->circumcenter();
        cit->info().circumcenter = circumcenter;

        // Step 5b: Compute scalar value at circumcenter via trilinear interpolation
        float scalar = trilinear_interpolate(circumcenter, grid);
        cit->info().circumcenter_scalar = scalar;

        // Step 5c-d: Set flag_positive based on comparison with isovalue
        cit->info().flag_positive = (scalar >= isovalue);

        if (cit->info().flag_positive) {
            positive_count++;
        } else {
            negative_count++;
        }
    }

    DEBUG_PRINT("[DEL-CELL] Cell classification: "
                << positive_count << " positive, "
                << negative_count << " negative");
}

// ============================================================================
// Step 6: Mark Isosurface Facets
// ============================================================================

void mark_isosurface_facets(Delaunay& dt) {
    DEBUG_PRINT("[DEL-FACET] Marking isosurface facets...");

    int isosurface_facet_count = 0;

    for (auto cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit) {
        for (int i = 0; i < 4; ++i) {
            // Get neighbor cell through facet i
            Cell_handle neighbor = cit->neighbor(i);

            // Skip if neighbor is infinite (boundary facet)
            if (dt.is_infinite(neighbor)) {
                cit->info().facet_is_isosurface[i] = false;
                continue;
            }

            // Step 6a: Check if this cell is positive and neighbor is negative
            bool cellA_positive = cit->info().flag_positive;
            bool cellB_positive = neighbor->info().flag_positive;

            // Mark as isosurface only from positive side to avoid double processing
            // This ensures consistent orientation: outward normal points from positive to negative
            if (cellA_positive && !cellB_positive) {
                cit->info().facet_is_isosurface[i] = true;
                isosurface_facet_count++;
            } else {
                cit->info().facet_is_isosurface[i] = false;
            }
        }
    }

    DEBUG_PRINT("[DEL-FACET] Marked " << isosurface_facet_count << " isosurface facets");
}

// ============================================================================
// Step 7: Mark Active Vertices
// ============================================================================

void mark_active_vertices(Delaunay& dt) {
    DEBUG_PRINT("[DEL-VERT] Marking active vertices...");

    // Initialize all vertices as inactive
    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
        vit->info().active = false;
    }

    int active_count = 0;

    // Check each isosurface facet and mark incident vertices as active
    for (auto cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit) {
        for (int i = 0; i < 4; ++i) {
            if (cit->info().facet_is_isosurface[i]) {
                // Facet i is opposite to vertex i
                // So the vertices of facet i are vertices j where j != i
                for (int j = 0; j < 4; ++j) {
                    if (j != i) {
                        Vertex_handle v = cit->vertex(j);
                        if (!v->info().active) {
                            v->info().active = true;
                            active_count++;
                        }
                    }
                }
            }
        }
    }

    DEBUG_PRINT("[DEL-VERT] Marked " << active_count << " active vertices");
}
