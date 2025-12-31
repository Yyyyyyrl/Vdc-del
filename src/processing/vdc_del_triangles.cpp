//! @file vdc_del_triangles.cpp
//! @brief Implementation of triangle generation for Delaunay-based isosurface extraction.
//!
//! This file implements Step 10 of the Delaunay-based VDC algorithm:
//! - Generate isosurface triangles from marked isosurface facets

#include "processing/vdc_del_isosurface.h"
#include "core/vdc_debug.h"

// ============================================================================
// Step 10: Generate Isosurface Triangles
// ============================================================================

void generate_isosurface_triangles(
    const Delaunay& dt,
    DelaunayIsosurface& iso_surface
) {
    DEBUG_PRINT("[DEL-TRI] Generating isosurface triangles...");

    iso_surface.clear();

    int isovertex_idx = 0;

    iso_surface.isovertices.reserve(dt.number_of_vertices());
    iso_surface.isovertex_delaunay_vertex.reserve(dt.number_of_vertices());
    iso_surface.isovertex_cycle_index.reserve(dt.number_of_vertices());
    iso_surface.vertex_cycle_to_isovertex.reserve(dt.number_of_vertices());

    // Collect all isovertices from active vertices
    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
        if (!vit->info().active) continue;

        int v_idx = vit->info().index;
        const auto& cycle_isovertices = vit->info().cycle_isovertices;

        for (size_t c = 0; c < cycle_isovertices.size(); ++c) {
            iso_surface.isovertices.push_back(cycle_isovertices[c]);
            iso_surface.isovertex_delaunay_vertex.push_back(v_idx);
            iso_surface.isovertex_cycle_index.push_back(static_cast<int>(c));
            iso_surface.vertex_cycle_to_isovertex[{v_idx, static_cast<int>(c)}] = isovertex_idx++;
        }
    }

    DEBUG_PRINT("[DEL-TRI] Collected " << iso_surface.isovertices.size() << " isovertices");
    iso_surface.triangles.reserve(iso_surface.isovertices.size() * 2);

    // Generate triangles from isosurface facets
    int degenerate_count = 0;

    for (auto cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit) {
        for (int i = 0; i < 4; ++i) {
            if (!cit->info().facet_is_isosurface[i]) continue;

            // Get the 3 vertices of facet i (opposite vertex i)
            // Using CGAL's convention: facet i has vertices at positions (i+1)%4, (i+2)%4, (i+3)%4
            std::array<Vertex_handle, 3> facet_verts;
            int k = 0;
            for (int j = 0; j < 4; ++j) {
                if (j != i) {
                    facet_verts[k++] = cit->vertex(j);
                }
            }

            // For each vertex, find the cycle containing this facet
            // and look up the corresponding isovertex
            std::array<int, 3> tri_isovertices;
            bool valid_triangle = true;

            int cell_idx = cit->info().index;

            for (int j = 0; j < 3; ++j) {
                Vertex_handle v = facet_verts[j];
                int v_idx = v->info().index;

                // Find which cycle this facet belongs to for this vertex
                int cycle_idx = find_cycle_containing_facet(v, cell_idx, i);

                if (cycle_idx < 0) {
                    // This shouldn't happen if facet marking and cycle detection are correct
                    DEBUG_PRINT("[DEL-TRI] Warning: Facet (" << cell_idx << ", " << i
                                << ") not found in any cycle of vertex " << v_idx);
                    valid_triangle = false;
                    break;
                }

                // Look up the isovertex index
                auto it = iso_surface.vertex_cycle_to_isovertex.find({v_idx, cycle_idx});
                if (it == iso_surface.vertex_cycle_to_isovertex.end()) {
                    DEBUG_PRINT("[DEL-TRI] Warning: No isovertex for vertex "
                                << v_idx << " cycle " << cycle_idx);
                    valid_triangle = false;
                    break;
                }

                tri_isovertices[j] = it->second;
            }

            if (!valid_triangle) {
                degenerate_count++;
                continue;
            }

            // Check for degenerate triangle (all same vertex)
            if (tri_isovertices[0] == tri_isovertices[1] ||
                tri_isovertices[1] == tri_isovertices[2] ||
                tri_isovertices[0] == tri_isovertices[2]) {
                degenerate_count++;
                continue;
            }

            // Ensure correct orientation
            // The facet is from the positive cell, so the orientation should give
            // outward normal pointing from positive to negative region.
            //
            // CGAL's facet ordering: for facet i, vertices are at (i+1, i+2, i+3) mod 4
            // The outward normal (away from vertex i) follows the right-hand rule.
            //
            // Need to check if i is even or odd for orientation:
            // - Even i: vertices in order give outward normal
            // - Odd i: vertices reversed give outward normal
            //
            // Using standard CGAL convention:
            if (i % 2 == 0) {
                // Even: (i+1, i+2, i+3) gives outward normal
                iso_surface.triangles.push_back({
                    tri_isovertices[0],
                    tri_isovertices[1],
                    tri_isovertices[2]
                });
            } else {
                // Odd: reverse order for outward normal
                iso_surface.triangles.push_back({
                    tri_isovertices[0],
                    tri_isovertices[2],
                    tri_isovertices[1]
                });
            }
        }
    }

    DEBUG_PRINT("[DEL-TRI] Generated " << iso_surface.triangles.size() << " triangles"
                << (degenerate_count > 0 ? " (" + std::to_string(degenerate_count) + " degenerate skipped)" : ""));
}
