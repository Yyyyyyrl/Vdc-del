//! @file vdc_del_cells.cpp
//! @brief Implementation of cell-level processing for Delaunay-based isosurface extraction.
//!
//! This file implements Steps 5-7 of the Delaunay-based VDC algorithm:
//! - Step 5: Compute circumcenters and scalar values for Delaunay cells
//! - Step 6: Mark isosurface facets based on cell signs
//! - Step 7: Mark active vertices incident on isosurface facets

#include "processing/vdc_del_isosurface.h"
#include "core/vdc_debug.h"
#include <array>
#include <cmath>
#include <limits>
#include <vector>

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
// Optional: Flip cell signs to avoid small isofacet dihedrals
// ============================================================================

static double vec_dot(const Vector3& a, const Vector3& b) {
    return a.x() * b.x() + a.y() * b.y() + a.z() * b.z();
}

static Vector3 vec_cross(const Vector3& a, const Vector3& b) {
    return Vector3(
        a.y() * b.z() - a.z() * b.y(),
        a.z() * b.x() - a.x() * b.z(),
        a.x() * b.y() - a.y() * b.x());
}

static double vec_norm(const Vector3& v0) {
    return std::sqrt(v0.squared_length());
}

static double clamp_unit(double x) {
    if (x < -1.0) return -1.0;
    if (x > 1.0) return 1.0;
    return x;
}

static bool cell_has_infinite_or_dummy_vertices(const Delaunay& dt, const Cell_handle& cell) {
    for (int i = 0; i < 4; ++i) {
        Vertex_handle vh = cell->vertex(i);
        if (vh == Vertex_handle() || dt.is_infinite(vh) || vh->info().is_dummy) {
            return true;
        }
    }
    return false;
}

int flip_cell_signs_for_small_isofacet_dihedral(
    Delaunay& dt,
    double cos_dihedral_angle_threshold
) {
    // Clamp threshold to valid range.
    cos_dihedral_angle_threshold = clamp_unit(cos_dihedral_angle_threshold);

    std::vector<Cell_handle> to_flip;
    to_flip.reserve(1024);

    int candidate_cells = 0;

    for (auto cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit) {
        Cell_handle cell = cit;
        if (cell == Cell_handle() || dt.is_infinite(cell)) {
            continue;
        }
        if (cell_has_infinite_or_dummy_vertices(dt, cell)) {
            continue;
        }

        std::array<int, 4> isofaces{};
        int isoface_count = 0;
        bool has_infinite_neighbor = false;

        const bool cell_pos = cell->info().flag_positive;
        for (int fi = 0; fi < 4; ++fi) {
            Cell_handle neigh = cell->neighbor(fi);
            if (neigh == Cell_handle() || dt.is_infinite(neigh)) {
                has_infinite_neighbor = true;
                break;
            }
            if (neigh->info().flag_positive != cell_pos) {
                isofaces[isoface_count++] = fi;
            }
        }

        // Only apply to interior cells (all neighbors finite) with 3 or 4 isosurface facets.
        if (has_infinite_neighbor) {
            continue;
        }
        if (isoface_count < 3) {
            continue;
        }

        ++candidate_cells;

        // Gather tetrahedron vertex positions.
        Point p[4];
        for (int i = 0; i < 4; ++i) {
            p[i] = cell->vertex(i)->point();
        }

        // Compute outward normals for each face (face index is opposite vertex index).
        Vector3 face_normal[4];
        bool degenerate = false;
        for (int fi = 0; fi < 4; ++fi) {
            const int a = CellInfo::FacetVertexIndex(fi, 0);
            const int b = CellInfo::FacetVertexIndex(fi, 1);
            const int c = CellInfo::FacetVertexIndex(fi, 2);

            Vector3 n = vec_cross(Vector3(p[a], p[b]), Vector3(p[a], p[c]));
            if (vec_norm(n) < 1e-18) {
                degenerate = true;
                break;
            }

            // Orient outward: ensure normal points away from the opposite vertex p[fi].
            if (vec_dot(n, Vector3(p[a], p[fi])) > 0.0) {
                n = -n;
            }
            face_normal[fi] = n;
        }

        if (degenerate) {
            continue;
        }

        bool trigger_flip = false;
        for (int a = 0; a < isoface_count && !trigger_flip; ++a) {
            for (int b = a + 1; b < isoface_count; ++b) {
                const int f0 = isofaces[a];
                const int f1 = isofaces[b];

                const double n0 = vec_norm(face_normal[f0]);
                const double n1 = vec_norm(face_normal[f1]);
                if (n0 < 1e-18 || n1 < 1e-18) {
                    continue;
                }

                // Internal dihedral angle = pi - angle_between(outward_normals).
                // So cos(dihedral) = -cos(angle_between(outward_normals)).
                const double cos_angle = clamp_unit(vec_dot(face_normal[f0], face_normal[f1]) / (n0 * n1));
                const double cos_dihedral = -cos_angle;

                if (cos_dihedral > cos_dihedral_angle_threshold) {
                    trigger_flip = true;
                    break;
                }
            }
        }

        if (trigger_flip) {
            to_flip.push_back(cell);
        }
    }

    for (Cell_handle cell : to_flip) {
        cell->info().flag_positive = !cell->info().flag_positive;
    }

    DEBUG_PRINT("[DEL-CELL] Small-dihedral flip pass: candidates=" << candidate_cells
                << ", flipped=" << to_flip.size());

    return static_cast<int>(to_flip.size());
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
