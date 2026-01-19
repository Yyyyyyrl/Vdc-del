//! @file vdc_del_isosurface.h
//! @brief Data structures and high-level function declarations for Delaunay-based isosurface extraction.
//!
//! This header provides the main API for isosurface extraction:
//! - DelaunayIsosurface: Output data structure
//! - High-level processing functions (compute_cell_circumcenters_and_scalars, etc.)
//!
//! For detailed cycle computation types and helpers, see vdc_del_cycles.h.

#ifndef VDC_DEL_ISOSURFACE_H
#define VDC_DEL_ISOSURFACE_H

#include "core/vdc_type.h"
#include "processing/vdc_grid.h"
#include "processing/vdc_del_cycles.h"  // For shared types (FacetKey, ResolutionResult, etc.)
#include <unordered_map>
#include <array>

// ============================================================================
// Output Data Structure
// ============================================================================

/**
 * @brief Structure to hold the extracted isosurface from Delaunay-based algorithm.
 *
 * Contains isosurface vertices (one or more per active Delaunay vertex),
 * triangles, and mappings between Delaunay structures and isosurface elements.
 */
struct DelaunayIsosurface {
    /**
     * @brief Hash functor for (vertex_index, cycle_index) pairs.
     */
    struct VertexCycleKeyHash {
        std::size_t operator()(const std::pair<int, int>& key) const noexcept {
            const std::size_t h1 = std::hash<int>{}(key.first);
            const std::size_t h2 = std::hash<int>{}(key.second);
            return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
        }
    };

    std::vector<Point> isovertices;                      ///< Isosurface vertex positions
    std::vector<int> isovertex_delaunay_vertex;          ///< Delaunay vertex index for each isovertex
    std::vector<int> isovertex_cycle_index;              ///< Cycle index for each isovertex
    std::vector<std::array<int, 3>> triangles;           ///< Output triangles (indices into isovertices)
    std::array<double, 3> vertex_scale{1.0, 1.0, 1.0};   ///< Per-axis scale to convert from grid units to physical space

    //! @brief Mapping from (delaunay_vertex_index, cycle_index) -> isovertex index
    std::unordered_map<std::pair<int, int>, int, VertexCycleKeyHash> vertex_cycle_to_isovertex;

    //! @brief Clear all data
    void clear() {
        isovertices.clear();
        isovertex_delaunay_vertex.clear();
        isovertex_cycle_index.clear();
        triangles.clear();
        vertex_cycle_to_isovertex.clear();
    }

    //! @brief Get the number of isosurface vertices
    size_t num_vertices() const { return isovertices.size(); }

    //! @brief Get the number of triangles
    size_t num_triangles() const { return triangles.size(); }
};

// ============================================================================
// High-Level Processing Functions
// ============================================================================

/**
 * @brief Compute circumcenters and scalar values for all finite Delaunay cells.
 *
 * For each finite Delaunay cell:
 * - Computes the circumcenter (center of circumscribed sphere)
 * - Interpolates the scalar value at the circumcenter from the grid
 * - Sets flag_positive based on whether the scalar >= isovalue
 *
 * @param dt The Delaunay triangulation.
 * @param grid The scalar field grid.
 * @param isovalue The isovalue threshold.
 */
void compute_cell_circumcenters_and_scalars(
    Delaunay& dt,
    const UnifiedGrid& grid,
    float isovalue
);

/**
 * @brief Flip cell signs to avoid very small dihedral angles between isosurface facets.
 *
 * For each finite Delaunay cell (tetrahedron), count the number of facets whose
 * neighboring cell has opposite sign. If that count is 3 or 4, compute internal
 * dihedral angles between the corresponding cell facets. If any dihedral angle
 * is smaller than the threshold (equivalently, cosine is greater than the given
 * cosine threshold), flip the sign of the cell's Voronoi vertex (circumcenter).
 *
 * Notes:
 * - The flip decision is computed from the current cell signs and then applied
 *   simultaneously to avoid order dependence.
 * - Cells with dummy vertices or infinite neighbors are skipped.
 *
 * @param dt The Delaunay triangulation (must have cell signs computed first).
 * @param cos_dihedral_angle_threshold Cosine threshold (default is cos(5 deg)).
 * @return Number of flipped cells.
 */
int flip_cell_signs_for_small_isofacet_dihedral(
    Delaunay& dt,
    double cos_dihedral_angle_threshold
);

/**
 * @brief Mark Delaunay facets that cross the isosurface.
 *
 * A facet is marked as isosurface if it separates a positive cell from a negative cell.
 * Only marks from the positive side to avoid double processing.
 *
 * @param dt The Delaunay triangulation (must have circumcenters computed first).
 */
void mark_isosurface_facets(Delaunay& dt);

/**
 * @brief Mark active Delaunay vertices.
 *
 * A vertex is active if it is incident on any isosurface facet.
 *
 * @param dt The Delaunay triangulation (must have isosurface facets marked first).
 */
void mark_active_vertices(Delaunay& dt);

// ============================================================================
// Cycle Modification
// ============================================================================

/**
 * @brief Result structure for modify_cycles_pass function.
 */
struct ModifyCyclesResult {
    int total_flips = 0;       ///< Total number of edge matching flips
    int problematic_edges = 0; ///< Number of problematic edges found
    int iterations = 0;        ///< Number of fix iterations performed
};

/**
 * @brief Modify cycles to fix non-manifold configurations from problematic matchings.
 *
 * Implements the mod_cyc algorithm:
 * - Detects edges where multiple facets share the same cycle component pair
 * - Flips the bipolar matching method (SEP_NEG <-> SEP_POS) to resolve conflicts
 * - Recomputes cycles for affected vertices
 *
 * @param dt The Delaunay triangulation (must have cycles computed first).
 * @return ModifyCyclesResult containing statistics about the modifications.
 */
ModifyCyclesResult modify_cycles_pass(Delaunay& dt);

// ============================================================================
// Triangle Generation
// ============================================================================

/**
 * @brief Generate isosurface triangles from marked facets.
 *
 * For each isosurface facet, creates a triangle by looking up the isosurface
 * vertex for each of the facet's three Delaunay vertices.
 *
 * @param dt The Delaunay triangulation (must have cycle isovertices computed).
 * @param iso_surface Output structure to receive vertices and triangles.
 */
void generate_isosurface_triangles(
    const Delaunay& dt,
    DelaunayIsosurface& iso_surface
);

#endif // VDC_DEL_ISOSURFACE_H
