//! @file vdc_del_cycles.h
//! @brief Header for cycle detection and isovertex computation in Delaunay-based isosurface extraction.
//!
//! This header provides data structures and function declarations for:
//! - Step 8: Computing cycles of isosurface facets around active Delaunay vertices
//! - Step 9: Computing isosurface vertex positions for each cycle
//! - Self-intersection detection and resolution for multi-cycle vertices

#ifndef VDC_DEL_CYCLES_H
#define VDC_DEL_CYCLES_H

#include "core/vdc_type.h"
#include "processing/vdc_grid.h"
#include <vector>
#include <cstdint>

// ============================================================================
// Facet and Edge Key Types
// ============================================================================

/**
 * @brief Identifies a Delaunay facet by its positive cell index and local facet index.
 *
 * Isosurface facets are stored using the convention that the facet is represented
 * from the positive cell side (where scalar >= isovalue). This ensures consistent
 * orientation across the algorithm.
 */
struct FacetKey {
    int cell_index = -1;   ///< Index of the positive cell containing this facet
    int facet_index = -1;  ///< Local facet index within the cell (0-3)

    /**
     * @brief Equality comparison operator.
     * @param other The FacetKey to compare against.
     * @return true if both cell_index and facet_index match.
     */
    bool operator==(const FacetKey& other) const {
        return cell_index == other.cell_index && facet_index == other.facet_index;
    }
};

/**
 * @brief Hash functor for FacetKey, enabling use in unordered containers.
 */
struct FacetKeyHash {
    /**
     * @brief Computes hash value for a FacetKey.
     * @param key The FacetKey to hash.
     * @return Hash value combining cell_index and facet_index.
     */
    size_t operator()(const FacetKey& key) const noexcept {
        return (static_cast<size_t>(static_cast<uint32_t>(key.cell_index)) << 32) ^
               static_cast<size_t>(static_cast<uint32_t>(key.facet_index));
    }
};

/**
 * @brief Identifies a Delaunay edge by its two vertex indices (ordered v0 < v1).
 *
 * This provides a canonical representation of edges for use as map keys,
 * ensuring the same edge is identified regardless of vertex order.
 */
struct EdgeKey {
    int v0 = -1;  ///< Lower vertex index
    int v1 = -1;  ///< Higher vertex index

    /**
     * @brief Factory method to create an EdgeKey from two vertex indices.
     *
     * Automatically orders vertices so v0 < v1.
     *
     * @param a First vertex index.
     * @param b Second vertex index.
     * @return EdgeKey with vertices in canonical order.
     */
    static EdgeKey FromVertexIndices(int a, int b) {
        if (a < b) {
            return EdgeKey{a, b};
        }
        return EdgeKey{b, a};
    }

    /**
     * @brief Equality comparison operator.
     * @param other The EdgeKey to compare against.
     * @return true if both vertex indices match.
     */
    bool operator==(const EdgeKey& other) const {
        return v0 == other.v0 && v1 == other.v1;
    }
};

/**
 * @brief Hash functor for EdgeKey, enabling use in unordered containers.
 */
struct EdgeKeyHash {
    /**
     * @brief Computes hash value for an EdgeKey.
     * @param key The EdgeKey to hash.
     * @return Hash value combining v0 and v1.
     */
    size_t operator()(const EdgeKey& key) const noexcept {
        return (static_cast<size_t>(static_cast<uint32_t>(key.v0)) << 32) ^
               static_cast<size_t>(static_cast<uint32_t>(key.v1));
    }
};

// ============================================================================
// Facet Cycle Matching Data
// ============================================================================

/**
 * @brief Stores per-vertex-slot matching data for an isosurface facet.
 *
 * Each facet has 3 vertices. For each vertex slot (0, 1, 2), this structure
 * records which facet is matched to form cycles around that vertex.
 */
struct FacetMatchingData {
    FacetKey cycleMatchingFacet[3];                ///< Matched facet for each vertex slot
    bool hasMatch[3] = {false, false, false};      ///< Whether a match exists for each slot
};

// ============================================================================
// Self-Intersection Resolution Types
// ============================================================================

/**
 * @brief Outcome status of self-intersection resolution.
 */
enum class ResolutionStatus {
    NOT_NEEDED,       ///< No self-intersection detected
    RESOLVED_SPHERE,  ///< Resolved by sphere-based projection
    RESOLVED_FALLBACK,///< Resolved by fallback distribution strategy
    UNRESOLVED        ///< Could not resolve self-intersection
};

/**
 * @brief Strategy used to resolve self-intersection.
 */
enum class ResolutionStrategy {
    NONE,                  ///< No strategy applied
    POSITION_MULTI_ISOV,   ///< PositionMultiIsov separation tests
    TWO_CYCLE_DIAMETRIC,   ///< Two cycles placed diametrically opposite
    FIBONACCI_FALLBACK     ///< Fibonacci spherical distribution fallback
};

/**
 * @brief Result of a self-intersection resolution attempt.
 */
struct ResolutionResult {
    ResolutionStatus status = ResolutionStatus::NOT_NEEDED;    ///< Resolution outcome
    ResolutionStrategy strategy = ResolutionStrategy::NONE;    ///< Strategy that succeeded
};

/**
 * @brief Stores matching method state for a Delaunay edge, used in cycle modification.
 *
 * The modify-cycles pass can flip the bipolar matching method (SEP_NEG <-> SEP_POS)
 * to resolve non-manifold configurations.
 */
struct EdgeMatchingState {
    BIPOLAR_MATCH_METHOD method = BIPOLAR_MATCH_METHOD::SEP_NEG;  ///< Current matching method
    int unconstrained_offset = 0;  ///< Offset for unconstrained matching
    bool flipped = false;          ///< Whether matching was flipped this pass
    int flip_count = 0;            ///< Total flip count to prevent oscillation
};

// ============================================================================
// Statistics Tracking
// ============================================================================

/**
 * @brief Statistics collected during isovertex computation.
 *
 * Tracks counts of single/multi-cycle vertices and self-intersection
 * resolution outcomes for reporting.
 */
struct IsovertexComputationStats {
    int single_cycle_count = 0;              ///< Number of single-cycle vertices
    int multi_cycle_count = 0;               ///< Number of multi-cycle vertices
    int self_intersection_detected = 0;      ///< Total self-intersections detected
    int self_intersection_resolved = 0;      ///< Resolved by sphere projection
    int self_intersection_fallback = 0;      ///< Resolved by fallback strategy
    int self_intersection_unresolved = 0;    ///< Unresolved self-intersections
    bool had_unresolved = false;             ///< Whether any case remained unresolved

    // Per-strategy counts
    int64_t strat_position_multi_isov = 0;   ///< Uses of POSITION_MULTI_ISOV strategy
    int64_t strat_two_cycle_diametric = 0;   ///< Uses of TWO_CYCLE_DIAMETRIC strategy
    int64_t strat_fibonacci_fallback = 0;    ///< Uses of FIBONACCI_FALLBACK strategy
};

// ============================================================================
// Core Cycle Functions
// ============================================================================

/**
 * @brief Compute cycles of isosurface facets around each active vertex.
 *
 * For each active Delaunay vertex, finds all isosurface facets incident on it
 * and groups them into cycles (connected components) using bipolar edge matching.
 *
 * @param dt The Delaunay triangulation (must have active vertices marked).
 */
void compute_facet_cycles(Delaunay& dt);

/**
 * @brief Compute isosurface vertex positions for each cycle around each active vertex.
 *
 * For single-cycle vertices: positions at the site-associated isosurface sample.
 * For multi-cycle vertices: computes centroids of Voronoi-edge/isovalue intersections
 * per cycle and projects them onto a sphere around the Delaunay site.
 *
 * @param dt The Delaunay triangulation (must have cycles computed).
 * @param grid The scalar field grid.
 * @param isovalue The isovalue threshold.
 * @param position_on_isov If true, Delaunay sites are on isosurface samples.
 */
void compute_cycle_isovertices(
    Delaunay& dt,
    const UnifiedGrid& grid,
    float isovalue,
    bool position_on_isov
);

// ============================================================================
// Helper Functions
// ============================================================================

/**
 * @brief Find which cycle contains a given facet for a vertex.
 *
 * @param v_handle Handle to the Delaunay vertex.
 * @param cell_index Index of the cell containing the facet.
 * @param facet_index Index of the facet within the cell (0-3).
 * @return Cycle index (0-based), or -1 if not found.
 */
int find_cycle_containing_facet(
    Vertex_handle v_handle,
    int cell_index,
    int facet_index
);

#endif // VDC_DEL_CYCLES_H
