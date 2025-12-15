# Plan: Self-Intersection Handling for vdc-del (Step 9.b-c)

## Overview

Implement self-intersection detection and resolution for multi-cycle vertices in the Delaunay-based VDC algorithm. This corresponds to Steps 9.b and 9.c in `vdc-DelaunayBased.txt`.

## Problem Statement

When a Delaunay vertex has multiple cycles of isosurface facets (representing multiple isosurface sheets passing through the same Voronoi cell), the computed isovertex positions for each cycle may cause the resulting triangles to self-intersect. This happens when the projected centroids are positioned such that triangles from different cycles overlap.

## Algorithm from vdc-DelaunayBased.txt

```
9. For each active Delaunay vertex v associated with 2 or more isosurface vertices:
   a. For each cycle C of Delaunay facets around v:
      i.   w = C.isov
      ii.  Compute centroid of Voronoi edge/isosurface intersections
      iii. Project centroid onto sphere around v
      iv.  Set initial position of w to projected centroid
   b. Check whether positions of isosurface vertices w creates self-intersection of cycle facets.
   c. If projection creates self-intersection, compute new positions for isosurface vertices w
      that do not create self-intersection.
```

## Current Implementation Status

**File:** `~/Documents/Vdc-del/src/processing/vdc_del_cycles.cpp`

Currently implemented:
- Step 9.a.i-ii: Centroid computation via `compute_centroid_of_voronoi_edge_and_isosurface()`
- Step 9.a.iii-iv: `project_to_sphere()` exists but not used; centroids used directly

**NOT implemented:**
- Step 9.b: Self-intersection detection
- Step 9.c: Position adjustment to resolve self-intersections

---

## Implementation Plan

### Phase 1: Self-Intersection Detection (Step 9.b)

#### 1.1 Define What Constitutes Self-Intersection

For a multi-cycle vertex v with cycles C1, C2, ..., Cn:
- Each cycle Ci has an isovertex position wi
- Each isosurface facet f in cycle Ci generates a triangle using wi and isovertices from adjacent vertices
- Self-intersection occurs when triangles from different cycles intersect

**Detection approach:** Check if any two triangles from different cycles of the same vertex intersect.

#### 1.2 Collect Triangles Per Cycle

For each cycle C around vertex v:
1. Get all isosurface facets in the cycle
2. For each facet f, determine the triangle vertices:
   - One vertex is the cycle's isovertex (wi)
   - Other two vertices come from the other two Delaunay vertices of facet f
3. Store triangles grouped by cycle

#### 1.3 Triangle-Triangle Intersection Test

Implement or use CGAL's triangle intersection test:
```cpp
bool triangles_intersect(
    const Point& a1, const Point& a2, const Point& a3,  // Triangle A
    const Point& b1, const Point& b2, const Point& b3   // Triangle B
);
```

Use `CGAL::do_intersect()` for Triangle_3 objects.

#### 1.4 Check All Cycle Pairs

```cpp
bool check_self_intersection(
    Vertex_handle v,
    const std::vector<Point>& cycle_isovertices,
    const Delaunay& dt
) {
    // For each pair of cycles (i, j) where i < j
    // Check if any triangle from cycle i intersects any triangle from cycle j
    // Return true if any intersection found
}
```

---

### Phase 2: Self-Intersection Resolution (Step 9.c)

#### 2.1 Resolution Strategy Options

**Option A: Sphere Projection with Separation**
- Project all cycle centroids onto a sphere around v
- Ensure minimum angular separation between projected points
- If too close, push them apart along the sphere surface

**Option B: Iterative Adjustment**
- Start with computed centroids
- If self-intersection detected, move isovertices toward v (shrink)
- Or move them apart from each other
- Iterate until no self-intersection

**Option C: Convex Hull Partitioning**
- Compute convex hull of all cycle centroids
- Assign each isovertex to a vertex of the hull
- Ensures maximum separation

**Recommended: Option A (Sphere Projection with Separation)**
- Aligns with Step 9.a.iii in the algorithm
- Provides geometric guarantee of separation
- Simpler to implement

#### 2.2 Sphere Projection Algorithm

```cpp
void resolve_self_intersection(
    Vertex_handle v,
    std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    double sphere_radius
) {
    Point center = v->point();
    int n = cycle_isovertices.size();

    // 1. Project all centroids onto sphere
    for (int i = 0; i < n; ++i) {
        cycle_isovertices[i] = project_to_sphere(cycle_isovertices[i], center, sphere_radius);
    }

    // 2. Check minimum angular separation
    double min_angle = compute_min_angular_separation(cycle_isovertices, center);
    double required_angle = 2.0 * M_PI / (3.0 * n);  // Heuristic threshold

    // 3. If too close, redistribute on sphere
    if (min_angle < required_angle) {
        redistribute_on_sphere(cycle_isovertices, center, sphere_radius);
    }

    // 4. Verify no self-intersection; if still intersecting, fall back to vertex position
    if (check_self_intersection(v, cycle_isovertices, dt)) {
        // Fallback: place all isovertices at vertex center
        for (int i = 0; i < n; ++i) {
            cycle_isovertices[i] = center;
        }
    }
}
```

#### 2.3 Sphere Radius Selection

The sphere radius should be small enough to keep isovertices near the Delaunay vertex but large enough to provide separation:
- Use a fraction of the minimum edge length incident on v
- Or use a fraction of the grid spacing
- Typical value: `radius = 0.1 * min_incident_edge_length`

---

### Phase 3: Integration

#### 3.1 Modify `compute_cycle_isovertices()`

Update the multi-cycle branch in `vdc_del_cycles.cpp`:

```cpp
} else {
    // Multiple cycles: compute centroids for each cycle
    for (size_t c = 0; c < cycles.size(); ++c) {
        Point centroid = compute_centroid_of_voronoi_edge_and_isosurface(
            dt, grid, isovalue, vit, cycles[c]);
        isovertices[c] = centroid;
    }

    // Step 9.b: Check for self-intersection
    if (check_self_intersection(vit, isovertices, dt)) {
        // Step 9.c: Resolve self-intersection
        double radius = compute_sphere_radius(vit, dt, grid);
        resolve_self_intersection(vit, isovertices, dt, radius);
    }

    multi_cycle_count++;
}
```

#### 3.2 Add New Functions to Header

Update `include/processing/vdc_del_isosurface.h`:

```cpp
//! @brief Check if cycle isovertex positions create self-intersecting triangles
bool check_self_intersection(
    Vertex_handle v,
    const std::vector<Point>& cycle_isovertices,
    const Delaunay& dt
);

//! @brief Resolve self-intersection by adjusting isovertex positions
void resolve_self_intersection(
    Vertex_handle v,
    std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    double sphere_radius
);

//! @brief Compute appropriate sphere radius for a vertex
double compute_sphere_radius(
    Vertex_handle v,
    const Delaunay& dt,
    const UnifiedGrid& grid
);
```

---

## Files to Modify

| File | Changes |
|------|---------|
| `src/processing/vdc_del_cycles.cpp` | Add self-intersection check and resolution functions; integrate into `compute_cycle_isovertices()` |
| `include/processing/vdc_del_isosurface.h` | Add function declarations for new self-intersection handling functions |

---

## Implementation Steps

1. **Add helper function to collect triangles for a cycle**
   - Given a cycle and its isovertex, enumerate all triangles that would be generated

2. **Implement `check_self_intersection()`**
   - Use CGAL's `do_intersect()` for Triangle_3
   - Check all triangle pairs from different cycles

3. **Implement `compute_sphere_radius()`**
   - Compute minimum incident edge length
   - Return a fraction (e.g., 0.1) of that length

4. **Implement `resolve_self_intersection()`**
   - Project centroids to sphere
   - Check angular separation
   - Redistribute if needed
   - Fallback to vertex position if still intersecting

5. **Integrate into `compute_cycle_isovertices()`**
   - Add self-intersection check after computing centroids
   - Call resolution if needed

6. **Add debug output**
   - Report number of vertices with self-intersections detected
   - Report number resolved vs fallback

7. **Test with datasets known to have multi-cycle vertices**
   - Use fuel.nhdr or other complex datasets
   - Verify mesh quality improvement

---

## Testing Strategy

1. **Unit test self-intersection detection**
   - Create known intersecting triangle configurations
   - Verify detection returns true

2. **Integration test**
   - Run on datasets with multi-cycle vertices
   - Compare output mesh quality before/after

3. **Visual inspection**
   - Export meshes and view in mesh viewer
   - Check for visible self-intersections

---

## Estimated Complexity

- Self-intersection detection: ~100 lines
- Resolution algorithm: ~80 lines
- Integration and helpers: ~50 lines
- **Total: ~230 lines of new code**
