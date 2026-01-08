# Delaunay-Based Dual Contouring (VDC-del) - 3D Isosurface Extraction

## Overview

VDC-del is a C++ implementation of a Delaunay-based dual contouring method for high-quality isosurface extraction from 3D scalar fields. Unlike standard dual contouring which uses a grid-based approach, this algorithm works directly in the Delaunay domain to generate high-quality triangle meshes.

## Installation

### Prerequisites
- CMake (3.12+)
- CGAL (5.x) with Core component
- Zlib
- Teem (for NRRD/NHDR IO)

### Building from Source

```bash
git clone <repository-url>
cd Vdc-del
mkdir build && cd build
cmake .. -DTEEM_ROOT=/path/to/teem  # Optional if Teem is in default path
make -j
```

## Usage

### Basic Command Line Syntax
```bash
./vdc-del [OPTIONS] <isovalue> <(nhdr/nrrd) raw data file path>
```

### Options
- `-o {output_filename}`: Specify output filename (default: derived from input filename)
- `-off`: Generate output in .off format (default)
- `-ply`: Generate output in .ply format
- `-sep_dist {D}`: Separation distance (D) in cubes/refined subcubes (default 1 = off)
- `-sep_split {K}`: Number of splits (K) per axis, using a refined grid of factor K+1 (default 0 = no split)
- `-supersample {factor}`: Supersample the input data by the given factor
- `-position_delv_on_isov`: Position Delaunay vertices exactly on the isosurface crossing point
- `-refine_small_angles`: Insert extra Delaunay sites at circumsphere centers to improve small angles near the isosurface
- `-min_angle {deg}`: Trigger refinement if any isosurface-facet triangle angle (per-site iso-sample) is below this threshold (default: 20; also used when `-refine_small_angles` is set without `-min_angle`/`-max_angle`)
- `-max_angle {deg}`: Trigger refinement if any isosurface-facet triangle angle (per-site iso-sample) is above this threshold (default: off)
- `-terse`: Print only the number of vertices/triangles and the output file
- `-multi_isov`: Use multi iso-vertices mode (default)
- `-single_isov`: Use single iso-vertices mode
- `-timing_stats`: Print timing statistics after the run
- `-debug`: Enable debug logging
- `-no_modcyc`: Disable modify-cycles pass (enabled by default)
- `-help`: Print help message

### Examples

Basic run (OFF output, multi-isov enabled by default):
```bash
./vdc-del 70.5 ./data/volvis/fuel.nhdr
```

PLY output with supersampling x2:
```bash
./vdc-del -ply -supersample 2 70.5 ./data/volvis/fuel.nhdr
```

Single iso-vertices mode with separation:
```bash
./vdc-del -single_isov -sep_dist 2 -sep_split 2 70.5 ./data/volvis/fuel.nhdr
```

Terse output (useful for batch runs):
```bash
./vdc-del -terse 70.5 ./data/volvis/fuel.nhdr
```

Angle-based refinement (adds extra Delaunay sites near the isosurface):
```bash
./vdc-del -refine_small_angles -min_angle 20 70.5 ./data/volvis/fuel.nhdr
```

## File Structure

```
.
├── CMakeLists.txt            # CMake build configuration
├── README.md                 # Project documentation
├── include/                  # Header files (15 total)
│   ├── core/                 # Core infrastructure (7 files)
│   │   ├── vdc.h             # Main umbrella header
│   │   ├── vdc_type.h        # Type definitions and CGAL types
│   │   ├── vdc_utilities.h   # Utility functions and helpers
│   │   ├── vdc_debug.h       # Debug output and logging
│   │   ├── vdc_commandline.h # CLI parsing
│   │   ├── vdc_timing.h      # Performance timing
│   │   └── vdc_stats.h       # Statistics collection and reporting
│   ├── processing/           # Algorithm headers (7 files)
│   │   ├── vdc_func.h        # High-level algorithm orchestration
│   │   ├── vdc_delaunay.h    # Delaunay triangulation structures
│   │   ├── vdc_del_isosurface.h  # Isosurface extraction structures
│   │   ├── vdc_del_cycles.h  # Cycle detection and self-intersection
│   │   ├── vdc_grid.h        # Grid data structures
│   │   ├── vdc_refinement.h  # Angle-based refinement
│   │   └── vdc_sep_isov.h    # Isosurface vertex separation methods
│   └── vdc_io.h              # I/O operations (NRRD, OFF, PLY)
├── src/                      # Source files (14 total)
│   ├── core/                 # Core implementations (5 files)
│   │   ├── vdc_commandline.cpp
│   │   ├── vdc_debug.cpp
│   │   ├── vdc_utilities.cpp
│   │   ├── vdc_timing.cpp
│   │   └── vdc_stats.cpp
│   ├── processing/           # Algorithm implementations (7 files)
│   │   ├── vdc_func_delaunay.cpp    # Delaunay triangulation construction
│   │   ├── vdc_del_cells.cpp        # Cell circumcenters and facet marking
│   │   ├── vdc_del_cycles.cpp       # Cycle detection and isovertex positioning
│   │   ├── vdc_del_triangles.cpp    # Triangle generation
│   │   ├── vdc_grid.cpp             # Grid operations
│   │   ├── vdc_sep_isov.cpp         # Separation implementations
│   │   └── vdc_refinement.cpp       # Angle refinement
│   ├── vdc_del_main.cpp      # Main application entry point
│   └── vdc_io.cpp            # I/O implementations
```

**Directory Organization:**
- **`core/`**: Infrastructure code (types, utilities, I/O, CLI, debugging, timing, stats)
- **`processing/`**: All computational geometry algorithms (Delaunay, grid, isosurface extraction)

## API Documentation

The project provides comprehensive API documentation through Doxygen-style comments. Key modules:

- **vdc_del_cycles.h**: Cycle detection, isovertex computation, self-intersection resolution
- **vdc_del_isosurface.h**: Isosurface data structures and extraction functions
- **vdc_grid.h**: Grid loading, supersampling, and interpolation

To generate documentation:
```bash
doxygen Doxyfile
```
