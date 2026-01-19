cmake_minimum_required(VERSION 3.12)
project(VDC-DEL
    VERSION 1.0
    DESCRIPTION "Delaunay-based dual contouring for isosurface extraction"
    LANGUAGES CXX
)

# ─── Compiler settings ────────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()
# Backfill optimization flags when CMake leaves them empty (seen on macOS/AppleClang with CMake 4.1).
if(CMAKE_CXX_FLAGS_RELEASE STREQUAL "")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
endif()
if(CMAKE_CXX_FLAGS_RELWITHDEBINFO STREQUAL "")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")
endif()
if(CMAKE_CXX_FLAGS_MINSIZEREL STREQUAL "")
  set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG")
endif()
if(CMAKE_CXX_FLAGS_DEBUG STREQUAL "")
  set(CMAKE_CXX_FLAGS_DEBUG "-g")
endif()

# ─── Find dependencies ────────────────────────────────────────────────────────
find_package(CGAL REQUIRED COMPONENTS Core)
find_package(ZLIB REQUIRED)

# Teem doesn't ship a CMake config, so try:
#  1) an optional TEEM_ROOT hint
#  2) fallback to find_path/find_library
option(TEEM_ROOT "Root directory of your Teem install" "")

if(TEEM_ROOT)
  set(TEEM_INCLUDE_DIR  "${TEEM_ROOT}/include")
  set(TEEM_LIBRARY      "${TEEM_ROOT}/lib/libteem${CMAKE_SHARED_LIBRARY_SUFFIX}")
else()
  find_path(TEEM_INCLUDE_DIR NAMES teem/ten.h
            DOC "Where to find teemplate headers (e.g. ten.h)")
  find_library(TEEM_LIBRARY NAMES teem
               DOC "Path to libteem")
endif()

if(NOT TEEM_INCLUDE_DIR OR NOT TEEM_LIBRARY)
  message(FATAL_ERROR
    "Could not locate Teem!
     Please install Teem or set -DTEEM_ROOT=/path/to/teem-prefix")
endif()

message(STATUS "Found Teem include: ${TEEM_INCLUDE_DIR}")
message(STATUS "Found Teem lib:     ${TEEM_LIBRARY}")

# ─── Include & link directories ──────────────────────────────────────────────
include_directories(${TEEM_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/include)

# ─── Core source files ───────────────────────────────────────────────────────
set(CORE_SOURCES
    src/core/vdc_commandline.cpp
    src/core/vdc_utilities.cpp
    src/core/vdc_stats.cpp
    src/core/vdc_timing.cpp
    src/vdc_debug.cpp
)

# ─── Processing source files ─────────────────────────────────────────────────
set(PROCESSING_SOURCES
    src/processing/vdc_grid.cpp
    src/processing/vdc_func_delaunay.cpp
    src/processing/vdc_sep_isov.cpp
    src/processing/vdc_refinement.cpp
    src/processing/vdc_del_cells.cpp
    src/processing/vdc_del_cycles.cpp
    src/processing/vdc_del_triangles.cpp
)

# ─── I/O source files ────────────────────────────────────────────────────────
set(IO_SOURCES
    src/vdc_io.cpp
)

# ─── Main executable ─────────────────────────────────────────────────────────
add_executable(vdc-del
    src/vdc_del_main.cpp
    ${CORE_SOURCES}
    ${PROCESSING_SOURCES}
    ${IO_SOURCES}
)

# ─── Link libraries ──────────────────────────────────────────────────────────
target_link_libraries(vdc-del
    PRIVATE
      CGAL::CGAL
      ZLIB::ZLIB
      ${TEEM_LIBRARY}
)

# ─── Installation (optional) ─────────────────────────────────────────────────
install(TARGETS vdc-del
        RUNTIME DESTINATION bin
)

# ─── Usage info ──────────────────────────────────────────────────────────────
message(STATUS "To build:")
message(STATUS "  mkdir build && cd build")
message(STATUS "  cmake ..")
message(STATUS "  make")
message(STATUS "  ./vdc-del [your options]")
