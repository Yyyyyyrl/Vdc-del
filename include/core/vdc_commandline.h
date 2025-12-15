//! @file vdc_commandline.h
//! @brief Header file for parsing command-line arguments and providing help messages.

#ifndef VDC_COMMANDLINE_H
#define VDC_COMMANDLINE_H

#include "core/vdc_type.h"
#include "core/vdc_utilities.h"

//! @brief Structure to hold top-level parameters parsed from command-line arguments.
/*!
 * This structure consolidates all configurable parameters for the program,
 * replacing the need for multiple global variables.
 */
struct VdcParam {
    std::string file_path;         //!< Path to the input raw data file (nhdr/nrrd format).
    float isovalue;                //!< The isovalue used for isosurface extraction.
    std::string output_format;     //!< The format of the output file ("off" or "ply").
    std::string output_filename;   //!< The name of the output file.
    std::string out_csv_name;      //!< The name of the CSV file for Voronoi diagram export.

    bool out_csv;                  //!< Flag to enable exporting Voronoi diagram to CSV.
    bool sep;                      //!< Flag to enable separation of active cubes.
    bool multi_isov;               //!< Flag to enable multi-isosurface mode.
    bool noclip;                   //!< Flag to disable centroid clipping in multi-cycle cases.
    bool supersample;              //!< Flag to enable supersampling of the input data.
    bool add_bounding_cells;       //!< Flag to include bounding cells in the Voronoi diagram.
    bool convex_hull;              //!< Flag to enable convex hull computation in building voronoi cells
    bool test_vor = false;         //!< Flag for testing the Voronoi diagram construction
    bool position_delv_on_isov = false; //!< Flag to position Delaunay vertices on isosurface vertices
    bool mod_cyc = false;          //!< Guard: run modify-cycles pass (facet rematching + cycle recompute)
    bool summary_stats = false;    //!< Guard: print summary statistics at the end of the run
    bool timing_stats = false;     //!< Guard: print timing statistics at the end of the run
    bool check_bipolar_max = false; //!< Guard: check and report maximum bipolar matches per facet
    bool refine_small_angles = false; //!< Guard: enable facet-centric surface refinement
    bool refine_min_angle_enabled = false; //!< Guard: enable min-angle-driven refinement
    bool refine_max_angle_enabled = false; //!< Guard: enable max-angle-driven refinement
    bool no_check = false;             //!< Guard: skip VoronoiDiagram validity checks (faster but less safe)

    int supersample_r;             //!< Factor by which the input data is supersampled.
    double collapse_eps;           //!< Absolute collapse threshold in world units (optional).
    int sep_dist;                  //!< Separation clearance measured in subcubes.
    int sep_split;                 //!< Number of splits (K) for refined subgrid (factor = K+1).
    double refine_min_surface_angle_deg; //!< Target min surface angle (deg) for optional refinement
    double refine_max_surface_angle_deg; //!< Target max surface angle (deg) triggering refinement when large
    int refine_insert_resolution;        //!< 1=cube center, 2=2x2x2, 3=3x3x3 subcell
    double refine_min_spacing;           //!< Min distance to existing vertices (world units); <0 auto derives
    int refine_max_iterations;           //!< Max refinement iterations
    bool refine_snap_to_grid;            //!< Use discrete subcells instead of continuous iso bisection

    //! @brief Constructor to initialize default parameter values.
    VdcParam()
        : file_path(""),
          isovalue(0.0f),
          output_format("off"),
          output_filename(""),
          out_csv_name("voronoi.csv"),
          out_csv(false),
          sep(false),
          multi_isov(true),
          noclip(false),
          supersample(false),
          add_bounding_cells(false),
          convex_hull(false),
          position_delv_on_isov(false),
          supersample_r(1),
          collapse_eps(-1.0),
          sep_dist(1),
          sep_split(0),
          mod_cyc(true),
          summary_stats(false),
          timing_stats(false),
          check_bipolar_max(false),
          refine_small_angles(false),
          refine_min_angle_enabled(false),
          refine_max_angle_enabled(false),
          refine_min_surface_angle_deg(20.0),
          refine_max_surface_angle_deg(-1.0),
          refine_insert_resolution(3),
          refine_min_spacing(-1.0),
          refine_max_iterations(1),
          refine_snap_to_grid(true)
    {}

    //! @brief Print VDC parameters for debugging
    template <typename OSTREAM_TYPE>
    void Print(OSTREAM_TYPE & out) const {
        out << "VdcParam:\n";
        out << "  File path: " << file_path << "\n";
        out << "  Isovalue: " << isovalue << "\n";
        out << "  Output format: " << output_format << "\n";
        out << "  Output filename: " << output_filename << "\n";
        out << "  Out CSV name: " << out_csv_name << "\n";
        out << "  Out CSV: " << (out_csv ? "true" : "false") << "\n";
        out << "  Separation: " << (sep ? "true" : "false") << " (dist=" << sep_dist << ", split=" << sep_split << ")\n";
        out << "  Multi isov: " << (multi_isov ? "true" : "false") << "\n";
        out << "  No clip: " << (noclip ? "true" : "false") << "\n";
        out << "  Supersample: " << (supersample ? "true" : "false") << "\n";
        out << "  Add bounding cells: " << (add_bounding_cells ? "true" : "false") << "\n";
        out << "  Convex hull: " << (convex_hull ? "true" : "false") << "\n";
        out << "  Test vor: " << (test_vor ? "true" : "false") << "\n";
        out << "  Mod cyc: " << (mod_cyc ? "true" : "false") << "\n";
        out << "  Position DelV on IsoV: " << (position_delv_on_isov ? "true" : "false") << "\n";
        out << "  Summary stats: " << (summary_stats ? "true" : "false") << "\n";
        out << "  Timing stats: " << (timing_stats ? "true" : "false") << "\n";
        out << "  Check bipolar max: " << (check_bipolar_max ? "true" : "false") << "\n";
        out << "  No check: " << (no_check ? "true" : "false") << "\n";
        out << "  Supersample r: " << supersample_r << "\n";
        out << "  Collapse eps: " << collapse_eps << "\n";
        out << "  Refine small angles: " << (refine_small_angles ? "true" : "false") << "\n";
        out << "  Refine min-angle enabled: " << (refine_min_angle_enabled ? "true" : "false") << "\n";
        out << "  Refine min surface angle (deg): " << refine_min_surface_angle_deg << "\n";
        out << "  Refine max-angle enabled: " << (refine_max_angle_enabled ? "true" : "false") << "\n";
        out << "  Refine max surface angle (deg): " << refine_max_surface_angle_deg << "\n";
        out << "  Refine insert resolution: " << refine_insert_resolution << "\n";
        out << "  Refine min spacing: " << refine_min_spacing << "\n";
        out << "  Refine max iterations: " << refine_max_iterations << "\n";
        out << "  Refine snap-to-grid: " << (refine_snap_to_grid ? "true" : "false") << "\n";
    }
};

//! @brief Prints the help message to the console.
/*!
 * This function outputs usage information and available options for the program,
 * including details about input/output configurations and processing modes.
 */
void print_help();

//! @brief Parses the command-line arguments to populate program parameters.
/*!
 * This function processes command-line arguments to configure the program's behavior.
 * It sets parameters such as output format, supersampling factor, and mode selection.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line argument strings.
 * @param vp A reference to a `VdcParam` object where parsed parameters are stored.
 */
void parse_arguments(int argc, char *argv[], VdcParam &vp);

#endif // VDC_COMMANDLINE_H
