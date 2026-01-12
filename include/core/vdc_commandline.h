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

    bool sep;                      //!< Flag to enable separation of active cubes.
    bool multi_isov;               //!< Flag to enable multi-isosurface mode.
    bool supersample;              //!< Flag to enable supersampling of the input data.
    bool position_delv_on_isov = false; //!< Flag to position Delaunay vertices on isosurface vertices
    bool position_multi_isov_on_delv = false; //!< Debug: place all multi-cycle isovertices at the Delaunay vertex
    bool reposition_multi_isovA = false; //!< Reposition multi isovertices using only hyperplane separation and reflection
    bool reposition_multi_isovA_trace = false; //!< Like -reposition_multi_isovA, but also dumps A trace (local/ + final/) as OFF/TXT files
    bool terse = false;            //!< Guard: print only vertices/triangles and output file
    bool timing_stats = false;     //!< Guard: print timing statistics at the end of the run
    bool refine_small_angles = false; //!< Guard: enable facet-centric surface refinement
    bool refine_min_angle_enabled = false; //!< Guard: enable min-angle-driven refinement
    bool refine_max_angle_enabled = false; //!< Guard: enable max-angle-driven refinement
    bool mod_cyc = true;                   //!< Guard: run modify-cycles pass to fix non-manifolds (default: enabled)

    int supersample_r;             //!< Factor by which the input data is supersampled.
    int sep_dist;                  //!< Separation clearance measured in subcubes.
    int sep_split;                 //!< Number of splits (K) for refined subgrid (factor = K+1).
    double refine_min_surface_angle_deg; //!< Target min surface angle (deg) for optional refinement
    double refine_max_surface_angle_deg; //!< Target max surface angle (deg) triggering refinement when large
    int refine_insert_resolution;        //!< 1=cube center, 2=2x2x2, 3=3x3x3 subcell

    //! @brief Constructor to initialize default parameter values.
    VdcParam()
        : file_path(""),
          isovalue(0.0f),
          output_format("off"),
          output_filename(""),
          sep(false),
          multi_isov(true),
          supersample(false),
          position_delv_on_isov(false),
          position_multi_isov_on_delv(false),
          reposition_multi_isovA(false),
          reposition_multi_isovA_trace(false),
          terse(false),
          timing_stats(false),
          refine_small_angles(false),
          refine_min_angle_enabled(false),
          refine_max_angle_enabled(false),
          mod_cyc(true),
          supersample_r(1),
          sep_dist(1),
          sep_split(0),
          refine_min_surface_angle_deg(20.0),
          refine_max_surface_angle_deg(-1.0),
          refine_insert_resolution(3)
    {}

    //! @brief Print VDC parameters for debugging
    template <typename OSTREAM_TYPE>
    void Print(OSTREAM_TYPE & out) const {
        out << "VdcParam:\n";
        out << "  File path: " << file_path << "\n";
        out << "  Isovalue: " << isovalue << "\n";
        out << "  Output format: " << output_format << "\n";
        out << "  Output filename: " << output_filename << "\n";
        out << "  Separation: " << (sep ? "true" : "false") << " (dist=" << sep_dist << ", split=" << sep_split << ")\n";
        out << "  Multi isov: " << (multi_isov ? "true" : "false") << "\n";
        out << "  Supersample: " << (supersample ? "true" : "false") << "\n";
        out << "  Position DelV on IsoV: " << (position_delv_on_isov ? "true" : "false") << "\n";
        out << "  Position multi IsoV on DelV: " << (position_multi_isov_on_delv ? "true" : "false") << "\n";
        out << "  Reposition multi IsoV A: " << (reposition_multi_isovA ? "true" : "false") << "\n";
        out << "  Reposition multi IsoV A trace: " << (reposition_multi_isovA_trace ? "true" : "false") << "\n";
        out << "  Terse: " << (terse ? "true" : "false") << "\n";
        out << "  Timing stats: " << (timing_stats ? "true" : "false") << "\n";
        out << "  Supersample r: " << supersample_r << "\n";
        out << "  Refine small angles: " << (refine_small_angles ? "true" : "false") << "\n";
        out << "  Refine min surface angle (deg): " << refine_min_surface_angle_deg << "\n";
        out << "  Refine max surface angle (deg): " << refine_max_surface_angle_deg << "\n";
        out << "  Refine insert resolution: " << refine_insert_resolution << "\n";
        out << "  Mod cyc: " << (mod_cyc ? "true" : "false") << "\n";
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
