#include "core/vdc_commandline.h"
#include "core/vdc_debug.h"

#include <cstdlib>
#include <filesystem>

//! Prints the help message for the program.
void print_help()
{
    std::cout << "Usage: vdc-del [OPTIONS] <isovalue> <(nhdr/nrrd) raw data file path>\n\n";
    std::cout << "OPTIONS:\n";
    std::cout << "  -o {output_filename}        : Specify output filename (default: derived from input filename).\n";
    std::cout << "  -off                        : Generate output in .off format (default).\n";
    std::cout << "  -ply                        : Generate output in .ply format.\n";
    std::cout << "  -sep_dist {D}               : Separation distance in refined subcubes (default: 0 = off).\n";
    std::cout << "  -sep_split {K}              : Number of splits per axis (refined factor K+1, default: 0).\n";
    std::cout << "  -supersample {factor}       : Supersample the input data by the given factor.\n";
    std::cout << "  -multi_isov                 : Use multi iso-vertices mode (default).\n";
    std::cout << "  -single_isov                : Use single iso-vertices mode.\n";
    std::cout << "  -position_delv_on_isov      : Position Delaunay vertices on isosurface vertices instead of cube iso-crossings.\n";
    std::cout << "  -refine_small_angles        : Insert Delaunay sites at circumsphere centers to improve small angles near the isosurface.\n";
    std::cout << "  -min_angle {deg}            : Trigger refinement if any isosurface-facet triangle angle (per-site iso-sample) is below this threshold.\n";
    std::cout << "  -max_angle {deg}            : Trigger refinement if any isosurface-facet triangle angle (per-site iso-sample) is above this threshold.\n";
    std::cout << "  -terse                      : Print only the number of vertices/triangles and the output file.\n";
    std::cout << "  -timing_stats               : Print timing statistics after the run.\n";
    std::cout << "  -debug                      : Enable debug logging ([DEBUG]/[ISO]/[ISO-MATCH]/[CYC-MOD]).\n";
    std::cout << "  -no_modcyc                  : Disable modify-cycles pass (enabled by default).\n";
    std::cout << "  -help                       : Print this help message.\n";
}

//! Parses command-line arguments and configures program settings.
void parse_arguments(int argc, char *argv[], VdcParam &vp)
{
    // Print help and exit if there are insufficient arguments.
    if (argc < 3)
    {
        print_help();
        exit(EXIT_FAILURE);
    }

    // Parse optional arguments (those starting with '-').
    int i = 1;
    bool sep_requested = false;

    while (i < argc && argv[i][0] == '-')
    {
        std::string arg = argv[i];

        if (arg == "-o" && i + 1 < argc)
        {
            vp.output_filename = argv[++i]; // Set custom output filename.
        }
        else if (arg == "-off")
        {
            vp.output_format = "off"; // Set output format to .off.
        }
        else if (arg == "-ply")
        {
            vp.output_format = "ply"; // Set output format to .ply.
        }
        else if (arg == "-sep_dist" && i + 1 < argc)
        {
            vp.sep_dist = std::atoi(argv[++i]); // Separation distance in refined subcubes.
            sep_requested = true;
        }
        else if (arg == "-sep_split" && i + 1 < argc)
        {
            vp.sep_split = std::atoi(argv[++i]); // Number of splits (K).
            sep_requested = true;
        }
        else if (arg == "-supersample" && i + 1 < argc)
        {
            vp.supersample = true;                     // Enable supersampling.
            vp.supersample_r = std::atoi(argv[++i]);   // Set supersampling factor.
        }
        else if (arg == "-multi_isov")
        {
            vp.multi_isov = true; // Enable multi-isovertex mode.
        }
        else if (arg == "-single_isov")
        {
            vp.multi_isov = false; // Enable single-isovertex mode.
        }
        else if (arg == "-position_delv_on_isov")
        {
            vp.position_delv_on_isov = true; // Place Delaunay vertices at isosurface vertex locations.
        }
        else if (arg == "-refine_small_angles")
        {
            vp.refine_small_angles = true;
        }
        else if (arg == "-min_angle" && i + 1 < argc)
        {
            vp.refine_min_surface_angle_deg = std::atof(argv[++i]);
            vp.refine_min_angle_enabled = true;
            vp.refine_small_angles = true;
        }
        else if (arg == "-max_angle" && i + 1 < argc)
        {
            vp.refine_max_surface_angle_deg = std::atof(argv[++i]);
            vp.refine_max_angle_enabled = true;
            vp.refine_small_angles = true;
        }
        else if (arg == "-terse")
        {
            vp.terse = true; // Reduce output: only print vertex/triangle counts and output file.
        }
        else if (arg == "-help")
        {
            print_help();
            exit(EXIT_SUCCESS);
        }
        else if (arg == "-timing_stats")
        {
            vp.timing_stats = true; // Enable timing statistics report
        }
        else if (arg == "-debug")
        {
            debug = true; // Enable global debug logging
        }
        else if (arg == "-no_modcyc")
        {
            vp.mod_cyc = false; // Disable modify-cycles pass
        }
        else
        {
            // Handle unknown options.
            std::cerr << "Unknown option: " << arg << "\n";
            print_help();
            exit(EXIT_FAILURE);
        }
        ++i;
    }

    // Parse required arguments: isovalue and file path.
    if (i + 2 > argc)
    {
        std::cerr << "Error: Missing required arguments.\n";
        print_help();
        exit(EXIT_FAILURE);
    }

    vp.isovalue = std::atof(argv[i++]); // Parse isovalue as a floating-point number.
    vp.file_path = argv[i++];           // Parse the raw data file path.

    // Normalize separation parameters
    if (vp.sep_dist < 0) vp.sep_dist = 0;
    if (vp.sep_split < 0) vp.sep_split = 0;
    vp.sep = sep_requested || (vp.sep_dist > 1 || vp.sep_split > 0);

    // Generate default output filename if not specified.
    if (vp.output_filename.empty())
    {
        // Extract base name from the file path (strip directory and extension).
        const std::filesystem::path input_path(vp.file_path);
        std::string base_name = input_path.stem().string();
        if (base_name.empty())
        {
            base_name = input_path.filename().string();
        }
        if (base_name.empty())
        {
            base_name = "vdc_output";
        }
        vp.output_filename = base_name;

        // Append processing details to the filename.
        if (!vp.multi_isov)
        {
            vp.output_filename += "_single-isov";
        }

        if (vp.supersample)
        {
            vp.output_filename += "_sup" + std::to_string(vp.supersample_r);
        }

        if (vp.sep)
        {
            vp.output_filename += "_sep-isov-D" + std::to_string(vp.sep_dist)
                                  + "-S" + std::to_string(vp.sep_split);
        }

        if (!vp.mod_cyc)
        {
            vp.output_filename += "_no-mod-cyc";
        }

        if (vp.refine_small_angles)
        {
            vp.output_filename += "_refine";
            if (vp.refine_min_angle_enabled)
            {
                vp.output_filename += "-min" + std::to_string(static_cast<int>(vp.refine_min_surface_angle_deg));
            }
            if (vp.refine_max_angle_enabled)
            {
                vp.output_filename += "-max" + std::to_string(static_cast<int>(vp.refine_max_surface_angle_deg));
            }
        }

        // Add file format extension.
        vp.output_filename += "." + vp.output_format;
    }
}
