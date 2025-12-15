#include "core/vdc_debug.h"
#include <iostream>
#include <fstream>


/* Helper Functions */

//! @brief Global variable to enable or disable debug output.
/*!
 * When set to `true`, debug-related print commands are executed.
 */
bool debug = false;

//! @brief Global variable to enable or disable progress indicators.
/*!
 * When set to `true`, progress-related print commands are executed.
 */
bool indicator = true;

//! @brief Prints information about a Delaunay cell.
/*!
 * Outputs the vertices of the given Delaunay cell to the standard error stream.
 * 
 * @param c The cell whose information is to be printed.
 */
void print_cell(Delaunay::Cell c)
{
    using namespace std;

    cerr << "Cell: [";
    for (int i = 0; i < 4; i++)
    {
        // Print the coordinates of the vertex at index `i`.
        cerr << "(" << c.vertex(i)->point() << ")";
        if (i < 3)
        {
            cerr << ","; // Add a comma between vertices, except after the last one.
        }
    }
    cerr << "]" << endl;
}

//! @brief Prints information about a Delaunay facet.
/*!
 * Outputs the vertices of the given Delaunay facet to the console, along with
 * the facet's index within its associated cell.
 * 
 * @param f The facet to be printed. A facet is represented as a pair of a cell handle
 *          and an index of the facet within the cell.
 */
void print_facet(Facet f)
{
    // Extract the index of the facet within the cell.
    int iFacet = f.second;

    // Determine the indices of the three vertices forming the facet.
    int d1, d2, d3;
    d1 = (iFacet + 1) % 4;
    d2 = (iFacet + 2) % 4;
    d3 = (iFacet + 3) % 4;

    // Retrieve the cell containing the facet.
    Cell_handle c = f.first;

    // Print the vertices of the facet.
    std::cout << "Facet: " << c->vertex(d1)->point() << ", "
              << c->vertex(d2)->point() << ", "
              << c->vertex(d3)->point() << std::endl;

    // Print the facet's index within the cell.
    std::cout << "ifacet: " << iFacet << std::endl;
}

//! @brief Writes dummy points to a CSV file for debugging purposes.
/*!
 * Outputs the grid dimensions, spacing, and the coordinates of dummy points to a CSV file
 * named `dummy_points.csv`. This function can be useful for visualizing dummy points during debugging.
 * 
 * @param grid The grid containing the dimensions and spacing information.
 * @param dummy_points A vector of points to be written to the CSV file.
 */
void write_dummy_points(UnifiedGrid &grid, std::vector<Point> dummy_points)
{
    // Temporary method for writing dummy points to a CSV file.
    // This is always executed since the `if (true)` condition is hardcoded.

    // Open the output file.
    std::ofstream ofs("dummy_points.csv");

    // Write grid metadata (dimensions and spacing) as the first line.
    ofs << grid.num_cells[0] << "," // Grid size along the x-axis.
        << grid.num_cells[1] << "," // Grid size along the y-axis.
        << grid.num_cells[2] << "," // Grid size along the z-axis.
        << grid.physical_spacing[0] << "," // Grid spacing along the x-axis.
        << grid.physical_spacing[1] << "," // Grid spacing along the y-axis.
        << grid.physical_spacing[2] << "\n"; // Grid spacing along the z-axis.

    // Write the column headers for the dummy points.
    ofs << "x,y,z\n";

    // Write the coordinates of each dummy point.
    for (const auto &p : dummy_points)
    {
        ofs << p.x() << "," << p.y() << "," << p.z() << "\n";
    }

    // Close the output file.
    ofs.close();
}
