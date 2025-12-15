#include "core/vdc_utilities.h"

//! Converts a CGAL::Object to a string for debugging purposes.
std::string objectToString(const Object &obj)
{
    std::ostringstream stream;
    Segment3 seg;
    Ray3 ray;
    Line3 line;

    if (CGAL::assign(seg, obj))
    {
        stream << "Segment: " << seg;
    }
    else if (CGAL::assign(ray, obj))
    {
        stream << "Ray: " << ray;
    }
    else if (CGAL::assign(line, obj))
    {
        stream << "Line: " << line;
    }
    else {
        stream << "unkown edge type";
    }
    return stream.str();
}

//! Checks if a CGAL::Object is degenerate.
bool is_degenerate(const Object &obj)
{
    Segment3 seg;
    if (CGAL::assign(seg, obj))
    {
        return seg.source() == seg.target(); // Degenerate if start == end.
    }
    return false; // Rays and lines are not degenerate in the same way.
}

//! Checks if a Delaunay cell is degenerate.
bool is_degenerate(Delaunay::Cell_handle cell)
{
    Point p0 = cell->vertex(0)->point();
    Point p1 = cell->vertex(1)->point();
    Point p2 = cell->vertex(2)->point();
    Point p3 = cell->vertex(3)->point();

    // Compute tetrahedron volume using CGAL.
    K::FT vol = CGAL::volume(p0, p1, p2, p3);

    // Consider degenerate if volume is near zero.
    return CGAL::abs(vol) < 1e-4;
}

//! Checks if two scalar values are bipolar.
bool is_bipolar(float val1, float val2, float isovalue)
{
    return ((val1 < isovalue) && (val2 >= isovalue)) || ((val1 >= isovalue) && (val2 < isovalue));
}

//! Computes the centroid of a set of points (with optional supersampling).
Point compute_centroid(const std::vector<Point> &points, bool supersample, int ratio)
{
    float sumX = 0, sumY = 0, sumZ = 0;

    for (const auto &pt : points)
    {
        sumX += pt.x();
        sumY += pt.y();
        sumZ += pt.z();
    }

    int n = points.size();
    float x = sumX / n;
    float y = sumY / n;
    float z = sumZ / n;

    return Point(x, y, z);
}

//! Computes the centroid of a set of points using CGAL.
Point compute_centroid(const std::vector<Point> &points)
{
    return CGAL::centroid(points.begin(), points.end());
}

//! Gets the corners of a cube given its center and side length.
std::array<Point, 8> get_cube_corners(const Point &center, float side_length)
{
    float half_side = side_length / 2.0;

    return {{
        Point(center.x() - half_side, center.y() - half_side, center.z() - half_side), // 0
        Point(center.x() + half_side, center.y() - half_side, center.z() - half_side), // 1
        Point(center.x() + half_side, center.y() + half_side, center.z() - half_side), // 2
        Point(center.x() - half_side, center.y() + half_side, center.z() - half_side), // 3
        Point(center.x() - half_side, center.y() - half_side, center.z() + half_side), // 4
        Point(center.x() + half_side, center.y() - half_side, center.z() + half_side), // 5
        Point(center.x() + half_side, center.y() + half_side, center.z() + half_side), // 6
        Point(center.x() - half_side, center.y() + half_side, center.z() + half_side)  // 7
    }};
}

//! Determines the orientation of a facet.
int get_orientation(const int iFacet, const Point v1, const Point v2, const float f1, const float f2)
{
    bool flag_v1_positive = (f1 >= f2);

    if ((iFacet % 2) == 0)
    {
        return flag_v1_positive ? 1 : -1;
    }
    else
    {
        return flag_v1_positive ? -1 : 1;
    }
}

//! Checks if a scalar value is positive or zero.
bool is_positive(double value, float isovalue)
{
    return value >= isovalue; // Compare to global isovalue.
}

bool read_points_from_file(const std::string &filename, std::vector<Point> &points)
{
    std::ifstream file(filename);
    if (!file)
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        double x, y, z;
        char comma1, comma2;

        if (!(ss >> x >> comma1 >> y >> comma2 >> z) || comma1 != ',' || comma2 != ',')
        {
            std::cerr << "Warning: Skipping invalid line: " << line << std::endl;
            continue;
        }

        points.emplace_back(x, y, z);
    }

    file.close();
    return true;
}

void write_triangulation(Delaunay dt, std::vector<Point> &points, std::string &input_filename)
{
    // Save points and edges to file
    std::ofstream file("triangulation.txt");
    if (!file)
    {
        std::cerr << "Error opening output file.\n";
        exit(EXIT_FAILURE);
    }

    // Write points
    file << "POINTS\n";
    for (const auto &p : points)
    {
        file << p.x() << " " << p.y() << " " << p.z() << "\n";
    }

    // Write edges
    file << "EDGES\n";
    for (auto e = dt.finite_edges_begin(); e != dt.finite_edges_end(); ++e)
    {
        auto v1 = e->first->vertex(e->second);
        auto v2 = e->first->vertex(e->third);
        file << v1->point().x() << " " << v1->point().y() << " " << v1->point().z() << " ";
        file << v2->point().x() << " " << v2->point().y() << " " << v2->point().z() << "\n";
    }

    file.close();
    std::cout << "Triangulation saved to triangulation.txt\n";
}
