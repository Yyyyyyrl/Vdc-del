#include "vdc_debug.h"

#include "core/vdc_commandline.h"
#include "processing/vdc_del_cycles.h"
#include "processing/vdc_del_isosurface.h"
#include "processing/vdc_func.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

bool debug = false;
bool indicator = true;

void print_cell(Delaunay::Cell c)
{
    using namespace std;

    cerr << "Cell: [";
    for (int i = 0; i < 4; i++)
    {
        cerr << "(" << c.vertex(i)->point() << ")";
        if (i < 3)
        {
            cerr << ",";
        }
    }
    cerr << "]" << endl;
}

void print_facet(Facet f)
{
    int iFacet = f.second;
    int d1, d2, d3;
    d1 = (iFacet + 1) % 4;
    d2 = (iFacet + 2) % 4;
    d3 = (iFacet + 3) % 4;

    Cell_handle c = f.first;

    std::cout << "Facet: " << c->vertex(d1)->point() << ", "
              << c->vertex(d2)->point() << ", "
              << c->vertex(d3)->point() << std::endl;
    std::cout << "ifacet: " << iFacet << std::endl;
}

void write_dummy_points(UnifiedGrid &grid, std::vector<Point> dummy_points)
{
    std::ofstream ofs("dummy_points.csv");
    ofs << grid.num_cells[0] << ","
        << grid.num_cells[1] << ","
        << grid.num_cells[2] << ","
        << grid.physical_spacing[0] << ","
        << grid.physical_spacing[1] << ","
        << grid.physical_spacing[2] << "\n";
    ofs << "x,y,z\n";
    for (const auto &p : dummy_points)
    {
        ofs << p.x() << "," << p.y() << "," << p.z() << "\n";
    }
    ofs.close();
}

namespace {

std::string format_float_for_path(float v) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(6) << static_cast<double>(v);
    std::string s = ss.str();
    while (!s.empty() && s.back() == '0') {
        s.pop_back();
    }
    if (!s.empty() && s.back() == '.') {
        s.pop_back();
    }
    if (s.empty()) {
        s = "0";
    }
    return s;
}

std::string trace_config_tag(const VdcParam& param) {
    std::vector<std::string> parts;

    if (param.supersample) {
        parts.push_back("sup" + std::to_string(param.supersample_r));
    }

    if (param.sep) {
        std::ostringstream ss;
        ss << "sep" << param.sep_split << "_dist" << param.sep_dist;
        parts.push_back(ss.str());
    }

    if (!param.multi_isov) {
        parts.push_back("single_isov");
    }
    if (param.position_delv_on_isov) {
        parts.push_back("delv_on_isov");
    }
    if (!param.mod_cyc) {
        parts.push_back("no_modcyc");
    }

    if (parts.empty()) {
        return "default";
    }

    std::ostringstream out;
    for (size_t i = 0; i < parts.size(); ++i) {
        if (i > 0) {
            out << "_";
        }
        out << parts[i];
    }
    return out.str();
}

void json_write_string(std::ostream& out, const std::string& s) {
    out << '"';
    for (const char ch : s) {
        switch (ch) {
            case '\\': out << "\\\\"; break;
            case '"': out << "\\\""; break;
            case '\n': out << "\\n"; break;
            case '\r': out << "\\r"; break;
            case '\t': out << "\\t"; break;
            default: out << ch; break;
        }
    }
    out << '"';
}

void json_write_point(std::ostream& out, const Point& p) {
    out << "[" << p.x() << "," << p.y() << "," << p.z() << "]";
}

void json_write_int3(std::ostream& out, const int v[3]) {
    out << "[" << v[0] << "," << v[1] << "," << v[2] << "]";
}

void ensure_parent_dir_exists(const std::filesystem::path& path) {
    const std::filesystem::path parent = path.parent_path();
    if (parent.empty()) {
        return;
    }
    std::error_code ec;
    std::filesystem::create_directories(parent, ec);
}

bool cell_vertices_include_dummy(const Delaunay& dt, const Cell_handle& cell) {
    if (cell == Cell_handle() || dt.is_infinite(cell)) {
        return true;
    }
    for (int i = 0; i < 4; ++i) {
        Vertex_handle vh = cell->vertex(i);
        if (dt.is_infinite(vh) || vh->info().is_dummy) {
            return true;
        }
    }
    return false;
}

template <typename CellPtrLike>
bool cell_ptr_like_vertices_include_dummy(const Delaunay& dt, const CellPtrLike& cell) {
    for (int i = 0; i < 4; ++i) {
        Vertex_handle vh = cell->vertex(i);
        if (dt.is_infinite(vh) || vh->info().is_dummy) {
            return true;
        }
    }
    return false;
}

std::vector<Cell_handle> build_cell_index_vector_local(const Delaunay& dt) {
    std::vector<Cell_handle> cells(static_cast<size_t>(dt.number_of_finite_cells()));
    for (auto cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit) {
        const int idx = cit->info().index;
        if (idx >= 0 && static_cast<size_t>(idx) < cells.size()) {
            cells[static_cast<size_t>(idx)] = cit;
        }
    }
    return cells;
}

Cell_handle lookup_cell_local(const std::vector<Cell_handle>& cells, int idx) {
    if (idx < 0 || static_cast<size_t>(idx) >= cells.size()) {
        return Cell_handle();
    }
    return cells[static_cast<size_t>(idx)];
}

Vertex_handle find_vertex_by_index(const Delaunay& dt, int idx) {
    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
        if (vit->info().index == idx) {
            return vit;
        }
    }
    return Vertex_handle();
}

} // namespace

namespace vdc_debug {

void dump_site_selection_json(
    const std::filesystem::path& path,
    const UnifiedGrid& grid,
    float isovalue,
    int sep_dist,
    int sep_split,
    const std::vector<Cube>& pre_sep,
    const std::vector<Cube>& post_sep
) {
    ensure_parent_dir_exists(path);
    std::ofstream out(path);
    if (!out) {
        std::cerr << "Warning: failed to write JSON: " << path << "\n";
        return;
    }

    out << std::setprecision(17);
    out << "{\n";
    out << "  \"isovalue\": " << static_cast<double>(isovalue) << ",\n";
    out << "  \"grid_num_cells\": [" << grid.num_cells[0] << "," << grid.num_cells[1] << "," << grid.num_cells[2] << "],\n";
    out << "  \"grid_spacing\": [" << grid.spacing[0] << "," << grid.spacing[1] << "," << grid.spacing[2] << "],\n";
    out << "  \"sep_dist\": " << sep_dist << ",\n";
    out << "  \"sep_split\": " << sep_split << ",\n";

    auto write_cubes = [&](const char* key, const std::vector<Cube>& cubes) {
        out << "  \"" << key << "\": [\n";
        for (size_t i = 0; i < cubes.size(); ++i) {
            const Cube& c = cubes[i];
            out << "    {\n";
            out << "      \"indices\": ";
            json_write_int3(out, c.indices);
            out << ",\n";
            out << "      \"rep_vertex\": ";
            json_write_point(out, c.repVertex);
            out << ",\n";
            out << "      \"cube_center\": ";
            json_write_point(out, c.cubeCenter);
            out << ",\n";
            out << "      \"iso_crossing\": ";
            json_write_point(out, c.accurateIsoCrossing);
            out << "\n";
            out << "    }" << (i + 1 < cubes.size() ? "," : "") << "\n";
        }
        out << "  ]";
    };

    write_cubes("active_cubes_pre_sep", pre_sep);
    out << ",\n";
    write_cubes("active_cubes_post_sep", post_sep);
    out << "\n";
    out << "}\n";
}

void dump_isosurface_facet_example_json(
    const std::filesystem::path& path,
    const VdcParam& param,
    const Delaunay& dt
) {
    ensure_parent_dir_exists(path);
    std::ofstream out(path);
    if (!out) {
        std::cerr << "Warning: failed to write JSON: " << path << "\n";
        return;
    }

    // Pick the first isosurface facet whose two incident tetrahedra do not use dummy vertices.
    auto pos_cell_it = dt.finite_cells_end();
    Cell_handle neg_cell;
    int pos_facet_idx = -1;
    int neg_facet_idx = -1;

    for (auto cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit) {
        for (int i = 0; i < 4; ++i) {
            if (!cit->info().facet_is_isosurface[i]) {
                continue;
            }
            Cell_handle neighbor = cit->neighbor(i);
            if (dt.is_infinite(neighbor)) {
                continue;
            }
            if (cell_ptr_like_vertices_include_dummy(dt, cit) || cell_vertices_include_dummy(dt, neighbor)) {
                continue;
            }

            // Find the matching facet index in neighbor such that neighbor->neighbor(j) == cit.
            int j_found = -1;
            for (int j = 0; j < 4; ++j) {
                Cell_handle nbr_nbr = neighbor->neighbor(j);
                if (dt.is_infinite(nbr_nbr)) {
                    continue;
                }
                if (nbr_nbr->info().index == cit->info().index) {
                    j_found = j;
                    break;
                }
            }
            if (j_found < 0) {
                continue;
            }

            pos_cell_it = cit;
            neg_cell = neighbor;
            pos_facet_idx = i;
            neg_facet_idx = j_found;
            break;
        }
        if (pos_cell_it != dt.finite_cells_end()) {
            break;
        }
    }

    out << std::setprecision(17);
    out << "{\n";
    out << "  \"dataset\": ";
    json_write_string(out, param.file_path);
    out << ",\n";
    out << "  \"isovalue\": " << static_cast<double>(param.isovalue) << ",\n";

    if (pos_cell_it == dt.finite_cells_end() || neg_cell == Cell_handle() || pos_facet_idx < 0 || neg_facet_idx < 0) {
        out << "  \"error\": \"No suitable interior isosurface facet found (dummy-free).\"\n";
        out << "}\n";
        return;
    }

    auto write_tet_vertices_pos = [&]() {
        out << "    [\n";
        for (int vi = 0; vi < 4; ++vi) {
            Vertex_handle vh = pos_cell_it->vertex(vi);
            out << "      {\"index\": " << vh->info().index << ", \"pos\": ";
            json_write_point(out, vh->point());
            out << "}" << (vi < 3 ? "," : "") << "\n";
        }
        out << "    ]\n";
    };
    auto write_tet_vertices_neg = [&]() {
        out << "    [\n";
        for (int vi = 0; vi < 4; ++vi) {
            Vertex_handle vh = neg_cell->vertex(vi);
            out << "      {\"index\": " << vh->info().index << ", \"pos\": ";
            json_write_point(out, vh->point());
            out << "}" << (vi < 3 ? "," : "") << "\n";
        }
        out << "    ]\n";
    };

    // Shared facet vertices (3 points) from the positive cell side.
    Point facet_pts[3];
    int facet_pt_idx = 0;
    for (int vi = 0; vi < 4; ++vi) {
        if (vi == pos_facet_idx) {
            continue;
        }
        facet_pts[facet_pt_idx++] = pos_cell_it->vertex(vi)->point();
    }

    const Point cA = pos_cell_it->info().circumcenter;
    const Point cB = neg_cell->info().circumcenter;
    const double sA = static_cast<double>(pos_cell_it->info().circumcenter_scalar);
    const double sB = static_cast<double>(neg_cell->info().circumcenter_scalar);
    const double sigma = static_cast<double>(param.isovalue);

    // Interpolate along Voronoi edge to isovalue.
    Point q = cA;
    double wA = 0.5;
    double wB = 0.5;
    if (sA != sB) {
        wA = (sigma - sB) / (sA - sB);
        wB = (sA - sigma) / (sA - sB);
        q = Point(
            wA * cA.x() + wB * cB.x(),
            wA * cA.y() + wB * cB.y(),
            wA * cA.z() + wB * cB.z());
    }

    out << "  \"facet\": {\n";
    out << "    \"positive_cell_index\": " << pos_cell_it->info().index << ",\n";
    out << "    \"positive_facet_index\": " << pos_facet_idx << ",\n";
    out << "    \"negative_cell_index\": " << neg_cell->info().index << ",\n";
    out << "    \"negative_facet_index\": " << neg_facet_idx << ",\n";
    out << "    \"positive_tet_vertices\":\n";
    write_tet_vertices_pos();
    out << "    ,\"negative_tet_vertices\":\n";
    write_tet_vertices_neg();
    out << "    ,\"shared_facet_vertices\": [";
    for (int i = 0; i < 3; ++i) {
        json_write_point(out, facet_pts[i]);
        if (i < 2) out << ",";
    }
    out << "]\n";
    out << "  },\n";

    out << "  \"circumcenters\": {\n";
    out << "    \"positive\": {\"pos\": ";
    json_write_point(out, cA);
    out << ", \"scalar\": " << sA << "},\n";
    out << "    \"negative\": {\"pos\": ";
    json_write_point(out, cB);
    out << ", \"scalar\": " << sB << "}\n";
    out << "  },\n";

    out << "  \"voronoi_edge\": {\n";
    out << "    \"a\": ";
    json_write_point(out, cA);
    out << ",\n";
    out << "    \"b\": ";
    json_write_point(out, cB);
    out << ",\n";
    out << "    \"q\": ";
    json_write_point(out, q);
    out << ",\n";
    out << "    \"weights\": {\"wA\": " << wA << ", \"wB\": " << wB << "}\n";
    out << "  }\n";
    out << "}\n";
}

void dump_multicycle_vertex_example_json(
    const std::filesystem::path& path,
    const VdcParam& param,
    const Delaunay& dt
) {
    ensure_parent_dir_exists(path);
    std::ofstream out(path);
    if (!out) {
        std::cerr << "Warning: failed to write JSON: " << path << "\n";
        return;
    }

    const std::vector<Cell_handle> cell_by_index = build_cell_index_vector_local(dt);

    Vertex_handle v;
    if (param.dump_multicycle_vertex >= 0) {
        v = find_vertex_by_index(dt, param.dump_multicycle_vertex);
    }
    if (v == Vertex_handle()) {
        // Pick a small, clean example: prefer exactly 2 cycles and fewer facets.
        int best_total_facets = std::numeric_limits<int>::max();
        for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
            if (!vit->info().active) continue;
            if (vit->info().is_dummy) continue;
            const auto& cycles = vit->info().facet_cycles;
            if (cycles.size() != 2) continue;

            const int total = static_cast<int>(cycles[0].size() + cycles[1].size());
            if (total < best_total_facets && total >= 6) {
                best_total_facets = total;
                v = vit;
            }
        }

        if (v == Vertex_handle()) {
            for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
                if (!vit->info().active) continue;
                if (vit->info().is_dummy) continue;
                if (vit->info().facet_cycles.size() >= 2) {
                    v = vit;
                    break;
                }
            }
        }
    }

    out << std::setprecision(17);
    out << "{\n";
    out << "  \"dataset\": ";
    json_write_string(out, param.file_path);
    out << ",\n";
    out << "  \"isovalue\": " << static_cast<double>(param.isovalue) << ",\n";

    if (v == Vertex_handle()) {
        out << "  \"error\": \"No multi-cycle vertex found.\"\n";
        out << "}\n";
        return;
    }

    const int vidx = v->info().index;
    const auto& cycles = v->info().facet_cycles;

    out << "  \"vertex_index\": " << vidx << ",\n";
    out << "  \"vertex_pos\": ";
    json_write_point(out, v->point());
    out << ",\n";
    out << "  \"num_cycles\": " << cycles.size() << ",\n";
    out << "  \"has_isov_sample\": " << (v->info().has_isov_sample ? "true" : "false") << ",\n";
    out << "  \"isov_sample\": ";
    json_write_point(out, v->info().isov);
    out << ",\n";

    // Delaunay neighborhood quality metrics (helps diagnose slivers / far circumcenters).
    struct WorstCell {
        int cell_index = -1;
        double min_dihedral_deg = std::numeric_limits<double>::infinity();
        double circumradius = 0.0;
        struct Vtx { int index = -1; Point pos; };
        Vtx vtx[4];
    };

    auto vec_dot = [](const Vector3& a, const Vector3& b) -> double {
        return a.x() * b.x() + a.y() * b.y() + a.z() * b.z();
    };
    auto vec_cross = [](const Vector3& a, const Vector3& b) -> Vector3 {
        return Vector3(
            a.y() * b.z() - a.z() * b.y(),
            a.z() * b.x() - a.x() * b.z(),
            a.x() * b.y() - a.y() * b.x());
    };
    auto vec_norm = [](const Vector3& v0) -> double {
        return std::sqrt(v0.squared_length());
    };
    auto clamp_unit = [](double x) -> double {
        if (x < -1.0) return -1.0;
        if (x > 1.0) return 1.0;
        return x;
    };
    auto angle_between = [&](const Vector3& a, const Vector3& b) -> double {
        const double na = vec_norm(a);
        const double nb = vec_norm(b);
        if (na < 1e-18 || nb < 1e-18) {
            return 0.0;
        }
        const double c = clamp_unit(vec_dot(a, b) / (na * nb));
        return std::acos(c);
    };
    auto outward_normal = [&](const Point& pi, const Point& pj, const Point& pk, const Point& pl) -> Vector3 {
        Vector3 n = vec_cross(Vector3(pi, pj), Vector3(pi, pk));
        // Ensure the normal points away from the opposite vertex (pl).
        if (vec_dot(n, Vector3(pi, pl)) > 0.0) {
            n = -n;
        }
        return n;
    };
    auto internal_dihedral_rad = [&](const Point& pi, const Point& pj, const Point& pk, const Point& pl) -> double {
        // Dihedral angle along edge (pi,pj) between faces (pi,pj,pk) and (pi,pj,pl).
        const Vector3 n1 = outward_normal(pi, pj, pk, pl);
        const Vector3 n2 = outward_normal(pi, pj, pl, pk);
        const double ang = angle_between(n1, n2); // angle between outward face normals
        constexpr double kPi = 3.141592653589793238462643383279502884;
        return kPi - ang;
    };
    auto min_cell_dihedral_deg = [&](const Point& p0, const Point& p1, const Point& p2, const Point& p3) -> double {
        constexpr double kRadToDeg = 57.295779513082320876798154814105170332;
        double m = std::numeric_limits<double>::infinity();
        m = std::min(m, internal_dihedral_rad(p0, p1, p2, p3));
        m = std::min(m, internal_dihedral_rad(p0, p2, p1, p3));
        m = std::min(m, internal_dihedral_rad(p0, p3, p1, p2));
        m = std::min(m, internal_dihedral_rad(p1, p2, p0, p3));
        m = std::min(m, internal_dihedral_rad(p1, p3, p0, p2));
        m = std::min(m, internal_dihedral_rad(p2, p3, p0, p1));
        return m * kRadToDeg;
    };

    // Compute min incident edge length (excluding dummy/infinite neighbors).
    double min_edge_length = std::numeric_limits<double>::infinity();
    {
        std::vector<Cell_handle> incident_cells;
        incident_cells.reserve(256);
        dt.incident_cells(v, std::back_inserter(incident_cells));

        const Point p = v->point();
        for (Cell_handle ch : incident_cells) {
            if (ch == Cell_handle() || dt.is_infinite(ch)) {
                continue;
            }
            for (int i = 0; i < 4; ++i) {
                Vertex_handle vh = ch->vertex(i);
                if (vh == v) {
                    continue;
                }
                if (dt.is_infinite(vh) || vh->info().is_dummy) {
                    continue;
                }
                const Point q = vh->point();
                const double dx = p.x() - q.x();
                const double dy = p.y() - q.y();
                const double dz = p.z() - q.z();
                const double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                min_edge_length = std::min(min_edge_length, dist);
            }
        }
    }

    // Compute per-cell quality stats (finite, non-dummy).
    int incident_cells_total = 0;
    int incident_cells_valid = 0;
    double min_dihedral_deg = std::numeric_limits<double>::infinity();
    double max_dihedral_deg = 0.0;
    double min_circumradius = std::numeric_limits<double>::infinity();
    double max_circumradius = 0.0;
    WorstCell worst;

    {
        std::vector<Cell_handle> incident_cells;
        incident_cells.reserve(256);
        dt.incident_cells(v, std::back_inserter(incident_cells));
        incident_cells_total = static_cast<int>(incident_cells.size());

        for (Cell_handle ch : incident_cells) {
            if (ch == Cell_handle() || dt.is_infinite(ch)) {
                continue;
            }
            if (cell_vertices_include_dummy(dt, ch)) {
                continue;
            }

            ++incident_cells_valid;

            Point p[4];
            int vid[4] = {-1, -1, -1, -1};
            for (int i = 0; i < 4; ++i) {
                Vertex_handle vh = ch->vertex(i);
                p[i] = vh->point();
                vid[i] = vh->info().index;
            }

            const double d01 = min_cell_dihedral_deg(p[0], p[1], p[2], p[3]); // includes all edges
            // min_cell_dihedral_deg returns the min over the 6 internal dihedrals;
            // track global min and (conservatively) the max of these minima.
            min_dihedral_deg = std::min(min_dihedral_deg, d01);
            max_dihedral_deg = std::max(max_dihedral_deg, d01);

            const Point cc = ch->info().circumcenter;
            const double dx = cc.x() - p[0].x();
            const double dy = cc.y() - p[0].y();
            const double dz = cc.z() - p[0].z();
            const double cr = std::sqrt(dx * dx + dy * dy + dz * dz);
            min_circumradius = std::min(min_circumradius, cr);
            max_circumradius = std::max(max_circumradius, cr);

            if (d01 < worst.min_dihedral_deg) {
                worst.cell_index = ch->info().index;
                worst.min_dihedral_deg = d01;
                worst.circumradius = cr;
                for (int i = 0; i < 4; ++i) {
                    worst.vtx[i].index = vid[i];
                    worst.vtx[i].pos = p[i];
                }
            }
        }
    }

    if (!std::isfinite(min_edge_length)) {
        // Fallback: should be rare (isolated vertex); use a conservative scale.
        min_edge_length = 1.0;
    }
    const double sphere_radius_base = 0.1 * min_edge_length;
    const double sphere_radius_scale = 1.0;
    const double sphere_radius = sphere_radius_base;
    out << "  \"delaunay_neighborhood\": {\n";
    out << "    \"incident_cells_total\": " << incident_cells_total << ",\n";
    out << "    \"incident_cells_valid\": " << incident_cells_valid << ",\n";
    out << "    \"min_incident_edge_length\": " << min_edge_length << ",\n";
    out << "    \"sphere_radius_from_min_edge\": " << sphere_radius_base << ",\n";
    out << "    \"sphere_radius_sliver_scale\": " << sphere_radius_scale << ",\n";
    out << "    \"sphere_radius_sliver_scaled\": " << sphere_radius << ",\n";
    out << "    \"min_cell_min_dihedral_deg\": " << (std::isfinite(min_dihedral_deg) ? min_dihedral_deg : -1.0) << ",\n";
    out << "    \"max_cell_min_dihedral_deg\": " << max_dihedral_deg << ",\n";
    out << "    \"min_cell_circumradius\": " << (std::isfinite(min_circumradius) ? min_circumradius : -1.0) << ",\n";
    out << "    \"max_cell_circumradius\": " << max_circumradius << ",\n";
    out << "    \"max_circumradius_over_min_edge\": " << (min_edge_length > 1e-18 ? (max_circumradius / min_edge_length) : -1.0) << ",\n";
    out << "    \"worst_cell_by_min_dihedral\": ";
    if (worst.cell_index < 0 || !std::isfinite(worst.min_dihedral_deg)) {
        out << "null\n";
    } else {
        out << "{\n";
        out << "      \"cell_index\": " << worst.cell_index << ",\n";
        out << "      \"min_dihedral_deg\": " << worst.min_dihedral_deg << ",\n";
        out << "      \"circumradius\": " << worst.circumradius << ",\n";
        out << "      \"vertices\": [\n";
        for (int i = 0; i < 4; ++i) {
            out << "        {\"index\": " << worst.vtx[i].index << ", \"pos\": ";
            json_write_point(out, worst.vtx[i].pos);
            out << "}" << (i + 1 < 4 ? "," : "") << "\n";
        }
        out << "      ]\n";
        out << "    }\n";
    }
    out << "  },\n";

    // Voronoi cell vertices are (a subset of) circumcenters of Delaunay cells incident to v.
    {
        std::vector<Cell_handle> incident_cells;
        incident_cells.reserve(256);
        dt.incident_cells(v, std::back_inserter(incident_cells));

        std::vector<Point> circumcenters;
        circumcenters.reserve(incident_cells.size());
        for (Cell_handle ch : incident_cells) {
            if (ch == Cell_handle() || dt.is_infinite(ch)) {
                continue;
            }
            circumcenters.push_back(ch->info().circumcenter);
        }

        // Deduplicate exact duplicates (rare but possible in degenerate cases).
        std::sort(circumcenters.begin(), circumcenters.end(), [](const Point& a, const Point& b) {
            if (a.x() != b.x()) return a.x() < b.x();
            if (a.y() != b.y()) return a.y() < b.y();
            return a.z() < b.z();
        });
        circumcenters.erase(std::unique(circumcenters.begin(), circumcenters.end(), [](const Point& a, const Point& b) {
            return a.x() == b.x() && a.y() == b.y() && a.z() == b.z();
        }), circumcenters.end());

        out << "  \"voronoi_cell_vertices\": [\n";
        for (size_t i = 0; i < circumcenters.size(); ++i) {
            out << "    ";
            json_write_point(out, circumcenters[i]);
            out << (i + 1 < circumcenters.size() ? "," : "") << "\n";
        }
        out << "  ],\n";
    }

    out << "  \"cycles\": [\n";
    for (size_t c = 0; c < cycles.size(); ++c) {
        const auto& facets = cycles[c];
        const Point isov_p = (c < v->info().cycle_isovertices.size())
                                 ? v->info().cycle_isovertices[c]
                                 : v->point();

        // Collect Voronoi-edge/isovalue intersections q_f for this cycle.
        std::vector<Point> q_points;
        q_points.reserve(facets.size());
        struct EdgeSample { Point a; double sa; Point b; double sb; Point q; };
        std::vector<EdgeSample> edges;
        edges.reserve(facets.size());

        for (const auto& facet_key : facets) {
            const int cell_idx = facet_key.first;
            const int facet_idx = facet_key.second;
            Cell_handle cell = lookup_cell_local(cell_by_index, cell_idx);
            if (cell == Cell_handle() || dt.is_infinite(cell)) {
                continue;
            }
            if (facet_idx < 0 || facet_idx >= 4) {
                continue;
            }

            Cell_handle nbr = cell->neighbor(facet_idx);
            if (nbr == Cell_handle() || dt.is_infinite(nbr)) {
                continue;
            }

            const Point a = cell->info().circumcenter;
            const Point b = nbr->info().circumcenter;
            const double sa = static_cast<double>(cell->info().circumcenter_scalar);
            const double sb = static_cast<double>(nbr->info().circumcenter_scalar);
            const double sigma = static_cast<double>(param.isovalue);
            if (sa == sb) {
                continue;
            }

            const double wA = (sigma - sb) / (sa - sb);
            const double wB = (sa - sigma) / (sa - sb);
            const Point q(
                wA * a.x() + wB * b.x(),
                wA * a.y() + wB * b.y(),
                wA * a.z() + wB * b.z());

            q_points.push_back(q);
            edges.push_back(EdgeSample{a, sa, b, sb, q});
        }

        Point centroid = v->point();
        if (!q_points.empty()) {
            double cx = 0.0, cy = 0.0, cz = 0.0;
            for (const Point& q : q_points) {
                cx += q.x();
                cy += q.y();
                cz += q.z();
            }
            const double inv = 1.0 / static_cast<double>(q_points.size());
            centroid = Point(cx * inv, cy * inv, cz * inv);
        }

        // Collect cycle triangle fan geometry (same connectivity as output mesh).
        struct Tri { Point p0; Point p1; Point p2; };
        std::vector<Tri> tris;
        tris.reserve(facets.size());

        for (const auto& facet_key : facets) {
            const int cell_idx = facet_key.first;
            const int facet_idx = facet_key.second;
            Cell_handle cell = lookup_cell_local(cell_by_index, cell_idx);
            if (cell == Cell_handle() || dt.is_infinite(cell)) {
                continue;
            }
            if (facet_idx < 0 || facet_idx >= 4) {
                continue;
            }

            // Find the two other vertices on this facet besides v.
            Vertex_handle others[2] = {Vertex_handle(), Vertex_handle()};
            int other_count = 0;
            for (int t = 0; t < 3; ++t) {
                const int cv = (facet_idx + 1 + t) % 4;
                Vertex_handle vh = cell->vertex(cv);
                if (vh == v) {
                    continue;
                }
                if (other_count < 2) {
                    others[other_count] = vh;
                }
                ++other_count;
            }
            if (other_count != 2 || others[0] == Vertex_handle() || others[1] == Vertex_handle()) {
                continue;
            }

            const int c1 = find_cycle_containing_facet(others[0], cell_idx, facet_idx);
            const int c2 = find_cycle_containing_facet(others[1], cell_idx, facet_idx);
            if (c1 < 0 || c2 < 0) {
                continue;
            }
            if (c1 >= static_cast<int>(others[0]->info().cycle_isovertices.size()) ||
                c2 >= static_cast<int>(others[1]->info().cycle_isovertices.size())) {
                continue;
            }

            const Point p1 = others[0]->info().cycle_isovertices[static_cast<size_t>(c1)];
            const Point p2 = others[1]->info().cycle_isovertices[static_cast<size_t>(c2)];
            tris.push_back(Tri{isov_p, p1, p2});
        }

        out << "    {\n";
        out << "      \"cycle_index\": " << c << ",\n";
        out << "      \"facet_count\": " << facets.size() << ",\n";
        out << "      \"iso_vertex\": ";
        json_write_point(out, isov_p);
        out << ",\n";
        out << "      \"centroid\": ";
        json_write_point(out, centroid);
        out << ",\n";

        // Facet keys (cell_index, facet_index)
        out << "      \"facets\": [";
        for (size_t fi = 0; fi < facets.size(); ++fi) {
            out << "[" << facets[fi].first << "," << facets[fi].second << "]";
            if (fi + 1 < facets.size()) out << ",";
        }
        out << "],\n";

        out << "      \"voronoi_edges\": [\n";
        for (size_t ei = 0; ei < edges.size(); ++ei) {
            const auto& e = edges[ei];
            out << "        {\"a\": ";
            json_write_point(out, e.a);
            out << ", \"sa\": " << e.sa << ", \"b\": ";
            json_write_point(out, e.b);
            out << ", \"sb\": " << e.sb << ", \"q\": ";
            json_write_point(out, e.q);
            out << "}" << (ei + 1 < edges.size() ? "," : "") << "\n";
        }
        out << "      ],\n";

        out << "      \"triangles\": [\n";
        for (size_t ti = 0; ti < tris.size(); ++ti) {
            const auto& t = tris[ti];
            out << "        [";
            json_write_point(out, t.p0);
            out << ",";
            json_write_point(out, t.p1);
            out << ",";
            json_write_point(out, t.p2);
            out << "]" << (ti + 1 < tris.size() ? "," : "") << "\n";
        }
        out << "      ]\n";

        out << "    }" << (c + 1 < cycles.size() ? "," : "") << "\n";
    }
    out << "  ]\n";
    out << "}\n";
}

void write_delv_off(
    const Delaunay& dt,
    const std::string& filename,
    bool has_bbox,
    const double* bbox_min,
    const double* bbox_max
) {
    auto in_bbox = [&](const Point& p) -> bool {
        if (!has_bbox || !bbox_min || !bbox_max) return true;
        return p.x() >= bbox_min[0] && p.x() <= bbox_max[0] &&
               p.y() >= bbox_min[1] && p.y() <= bbox_max[1] &&
               p.z() >= bbox_min[2] && p.z() <= bbox_max[2];
    };

    std::unordered_map<Vertex_handle, int> vertex_index;
    std::vector<Point> vertices;
    int idx = 0;
    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
        if (in_bbox(vit->point())) {
            vertex_index[vit] = idx++;
            vertices.push_back(vit->point());
        }
    }

    std::vector<std::array<int, 3>> facets;
    for (auto fit = dt.finite_facets_begin(); fit != dt.finite_facets_end(); ++fit) {
        Cell_handle c = fit->first;
        int i = fit->second;

        std::array<int, 3> tri;
        int ti = 0;
        bool all_in_bbox = true;
        for (int j = 0; j < 4; ++j) {
            if (j == i) continue;
            Vertex_handle vh = c->vertex(j);
            if (dt.is_infinite(vh)) goto skip_facet;
            auto it = vertex_index.find(vh);
            if (it == vertex_index.end()) {
                all_in_bbox = false;
                break;
            }
            tri[ti++] = it->second;
        }
        if (!all_in_bbox) goto skip_facet;
        if (i % 2 == 0) {
            std::swap(tri[0], tri[1]);
        }
        facets.push_back(tri);
        skip_facet:;
    }

    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: Failed to open " << filename << " for writing.\n";
        return;
    }

    out << "OFF\n";
    out << vertices.size() << " " << facets.size() << " 0\n";
    out << std::setprecision(17);
    for (const auto& p : vertices) {
        out << p.x() << " " << p.y() << " " << p.z() << "\n";
    }
    for (const auto& f : facets) {
        out << "3 " << f[0] << " " << f[1] << " " << f[2] << "\n";
    }

    std::cout << "Wrote Delaunay triangulation to: " << filename << " ("
              << vertices.size() << " vertices, " << facets.size() << " facets)";
    if (has_bbox) {
        std::cout << " [cropped to bbox]";
    }
    std::cout << "\n";
}

void dump_duplicate_isovertex_positions(
    const std::filesystem::path& trace_root,
    const DelaunayIsosurface& iso_surface
) {
    if (trace_root.empty()) {
        return;
    }

    struct Key {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;

        bool operator==(const Key& other) const noexcept {
            return x == other.x && y == other.y && z == other.z;
        }
    };

    struct KeyHash {
        size_t operator()(const Key& k) const noexcept {
            const size_t h1 = std::hash<double>{}(k.x);
            const size_t h2 = std::hash<double>{}(k.y);
            const size_t h3 = std::hash<double>{}(k.z);
            size_t h = h1;
            h ^= (h2 + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
            h ^= (h3 + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
            return h;
        }
    };

    std::unordered_map<Key, std::vector<int>, KeyHash> groups;
    groups.reserve(iso_surface.isovertices.size());

    constexpr double QUANT = 1e8;

    const auto& s = iso_surface.vertex_scale;
    for (size_t i = 0; i < iso_surface.isovertices.size(); ++i) {
        const Point& p = iso_surface.isovertices[i];
        const Key k{
            std::round(p.x() * s[0] * QUANT) / QUANT,
            std::round(p.y() * s[1] * QUANT) / QUANT,
            std::round(p.z() * s[2] * QUANT) / QUANT,
        };
        groups[k].push_back(static_cast<int>(i));
    }

    std::vector<std::pair<Key, std::vector<int>>> dups;
    dups.reserve(32);
    for (auto& kv : groups) {
        if (kv.second.size() > 1) {
            dups.emplace_back(kv.first, std::move(kv.second));
        }
    }

    std::sort(dups.begin(), dups.end(), [](const auto& a, const auto& b) {
        return a.second.size() > b.second.size();
    });

    std::error_code ec;
    std::filesystem::create_directories(trace_root, ec);

    std::ofstream out(trace_root / "global_duplicate_isovertex_positions.txt");
    if (!out) {
        return;
    }

    out << std::setprecision(17);
    out << "count_groups " << dups.size() << "\n";

    if (dups.empty()) {
        return;
    }

    for (const auto& g : dups) {
        out << "pos " << g.first.x << " " << g.first.y << " " << g.first.z
            << " count " << g.second.size() << "\n";
        for (const int isov_idx : g.second) {
            const int delv = (isov_idx >= 0 && static_cast<size_t>(isov_idx) < iso_surface.isovertex_delaunay_vertex.size())
                                 ? iso_surface.isovertex_delaunay_vertex[static_cast<size_t>(isov_idx)]
                                 : -1;
            const int cyc = (isov_idx >= 0 && static_cast<size_t>(isov_idx) < iso_surface.isovertex_cycle_index.size())
                                 ? iso_surface.isovertex_cycle_index[static_cast<size_t>(isov_idx)]
                                 : -1;
            out << "  isov " << isov_idx << " delv " << delv << " cycle " << cyc << "\n";
        }
    }
}

void dump_isovertex_map_txt(
    const std::filesystem::path& path,
    const DelaunayIsosurface& iso_surface
) {
    if (path.empty()) {
        return;
    }

    ensure_parent_dir_exists(path);
    std::ofstream out(path);
    if (!out) {
        std::cerr << "Warning: failed to write isovertex map: " << path << "\n";
        return;
    }

    out << std::setprecision(17);
    out << "# isovertex_index delv_index cycle_index x y z\n";
    const auto& s = iso_surface.vertex_scale;
    for (size_t i = 0; i < iso_surface.isovertices.size(); ++i) {
        const Point& p = iso_surface.isovertices[i];
        const int delv = (i < iso_surface.isovertex_delaunay_vertex.size())
                             ? iso_surface.isovertex_delaunay_vertex[i]
                             : -1;
        const int cyc = (i < iso_surface.isovertex_cycle_index.size())
                             ? iso_surface.isovertex_cycle_index[i]
                             : -1;

        out << i << " " << delv << " " << cyc << " "
            << (p.x() * s[0]) << " " << (p.y() * s[1]) << " " << (p.z() * s[2]) << "\n";
    }
}

std::string prepare_multi_isov_trace_dir(const char* argv0, const VdcParam& param) {
    std::filesystem::path repo_root;
    {
        std::error_code ec;
        std::filesystem::path exe_path = std::filesystem::weakly_canonical(argv0, ec);
        if (ec || exe_path.empty()) {
            ec.clear();
            exe_path = std::filesystem::absolute(argv0, ec);
        }

        if (!exe_path.empty()) {
            repo_root = exe_path.parent_path().parent_path();
            if (!std::filesystem::exists(repo_root / "CMakeLists.txt")) {
                repo_root.clear();
            }
        }
        if (repo_root.empty()) {
            repo_root = std::filesystem::current_path();
        }
    }

    const std::filesystem::path base = repo_root / "simple_multi_failures_trace";
    const std::filesystem::path input_path(param.file_path);
    const std::string dataset = input_path.stem().string().empty() ? "dataset" : input_path.stem().string();
    const std::string iso = format_float_for_path(param.isovalue);
    const std::string cfg = trace_config_tag(param);
    const std::filesystem::path dir = base / dataset / ("iso" + iso + "_" + cfg);

    {
        std::error_code rm_ec;
        std::filesystem::remove_all(dir / "local", rm_ec);
        rm_ec.clear();
        std::filesystem::remove_all(dir / "final", rm_ec);
    }

    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    if (ec) {
        std::cerr << "[DEL-SELFI-TRACE] Warning: failed to create directory: " << dir
                  << " (" << ec.message() << ")\n";
        return "";
    }

    std::cerr << "[DEL-SELFI-TRACE] Dumping multi-isov trace under: " << dir.string()
              << " (subdirs: local/, final/)\n";
    return dir.string();
}

} // namespace vdc_debug

// -----------------------------------------------------------------------------
// Self-intersection tracing + dump utilities (moved from vdc_del_cycles.cpp)
// -----------------------------------------------------------------------------

// Declared in src/processing/vdc_del_cycles.cpp (algorithm code).
ResolutionResult try_resolve_multicycle_by_cycle_separation_tests(
    Vertex_handle v,
    std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    std::vector<Point>* geometric_attempt_positions_out);

// Declared in src/processing/vdc_del_cycles.cpp (algorithm code).
bool check_self_intersection(
    Vertex_handle v,
    const std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index);

namespace {

using Triangle_3 = K::Triangle_3;

const std::unordered_set<int>& trace_selfi_vertices_internal() {
    static std::unordered_set<int> traced;
    static bool initialized = false;
    if (initialized) {
        return traced;
    }
    initialized = true;

    auto insert_int = [&](const std::string& token) {
        try {
            size_t pos = 0;
            const int v = std::stoi(token, &pos);
            if (pos == token.size()) {
                traced.insert(v);
            }
        } catch (...) {
            // Ignore malformed tokens.
        }
    };

    if (const char* env = std::getenv("VDC_TRACE_SELFI_VERTICES")) {
        std::stringstream ss(env);
        std::string tok;
        while (std::getline(ss, tok, ',')) {
            if (!tok.empty()) {
                insert_int(tok);
            }
        }
        return traced;
    }

    if (const char* env = std::getenv("VDC_TRACE_SELFI_VERTEX")) {
        insert_int(env);
    }

    return traced;
}

std::filesystem::path selfi_dump_path_for_vertex(
    const int vertex_index,
    const char* stage,
    const char* suffix
) {
    const std::filesystem::path base(vdc_debug::trace_selfi_dump_dir());
    std::ostringstream name;
    name << "selfi_v" << vertex_index << "_" << stage << suffix;
    return base / name.str();
}

struct TrianglePoints {
    Point p0;
    Point p1;
    Point p2;
};

std::vector<TrianglePoints> collect_cycle_triangles_points(
    const Delaunay& dt,
    Vertex_handle v,
    int cycle_idx,
    const Point& cycle_isovertex,
    const std::vector<Cell_handle>& cell_by_index
) {
    std::vector<TrianglePoints> triangles;
    const auto& cycles = v->info().facet_cycles;

    if (cycle_idx < 0 || cycle_idx >= static_cast<int>(cycles.size())) {
        return triangles;
    }

    const auto& cycle_facets = cycles[static_cast<size_t>(cycle_idx)];
    triangles.reserve(cycle_facets.size());

    for (const auto& [cell_idx, facet_idx] : cycle_facets) {
        Cell_handle ch = lookup_cell_local(cell_by_index, cell_idx);
        if (ch == Cell_handle()) {
            continue;
        }

        Vertex_handle other_verts[2] = {Vertex_handle(), Vertex_handle()};
        int other_slots[2] = {-1, -1};
        int other_count = 0;
        for (int t = 0; t < 3; ++t) {
            const int cell_vertex_idx = (facet_idx + 1 + t) % 4;
            Vertex_handle vh = ch->vertex(cell_vertex_idx);
            if (vh == v) {
                continue;
            }
            if (other_count < 2) {
                other_verts[other_count] = vh;
                other_slots[other_count] = t;
            }
            ++other_count;
        }

        if (other_count != 2) {
            continue;
        }

        Vertex_handle v1 = other_verts[0];
        Vertex_handle v2 = other_verts[1];

        int c1 = -1;
        if (facet_idx >= 0 && facet_idx < 4 &&
            other_slots[0] >= 0 && other_slots[0] < 3) {
            c1 = ch->info().facet_info[facet_idx].dualCellEdgeIndex[other_slots[0]];
        }
        if (c1 < 0) {
            c1 = find_cycle_containing_facet(v1, cell_idx, facet_idx);
        }

        int c2 = -1;
        if (facet_idx >= 0 && facet_idx < 4 &&
            other_slots[1] >= 0 && other_slots[1] < 3) {
            c2 = ch->info().facet_info[facet_idx].dualCellEdgeIndex[other_slots[1]];
        }
        if (c2 < 0) {
            c2 = find_cycle_containing_facet(v2, cell_idx, facet_idx);
        }

        if (c1 < 0 || c2 < 0) {
            continue;
        }
        if (c1 >= static_cast<int>(v1->info().cycle_isovertices.size()) ||
            c2 >= static_cast<int>(v2->info().cycle_isovertices.size())) {
            continue;
        }

        triangles.push_back(TrianglePoints{
            cycle_isovertex,
            v1->info().cycle_isovertices[static_cast<size_t>(c1)],
            v2->info().cycle_isovertices[static_cast<size_t>(c2)],
        });
    }

    return triangles;
}

bool write_vtk_selfi_local_triangles(
    const std::filesystem::path& path,
    const std::vector<std::vector<TrianglePoints>>& cycle_triangles
) {
    size_t num_polys = 0;
    for (const auto& tris : cycle_triangles) {
        num_polys += tris.size();
    }
    const size_t num_points = num_polys * 3;

    std::ofstream out(path);
    if (!out) {
        return false;
    }

    out << "# vtk DataFile Version 3.0\n";
    out << "vdc-del self-intersection local fan\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";
    out << "POINTS " << num_points << " float\n";

    for (const auto& tris : cycle_triangles) {
        for (const auto& tri : tris) {
            const Point pts[3] = {tri.p0, tri.p1, tri.p2};
            for (int k = 0; k < 3; ++k) {
                const Point p = pts[k];
                out << static_cast<float>(p.x()) << " "
                    << static_cast<float>(p.y()) << " "
                    << static_cast<float>(p.z()) << "\n";
            }
        }
    }

    out << "POLYGONS " << num_polys << " " << (num_polys * 4) << "\n";
    size_t poly_idx = 0;
    for (const auto& tris : cycle_triangles) {
        for (size_t t = 0; t < tris.size(); ++t) {
            const size_t base = poly_idx * 3;
            out << "3 " << base << " " << (base + 1) << " " << (base + 2) << "\n";
            ++poly_idx;
        }
    }

    out << "CELL_DATA " << num_polys << "\n";
    out << "SCALARS cycle int 1\n";
    out << "LOOKUP_TABLE default\n";
    for (size_t c = 0; c < cycle_triangles.size(); ++c) {
        for (size_t t = 0; t < cycle_triangles[c].size(); ++t) {
            out << static_cast<int>(c) << "\n";
        }
    }

    return true;
}

bool write_vtk_selfi_points(
    const std::filesystem::path& path,
    Vertex_handle v,
    const std::vector<Point>& cycle_isovertices
) {
    const int num_cycles = static_cast<int>(cycle_isovertices.size());
    const int num_points = num_cycles + 1; // includes center

    std::ofstream out(path);
    if (!out) {
        return false;
    }

    out << "# vtk DataFile Version 3.0\n";
    out << "vdc-del self-intersection points\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";
    out << "POINTS " << num_points << " float\n";

    const Point center = v->point();
    out << static_cast<float>(center.x()) << " "
        << static_cast<float>(center.y()) << " "
        << static_cast<float>(center.z()) << "\n";

    for (int c = 0; c < num_cycles; ++c) {
        const Point p = cycle_isovertices[static_cast<size_t>(c)];
        out << static_cast<float>(p.x()) << " "
            << static_cast<float>(p.y()) << " "
            << static_cast<float>(p.z()) << "\n";
    }

    out << "VERTICES " << num_points << " " << (num_points * 2) << "\n";
    for (int i = 0; i < num_points; ++i) {
        out << "1 " << i << "\n";
    }

    out << "POINT_DATA " << num_points << "\n";
    out << "SCALARS cycle int 1\n";
    out << "LOOKUP_TABLE default\n";
    out << "-1\n"; // center
    for (int c = 0; c < num_cycles; ++c) {
        out << c << "\n";
    }

    return true;
}

// Intersection helpers for simple_multi_failures summaries (debug-only).

enum class TriIntersectionKind {
    NONE = 0,
    DISJOINT_TRIANGLE_TRIANGLE,
    SHARED_VERTEX_T1_OPPOSITE_SEGMENT,
    SHARED_VERTEX_T2_OPPOSITE_SEGMENT,
    SHARED_VERTEX_INTERIOR_SEGMENT_T1_0,
    SHARED_VERTEX_INTERIOR_SEGMENT_T1_1,
    SHARED_VERTEX_INTERIOR_SEGMENT_T2_0,
    SHARED_VERTEX_INTERIOR_SEGMENT_T2_1,
};

struct TriIntersectionTrace {
    int shared_count = 0;
    Point shared_point = Point(0, 0, 0);
    TriIntersectionKind kind = TriIntersectionKind::NONE;
};

bool triangles_have_nontrivial_intersection(
    const Triangle_3& t1,
    const Triangle_3& t2,
    TriIntersectionTrace* trace_out
) {
    TriIntersectionTrace trace;
    auto segment_interior_intersects_triangle = [](
        const Point& shared,
        const Point& other,
        const Triangle_3& tri
    ) -> bool {
        constexpr double t = 1e-6;
        const double dx = other.x() - shared.x();
        const double dy = other.y() - shared.y();
        const double dz = other.z() - shared.z();
        const double len2 = dx * dx + dy * dy + dz * dz;
        if (len2 < 1e-24) {
            return false;
        }

        const Point nudged(
            shared.x() + t * dx,
            shared.y() + t * dy,
            shared.z() + t * dz);
        if (nudged == other) {
            return false;
        }

        return CGAL::do_intersect(Segment3(nudged, other), tri);
    };

    const std::array<Point, 3> t1v = {t1.vertex(0), t1.vertex(1), t1.vertex(2)};
    const std::array<Point, 3> t2v = {t2.vertex(0), t2.vertex(1), t2.vertex(2)};
    int shared_count = 0;
    Point shared_point = t1v[0];
    for (const Point& a : t1v) {
        for (const Point& b : t2v) {
            if (a == b) {
                if (shared_count == 0) {
                    shared_point = a;
                }
                ++shared_count;
                break;
            }
        }
        if (shared_count >= 2) {
            break;
        }
    }
    trace.shared_count = shared_count;
    trace.shared_point = shared_point;

    if (shared_count == 0) {
        if (CGAL::do_intersect(t1, t2)) {
            trace.kind = TriIntersectionKind::DISJOINT_TRIANGLE_TRIANGLE;
            if (trace_out) {
                *trace_out = trace;
            }
            return true;
        }
        if (trace_out) {
            *trace_out = trace;
        }
        return false;
    }

    if (shared_count >= 2) {
        if (trace_out) {
            *trace_out = trace;
        }
        return false;
    }

    const Point& s = shared_point;

    std::array<Point, 2> t1_other;
    std::array<Point, 2> t2_other;
    int t1_count = 0;
    int t2_count = 0;
    for (const Point& p : t1v) {
        if (p != s && t1_count < 2) {
            t1_other[static_cast<size_t>(t1_count++)] = p;
        }
    }
    for (const Point& p : t2v) {
        if (p != s && t2_count < 2) {
            t2_other[static_cast<size_t>(t2_count++)] = p;
        }
    }
    if (t1_count != 2 || t2_count != 2) {
        if (trace_out) {
            *trace_out = trace;
        }
        return false;
    }

    if (segment_interior_intersects_triangle(s, t1_other[0], t2)) {
        trace.kind = TriIntersectionKind::SHARED_VERTEX_INTERIOR_SEGMENT_T1_0;
        if (trace_out) {
            *trace_out = trace;
        }
        return true;
    }
    if (segment_interior_intersects_triangle(s, t1_other[1], t2)) {
        trace.kind = TriIntersectionKind::SHARED_VERTEX_INTERIOR_SEGMENT_T1_1;
        if (trace_out) {
            *trace_out = trace;
        }
        return true;
    }
    if (segment_interior_intersects_triangle(s, t2_other[0], t1)) {
        trace.kind = TriIntersectionKind::SHARED_VERTEX_INTERIOR_SEGMENT_T2_0;
        if (trace_out) {
            *trace_out = trace;
        }
        return true;
    }
    if (segment_interior_intersects_triangle(s, t2_other[1], t1)) {
        trace.kind = TriIntersectionKind::SHARED_VERTEX_INTERIOR_SEGMENT_T2_1;
        if (trace_out) {
            *trace_out = trace;
        }
        return true;
    }

    if (trace_out) {
        *trace_out = trace;
    }
    return false;
}

bool triangles_have_nontrivial_intersection(const Triangle_3& t1, const Triangle_3& t2) {
    return triangles_have_nontrivial_intersection(t1, t2, nullptr);
}

struct IsovertexKey {
    int delv = -1;
    int cycle = -1;

    bool operator==(const IsovertexKey& o) const noexcept {
        return delv == o.delv && cycle == o.cycle;
    }
};

struct IsovertexKeyHash {
    size_t operator()(const IsovertexKey& k) const noexcept {
        const uint64_t a = static_cast<uint64_t>(static_cast<uint32_t>(k.delv));
        const uint64_t b = static_cast<uint64_t>(static_cast<uint32_t>(k.cycle));
        const uint64_t x = (a << 32) ^ b;
        return static_cast<size_t>(x ^ (x >> 33));
    }
};

struct KeyedTriangle {
    std::array<IsovertexKey, 3> key;
    std::array<Point, 3> p;
};

std::vector<KeyedTriangle> collect_cycle_triangles_keyed(
    const Delaunay& dt,
    Vertex_handle v,
    int cycle_idx,
    const Point& cycle_isovertex,
    const std::vector<Cell_handle>& cell_by_index
) {
    std::vector<KeyedTriangle> triangles;
    const auto& cycles = v->info().facet_cycles;

    if (cycle_idx < 0 || cycle_idx >= static_cast<int>(cycles.size())) {
        return triangles;
    }

    const int v_idx = v->info().index;
    const auto& cycle_facets = cycles[static_cast<size_t>(cycle_idx)];
    triangles.reserve(cycle_facets.size());

    for (const auto& [cell_idx, facet_idx] : cycle_facets) {
        Cell_handle ch = lookup_cell_local(cell_by_index, cell_idx);
        if (ch == Cell_handle()) {
            continue;
        }

        Vertex_handle other_verts[2] = {Vertex_handle(), Vertex_handle()};
        int other_slots[2] = {-1, -1};
        int other_count = 0;
        for (int t = 0; t < 3; ++t) {
            const int cell_vertex_idx = (facet_idx + 1 + t) % 4;
            Vertex_handle vh = ch->vertex(cell_vertex_idx);
            if (vh == v) {
                continue;
            }
            if (other_count < 2) {
                other_verts[other_count] = vh;
                other_slots[other_count] = t;
            }
            ++other_count;
        }

        if (other_count != 2) {
            continue;
        }

        Vertex_handle v1 = other_verts[0];
        Vertex_handle v2 = other_verts[1];

        int c1 = -1;
        if (facet_idx >= 0 && facet_idx < 4 &&
            other_slots[0] >= 0 && other_slots[0] < 3) {
            c1 = ch->info().facet_info[facet_idx].dualCellEdgeIndex[other_slots[0]];
        }
        if (c1 < 0) {
            c1 = find_cycle_containing_facet(v1, cell_idx, facet_idx);
        }

        int c2 = -1;
        if (facet_idx >= 0 && facet_idx < 4 &&
            other_slots[1] >= 0 && other_slots[1] < 3) {
            c2 = ch->info().facet_info[facet_idx].dualCellEdgeIndex[other_slots[1]];
        }
        if (c2 < 0) {
            c2 = find_cycle_containing_facet(v2, cell_idx, facet_idx);
        }

        if (c1 < 0 || c2 < 0) {
            continue;
        }
        if (c1 >= static_cast<int>(v1->info().cycle_isovertices.size()) ||
            c2 >= static_cast<int>(v2->info().cycle_isovertices.size())) {
            continue;
        }

        const int v1_idx = v1->info().index;
        const int v2_idx = v2->info().index;
        const Point p2 = v1->info().cycle_isovertices[static_cast<size_t>(c1)];
        const Point p3 = v2->info().cycle_isovertices[static_cast<size_t>(c2)];

        KeyedTriangle tri;
        tri.key = {
            IsovertexKey{v_idx, cycle_idx},
            IsovertexKey{v1_idx, c1},
            IsovertexKey{v2_idx, c2},
        };
        tri.p = {cycle_isovertex, p2, p3};
        triangles.push_back(std::move(tri));
    }

    return triangles;
}

std::vector<std::vector<KeyedTriangle>> collect_selfi_cycle_triangles_keyed(
    const Delaunay& dt,
    Vertex_handle v,
    const std::vector<Point>& cycle_isovertices,
    const std::vector<Cell_handle>& cell_by_index
) {
    const auto& cycles = v->info().facet_cycles;
    std::vector<std::vector<KeyedTriangle>> cycle_triangles(cycles.size());
    for (size_t c = 0; c < cycles.size(); ++c) {
        if (c >= cycle_isovertices.size()) {
            break;
        }
        cycle_triangles[c] = collect_cycle_triangles_keyed(
            dt, v, static_cast<int>(c), cycle_isovertices[c], cell_by_index);
    }
    return cycle_triangles;
}

bool write_off_keyed_triangles(
    const std::filesystem::path& path,
    const std::vector<KeyedTriangle>& triangles
) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }

    std::vector<Point> vertices;
    vertices.reserve(triangles.size() * 3);

    std::unordered_map<IsovertexKey, size_t, IsovertexKeyHash> vtx_index;
    vtx_index.reserve(triangles.size() * 3);

    std::vector<std::array<size_t, 3>> faces;
    faces.reserve(triangles.size());

    for (const auto& tri : triangles) {
        std::array<size_t, 3> f{};
        for (int k = 0; k < 3; ++k) {
            const IsovertexKey key = tri.key[static_cast<size_t>(k)];
            const auto it = vtx_index.find(key);
            if (it != vtx_index.end()) {
                f[static_cast<size_t>(k)] = it->second;
                continue;
            }

            const size_t idx = vertices.size();
            vertices.push_back(tri.p[static_cast<size_t>(k)]);
            vtx_index.emplace(key, idx);
            f[static_cast<size_t>(k)] = idx;
        }
        faces.push_back(f);
    }

    out << "OFF\n";
    out << vertices.size() << " " << faces.size() << " 0\n";
    out << std::setprecision(17);

    for (const Point p : vertices) {
        out << p.x() << " " << p.y() << " " << p.z() << "\n";
    }

    for (const auto& f : faces) {
        out << "3 " << f[0] << " " << f[1] << " " << f[2] << "\n";
    }

    return true;
}

bool cycle_triangles_have_within_cycle_intersection(
    const std::vector<std::vector<Triangle_3>>& cycle_triangles
) {
    for (const auto& tris : cycle_triangles) {
        for (size_t i = 0; i < tris.size(); ++i) {
            for (size_t j = i + 1; j < tris.size(); ++j) {
                if (triangles_have_nontrivial_intersection(tris[i], tris[j])) {
                    return true;
                }
            }
        }
    }
    return false;
}

std::optional<std::pair<int, int>> first_between_cycle_intersection_pair(
    const std::vector<std::vector<Triangle_3>>& cycle_triangles
) {
    for (size_t i = 0; i < cycle_triangles.size(); ++i) {
        for (size_t j = i + 1; j < cycle_triangles.size(); ++j) {
            for (const auto& t1 : cycle_triangles[i]) {
                for (const auto& t2 : cycle_triangles[j]) {
                    if (triangles_have_nontrivial_intersection(t1, t2)) {
                        return std::pair<int, int>(static_cast<int>(i), static_cast<int>(j));
                    }
                }
            }
        }
    }
    return std::nullopt;
}

struct SelfIntersectionWitnessTrace {
    int cycle0 = -1;
    int cycle1 = -1;
    size_t tri0 = 0;
    size_t tri1 = 0;
    TriIntersectionTrace trace;
};

std::optional<SelfIntersectionWitnessTrace> first_within_cycle_intersection_trace(
    const std::vector<std::vector<Triangle_3>>& cycle_triangles
) {
    for (size_t c = 0; c < cycle_triangles.size(); ++c) {
        const auto& tris = cycle_triangles[c];
        for (size_t i = 0; i < tris.size(); ++i) {
            for (size_t j = i + 1; j < tris.size(); ++j) {
                TriIntersectionTrace trace;
                if (triangles_have_nontrivial_intersection(tris[i], tris[j], &trace)) {
                    return SelfIntersectionWitnessTrace{
                        static_cast<int>(c), static_cast<int>(c), i, j, trace};
                }
            }
        }
    }
    return std::nullopt;
}

std::optional<SelfIntersectionWitnessTrace> first_between_cycle_intersection_trace(
    const std::vector<std::vector<Triangle_3>>& cycle_triangles
) {
    for (size_t c0 = 0; c0 < cycle_triangles.size(); ++c0) {
        for (size_t c1 = c0 + 1; c1 < cycle_triangles.size(); ++c1) {
            for (size_t i = 0; i < cycle_triangles[c0].size(); ++i) {
                for (size_t j = 0; j < cycle_triangles[c1].size(); ++j) {
                    TriIntersectionTrace trace;
                    if (triangles_have_nontrivial_intersection(
                            cycle_triangles[c0][i], cycle_triangles[c1][j], &trace)) {
                        return SelfIntersectionWitnessTrace{
                            static_cast<int>(c0), static_cast<int>(c1), i, j, trace};
                    }
                }
            }
        }
    }
    return std::nullopt;
}

bool write_simple_multi_failure_summary_txt(
    const std::filesystem::path& path,
    Vertex_handle v,
    const std::vector<Point>& cycle_isovertices,
    const std::vector<std::vector<Triangle_3>>& cycle_triangles
) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }

    const int vertex_index = v->info().index;
    const Point center = v->point();
    const auto& cycles = v->info().facet_cycles;

    out << std::setprecision(17);
    out << "vertex_index " << vertex_index << "\n";
    out << "center " << center.x() << " " << center.y() << " " << center.z() << "\n";
    out << "num_cycles " << cycles.size() << "\n";

    for (size_t c = 0; c < cycles.size(); ++c) {
        const size_t facet_count = cycles[c].size();
        const Point p = (c < cycle_isovertices.size()) ? cycle_isovertices[c] : center;
        out << "cycle " << c << " facet_count " << facet_count << " isov "
            << p.x() << " " << p.y() << " " << p.z() << "\n";
    }

    auto kind_name = [](TriIntersectionKind k) -> const char* {
        switch (k) {
            case TriIntersectionKind::NONE: return "NONE";
            case TriIntersectionKind::DISJOINT_TRIANGLE_TRIANGLE: return "DISJOINT_TRIANGLE_TRIANGLE";
            case TriIntersectionKind::SHARED_VERTEX_T1_OPPOSITE_SEGMENT: return "SHARED_VERTEX_T1_OPPOSITE_SEGMENT";
            case TriIntersectionKind::SHARED_VERTEX_T2_OPPOSITE_SEGMENT: return "SHARED_VERTEX_T2_OPPOSITE_SEGMENT";
            case TriIntersectionKind::SHARED_VERTEX_INTERIOR_SEGMENT_T1_0: return "SHARED_VERTEX_INTERIOR_SEGMENT_T1_0";
            case TriIntersectionKind::SHARED_VERTEX_INTERIOR_SEGMENT_T1_1: return "SHARED_VERTEX_INTERIOR_SEGMENT_T1_1";
            case TriIntersectionKind::SHARED_VERTEX_INTERIOR_SEGMENT_T2_0: return "SHARED_VERTEX_INTERIOR_SEGMENT_T2_0";
            case TriIntersectionKind::SHARED_VERTEX_INTERIOR_SEGMENT_T2_1: return "SHARED_VERTEX_INTERIOR_SEGMENT_T2_1";
            default: return "UNKNOWN";
        }
    };

    auto dump_tri = [&](const char* label, const Triangle_3& t) {
        out << label << " "
            << t.vertex(0).x() << " " << t.vertex(0).y() << " " << t.vertex(0).z() << "  "
            << t.vertex(1).x() << " " << t.vertex(1).y() << " " << t.vertex(1).z() << "  "
            << t.vertex(2).x() << " " << t.vertex(2).y() << " " << t.vertex(2).z() << "\n";
    };

    const bool within = cycle_triangles_have_within_cycle_intersection(cycle_triangles);
    const auto between_pair = first_between_cycle_intersection_pair(cycle_triangles);

    std::vector<std::pair<int, int>> coincident_cycles;
    {
        const size_t n = std::min(cycles.size(), cycle_isovertices.size());
        for (size_t a = 0; a < n; ++a) {
            for (size_t b = a + 1; b < n; ++b) {
                if (cycle_isovertices[a] == cycle_isovertices[b]) {
                    coincident_cycles.emplace_back(static_cast<int>(a), static_cast<int>(b));
                }
            }
        }
    }

    const bool any = within || between_pair.has_value() || !coincident_cycles.empty();

    out << "self_intersection_any " << (any ? 1 : 0) << "\n";
    if (!coincident_cycles.empty()) {
        out << "self_intersection_coincident_isovertices\n";
        for (const auto& p : coincident_cycles) {
            out << "coincident_cycles " << p.first << "-" << p.second << "\n";
        }
    }

    const auto within_witness = first_within_cycle_intersection_trace(cycle_triangles);
    if (within) {
        out << "self_intersection_within_cycles\n";
    }
    if (within_witness) {
        out << "within_witness cycle " << within_witness->cycle0
            << " tri(" << within_witness->tri0 << "," << within_witness->tri1 << ")"
            << " kind=" << kind_name(within_witness->trace.kind)
            << " shared_count=" << within_witness->trace.shared_count
            << " shared_point=(" << within_witness->trace.shared_point.x()
            << "," << within_witness->trace.shared_point.y()
            << "," << within_witness->trace.shared_point.z() << ")\n";
        const Triangle_3& t0 = cycle_triangles[static_cast<size_t>(within_witness->cycle0)][within_witness->tri0];
        const Triangle_3& t1 = cycle_triangles[static_cast<size_t>(within_witness->cycle1)][within_witness->tri1];
        dump_tri("within_tri0", t0);
        dump_tri("within_tri1", t1);
    }

    const auto between_witness = first_between_cycle_intersection_trace(cycle_triangles);
    if (between_pair) {
        out << "self_intersection_between_cycles " << between_pair->first << "-" << between_pair->second << "\n";
    }
    if (between_witness) {
        out << "between_witness cycle(" << between_witness->cycle0 << "," << between_witness->cycle1 << ")"
            << " tri(" << between_witness->tri0 << "," << between_witness->tri1 << ")"
            << " kind=" << kind_name(between_witness->trace.kind)
            << " shared_count=" << between_witness->trace.shared_count
            << " shared_point=(" << between_witness->trace.shared_point.x()
            << "," << between_witness->trace.shared_point.y()
            << "," << between_witness->trace.shared_point.z() << ")\n";
        const Triangle_3& t0 = cycle_triangles[static_cast<size_t>(between_witness->cycle0)][between_witness->tri0];
        const Triangle_3& t1 = cycle_triangles[static_cast<size_t>(between_witness->cycle1)][between_witness->tri1];
        dump_tri("between_tri0", t0);
        dump_tri("between_tri1", t1);
    }

    return true;
}

Point reflect_through_center(const Point& center, const Point& p) {
    return Point(
        2.0 * center.x() - p.x(),
        2.0 * center.y() - p.y(),
        2.0 * center.z() - p.z());
}

} // namespace

namespace vdc_debug {

bool trace_selfi_enabled(const int vertex_index) {
    const auto& traced = trace_selfi_vertices_internal();
    return (!traced.empty()) && (traced.find(vertex_index) != traced.end());
}

const std::string& trace_selfi_dump_dir() {
    static std::string dir;
    static bool initialized = false;
    if (initialized) {
        return dir;
    }
    initialized = true;
    if (const char* env = std::getenv("VDC_TRACE_SELFI_DUMP_DIR")) {
        dir = env;
    }
    return dir;
}

bool log_unresolved_A_vertices() {
    static int cached = -1;
    if (cached >= 0) {
        return cached != 0;
    }
    const char* env = std::getenv("VDC_LOG_UNRESOLVED_A");
    if (!env || env[0] == '\0') {
        cached = 0;
        return false;
    }
    cached = (std::strcmp(env, "0") == 0) ? 0 : 1;
    return cached != 0;
}

bool log_selfi_within_cycle_vertices() {
    static int cached = -1;
    if (cached >= 0) {
        return cached != 0;
    }
    const char* env = std::getenv("VDC_LOG_SELFI_WITHIN");
    if (!env || env[0] == '\0') {
        cached = 0;
        return false;
    }
    cached = (std::strcmp(env, "0") == 0) ? 0 : 1;
    return cached != 0;
}

bool trace_selfi_intersection_details_enabled() {
    static int cached = -1;
    if (cached >= 0) {
        return cached != 0;
    }
    const char* env = std::getenv("VDC_TRACE_SELFI_INTERSECTION_DETAILS");
    if (!env || env[0] == '\0') {
        cached = 0;
        return false;
    }
    cached = (std::strcmp(env, "0") == 0) ? 0 : 1;
    return cached != 0;
}

void maybe_dump_selfi_stage(
    const Delaunay& dt,
    Vertex_handle v,
    const std::vector<Point>& cycle_isovertices,
    const std::vector<Cell_handle>& cell_by_index,
    const char* stage
) {
    const int idx = v->info().index;
    if (!trace_selfi_enabled(idx)) {
        return;
    }
    const std::string& dir = trace_selfi_dump_dir();
    if (dir.empty()) {
        return;
    }

    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    if (ec) {
        return;
    }

    const auto& cycles = v->info().facet_cycles;
    std::vector<std::vector<TrianglePoints>> cycle_triangles(cycles.size());
    for (size_t c = 0; c < cycles.size(); ++c) {
        if (c >= cycle_isovertices.size()) {
            break;
        }
        cycle_triangles[c] = collect_cycle_triangles_points(
            dt, v, static_cast<int>(c), cycle_isovertices[c], cell_by_index);
    }

    (void)write_vtk_selfi_local_triangles(
        selfi_dump_path_for_vertex(idx, stage, ".vtk"), cycle_triangles);
    (void)write_vtk_selfi_points(
        selfi_dump_path_for_vertex(idx, stage, "_points.vtk"), v, cycle_isovertices);
}

void maybe_dump_selfi_cycle_metadata(
    const Delaunay& dt,
    Vertex_handle v,
    const std::vector<Cell_handle>& cell_by_index
) {
    const int idx = v->info().index;
    if (!trace_selfi_enabled(idx)) {
        return;
    }
    const std::string& dir = trace_selfi_dump_dir();
    if (dir.empty()) {
        return;
    }

    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    if (ec) {
        return;
    }

    const std::filesystem::path path = selfi_dump_path_for_vertex(idx, "cycles", ".txt");
    std::ofstream out(path);
    if (!out) {
        return;
    }

    const Point center = v->point();
    out << "vertex_index: " << idx << "\n";
    out << "center: " << center.x() << " " << center.y() << " " << center.z() << "\n";
    out << "num_cycles: " << v->info().facet_cycles.size() << "\n";

    for (size_t c = 0; c < v->info().facet_cycles.size(); ++c) {
        const auto& cycle = v->info().facet_cycles[c];

        std::unordered_set<int> boundary_indices;
        std::vector<int> boundary_delv;
        boundary_indices.reserve(cycle.size() * 2);
        boundary_delv.reserve(cycle.size() * 2);

        for (const auto& [cell_idx, facet_idx] : cycle) {
            Cell_handle cell = lookup_cell_local(cell_by_index, cell_idx);
            if (cell == Cell_handle()) {
                continue;
            }
            if (dt.is_infinite(cell)) {
                continue;
            }
            if (facet_idx < 0 || facet_idx >= 4) {
                continue;
            }

            for (int j = 0; j < 4; ++j) {
                if (j == facet_idx) {
                    continue;
                }
                Vertex_handle vh = cell->vertex(j);
                if (vh == v) {
                    continue;
                }
                if (dt.is_infinite(vh) || vh->info().is_dummy) {
                    continue;
                }
                const int delv_idx = vh->info().index;
                if (boundary_indices.insert(delv_idx).second) {
                    boundary_delv.push_back(delv_idx);
                }
            }
        }

        out << "\ncycle " << c << ":\n";
        out << "  num_facets: " << cycle.size() << "\n";
        out << "  boundary_delv_count: " << boundary_delv.size() << "\n";
        out << "  boundary_delv_indices:";
        for (int delv_idx : boundary_delv) {
            out << " " << delv_idx;
        }
        out << "\n";
        out << "  facets (cell_index,facet_index):\n";
        for (const auto& [cell_idx, facet_idx] : cycle) {
            out << "    " << cell_idx << " " << facet_idx << "\n";
        }
    }
}

void dump_simple_multi_failure_stage(
    const std::filesystem::path& out_dir,
    Vertex_handle v,
    const std::vector<Point>& cycle_isovertices,
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index,
    const char* stage
) {
    const int vertex_index = v->info().index;
    const size_t num_cycles = v->info().facet_cycles.size();

    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);
    if (ec) {
        return;
    }

    std::ostringstream base;
    base << "selfi_v" << vertex_index << "_c" << num_cycles << "_" << stage;
    const std::string prefix = base.str();

    const auto cycle_tris_keyed = collect_selfi_cycle_triangles_keyed(dt, v, cycle_isovertices, cell_by_index);
    std::vector<std::vector<Triangle_3>> cycle_triangles(cycle_tris_keyed.size());
    for (size_t c = 0; c < cycle_tris_keyed.size(); ++c) {
        cycle_triangles[c].reserve(cycle_tris_keyed[c].size());
        for (const auto& tri : cycle_tris_keyed[c]) {
            cycle_triangles[c].push_back(Triangle_3(tri.p[0], tri.p[1], tri.p[2]));
        }
    }

    (void)write_simple_multi_failure_summary_txt(out_dir / (prefix + ".txt"), v, cycle_isovertices, cycle_triangles);

    std::vector<KeyedTriangle> all_tris;
    size_t total = 0;
    for (const auto& tris : cycle_tris_keyed) {
        total += tris.size();
    }
    all_tris.reserve(total);
    for (const auto& tris : cycle_tris_keyed) {
        all_tris.insert(all_tris.end(), tris.begin(), tris.end());
    }

    (void)write_off_keyed_triangles(out_dir / (prefix + "_all.off"), all_tris);

    for (size_t c = 0; c < cycle_tris_keyed.size(); ++c) {
        std::ostringstream name;
        name << prefix << "_cycle" << c << ".off";
        (void)write_off_keyed_triangles(out_dir / name.str(), cycle_tris_keyed[c]);
    }
}

void dump_multi_isov_trace_case(
    const std::filesystem::path& out_dir,
    Vertex_handle v,
    const std::vector<Point>& baseline_positions,
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index
) {
    dump_simple_multi_failure_stage(out_dir, v, baseline_positions, dt, cell_by_index, "baseline");

    const int num_cycles = static_cast<int>(baseline_positions.size());
    if (num_cycles < 2) {
        return;
    }

    const Point center = v->point();

    if (num_cycles == 2) {
        {
            std::vector<Point> cand = baseline_positions;
            cand[0] = reflect_through_center(center, baseline_positions[1]);
            dump_simple_multi_failure_stage(out_dir, v, cand, dt, cell_by_index, "reflect_A");
        }
        {
            std::vector<Point> cand = baseline_positions;
            cand[1] = reflect_through_center(center, baseline_positions[0]);
            dump_simple_multi_failure_stage(out_dir, v, cand, dt, cell_by_index, "reflect_B");
        }
        {
            std::vector<Point> cand = baseline_positions;
            cand[0] = reflect_through_center(center, baseline_positions[0]);
            dump_simple_multi_failure_stage(out_dir, v, cand, dt, cell_by_index, "reflect_self0");
        }
        {
            std::vector<Point> cand = baseline_positions;
            cand[1] = reflect_through_center(center, baseline_positions[1]);
            dump_simple_multi_failure_stage(out_dir, v, cand, dt, cell_by_index, "reflect_self1");
        }
        {
            std::vector<Point> cand = baseline_positions;
            cand[0] = reflect_through_center(center, baseline_positions[0]);
            cand[1] = reflect_through_center(center, baseline_positions[1]);
            dump_simple_multi_failure_stage(out_dir, v, cand, dt, cell_by_index, "reflect_self_both");
        }
        {
            std::vector<Point> cand = baseline_positions;
            cand[0] = reflect_through_center(center, baseline_positions[1]);
            cand[1] = reflect_through_center(center, baseline_positions[0]);
            dump_simple_multi_failure_stage(out_dir, v, cand, dt, cell_by_index, "reflect_cross_both");
        }
    } else {
        for (int c = 0; c < num_cycles; ++c) {
            std::vector<Point> cand = baseline_positions;
            cand[static_cast<size_t>(c)] =
                reflect_through_center(center, baseline_positions[static_cast<size_t>(c)]);
            std::ostringstream stage;
            stage << "reflect_cycle" << c;
            dump_simple_multi_failure_stage(out_dir, v, cand, dt, cell_by_index, stage.str().c_str());
        }

        for (int a = 0; a < num_cycles; ++a) {
            for (int b = 0; b < num_cycles; ++b) {
                if (a == b) continue;
                std::vector<Point> cand = baseline_positions;
                cand[static_cast<size_t>(a)] =
                    reflect_through_center(center, baseline_positions[static_cast<size_t>(b)]);
                std::ostringstream stage;
                stage << "reflect_cycle" << a << "_from" << b;
                dump_simple_multi_failure_stage(out_dir, v, cand, dt, cell_by_index, stage.str().c_str());
            }
        }
    }

    {
        std::vector<Point> attempt = baseline_positions;
        std::vector<Point> attempt_out;
        (void)try_resolve_multicycle_by_cycle_separation_tests(
            v, attempt, dt, cell_by_index, &attempt_out);
        if (!attempt_out.empty()) {
            dump_simple_multi_failure_stage(out_dir, v, attempt_out, dt, cell_by_index, "pairwise_attempt");
        }
    }
}

void finalize_multi_isov_trace(
    const CycleIsovertexOptions& options,
    const std::unordered_set<int>& local_unresolved_dumped_vertices,
    const std::vector<Vertex_handle>& multi_cycle_vertices,
    const Delaunay& dt,
    const std::vector<Cell_handle>& cell_by_index
) {
    if (!options.multi_isov_trace || options.multi_isov_trace_dir.empty()) {
        return;
    }

    const std::filesystem::path trace_root(options.multi_isov_trace_dir);
    const std::filesystem::path local_dir = trace_root / "local";
    const std::filesystem::path final_dir = trace_root / "final";

    {
        std::error_code ec;
        std::filesystem::create_directories(local_dir, ec);
        ec.clear();
        std::filesystem::create_directories(final_dir, ec);
    }

    {
        std::vector<int> local_unresolved;
        local_unresolved.reserve(local_unresolved_dumped_vertices.size());
        for (int idx : local_unresolved_dumped_vertices) {
            local_unresolved.push_back(idx);
        }
        std::sort(local_unresolved.begin(), local_unresolved.end());

        std::ofstream out(local_dir / "local_unresolved_vertices.txt");
        if (out) {
            out << "count " << local_unresolved.size() << "\n";
            out << "vertex_indices";
            for (int idx : local_unresolved) {
                out << " " << idx;
            }
            out << "\n";
        }
        std::cerr << "[DEL-SELFI-TRACE] Local unresolved vertices dumped: "
                  << local_unresolved.size() << "\n";
    }

    std::vector<int> final_unresolved;
    final_unresolved.reserve(32);

    for (Vertex_handle vh : multi_cycle_vertices) {
        if (!vh->info().active) continue;
        if (vh->info().is_dummy) continue;
        if (vh->info().facet_cycles.size() < 2) continue;

        const auto& isovertices = vh->info().cycle_isovertices;
        if (!check_self_intersection(vh, isovertices, dt, cell_by_index)) {
            continue;
        }

        final_unresolved.push_back(vh->info().index);
        dump_multi_isov_trace_case(final_dir, vh, isovertices, dt, cell_by_index);
    }

    {
        std::ofstream out(final_dir / "final_unresolved_vertices.txt");
        if (out) {
            out << "count " << final_unresolved.size() << "\n";
            out << "vertex_indices";
            for (int idx : final_unresolved) {
                out << " " << idx;
            }
            out << "\n";
        }
    }
    std::cerr << "[DEL-SELFI-TRACE] Final unresolved vertices: "
              << final_unresolved.size() << "\n";
}

} // namespace vdc_debug
