//
// Created by Charles Du on 4/27/21.
//


#include <CLI/CLI.hpp>
#include <string>
#include <nlohmann/json.hpp>

#include "Arrangement.h"

void parse_input_file(const std::string &input_file,
                      std::vector<Point> &verts, std::vector<Arc_Edge> &edges) {
    using json = nlohmann::json;
    std::ifstream fin(input_file.c_str());
    if (!fin) {
        throw std::runtime_error("Config file does not exist!");
    }
    json input;
    fin >> input;
    fin.close();

    assert(input.contains("verts"));
    verts.clear();
    for (const auto &p : input["verts"]) {
        verts.emplace_back(p[0], p[1]);
    }

    assert(input.contains("theta"));
    double theta = input["theta"];

    assert(input.contains("edges"));
    edges.clear();
    for (const auto &e : input["edges"]) {
        auto i = e[0].get<size_t>();
        auto j = e[1].get<size_t>();
        edges.emplace_back(Arc_Edge{i, j, Circular_Arc(verts[i], verts[j], theta)});
    }
}

void export_result(const std::string &output_file,
                   const std::vector<Point> &pts, const std::vector<SubArc_Edge> &pEdges,
                   const std::vector<bool> &is_intersection_point,
                   const std::vector<std::vector<size_t>> &cells,
                   const std::vector<int> &windings)
{
    using json = nlohmann::json;
    json output;

    // put data in json object
    output["pts"] = json::array();
    for (const auto & pt : pts) {
        output["pts"].push_back({pt.x(), pt.y()});

    }
    output["is_intersection"] = json::array();
    for (const bool b : is_intersection_point) {
        output["is_intersection"].push_back(b);
    }

    output["edges"] = json::array();
    output["angles"] = json::array();
    output["parent_edges"] = json::array();
    for (const auto &e : pEdges) {
        output["edges"].push_back({e.id1, e.id2});
        output["angles"].push_back(e.arc.get_arc_angle());
        output["parent_edges"].push_back(e.parent_id);
    }

    output["cells"] = json::array();
    for (const auto & cell : cells) {
        json js = json::array();
        for (unsigned long j : cell) {
            js.push_back(j);
        }
        output["cells"].push_back(js);
    }

    output["windings"] = json::array();
    for (int w : windings) {
        output["windings"].push_back(w);
    }

    // write JSON
    std::ofstream fout(output_file);
    fout << std::setw(4) << output << std::endl;
    fout.close();
}

int main(int argc, char **argv)
{
    struct {
        std::string input_file;
        std::string output_file;
    } args;

    CLI::App app{"Subdivide polyArc Test"};
    app.add_option("input_file", args.input_file, "input file")
            ->required();
    app.add_option("output_file", args.output_file, "output file")
            ->required();
    CLI11_PARSE(app, argc, argv);

    // load input
    std::vector<Point> verts;
    std::vector<Arc_Edge> edges;

    parse_input_file(args.input_file, verts, edges);

    // debug: check input
//    std::cout << "verts: " << std::endl;
//    for (const auto & p : verts) {
//        std::cout << "(" << p.x() << ", " << p.y() << ")" << std::endl;
//    }
//    std::cout << "edges: " << std::endl;
//    for (const auto & e : edges) {
//        std::cout << "(" << e.id1 << ", " << e.id2 << ")" << std::endl;
//    }

    // subdivide input arcs by their intersections
    std::vector<Point> pts;
    std::vector<SubArc_Edge> pEdges;
    std::vector<bool> is_intersection_point;
    std::vector<int>  arc1_of_intersection;
    std::vector<int>  arc2_of_intersection;
    Arrangement::subdivide_polyArc_by_intersection(verts, edges, pts, pEdges, is_intersection_point,
                                                   arc1_of_intersection, arc2_of_intersection);

    // debug: check output
//    std::cout << "pts: " << std::endl;
//    for (const auto & p : pts) {
//        std::cout << "(" << p.x() << ", " << p.y() << ")" << std::endl;
//    }
//    std::cout << "sub edges: " << std::endl;
//    for (const auto & e : pEdges) {
//        std::cout << "(" << e.id1 << ", " << e.id2 << ")" << std::endl;
//    }

    // decompose into cells
    std::vector<std::vector<size_t>> eIn;
    std::vector<std::vector<size_t>> eOut;
    std::vector<std::vector<size_t>> cells;
    Arrangement::decompose_into_cells(pts, pEdges, eIn, eOut, cells);

    // compute winding numbers for each cell
    std::vector<int> windings;
    Arrangement::compute_cell_windings(pEdges, eIn, eOut, cells, windings);


    // output result
    export_result(args.output_file, pts, pEdges, is_intersection_point, cells, windings);

    return 0;
}
