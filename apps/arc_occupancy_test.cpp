//
// Created by Charles Du on 4/28/21.
//

#include <CLI/CLI.hpp>
#include <string>
#include <nlohmann/json.hpp>

#include "ScopedTimer.h"
#include "Arc_Occupancy.h"

void parse_input_file(const std::string &input_file,
                      std::vector<Point> &verts, std::vector<std::pair<size_t,size_t>> &edges, double &theta) {
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
    theta = input["theta"];

    assert(input.contains("edges"));
    edges.clear();
    for (const auto &e : input["edges"]) {
        auto i = e[0].get<size_t>();
        auto j = e[1].get<size_t>();
        edges.emplace_back(i, j);
    }
}

void export_result(const std::string &output_file, double occupancy)
{
    using json = nlohmann::json;
    json output;

    // put data in json object
    output["occupancy"] = occupancy;

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
    std::vector<std::pair<size_t,size_t>> edges;
    double theta;

    parse_input_file(args.input_file, verts, edges, theta);

    // debug: check input
//    std::cout << "verts: " << std::endl;
//    for (const auto & p : verts) {
//        std::cout << "(" << p.x() << ", " << p.y() << ")" << std::endl;
//    }
//    std::cout << "edges: " << std::endl;
//    for (const auto & e : edges) {
//        std::cout << "(" << e.id1 << ", " << e.id2 << ")" << std::endl;
//    }

    // compute arc occupancy
    double occupancy;
    {
        ScopedTimer<> timer("energy");
        occupancy = Arc_Occupancy::compute_arc_occupancy(verts, edges, theta);
    }


    // debug: check output
//    std::cout << "pts: " << std::endl;
//    for (const auto & p : pts) {
//        std::cout << "(" << p.x() << ", " << p.y() << ")" << std::endl;
//    }
//    std::cout << "sub edges: " << std::endl;
//    for (const auto & e : pEdges) {
//        std::cout << "(" << e.id1 << ", " << e.id2 << ")" << std::endl;
//    }


    // output result
    export_result(args.output_file, occupancy);

    return 0;
}