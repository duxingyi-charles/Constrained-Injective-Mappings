//
// Created by Charles Du on 4/29/21.
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

void export_result(const std::string &output_file, double occupancy, const Eigen::Matrix2Xd &grad)
{
    using json = nlohmann::json;
    json output;

    // put data in json object
    output["occupancy"] = occupancy;

    output["grad"] = json::array();
    for (int i = 0; i < grad.cols(); ++i) {
        output["grad"].push_back({grad.col(i).x(), grad.col(i).y()});
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
    std::vector<std::pair<size_t,size_t>> edges;
    double theta;

    parse_input_file(args.input_file, verts, edges, theta);


    // compute arc occupancy
    double occupancy;
    Eigen::Matrix2Xd grad;
    {
        Arc_Occupancy arcOccupancy(theta);
        ScopedTimer<> timer("energy + gradient");
        occupancy = arcOccupancy.compute_arc_occupancy_with_gradient(verts, edges, grad);
    }


    // output result
    export_result(args.output_file, occupancy, grad);

    return 0;
}