//
// Created by Charles Du on 4/28/21.
//

#include <CLI/CLI.hpp>
#include <string>
#include <nlohmann/json.hpp>

#include "ScopedTimer.h"

#include "Arc_Occupancy.h"

void parse_input_file(const std::string &input_file,
                      std::vector<Point> &verts, std::vector<std::pair<size_t,size_t>> &edges, var &theta) {
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
        double x = p[0];
        double y = p[1];
        verts.emplace_back(x, y);
    }

    assert(input.contains("theta"));
    double t = input["theta"];
    theta = t;

    assert(input.contains("edges"));
    edges.clear();
    for (const auto &e : input["edges"]) {
        auto i = e[0].get<size_t>();
        auto j = e[1].get<size_t>();
        edges.emplace_back(i, j);
    }
}

void export_result(const std::string &output_file, double occupancy,
                   const Eigen::VectorXd &grad)
//                   ,const Eigen::MatrixXd &Hess)
{
    using json = nlohmann::json;
    json output;

    // put data in json object
    output["occupancy"] = occupancy;

    output["grad"] = json::array();
    for (int i = 0; i < grad.size(); ++i) {
        output["grad"].push_back(grad[i]);
    }

//    output["Hess"] = json::array();
//    for (int i = 0; i < Hess.rows(); ++i) {
//        json jrow = json::array();
//        for (int j = 0; j < Hess.cols(); ++j) {
//            jrow.push_back(Hess(i,j));
//        }
//        output["Hess"].push_back(jrow);
//    }

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
    var theta;

    parse_input_file(args.input_file, verts, edges, theta);


    //
    Eigen::VectorXvar flatten_verts(verts.size()*2);
    for (int i = 0; i < verts.size(); ++i) {
        flatten_verts[2*i] = verts[i].x();
        flatten_verts[2*i+1] = verts[i].y();
    }

    // compute arc occupancy
//    var occupancy = Arc_Occupancy::compute_arc_occupancy(verts, edges, theta);
    Arc_Occupancy obj(theta);
    var occupancy;
    {
        ScopedTimer<> timer("energy");
        occupancy = obj.compute_arc_occupancy(flatten_verts, edges);
    }
    double occu_val = val(occupancy);

    // derivative to scalar
//    auto [dOccu_dTheta] = autodiff::derivatives(occupancy, autodiff::wrt(theta));
//    std::cout << "dOccu_dTheta = " << dOccu_dTheta << std::endl;

    // gradient
//    for (auto & v : verts) {
//        auto [dx, dy] = autodiff::derivatives(occupancy, autodiff::wrt(v.x(), v.y()));
//        std::cout << "{ " << dx << ", " << dy << " }" << std::endl;
//    }

    Eigen::VectorXd grad;
    {
        ScopedTimer<> timer("gradient");
        grad = autodiff::gradient(occupancy, flatten_verts);
    }
//    std::cout << "grad = \n" << std::endl;
//    for (int i = 0; i < verts.size(); ++i) {
//        std::cout << "{ " << grad[2*i] << ", " << grad[2*i+1] << " }" << std::endl;
//    }

    // Hessian (and gradient)
//    Eigen::VectorXd g;
//    Eigen::MatrixXd H;
//    {
//        ScopedTimer<> timer("Hessian");
//        H = autodiff::hessian(occupancy, flatten_verts, g);
//    }
//    std::cout << "g = \n" << std::endl;
//    for (int i = 0; i < verts.size(); ++i) {
//        std::cout << "{ " << g[2*i] << ", " << g[2*i+1] << " }" << std::endl;
//    }
//    std::cout << "H = \n" << H << std::endl;



    // output result
    export_result(args.output_file, occu_val, grad);

    return 0;
}