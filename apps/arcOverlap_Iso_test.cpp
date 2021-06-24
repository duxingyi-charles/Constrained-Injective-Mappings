//
// Created by Charles Du on 6/20/21.
//

#include <CLI/CLI.hpp>
#include <string>
#include <nlohmann/json.hpp>

#include "ScopedTimer.h"
#include "Arc_Overlap_Iso_Formulation.h"

using namespace Eigen;



void parse_input_file(const std::string &input_file,
                      double &alpha, double &theta, std::string &form,
                      bool& scale_rest_mesh,
                      MatrixXd &rest_vertices, Matrix2Xd &init_vertices,
                      Matrix3Xi &faces, VectorXi &handles) {
    using json = nlohmann::json;
    std::ifstream fin(input_file.c_str());
    if (!fin) {
        throw std::runtime_error("Config file does not exist!");
    }
    json input;
    fin >> input;
    fin.close();

    assert(input.contains("alpha"));
    alpha = input["alpha"];

    assert(input.contains("theta"));
    theta = input["theta"];

    assert(input.contains("form"));
    form = input["form"];

    assert(input.contains("scale_rest_mesh"));
    scale_rest_mesh = input["scale_rest_mesh"];

    assert(input.contains("restV"));
    int nv = input["restV"].size();
    int vDim = input["restV"][0].size();
    rest_vertices.resize(vDim, nv);
    for (int i = 0; i < nv; ++i) {
        for (int j = 0; j < vDim; ++j) {
            rest_vertices(j, i) = input["restV"][i][j];
        }
    }

    assert(input.contains("initV"));
    nv = input["initV"].size();
    vDim = input["initV"][0].size();
    init_vertices.resize(vDim, nv);
    for (int i = 0; i < nv; ++i) {
        for (int j = 0; j < vDim; ++j) {
            init_vertices(j, i) = input["initV"][i][j];
        }
    }

    assert(input.contains("faces"));
    int nf = input["faces"].size();
    int simplex_size = input["faces"][0].size();
    faces.resize(simplex_size, nf);
    for (int i = 0; i < nf; ++i) {
        for (int j = 0; j < simplex_size; ++j) {
            faces(j, i) = input["faces"][i][j];
        }
    }

    assert(input.contains("handles"));
    int nHdls = input["handles"].size();
    handles.resize(nHdls);
    for (int i = 0; i < nHdls; ++i) {
        handles(i) = input["handles"][i];
    }

}

void export_result(const std::string &output_file, double energy, const Eigen::VectorXd &grad)
{
    using json = nlohmann::json;
    json output;

    // put data in json object
    output["energy"] = energy;

    output["grad"] = json::array();
    for (int i = 0; i < grad.size(); ++i) {
        output["grad"].push_back(grad(i));
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

    CLI::App app{"arcOverlap_Isometric Test"};
    app.add_option("input_file", args.input_file, "input file")
            ->required();
    app.add_option("output_file", args.output_file, "output file")
            ->required();
    CLI11_PARSE(app, argc, argv);

    // load input
    double alpha;
    double theta;
    std::string form;
    bool scale_rest_mesh;

    MatrixXd rest_vertices;
    Matrix2Xd init_vertices;
    Matrix3Xi faces;
    VectorXi handles;

    parse_input_file(args.input_file, alpha,theta, form, scale_rest_mesh,rest_vertices, init_vertices,
                     faces, handles);

    // debug: print out input
//    std::cout << "alphaRatio = " << alphaRatio << std::endl;
//    std::cout << "theta = " << theta << std::endl;
//    std::cout << "form = " << form << std::endl;
//    std::cout << "restV = \n" << rest_vertices << std::endl;
//    std::cout << "initV = \n" << init_vertices << std::endl;
//    std::cout << "faces = \n" << faces << std::endl;
//    std::cout << "handles = \n" << handles << std::endl;

    // initialize arc overlap formulation
    Arc_Overlap_Iso_Formulation formulation(rest_vertices,init_vertices,faces,handles,form,alpha,theta,scale_rest_mesh);

    // debug: print boundary edges
//    auto boundary_edges = formulation.get_boundary_edges();
//    std::cout << "boundary edges = " << std::endl;
//    for (const auto & e : boundary_edges) {
//        std::cout << e.first << ", " << e.second << std::endl;
//    }

//    auto freeI = formulation.get_freeI();
//    std::cout << "freeI = \n" << freeI << std::endl;

    // compute energy
//    double energy = formulation.compute_energy(formulation.get_x0());
//    std::cout << "energy (x0) = " << energy << std::endl;

    // compute energy and gradient
    VectorXd grad;
    double energy = formulation.compute_energy_with_gradient(formulation.get_x0(), grad);

    Eigen::VectorXd energyList;
    formulation.compute_energy(formulation.get_x0(), energyList);

    // energy breakdown
//    int nF = faces.cols();
//    double tlc_iso_energy = energyList.head(nF).sum();
//    double arc_occupancy = energyList(energyList.size()-1);
//    std::cout << "TLC_Iso energy = " << tlc_iso_energy << std::endl;

//    std::cout << "energy = " << energy << std::endl;
//    std::cout << "grad = \n" << grad << std::endl;


    // output result
    export_result(args.output_file, energy, grad);

    return 0;
}
