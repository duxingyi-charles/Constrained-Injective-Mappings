//
// Created by Charles Du on 3/29/22.
//

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include <chrono>

#include <nlopt.hpp>

#include "sTGC_Formulation.h"
#include "optimization_util.h"
#include "deformation_gradient_util.h"

using namespace Eigen;


// solver options
class SolverOptionManager
{
public:
    //default options
    SolverOptionManager() :
            form("Tutte"), alpha(1),
            lambda1(0.5), lambda2(0.), k(1.), // lambda1 + k * lambda2 = 1/2
            scale_rest_mesh(false), subtract_total_signed_area(true),
            aspect_ratio_threshold(0),
            ftol_abs(1e-8), ftol_rel(1e-8), xtol_abs(1e-8), xtol_rel(1e-8),
            maxeval(10000), algorithm("LBFGS"), stopCode("no_flip_degenerate"), record()
    {};
    //import options from file
    explicit SolverOptionManager(const char* filename) :
            form("Tutte"), alpha(1),
            lambda1(0.5), lambda2(0.), k(1.),
            scale_rest_mesh(false), subtract_total_signed_area(true),
            aspect_ratio_threshold(0),
            ftol_abs(1e-8), ftol_rel(1e-8), xtol_abs(1e-8), xtol_rel(1e-8),
            maxeval(10000), algorithm("LBFGS"), stopCode("no_flip_degenerate"), record()
    {
        if (!importOptions(filename))
        {
            std::cout << "SolverOptionManager Warn: default options are used." << std::endl;
        }
    };

    ~SolverOptionManager() = default;

    // energy formulation options
    std::string form;
    double alpha;
    bool scale_rest_mesh;
    double lambda1;
    double lambda2;
    double k;

    // whether to subtract total signed area
    bool subtract_total_signed_area;

    // rest triangle aspect ratio threshold
    // triangles above the threshold will be replaced by a better-shaped one with the same area
    // if the threshold <=1, no rest triangle will be changed
    double aspect_ratio_threshold;

    // optimization options
    double ftol_abs;
    double ftol_rel;
    double xtol_abs;
    double xtol_rel;
    int maxeval;
    std::string algorithm;
    std::string stopCode;
    std::vector<std::string> record;

    void printOptions()
    {
        std::cout << "form:\t" << form << "\n";
        std::cout << "alpha:\t" << alpha << "\n";
        std::cout << "lambda1:\t" << lambda1 << "\n";
        std::cout << "lambda2:\t" << lambda2 << "\n";
        std::cout << "k:\t" << k << "\n";
        std::cout << "scale_rest_mesh:\t" << scale_rest_mesh << "\n";
        std::cout << "subtract_total_signed_area:\t" << subtract_total_signed_area << "\n";
        std::cout << "aspect_ratio_threshold:\t" << aspect_ratio_threshold << "\n";
        std::cout << "ftol_abs:\t" << ftol_abs << "\n";
        std::cout << "ftol_rel:\t" << ftol_rel << "\n";
        std::cout << "xtol_abs:\t" << xtol_abs << "\n";
        std::cout << "xtol_rel:\t" << xtol_rel << "\n";
        std::cout << "maxeval:\t" << maxeval << "\n";
        std::cout << "algorithm:\t" << algorithm << "\n";
        std::cout << "stopCode:\t" << stopCode << "\n";
        std::cout << "record:  \t" << "{ ";
        for (auto & i : record)
        {
            std::cout << i << " ";
        }
        std::cout << "}" << std::endl;
    }

    bool importOptions(const char* filename)
    {
        //open the data file
        std::ifstream in_file(filename);

        if (!in_file.is_open())
        {
            std::cerr << "Failed to open solver option file:" << filename << "!" << std::endl;
            return false;
        }

        //read the file
        std::string abnormal = "normal";
        std::string optName;
        while (true)
        {
            in_file >> optName;
            if (optName != "form")
            {
                abnormal = "form";
                break;
            }
            in_file >> form;

            //            in_file >> optName;
            //            if (optName != "alphaRatio")
            //            {
            //                abnormal = "alphaRatio";
            //                break;
            //            }
            //            in_file >> alphaRatio;

            in_file >> optName;
            if (optName != "alpha")
            {
                abnormal = "alpha";
                break;
            }
            in_file >> alpha;

//            in_file >> optName;
//            if (optName != "theta") {
//                abnormal = "theta";
//                break;
//            }
//            in_file >> theta;

            in_file >> optName;
            if (optName != "lambda1")
            {
                abnormal = "lambda1";
                break;
            }
            in_file >> lambda1;

            in_file >> optName;
            if (optName != "lambda2")
            {
                abnormal = "lambda2";
                break;
            }
            in_file >> lambda2;

            in_file >> optName;
            if (optName != "k")
            {
                abnormal = "k";
                break;
            }
            in_file >> k;


            in_file >> optName;
            if (optName != "scale_rest_mesh") {
                abnormal = "scale_rest_mesh";
                break;
            }
            in_file >> scale_rest_mesh;

            in_file >> optName;
            if (optName != "subtract_total_signed_area") {
                abnormal = "subtract_total_signed_area";
                break;
            }
            in_file >> subtract_total_signed_area;

            in_file >> optName;
            if (optName != "aspect_ratio_threshold") {
                abnormal = "aspect_ratio_threshold";
                break;
            }
            in_file >> aspect_ratio_threshold;

            in_file >> optName;
            if (optName != "ftol_abs")
            {
                abnormal = "ftol_abs";
                break;
            }
            in_file >> ftol_abs;

            in_file >> optName;
            if (optName != "ftol_rel")
            {
                abnormal = "ftol_rel";
                break;
            }
            in_file >> ftol_rel;

            in_file >> optName;
            if (optName != "xtol_abs")
            {
                abnormal = "xtol_abs";
                break;
            }
            in_file >> xtol_abs;

            in_file >> optName;
            if (optName != "xtol_rel")
            {
                abnormal = "xtol_rel";
                break;
            }
            in_file >> xtol_rel;

            in_file >> optName;
            if (optName != "algorithm")
            {
                abnormal = "algorithm";
                break;
            }
            in_file >> algorithm;

            in_file >> optName;
            if (optName != "maxeval")
            {
                abnormal = "maxeval";
                break;
            }
            in_file >> maxeval;

            in_file >> optName;
            if (optName != "stopCode")
            {
                abnormal = "stopCode";
                break;
            }
            in_file >> stopCode;

            // record values
            in_file >> optName;
            if (optName != "record")
            {
                abnormal = "record";
                break;
            }

            in_file >> optName;
            if (optName != "vert") {
                abnormal = "vert";
                break;
            }
            int selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("vert");
            }

            in_file >> optName;
            if (optName != "energy") {
                abnormal = "energy";
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("energy");
            }

            in_file >> optName;
            if (optName != "minArea") {
                abnormal = "minArea";
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("minArea");
            }

            in_file >> optName;
            if (optName != "nbWindVert") {
                abnormal = "nbWindVert";
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("nbWindVert");
            }

            in_file >> optName;
            if (optName != "grad") {
                abnormal = "grad";
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("grad");
            }

            in_file >> optName;
            if (optName != "gradNorm") {
                abnormal = "gradNorm";
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("gradNorm");
            }

            in_file >> optName;
            if (optName != "initSingularValues") {
                abnormal = "initSingularValues";
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("initSingularValues");
            }

            in_file >> optName;
            if (optName != "resultSingularValues") {
                abnormal = "resultSingularValues";
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("resultSingularValues");
            }

            break;
        }

        in_file.close();

        if (abnormal != "normal")
        {
            std::cout << "Err:(" << abnormal << ") fail to import options from file. Check file format." << std::endl;
            return false;
        }

        return true;
    }
};




//
class Optimization_Data
{
public:
    Optimization_Data(const MatrixXd& restV,
                      const Matrix2Xd& initV,
                      const Matrix3Xi& restF,
                      const VectorXi& handles,
                      const std::string& form,
                      double alpha,
                      double lambda1, double lambda2, double k,
                      bool scale_rest_mesh,
                      bool subtract_total_signed_area,
                      double aspect_ratio_threshold) :
            F(restF), solutionFound(false), custom_criteria_met(false),
            locally_injective_Found(false),
            first_locally_injective_iteration(-1), iteration_count(0), max_iterations(1),
            first_locally_injective_V(initV),
            non_flip_Found(false), first_non_flip_iteration(-1),
            first_non_flip_V(initV),
            lastFunctionValue(HUGE_VAL), stopCode("none"),
            nb_feval(0), nb_geval(0),
            record_vert(false), record_energy(false), record_gradient(false),
            record_gradient_norm(false), record_minArea(false),
            record_nb_winded_interior_vertices(false),
            record_init_singular_values(false), record_result_singular_values(false),
            vertRecord(0), energyRecord(0), minAreaRecord(0), gradRecord(0),
            formulation(restV, initV, restF, handles, form, alpha,
                        lambda1, lambda2, k,
                        scale_rest_mesh, subtract_total_signed_area,
                        aspect_ratio_threshold),
            is_boundary_vertex(initV.cols(), false)
    {
        x0 = formulation.get_x0();

        // find boundary vertex
        const auto& boundary_edges = formulation.get_boundary_edges();
        for (const auto& e : boundary_edges) {
            is_boundary_vertex[e.first] = true;
            is_boundary_vertex[e.second] = true;
        }

        // initialize lastGradient
        lastGradient.setZero(x0.size());

        // singular values of initial mesh
        compute_tri_mesh_singular_values(formulation.get_scaled_rest_vertices(), initV, F,
                                         init_singular_values);
    };

    ~Optimization_Data() = default;

    sTGC_Formulation formulation;

    //    Matrix2Xd V;
    //    VectorXi freeI;
    Matrix3Xi F;
    //    Matrix3Xd restD;
    VectorXd x0;

    std::vector<bool> is_boundary_vertex;

    bool custom_criteria_met;
    bool solutionFound;
    double lastFunctionValue;

    VectorXd lastGradient;

    bool locally_injective_Found;
    int  iteration_count;
    int  max_iterations;
    int  first_locally_injective_iteration;
    Matrix2Xd first_locally_injective_V;
    bool non_flip_Found;
    int first_non_flip_iteration;
    Matrix2Xd first_non_flip_V;

    int nb_feval;
    int nb_geval;

    //from options
    std::string stopCode;
    //record data
    bool record_vert;
    bool record_energy;
    bool record_gradient;
    bool record_gradient_norm;
    bool record_minArea;
    bool record_nb_winded_interior_vertices;
    bool record_init_singular_values;
    bool record_result_singular_values;
    std::vector<Matrix2Xd> vertRecord;
    std::vector<double> minAreaRecord;
    std::vector<double> energyRecord;
    std::vector<VectorXd> gradRecord;
    std::vector<double> gradNormRecord;
    std::vector<int> nb_winded_interior_vertices_Record;
    Matrix2Xd init_singular_values;


    // record information we cared about
    void set_record_flags(const std::vector<std::string>& record)
    {
        for (const auto & i : record)
        {
            if (i == "vert")    record_vert = true;
            if (i == "energy")  record_energy = true;
            if (i == "minArea") record_minArea = true;
            if (i == "nbWindVert") record_nb_winded_interior_vertices = true;
            if (i == "grad")  record_gradient = true;
            if (i == "gradNorm") record_gradient_norm = true;
            if (i == "initSingularValues") record_init_singular_values = true;
            if (i == "resultSingularValues") record_result_singular_values = true;
        }
    }

    void record()
    {
        if (record_vert) vertRecord.push_back(formulation.get_V());
        if (record_gradient) gradRecord.push_back(lastGradient);
        if (record_gradient_norm) gradNormRecord.push_back(lastGradient.norm());
        //you need to make sure lastFunctionValue is up-to-date
        if (record_energy) energyRecord.push_back(lastFunctionValue);
        if (record_minArea)
        {
            double minA = compute_min_signed_mesh_area(formulation.get_V(), F);
            minAreaRecord.push_back(minA);
        }
        if (record_nb_winded_interior_vertices) {
            std::vector<size_t> winded_vertices;
            compute_winded_interior_vertices(formulation.get_V(), F, is_boundary_vertex, winded_vertices);
            nb_winded_interior_vertices_Record.push_back(winded_vertices.size());
        }
    }

    //custom stop criteria
    bool stopQ()
    {
        if (stopCode == "no_flip_degenerate") {
            return stopQ_no_flip_degenerate();
        }
        else if (stopCode == "locally_injective") {
            return stopQ_locally_injective();
        }
//        else if (stopCode == "globally_injective") {
//            return stopQ_globally_injective();
//        }
        else {
            //default
            return false;
        }
    }

    bool stopQ_no_flip_degenerate()
    {
        bool no_flip_degenerate = true;

        if (record_minArea && !minAreaRecord.empty()) {
            if (minAreaRecord.back() <= 0) no_flip_degenerate = false;
        }
        else {
            double minA = compute_min_signed_mesh_area(formulation.get_V(), F); //todo: only check free faces
            if (minA <= 0) no_flip_degenerate = false;
        }

        if (!non_flip_Found && no_flip_degenerate) {
            non_flip_Found = true;
            first_non_flip_iteration = iteration_count;
            first_non_flip_V = formulation.get_V();
        }

        return no_flip_degenerate;
    }

    bool stopQ_locally_injective() {
        if (!stopQ_no_flip_degenerate()) {
            return false;
        }
        else {
            bool is_locally_injective = false;
            // check if there is any over-winding interior vertices
            if (record_nb_winded_interior_vertices && !nb_winded_interior_vertices_Record.empty()) {
                if (nb_winded_interior_vertices_Record.back() == 0) {
                    is_locally_injective = true;
                }
            }
            else {
                std::vector<size_t> winded_vertices;
                // todo: only compute winding number if the one ring of the vertex is free
                compute_winded_interior_vertices(formulation.get_V(), F, is_boundary_vertex, winded_vertices);
                if (winded_vertices.empty()) {
                    is_locally_injective = true;
                }
            }

            if (is_locally_injective && !locally_injective_Found) {
                locally_injective_Found = true;
                first_locally_injective_iteration = iteration_count;
                first_locally_injective_V = formulation.get_V();
            }
            return is_locally_injective;
        }
    }

    //export result
    bool exportResult(const char* filename)
    {
        std::ofstream out_file(filename);
        if (!out_file.is_open()) {
            std::cerr << "Failed to open " << filename << "!" << std::endl;
            return false;
        }

        //precision of output
        typedef std::numeric_limits< double > dbl;
        out_file.precision(dbl::max_digits10);

        //write the file

        /// resV
        const auto& V = formulation.get_V();
        auto nv = V.cols();
        auto ndim = V.rows();
        out_file << "resV " << nv << " " << ndim << "\n";
        for (int i = 0; i < nv; ++i)
        {
            for (int j = 0; j < ndim; ++j)
            {
                out_file << V(j, i) << " ";
            }
        }
        out_file << std::endl;

        // first_non_flip_iter
        out_file << "first_non_flip_iter " << 1 << " " << 1 << "\n";
        out_file << first_non_flip_iteration << " ";
        out_file << std::endl;

        /// first_non_flip_V
        nv = first_non_flip_V.cols();
        ndim = first_non_flip_V.rows();
        out_file << "first_non_flip_V " << nv << " " << ndim << "\n";
        for (int i = 0; i < nv; ++i)
        {
            for (int j = 0; j < ndim; ++j)
            {
                out_file << first_non_flip_V(j, i) << " ";
            }
        }
        out_file << std::endl;

        // first_locally_injective_iter
        out_file << "first_locally_injective_iter " << 1 << " " << 1 << "\n";
        out_file << first_locally_injective_iteration << " ";
        out_file << std::endl;

        /// first_locally_injective_V
        nv = first_locally_injective_V.cols();
        ndim = first_locally_injective_V.rows();
        out_file << "first_locally_injective_V " << nv << " " << ndim << "\n";
        for (int i = 0; i < nv; ++i)
        {
            for (int j = 0; j < ndim; ++j)
            {
                out_file << first_locally_injective_V(j, i) << " ";
            }
        }
        out_file << std::endl;


        if (record_vert)
        {
            auto n_record = vertRecord.size();
            nv = vertRecord[0].cols();
            ndim = vertRecord[0].rows();
            out_file << "vert " << n_record << " " << nv << " " << ndim << std::endl;
            for (const auto& vert : vertRecord) {
                for (int i = 0; i < nv; ++i) {
                    for (int j = 0; j < ndim; ++j) {
                        out_file << vert(j, i) << " ";
                    }
                }
            }
            out_file << std::endl;
        }

        if (record_energy)
        {
            unsigned n_record = energyRecord.size();
            out_file << "energy " << n_record << std::endl;
            for (int i = 0; i < n_record; ++i)
            {
                out_file << energyRecord[i] << " ";
            }
            out_file << std::endl;
        }

        if (record_minArea)
        {
            unsigned n_record = minAreaRecord.size();
            out_file << "minArea " << n_record << std::endl;
            for (int i = 0; i < n_record; ++i)
            {
                out_file << minAreaRecord[i] << " ";
            }
            out_file << std::endl;
        }

        if (record_nb_winded_interior_vertices)
        {
            auto n_record = nb_winded_interior_vertices_Record.size();
            out_file << "nbWindVert " << n_record << std::endl;
            for (int i = 0; i < n_record; ++i)
            {
                out_file << nb_winded_interior_vertices_Record[i] << " ";
            }
            out_file << std::endl;
        }

        if (record_gradient)
        {
            auto n_record = gradRecord.size();
            ndim = gradRecord[0].size();
            out_file << "grad " << n_record << " " << " " << ndim << std::endl;
            for (const auto & grad : gradRecord) {
                for (int i = 0; i < ndim; ++i) {
                    out_file << grad(i) << " ";
                }
            }
            out_file << std::endl;
        }

        if (record_gradient_norm)
        {
            unsigned n_record = gradNormRecord.size();
            out_file << "gradNorm " << n_record << std::endl;
            for (int i = 0; i < n_record; ++i)
            {
                out_file << gradNormRecord[i] << " ";
            }
            out_file << std::endl;
        }

        if (record_init_singular_values)
        {
            auto ncol = init_singular_values.cols();
            auto nrow = init_singular_values.rows();
            out_file << "initSingularValues " << ncol << " " << nrow << "\n";
            for (int i = 0; i < ncol; ++i)
            {
                for (int j = 0; j < nrow; ++j)
                {
                    out_file << init_singular_values(j, i) << " ";
                }
            }
            out_file << std::endl;
        }

        if (record_result_singular_values)
        {
            Matrix2Xd result_singular_values;
            compute_tri_mesh_singular_values(formulation.get_scaled_rest_vertices(), formulation.get_V(), F,
                                             result_singular_values);
            auto ncol = result_singular_values.cols();
            auto nrow = result_singular_values.rows();
            out_file << "resultSingularValues " << ncol << " " << nrow << "\n";
            for (int i = 0; i < ncol; ++i)
            {
                for (int j = 0; j < nrow; ++j)
                {
                    out_file << result_singular_values(j, i) << " ";
                }
            }
            out_file << std::endl;
        }


        out_file.close();
        return true;
    }
};

// nan gradient exception
class null_Grad_Exception : public std::exception
{
    virtual const char* what() const throw()
    {
        return "nan gradient entry exception.";
    }
} null_grad_exception;

//
double objective_func(const std::vector<double>& x, std::vector<double>& grad, void* my_func_data)
{
    auto* data = (Optimization_Data*)my_func_data;

    if (data->solutionFound) {
        if (!grad.empty()) {
            for (double & i : grad) {
                i = 0;
            }
        }
        return data->lastFunctionValue;
    }

    // compute energy/gradient
    double energy;
    VectorXd g_vec;

    if (grad.empty()) {  // only energy is required
        energy = -1;
        data->iteration_count += 1;
        data->record();
        // custom stop criterion
        if (data->stopQ())
        {
            data->custom_criteria_met = true;
            data->solutionFound = true;
        }
        // max iter criterion
        if (data->iteration_count >= data->max_iterations) {
            data->solutionFound = true;
        }
        //test
//        data->nb_feval += 1;
    }
    else { // gradient is required
        // convert x to Eigen::VectorXd
        VectorXd x_vec(x.size());
        for (int i = 0; i < x.size(); ++i) {
            x_vec(i) = x[i];
        }
        //
        energy = data->formulation.compute_energy_with_gradient(x_vec, g_vec);
        for (int i = 0; i < g_vec.size(); ++i) {
            if (isnan(g_vec(i))) {
                std::cout << "g_vec(" << i << ") is nan." << std::endl;
                throw null_grad_exception;
            }
            grad[i] = g_vec(i);
        }
        //test
        data->nb_geval += 1;
        // record energy
        data->lastFunctionValue = energy;
        // record gradient
        data->lastGradient = g_vec;
    }

    return energy;
}


//
int main(int argc, char const* argv[])
{
    const char* dataFile = (argc > 1) ? argv[1] : "./input";
    const char* optFile = (argc > 2) ? argv[2] : "./solver_options";
    const char* resFile = (argc > 3) ? argv[3] : "./result";

    //import data
    MatrixXd restV;
    Matrix2Xd initV;
    Matrix3Xi F;
    VectorXi handles;

    bool succeed = importData(dataFile, restV, initV, F, handles);
    if (!succeed) {
        std::cout << "usage: sTGC_QN [inputFile] [solverOptionsFile] [resultFile]" << std::endl;
        return -1;
    }

    // normalize mesh to have unit content
    /*double init_total_signed_content = compute_total_signed_mesh_area(initV, F);
    if (init_total_signed_content > 0) {
        double scale = sqrt(1. / init_total_signed_content);
        initV *= scale;
        restV *= scale;
    }*/

    //import options
    SolverOptionManager options(optFile);
    //    std::cout << "--- options ---" << std::endl;
    //    options.printOptions();


    //init
    Optimization_Data data(restV, initV, F, handles, options.form,options.alpha,
                           options.lambda1, options.lambda2, options.k,
                           options.scale_rest_mesh, options.subtract_total_signed_area,
                           options.aspect_ratio_threshold);

    auto nv = restV.cols();
    auto nfree = nv - handles.size();
    auto embed_dim = initV.rows();
    unsigned problem_dim = nfree * embed_dim;

    //set algorithm
    nlopt::algorithm algo = nlopt::LD_LBFGS;
    if (options.algorithm == "LD_LBFGS")
    {
        algo = nlopt::LD_LBFGS;
    }

    nlopt::opt opt(algo, problem_dim);

    //set stop criteria
    opt.set_ftol_abs(options.ftol_abs);
    opt.set_ftol_rel(options.ftol_rel);
    opt.set_xtol_abs(options.xtol_abs);
    opt.set_xtol_rel(options.xtol_rel);
//    opt.set_maxeval(options.maxeval);
    // bug of NLopt: setting maxeval to -1 other than a positive integer will yield different optimization process
    //opt.set_maxeval(-1);
    // instead, we set maxeval to a number much larger than the user-specified max number of iterations
    opt.set_maxeval(options.maxeval * 100);

    //pass relevant options to Optimization_Data
    data.stopCode = options.stopCode;
    data.set_record_flags(options.record);
    data.max_iterations = options.maxeval;

    //
    opt.set_min_objective(objective_func, &data);

    std::vector<double> x(data.x0.size());
    for (int i = 0; i < data.x0.size(); ++i) {
        x[i] = data.x0(i);
    }
    double minf;

    //optimize
    try {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        nlopt::result result = opt.optimize(x, minf);
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Time difference: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " [microseconds]" << std::endl;
        std::cout << "result: ";
        switch (result) {
            case nlopt::SUCCESS:
                std::cout << "SUCCESS" << std::endl;
                break;
            case nlopt::STOPVAL_REACHED:
                std::cout << "STOPVAL_REACHED" << std::endl;
                break;
            case nlopt::FTOL_REACHED:
                std::cout << "FTOL_REACHED" << std::endl;
                break;
            case nlopt::XTOL_REACHED:
                std::cout << "XTOL_REACHED" << std::endl;
                break;
            case nlopt::MAXEVAL_REACHED:
                std::cout << "MAXEVAL_REACHED" << std::endl;
                break;
            case nlopt::MAXTIME_REACHED:
                std::cout << "MAXTIME_REACHED" << std::endl;
                break;
            default:
                std::cout << "unexpected return code!" << std::endl;
                break;
        }
        std::cout << "met custom stop criteria (" << options.stopCode << "): ";
        if (data.custom_criteria_met) std::cout << "yes" << std::endl;
        else std::cout << "no" << std::endl;
        //
        std::cout << data.iteration_count << " iterations" << std::endl;
//        std::cout << data.nb_feval << " pure function evaluations" << std::endl;
        std::cout << data.nb_geval << " function/gradient evaluations" << std::endl;
    }
    catch (std::exception& e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }

    //export result
    data.exportResult(resFile);

    return 0;
}
