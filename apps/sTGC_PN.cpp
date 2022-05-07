//
// Created by Charles Du on 3/30/22.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <chrono>

#include "sTGC_Formulation.h"
#include "optimization_util.h"
#include "deformation_gradient_util.h"

#include <Eigen/CholmodSupport>

using namespace Eigen;
typedef Eigen::CholmodSupernodalLLT<SpMat> CholmodSolver;

bool export_sparse_matrix(const char* filename, const SpMat& mat)
{
    std::ofstream out_file(filename);
    if (!out_file.is_open()) {
        std::cerr << "Failed to open " << filename << "!" << std::endl;
        return false;
    }

    //precision of output
    typedef std::numeric_limits< double > dbl;
    out_file.precision(dbl::max_digits10);

    out_file << mat.rows() << " " << mat.cols() << std::endl;
    for (int k = 0; k < mat.outerSize(); ++k) {
        for (SpMat::InnerIterator it(mat, k); it; ++it) {
            out_file << it.row() << " " << k << " " << it.value() << std::endl;
        }
    }

    out_file.close();
    return true;
}


// solver options
class SolverOptionManager
{
public:
    //default options
    SolverOptionManager() :
            form("Tutte"), alpha(1),
            lambda1(0.5), lambda2(0.), k(1.), // lambda1 + k * lambda2 = 1/2
            scale_rest_mesh(false), subtract_total_signed_area(true),
            ftol_abs(1e-8), ftol_rel(1e-8), xtol_abs(1e-8), xtol_rel(1e-8), gtol_abs(1e-8),
            maxeval(10000), line_search_gamma(1e-4), aspect_ratio_threshold(0),
            algorithm("Projected_Newton"), stopCode("no_flip_degenerate"),
            record(),
            save_vert(false) {};

    //import options from file
    explicit SolverOptionManager(const char* filename, const char *result_filename) :
            form("Tutte"), alpha(1),
            lambda1(0.5), lambda2(0.), k(1.),
            scale_rest_mesh(false), subtract_total_signed_area(true),
            ftol_abs(1e-8), ftol_rel(1e-8), xtol_abs(1e-8), xtol_rel(1e-8), gtol_abs(1e-8),
            maxeval(10000), line_search_gamma(1e-4), aspect_ratio_threshold(0),
            algorithm("Projected_Newton"), stopCode("no_flip_degenerate"),
            record(),
            save_vert(false), resFile(result_filename)
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
    double gtol_abs;
    int maxeval;
    double line_search_gamma;
    std::string algorithm;
    std::string stopCode;

    // record values
    std::vector<std::string> record;

    //save values
    std::string resFile;
    bool save_vert;

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
        std::cout << "line_search_gamma:\t" << line_search_gamma << "\n";
        std::cout << "algorithm:\t" << algorithm << "\n";
        std::cout << "stopCode:\t" << stopCode << "\n";
        std::cout << "record:  \t" << "{ ";
        for (auto & i : record)
        {
            std::cout << i << " ";
        }
        std::cout << "}" << std::endl;
        //
        std::cout << "result file: \t" << resFile << std::endl;
        std::cout << "save:  \t" << "{ ";
        if (save_vert) std::cout << "vert ";
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
            if (optName != "gtol_abs")
            {
                abnormal = "gtol_abs";
                break;
            }
            in_file >> gtol_abs;

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
            if (optName != "lineSearchGamma") {
                abnormal = "lineSearchGamma";
                break;
            }
            in_file >> line_search_gamma;

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
            if (optName != "searchDirection") {
                abnormal = "searchDirection";
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("searchDirection");
            }

            in_file >> optName;
            if (optName != "searchNorm") {
                abnormal = "searchNorm";
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("searchNorm");
            }

            in_file >> optName;
            if (optName != "stepSize") {
                abnormal = "stepSize";
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("stepSize");
            }

            in_file >> optName;
            if (optName != "stepNorm") {
                abnormal = "stepNorm";
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("stepNorm");
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

            // save values
            in_file >> optName;
            if (optName != "save") {
                abnormal = "save";
                break;
            }

            in_file >> optName;
            if (optName != "vert") {
                abnormal = "vert";
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                save_vert = true;
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
            first_locally_injective_iteration(-1), iteration_count(0),
            first_locally_injective_V(initV),
            non_flip_Found(false), first_non_flip_iteration(-1),
            first_non_flip_V(initV),
            lastFunctionValue(HUGE_VAL), stopCode("none"),
            record_vert(false), record_energy(false), record_gradient(false),
            record_gradient_norm(false), record_minArea(false),
            record_nb_winded_interior_vertices(false), record_searchDirection(false), record_searchNorm(false),
            record_init_singular_values(false), record_result_singular_values(false),
            record_stepNorm(false), record_stepSize(false),
            vertRecord(0), energyRecord(0), minAreaRecord(0), gradRecord(0),
            save_vert(false),
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

    bool locally_injective_Found;
    int  iteration_count;
    int  first_locally_injective_iteration;
    Matrix2Xd first_locally_injective_V;
    bool non_flip_Found;
    int first_non_flip_iteration;
    Matrix2Xd first_non_flip_V;

    //from options
    std::string stopCode;
    //record data
    bool record_vert;
    bool record_energy;
    bool record_gradient;
    bool record_gradient_norm;
    bool record_minArea;
    bool record_nb_winded_interior_vertices;
    bool record_searchDirection;
    bool record_searchNorm;
    bool record_stepSize;
    bool record_stepNorm;  // ||x_next - x||
    bool record_init_singular_values;
    bool record_result_singular_values;
    std::vector<Matrix2Xd> vertRecord;
    std::vector<double> minAreaRecord;
    std::vector<double> energyRecord;
    std::vector<int> nb_winded_interior_vertices_Record;
    std::vector<VectorXd> gradRecord;
    std::vector<double> gradNormRecord;
    std::vector<VectorXd> searchDirectionRecord;
    std::vector<double> searchNormRecord;
    std::vector<double> stepSizeRecord;
    std::vector<double> stepNormRecord;
    Matrix2Xd init_singular_values;


    // save data
    bool save_vert;

    // for what reason did the solver stop
    std::string stop_type;


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
            if (i == "searchDirection") record_searchDirection = true;
            if (i == "searchNorm") record_searchNorm = true;
            if (i == "stepSize") record_stepSize = true;
            if (i == "stepNorm") record_stepNorm = true;
            if (i == "initSingularValues") record_init_singular_values = true;
            if (i == "resultSingularValues") record_result_singular_values = true;
        }
    }

    void record()
    {
        if (record_vert) vertRecord.push_back(formulation.get_V());
//        if (record_gradient) gradRecord.push_back(lastGradient);
//        if (record_gradient_norm) gradNormRecord.push_back(lastGradient.norm());
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
        iteration_count += 1;
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
            //todo: only check free faces
            double minA = compute_min_signed_mesh_area(formulation.get_V(), F);
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

//        if (record_gradient)
//        {
//            auto n_record = gradRecord.size();
//            ndim = gradRecord[0].size();
//            out_file << "grad " << n_record << " " << " " << ndim << std::endl;
//            for (const auto & grad : gradRecord) {
//                for (int i = 0; i < ndim; ++i) {
//                    out_file << grad(i) << " ";
//                }
//            }
//            out_file << std::endl;
//        }
        if (record_gradient) {
            auto n_record = gradRecord.size();
            auto n_free = formulation.get_freeI().size();
            ndim = 2;
            out_file << "grad " << n_record << " " << n_free << " " << ndim << std::endl;
            for (auto i = 0; i < n_record; ++i) {
                for (auto j = 0; j < n_free; ++j) {
                    for (auto k = 0; k < ndim; ++k) {
                        out_file << gradRecord[i](j * ndim + k) << " ";
                    }
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

        if (record_searchDirection) {
            auto n_record = searchDirectionRecord.size();
            auto n_free = formulation.get_freeI().size();
            out_file << "searchDirection " << n_record << " " << n_free << " " << ndim << std::endl;
            for (auto i = 0; i < n_record; ++i) {
                for (auto j = 0; j < n_free; ++j) {
                    for (auto k = 0; k < ndim; ++k) {
                        out_file << searchDirectionRecord[i](j * ndim + k) << " ";
                    }
                }
            }
            out_file << std::endl;
        }

        if (record_searchNorm) {
            auto n_record = searchNormRecord.size();
            out_file << "searchNorm " << n_record << std::endl;
            for (auto i = 0; i < n_record; ++i) {
                out_file << searchNormRecord[i] << " ";
            }
            out_file << std::endl;
        }

        if (record_stepSize) {
            auto n_record = stepSizeRecord.size();
            out_file << "stepSize " << n_record << std::endl;
            for (auto i = 0; i < n_record; ++i) {
                out_file << stepSizeRecord[i] << " ";
            }
            out_file << std::endl;
        }

        if (record_stepNorm) {
            auto n_record = stepNormRecord.size();
            out_file << "stepNorm " << n_record << std::endl;
            for (auto i = 0; i < n_record; ++i) {
                out_file << stepNormRecord[i] << " ";
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

    void lineSearch(VectorXd &x, const VectorXd &p, double &step_size, const VectorXd &grad,
                    double energy, const VectorXd &energyList,
                    double &energy_next,
                    double shrink, double gamma) {
        double gp = gamma * grad.transpose() * p;
        VectorXd x_next = x + step_size * p;
        VectorXd energyList_next;
        energy_next = formulation.compute_energy(x_next, energyList_next);

        VectorXd energy_diff_list = energyList_next - energyList;
        double energy_diff = energy_diff_list.sum();

        while (energy_diff > step_size * gp) {
            step_size *= shrink;
            x_next = x + step_size * p;
            energy_next = formulation.compute_energy(x_next, energyList_next);

            energy_diff_list = energyList_next - energyList;
            energy_diff = energy_diff_list.sum();

//            std::cout << energy_diff << std::endl;
        }
        x = x_next;
    }
};


void projected_Newton(Optimization_Data &data, VectorXd &x, SolverOptionManager &options, double shrink = 0.7) {
    //handle options
    double ftol_rel = options.ftol_rel;
    double ftol_abs = options.ftol_abs;
    double xtol_rel = options.xtol_rel;
    double xtol_abs = options.xtol_abs;
    double gtol_abs = options.gtol_abs;
    int maxIter = options.maxeval;
    //
//    bool save_vert = options.save_vert;
    //handle options end

    //
    double energy;
    VectorXd energyList;
    VectorXd grad(x.size());
    SpMat mat(x.size(), x.size());

    double energy_next;

    //first iter: initialize solver
    energy = data.formulation.compute_energy_with_gradient_projectedHessian(x,energyList,grad,mat);

    // solver step monitor
    data.lastFunctionValue = energy;
    data.record();
    if (data.record_gradient) data.gradRecord.emplace_back(grad);
    if (data.record_gradient_norm) data.gradNormRecord.push_back(grad.norm());
    // custom stop criterion
    if (data.stopQ())
    {
        data.solutionFound = true;
        data.stop_type = "custom stop goal met";
        return;
    }
    // solver step monitor end

    //check gtol
    if (grad.norm() < gtol_abs) {
        data.stop_type = "gtol reached";
        return;
    }

    // initialize solver
    CholmodSolver solver;
    solver.analyzePattern(mat);
    // initialize solver end

    solver.factorize(mat);
    if (solver.info() != Success) {
        std::cout << "iter 0: decomposition failed" << std::endl;
        data.stop_type = "solver.factorize fail";
        return;
    }
    VectorXd p = solver.solve(-grad);
    if (solver.info() != Success) {
        std::cout << "iter 0: solving failed" << std::endl;
        data.stop_type = "solver.solve fail";
        return;
    }
    if (data.record_searchDirection) data.searchDirectionRecord.emplace_back(p);
    if (data.record_searchNorm) data.searchNormRecord.push_back(p.norm());

    // backtracking line search
    double step_size = 1.0;
    double x_norm = x.norm();
    data.lineSearch(x, p, step_size, grad, energy, energyList, energy_next, shrink, options.line_search_gamma);
    //
    if (data.record_stepSize) data.stepSizeRecord.push_back(step_size);
    double step_norm = p.norm() * step_size;
    if (data.record_stepNorm) data.stepNormRecord.push_back(step_norm);
    //check ftol
    if (fabs(energy_next - energy) < ftol_abs) {
        data.stop_type = "ftol_abs reached";
        return;
    }
    if (fabs((energy_next - energy) / energy) < ftol_rel) {
        data.stop_type = "ftol_rel reached";
        return;
    }
    // check xtol
    if (step_norm < xtol_abs) {
        data.stop_type = "xtol_abs reached";
        return;
    }
    if (step_norm / x_norm < xtol_rel) {
        data.stop_type = "xtol_rel reached";
        return;
    }


    for (int i = 1; i < maxIter; ++i) {
        energy = data.formulation.compute_energy_with_gradient_projectedHessian(x,energyList,grad,mat);

        // solver step monitor
        data.lastFunctionValue = energy;
        data.record();
        if (data.record_gradient) data.gradRecord.emplace_back(grad);
        if (data.record_gradient_norm) data.gradNormRecord.push_back(grad.norm());
        // custom stop criterion
        if (data.stopQ())
        {
            data.solutionFound = true;
            data.stop_type = "custom stop goal met";
            return;
        }
        // solver step monitor end

        //check gtol
        if (grad.norm() < gtol_abs) {
            data.stop_type = "gtol reached";
            return;
        }

        solver.factorize(mat);
        if (solver.info() != Success) {
            std::cout << "iter " << i << ": decomposition failed" << std::endl;
            data.stop_type = "solver.factorize fail";
            //export_sparse_matrix("D:\\research\\arcOverlap\\debug\\Hess.txt", mat);
            return;
        }
        p = solver.solve(-grad);
        if (solver.info() != Success) {
            std::cout << "iter " << i << ": solving failed" << std::endl;
            data.stop_type = "solver.solve fail";
            return;
        }
        if (data.record_searchDirection) data.searchDirectionRecord.emplace_back(p);
        if (data.record_searchNorm) data.searchNormRecord.push_back(p.norm());

        // backtracking line search
        step_size = 1.0;
        x_norm = x.norm();
        data.lineSearch(x, p, step_size, grad, energy, energyList, energy_next, shrink, options.line_search_gamma);
        //
        if (data.record_stepSize) data.stepSizeRecord.push_back(step_size);
        step_norm = p.norm() * step_size;
        if (data.record_stepNorm) data.stepNormRecord.push_back(step_norm);
        //check ftol
        if (fabs(energy_next - energy) < ftol_abs) {
            data.stop_type = "ftol_abs reached";
            return;
        }
        if (fabs((energy_next - energy) / energy) < ftol_rel) {
            data.stop_type = "ftol_rel reached";
            return;
        }
        // check xtol
        if (step_norm < xtol_abs) {
            data.stop_type = "xtol_abs reached";
            return;
        }
        if (step_norm / x_norm < xtol_rel) {
            data.stop_type = "xtol_rel reached";
            return;
        }
    }

    // reach max iterations
    data.stop_type = "max iteration reached";
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
        std::cout << "usage: sTGC_PN [inputFile] [solverOptionsFile] [resultFile]" << std::endl;
        return -1;
    }

    // normalize mesh to have unit content
    double init_total_signed_content = compute_total_signed_mesh_area(initV, F);
    if (init_total_signed_content > 0) {
        double scale = sqrt(1. / init_total_signed_content);
        initV *= scale;
        restV *= scale;
    }

    //import options
    SolverOptionManager options(optFile, resFile);
//    options.printOptions();

    //init
    Optimization_Data data(restV, initV, F, handles, options.form,options.alpha,
                           options.lambda1, options.lambda2, options.k,
                           options.scale_rest_mesh,
                           options.subtract_total_signed_area,
                           options.aspect_ratio_threshold);
    VectorXd x = data.x0;

    //pass relevant options to Optimization_Data
    data.stopCode = options.stopCode;
    data.set_record_flags(options.record);
    data.save_vert = options.save_vert;

    //test: consistency of different runs
    /*double energy1;
    VectorXd energyList1;
    VectorXd grad1(x.size());
    SpMat mat1(x.size(), x.size());
    energy1 = data.formulation.compute_energy_with_gradient_projectedHessian(x, energyList1, grad1, mat1);
    double energy2;
    VectorXd energyList2;
    VectorXd grad2(x.size());
    SpMat mat2(x.size(), x.size());
    energy2 = data.formulation.compute_energy_with_gradient_projectedHessian(x, energyList2, grad2, mat2);
    bool same_energy = true;
    for (size_t i = 0; i < energyList1.size(); i++)
    {
        if (energyList1(i) != energyList2(i)) {
            same_energy = false;
            break;
        }
    }
    std::cout << "same_energy: " << same_energy << std::endl;
    bool same_grad = true;
    for (size_t i = 0; i < grad1.size(); i++)
    {
        if (grad1(i) != grad2(i)) {
            same_grad = false;
            break;
        }
    }
    std::cout << "same_grad: " << same_grad << std::endl;
    bool same_Hess = true;
    for (int k = 0; k < mat1.outerSize(); ++k) {
        SpMat::InnerIterator it1(mat1, k);
        SpMat::InnerIterator it2(mat2, k);
        while (it1) {
            if (it1.row() != it2.row() || it1.value() != it2.value()) {
                same_Hess = false;
                break;
            }
            ++it1;
            ++it2;
        }
        if (!same_Hess) {
            break;
        }
    }
    std::cout << "same_Hess: " << same_Hess << std::endl;*/

    // test: export Hessian matrix
    /*double energy;
    VectorXd energyList;
    VectorXd grad(x.size());
    SpMat Hess(x.size(), x.size());
    energy = data.formulation.compute_energy_with_gradient_projectedHessian(x, energyList, grad, Hess);
    export_sparse_matrix("D:\\research\\arcOverlap\\debug\\Hess.txt", Hess);*/

    // test: compare gradient with finite difference
    /*VectorXd x0 = data.x0;
    VectorXd energy_list0;
    VectorXd g0;
    SpMat hess0;
    double energy0 = data.formulation.compute_energy_with_gradient_projectedHessian(x0, energy_list0, g0, hess0);
    double delta = 1e-6;
    VectorXd finite_diff(g0.size());
    for (int i = 0; i < finite_diff.size(); ++i) {
        VectorXd x1 = x0;
        x1(i) -= delta;
        VectorXd x2 = x0;
        x2(i) += delta;
        finite_diff(i) = (data.formulation.compute_energy(x2)-data.formulation.compute_energy(x1))/(2*delta);
    }
    std::cout << "gradient\t\t|\t finite diff (delta = " << delta << ")" << std::endl;
    for (int i = 0; i < g0.size(); ++i) {
        std::cout << i << ": " <<  g0(i) << "\t\t|\t" << finite_diff(i) << std::endl;
    }
    return 0;*/


    //projected newton
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    projected_Newton(data, x, options);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()
              << " [microseconds]" << std::endl;
    std::cout << "Stop: " << data.stop_type << std::endl;

    //export result
    data.exportResult(resFile);

    return 0;
}
