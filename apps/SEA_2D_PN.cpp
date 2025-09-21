//
// Created by Charles Du on 5/3/21.
// Smooth Excess Area (SEA) Projected-Newton solver for triangle meshes.
// Caveat: This is an experimental PN solver for SEA. The Hessian is approximated using TLC's Hessian.
// Experimental results show that this PN solver is only effective when the mesh boundary is fixed.
// Optimizing Global Injectivity for Constrained Parameterization, SIGGRAPH Asia 2021.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <algorithm>

#include <chrono>

#include "SEA_2D_Formulation.h"
#include "optimization_util.h"

#include <Eigen/CholmodSupport>

using namespace std;

typedef Eigen::CholmodSupernodalLLT<SpMat> CholmodSolver;

// export Eigen Matrix
bool exportMatrix(const std::string& filename, const MatrixXd &mat) {
    std::ofstream out_file(filename);
    if (!out_file.is_open()) {
        std::cerr << "Failed to open " << filename << "!" << std::endl;
        return false;
    }

    //precision of output
    typedef std::numeric_limits<double> dbl;
    out_file.precision(dbl::max_digits10);
    //
    for (auto i = 0; i < mat.cols(); ++i) {
        for (auto j = 0; j < mat.rows(); ++j) {
            out_file << mat(j, i) << " ";
        }
        out_file << std::endl;
    }
    //
    out_file.close();
    return true;
}

// solver options
class SolverOptionManager {
public:
    //default options
    SolverOptionManager() :
            form("Tutte"), alphaRatio(-1), alpha(1e-6), theta(0.1),
            ftol_abs(1e-8), ftol_rel(1e-8), xtol_abs(1e-8), xtol_rel(1e-8), gtol_abs(1e-8),
            maxeval(10000), line_search_gamma(1e-4),
            algorithm("ProjectedNewton"), stopCode("globally_injective"),
            record(),
            save_vert(false) {};

    //import options from file
    SolverOptionManager(const char *option_filename, const char *result_filename) :
            form("Tutte"), alphaRatio(-1), alpha(1e-6), theta(0.1),
            ftol_abs(1e-8), ftol_rel(1e-8), xtol_abs(1e-8), xtol_rel(1e-8), gtol_abs(1e-8),
            maxeval(10000), line_search_gamma(1e-4),
            algorithm("ProjectedNewton"), stopCode("globally_injective"),
            record(),
            save_vert(false), resFile(result_filename) {
        if (!importOptions(option_filename)) {
            std::cout << "SolverOptionManager Warn: default options are used." << std::endl;
        }
    };

    ~SolverOptionManager() = default;

    // energy formulation options
    std::string form;
    double alphaRatio;
    double alpha;
    double theta;

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

    //record values
    std::vector<std::string> record;

    //save values
    std::string resFile;
    bool save_vert;


    void printOptions() {
        std::cout << "form:\t" << form << "\n";
        std::cout << "alphaRatio:\t" << alphaRatio << "\n";
        std::cout << "alpha:\t" << alpha << "\n";
        std::cout << "theta:\t" << theta << "\n";
        std::cout << "ftol_abs:\t" << ftol_abs << "\n";
        std::cout << "ftol_rel:\t" << ftol_rel << "\n";
        std::cout << "xtol_abs:\t" << xtol_abs << "\n";
        std::cout << "xtol_rel:\t" << xtol_rel << "\n";
        std::cout << "gtol_abs:\t" << gtol_abs << "\n";
        std::cout << "maxeval:\t" << maxeval << "\n";
        std::cout << "line_search_gamma:\t" << line_search_gamma << "\n";
        std::cout << "algorithm:\t" << algorithm << "\n";
        std::cout << "stopCode:\t" << stopCode << "\n";
        std::cout << "record:  \t" <<  "{ ";
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

    bool importOptions(const char *filename) {
        //open the data file
        std::ifstream in_file(filename);

        if (!in_file.is_open()) {
            std::cerr << "Failed to open " << filename << "!" << std::endl;
            return false;
        }

        //read the file
        std::string normal = "normal";
        std::string optName;
        while (true) {
            in_file >> optName;
            if (optName != "form") {
                normal = "form";
                break;
            }
            in_file >> form;

            in_file >> optName;
            if (optName != "alphaRatio") {
                normal = "alphaRatio";
                break;
            }
            in_file >> alphaRatio;

            in_file >> optName;
            if (optName != "alpha") {
                normal = "alpha";
                break;
            }
            in_file >> alpha;

            in_file >> optName;
            if (optName != "theta") {
                normal = "theta";
                break;
            }
            in_file >> theta;

            in_file >> optName;
            if (optName != "ftol_abs") {
                normal = "ftol_abs";
                break;
            }
            in_file >> ftol_abs;

            in_file >> optName;
            if (optName != "ftol_rel") {
                normal = "ftol_rel";
                break;
            }
            in_file >> ftol_rel;

            in_file >> optName;
            if (optName != "xtol_abs") {
                normal = "xtol_abs";
                break;
            }
            in_file >> xtol_abs;

            in_file >> optName;
            if (optName != "xtol_rel") {
                normal = "xtol_rel";
                break;
            }
            in_file >> xtol_rel;

            in_file >> optName;
            if (optName != "gtol_abs") {
                normal = "gtol_abs";
                break;
            }
            in_file >> gtol_abs;

            in_file >> optName;
            if (optName != "algorithm") {
                normal = "algorithm";
                break;
            }
            in_file >> algorithm;

            in_file >> optName;
            if (optName != "maxeval") {
                normal = "maxeval";
                break;
            }
            in_file >> maxeval;

            in_file >> optName;
            if (optName != "lineSearchGamma") {
                normal = "lineSearchGamma";
                break;
            }
            in_file >> line_search_gamma;

            in_file >> optName;
            if (optName != "stopCode") {
                normal = "stopCode";
                break;
            }
            in_file >> stopCode;

            // record values
            in_file >> optName;
            if (optName != "record") {
                normal = "record";
                break;
            }

            in_file >> optName;
            if (optName != "vert") {
                normal = "vert";
                break;
            }
            int selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("vert");
            }

            in_file >> optName;
            if (optName != "energy") {
                normal = "energy";
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("energy");
            }

            in_file >> optName;
            if (optName != "minArea") {
                normal = "minArea";
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("minArea");
            }

            in_file >> optName;
            if (optName != "nbWindVert") {
                normal = "nbWindVert";
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("nbWindVert");
            }

            in_file >> optName;
            if (optName != "grad") {
                normal = "grad";
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("grad");
            }

            in_file >> optName;
            if (optName != "gradNorm") {
                normal = "gradNorm";
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("gradNorm");
            }

            in_file >> optName;
            if (optName != "searchDirection") {
                normal = "searchDirection";
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("searchDirection");
            }

            in_file >> optName;
            if (optName != "searchNorm") {
                normal = "searchNorm";
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("searchNorm");
            }

            in_file >> optName;
            if (optName != "stepSize") {
                normal = "stepSize";
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("stepSize");
            }

            in_file >> optName;
            if (optName != "stepNorm") {
                normal = "stepNorm";
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("stepNorm");
            }

            // save values
            in_file >> optName;
            if (optName != "save") {
                normal = "save";
                break;
            }

            in_file >> optName;
            if (optName != "vert") {
                normal = "vert";
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

        if (normal != "normal") {
            std::cout << "Err:(" << normal << ") fail to import options from file. Check file format." << std::endl;
            return false;
        }

        return true;
    }

};

class Optimization_Data
{
public:
    Optimization_Data(const MatrixXd  &restV,
                      const Matrix2Xd &initV,
                      const Matrix3Xi &restF,
                      const VectorXi &handles,
                      const std::string &form,
                      double alphaRatio,
                      double alpha,
                      double theta) :
            F(restF), solutionFound(false), locally_injective_Found(false),
            first_locally_injective_iteration(-1), iteration_count(0),
            first_locally_injective_V(initV),
            lastFunctionValue(HUGE_VAL), stopCode("none"),
            nb_feval(0), nb_geval(0), stop_type("unknown"),
            record_vert(false), record_energy(false), record_minArea(false),
            record_nb_winded_interior_vertices(false), record_gradient(false),
            record_gradient_norm(false), record_searchDirection(false), record_searchNorm(false),
            record_stepNorm(false), record_stepSize(false),
            vertRecord(0), energyRecord(0), minAreaRecord(0), nb_winded_interior_vertices_Record(),
            gradRecord(), gradNormRecord(), searchDirectionRecord(), searchNormRecord(),
            stepNormRecord(), stepSizeRecord(),
            save_vert(false),
            formulation(restV,initV,restF,handles,form,alphaRatio,alpha,theta),
            is_boundary_vertex(initV.cols(), false)
    {
        x0 = formulation.get_x0();

        // find boundary vertex
        const auto& boundary_edges = formulation.get_boundary_edges();
        for (const auto & e : boundary_edges) {
            is_boundary_vertex[e.first] = true;
            is_boundary_vertex[e.second] = true;
        }
    };

    ~ Optimization_Data() = default;

    SEA_2D_Formulation formulation;

//    Matrix2Xd V;
//    VectorXi freeI;
    Matrix3Xi F;
//    Matrix3Xd restD;
    VectorXd x0;

    std::vector<bool> is_boundary_vertex;

    bool solutionFound;
    double lastFunctionValue;

    bool locally_injective_Found;
    int  iteration_count;
    int  first_locally_injective_iteration;
    Matrix2Xd first_locally_injective_V;

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
    bool record_searchDirection;
    bool record_searchNorm;
    bool record_stepSize;
    bool record_stepNorm;  // ||x_next - x||
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

    // save data
    bool save_vert;

    // for what reason did the solver stop
    std::string stop_type;


    // record information we cared about
    void set_record_flags(const std::vector<std::string> &record)
    {
        for (const auto & i : record)
        {
            if (i == "vert")    record_vert = true;
            if (i == "energy")  record_energy = true;
            if (i == "minArea") record_minArea = true;
            if (i == "nbWindVert") record_nb_winded_interior_vertices = true;
            if (i == "grad") record_gradient = true;
            if (i == "gradNorm") record_gradient_norm = true;
            if (i == "searchDirection") record_searchDirection = true;
            if (i == "searchNorm") record_searchNorm = true;
            if (i == "stepSize") record_stepSize = true;
            if (i == "stepNorm") record_stepNorm = true;
        }
    }

    void record()
    {
        if (record_vert) vertRecord.push_back(formulation.get_V());
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
        } else if (stopCode == "locally_injective") {
            return stopQ_locally_injective();
        } else if (stopCode == "globally_injective") {
            return stopQ_globally_injective();
        } else {
            //default
            return false;
        }
    }

    bool stopQ_no_flip_degenerate()
    {
        bool no_flip_degenerate = true;

        if (record_minArea && !minAreaRecord.empty()) {
            if (minAreaRecord.back() <= 0) no_flip_degenerate = false;
        } else {
            double minA = compute_min_signed_mesh_area(formulation.get_V(), F);
            if (minA <= 0) no_flip_degenerate = false;
        }

        return no_flip_degenerate;
    }

    bool stopQ_locally_injective() {
        if (!stopQ_no_flip_degenerate()) {
            return false;
        } else {
            // check if there is any over-winding interior vertices
            if (record_nb_winded_interior_vertices && !nb_winded_interior_vertices_Record.empty()) {
                if (nb_winded_interior_vertices_Record.back() == 0) {
                    if (!locally_injective_Found) {
                        locally_injective_Found = true;
                        first_locally_injective_iteration = iteration_count;
                        first_locally_injective_V = formulation.get_V();
                    }
                    return true;
                }
                else return false;
            } else {
                std::vector<size_t> winded_vertices;
                compute_winded_interior_vertices(formulation.get_V(), F, is_boundary_vertex, winded_vertices);
                if (winded_vertices.empty()) {
                    if (!locally_injective_Found) {
                        locally_injective_Found = true;
                        first_locally_injective_iteration = iteration_count;
                        first_locally_injective_V = formulation.get_V();
                    }
                    return true;
                }
                else return false;
            }
        }
    }

    bool stopQ_globally_injective() {
        if (locally_injective_Found) {
            // locally injective mesh is already found
            if (!stopQ_no_flip_degenerate()) {
                return false;
            } else {
                return !(formulation.has_found_intersection);
            }
        } else {
            if (!stopQ_locally_injective()) {
                return false;
            } else {
                // record the first locally injective iteration
                locally_injective_Found = true;
                first_locally_injective_iteration = iteration_count;
                first_locally_injective_V = formulation.get_V();

                // check global injectivity:
                // if no arc-arc intersection is found in the last energy/gradient evaluation,
                // then global injectivity is surely achieved
                return !(formulation.has_found_intersection);
            }
        }
    }

    //export result
    bool exportResult(const char* filename)
    {
        std::ofstream out_file(filename);
        if (! out_file.is_open()) {
            std::cerr << "Failed to open " << filename << "!" << std::endl;
            return false;
        }

        //precision of output
        typedef std::numeric_limits< double > dbl;
        out_file.precision(dbl::max_digits10);

        //write the file

        /// resV
        const auto &V = formulation.get_V();
        auto nv = V.cols();
        auto ndim = V.rows();
        out_file << "resV " << nv << " " << ndim << "\n";
        for (int i = 0; i < nv; ++i)
        {
            for (int j = 0; j < ndim; ++j)
            {
                out_file << V(j,i) << " ";
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
                out_file << first_locally_injective_V(j,i) << " ";
            }
        }
        out_file << std::endl;


        if (record_vert)
        {
            auto n_record = vertRecord.size();
            nv = vertRecord[0].cols();
            ndim = vertRecord[0].rows();
            out_file << "vert " << n_record << " " << nv << " " << ndim << std::endl;
            for (const auto & vert : vertRecord) {
                for (int i = 0; i < nv; ++i) {
                    for (int j = 0; j < ndim; ++j) {
                        out_file << vert(j,i) << " ";
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

        if (record_gradient_norm) {
            auto n_record = gradNormRecord.size();
            out_file << "gradNorm " << n_record << std::endl;
            for (auto i = 0; i < n_record; ++i) {
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
    energy = data.formulation.compute_energy_with_gradient_approxProjectedHessian(x,energyList,grad,mat);

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
        cout << "iter 0: decomposition failed" << endl;
        data.stop_type = "solver.factorize fail";
        return;
    }
    VectorXd p = solver.solve(-grad);
    if (solver.info() != Success) {
        cout << "iter 0: solving failed" << endl;
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
    if (data.record_stepNorm) data.stepNormRecord.push_back(p.norm() * step_size);
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
        energy = data.formulation.compute_energy_with_gradient_approxProjectedHessian(x,energyList,grad,mat);

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
            cout << "iter " << i << ": decomposition failed" << endl;
            data.stop_type = "solver.factorize fail";
            return;
        }
        p = solver.solve(-grad);
        if (solver.info() != Success) {
            cout << "iter " << i << ": solving failed" << endl;
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


int main(int argc, char const *argv[]) {
    std::cout << "Experimental PN solver for SEA, not recommended for practical use" << std::endl;

    const char *dataFile = (argc > 1) ? argv[1] : "./input";
    const char *optFile = (argc > 2) ? argv[2] : "./solver_options";
    const char *resFile = (argc > 3) ? argv[3] : "./result";

    //import data
    MatrixXd restV;
    Matrix2Xd initV;
    Matrix3Xi F;
    VectorXi handles;


    bool succeed = importData(dataFile,restV,initV,F,handles);
    if (!succeed) {
        std::cout << "usage: SEA_2D_PN [inputFile] [solverOptionsFile] [resultFile]" << std::endl;
        return -1;
    }


    //import options
    SolverOptionManager options(optFile, resFile);

    // init
    Optimization_Data data(restV, initV, F, handles, options.form, options.alphaRatio, options.alpha, options.theta);
    VectorXd x = data.x0;

    //pass relevant options to Optimization_Data
    data.stopCode = options.stopCode;
    data.set_record_flags(options.record);
    data.save_vert = options.save_vert;

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