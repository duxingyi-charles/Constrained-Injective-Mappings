//
// Created by Charles Du on 4/30/21.
//

// 2D/3D lifted formulation
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include <chrono>

#include <nlopt.hpp>

#include "Arc_Overlap_Formulation.h"

using namespace Eigen;

bool importData(const char* filename,
                MatrixXd &restV,
                Matrix2Xd &initV,
                Matrix3Xi &F,
                VectorXi &handles
)
{
    std::ifstream in_file(filename);

    if (! in_file.is_open()) {
        std::cerr << "Failed to open " << filename << "!" << std::endl;
        return false;
    }

    //read the file
    size_t n, ndim;
    // restV
    in_file >> n >> ndim;
    restV.resize(ndim, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < ndim; ++j) {
            in_file >> restV(j,i);
        }
    }

    //initV
    in_file >> n >> ndim;
    initV.resize(ndim,n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < ndim; ++j) {
            in_file >> initV(j,i);
        }

    }

    //F
    size_t simplexSize;
    in_file >> n >> simplexSize;
    F.resize(simplexSize, n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < simplexSize; ++j)
        {
            in_file >> F(j,i);
        }
    }

    //handles
    in_file >> n;
    handles.resize(n);
    for (int i = 0; i < n; ++i)
    {
        in_file >> handles(i);
    }

    in_file.close();

    return true;
}

// solver options
class NloptOptionManager
{
public:
    //default options
    NloptOptionManager():
            form("Tutte"), alphaRatio(1e-6), alpha(-1), theta(0.1),
            ftol_abs(1e-8), ftol_rel(1e-8), xtol_abs(1e-8), xtol_rel(1e-8),
            maxeval(10000), algorithm("LBFGS"), stopCode("all_good"), record()
    {};
    //import options from file
    explicit NloptOptionManager(const char* filename):
            form("Tutte"), alphaRatio(1e-6), alpha(-1), theta(0.1),
            ftol_abs(1e-8), ftol_rel(1e-8), xtol_abs(1e-8), xtol_rel(1e-8),
            maxeval(10000), algorithm("LBFGS"), stopCode("all_good"), record()
    {
        if (!importOptions(filename))
        {
            std::cout << "NloptOptionManager Warn: default options are used." << std::endl;
        }
    };

    ~NloptOptionManager() = default;

    // energy formulation options
    std::string form;
    double alphaRatio;
    double alpha;
    double theta; //arc angle

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
        std::cout << "alphaRatio:\t" << alphaRatio << "\n";
        std::cout << "alpha:\t"    <<  alpha    << "\n";
        std::cout << "theta:\t"    <<  theta    << "\n";
        std::cout << "ftol_abs:\t" <<  ftol_abs << "\n";
        std::cout << "ftol_rel:\t" <<  ftol_rel << "\n";
        std::cout << "xtol_abs:\t" <<  xtol_abs << "\n";
        std::cout << "xtol_rel:\t" <<  xtol_rel << "\n";
        std::cout << "maxeval:\t" <<  maxeval << "\n";
        std::cout << "algorithm:\t" <<  algorithm << "\n";
        std::cout << "stopCode:\t" <<  stopCode << "\n";
        std::cout << "record:  \t" <<  "{ ";
        for (std::vector<std::string>::iterator i = record.begin(); i != record.end(); ++i)
        {
            std::cout << *i << " ";
        }
        std::cout << "}" << std::endl;
    }

    bool importOptions(const char* filename)
    {
        //open the data file
        std::ifstream in_file(filename);

        if (! in_file.is_open())
        {
            std::cerr << "Failed to open solver option file:" << filename << "!" << std::endl;
            return false;
        }

        //read the file
        unsigned abnormal = 0;
        std::string optName;
        while(true)
        {
            in_file >> optName;
            if (optName != "form")
            {
                abnormal = 1;
                break;
            }
            in_file >> form;

            in_file >> optName;
            if (optName != "alphaRatio")
            {
                abnormal = 2;
                break;
            }
            in_file >> alphaRatio;

            in_file >> optName;
            if (optName != "alpha")
            {
                abnormal = 3;
                break;
            }
            in_file >> alpha;

            in_file >> optName;
            if (optName != "theta") {
                abnormal = 15;
                break;
            }
            in_file >> theta;

            in_file >> optName;
            if (optName != "ftol_abs")
            {
                abnormal = 4;
                break;
            }
            in_file >> ftol_abs;

            in_file >> optName;
            if (optName != "ftol_rel")
            {
                abnormal = 5;
                break;
            }
            in_file >> ftol_rel;

            in_file >> optName;
            if (optName != "xtol_abs")
            {
                abnormal = 6;
                break;
            }
            in_file >> xtol_abs;

            in_file >> optName;
            if (optName != "xtol_rel")
            {
                abnormal = 7;
                break;
            }
            in_file >> xtol_rel;

            in_file >> optName;
            if (optName != "algorithm")
            {
                abnormal = 8;
                break;
            }
            in_file >> algorithm;

            in_file >> optName;
            if (optName != "maxeval")
            {
                abnormal = 9;
                break;
            }
            in_file >> maxeval;

            in_file >> optName;
            if (optName != "stopCode")
            {
                abnormal = 10;
                break;
            }
            in_file >> stopCode;

            // record values
            in_file >> optName;
            if (optName != "record")
            {
                abnormal = 11;
                break;
            }

            in_file >> optName;
            if (optName != "vert") {
                abnormal = 12;
                break;
            }
            int selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("vert");
            }

            in_file >> optName;
            if (optName != "energy") {
                abnormal = 13;
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("energy");
            }

            in_file >> optName;
            if (optName != "minArea") {
                abnormal = 14;
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("minArea");
            }

            break;
        }

        in_file.close();

        if (abnormal!=0)
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
    Optimization_Data(const MatrixXd  &restV,
                      const Matrix2Xd &initV,
                      const Matrix3Xi &restF,
                      const VectorXi &handles,
                      const std::string &form,
                      double alphaRatio,
                      double alpha,
                      double theta) :
                      F(restF), solutionFound(false),
            lastFunctionValue(HUGE_VAL), stopCode("none"),
            nb_feval(0),nb_geval(0),
            record_vert(false), record_energy(false), record_minArea(false),
            vertRecord(0),energyRecord(0), minAreaRecord(0),
            formulation(restV,initV,restF,handles,form,alphaRatio,alpha,theta)
    {
        x0 = formulation.get_x0();
    };

    ~ Optimization_Data() = default;

    Arc_Overlap_Formulation formulation;

//    Matrix2Xd V;
//    VectorXi freeI;
    Matrix3Xi F;
//    Matrix3Xd restD;
    VectorXd x0;

    bool solutionFound;
    double lastFunctionValue;

    int nb_feval;
    int nb_geval;

    //from options
    std::string stopCode;
    //record data
    bool record_vert;
    bool record_energy;
    bool record_minArea;
    std::vector<Matrix2Xd> vertRecord;
    std::vector<double> minAreaRecord;
    std::vector<double> energyRecord;


    // record information we cared about
    void set_record_flags(const std::vector<std::string> &record)
    {
        for (std::vector<std::string>::const_iterator i = record.begin(); i != record.end(); ++i)
        {
            if (*i == "vert")    record_vert = true;
            if (*i == "energy")  record_energy = true;
            if (*i == "minArea") record_minArea = true;
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
    }

    //custom stop criteria
    bool stopQ()
    {
        if (stopCode == "all_good") return stopQ_all_good();
        //default
        return false;
    }

    bool stopQ_all_good()
    {
        bool good = true;

        if (record_minArea && !minAreaRecord.empty()) {
            if (minAreaRecord.back() <= 0) good = false;
        } else {
            double minA = compute_min_signed_mesh_area(formulation.get_V(), F);
            if (minA <= 0) good = false;
        }

        return good;
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

        out_file.close();
        return true;
    }
};


//
double ojbective_func(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    Optimization_Data *data = (Optimization_Data *) my_func_data;

    if (data->solutionFound) {
        if (!grad.empty()) {
            for (int i = 0; i < grad.size(); ++i) {
                grad[i] = 0;
            }
        }
        return data->lastFunctionValue;
    }

    // convert x to Eigen::VectorXd
    VectorXd x_vec(x.size());
    for (int i = 0; i < x.size(); ++i) {
        x_vec(i) = x[i];
    }

    // compute energy/gradient
    double energy;
    VectorXd g_vec;
    if (grad.empty()) {  // only energy is required
        energy = data->formulation.compute_energy(x_vec);
        //test
        data->nb_feval += 1;
    } else { // gradient is required
        energy = data->formulation.compute_energy_with_gradient(x_vec, g_vec);
        for (int i = 0; i < g_vec.size(); ++i) {
            grad[i] = g_vec(i);
        }
        //test
        data->nb_geval += 1;
    }

    // custom stop criterion
    if (data->stopQ())
    {
        data->solutionFound = true;
    }

    //record information
    data->lastFunctionValue = energy;
    data->record();
    return energy;
}


//
int main(int argc, char const *argv[])
{
    const char* dataFile = (argc > 1) ? argv[1] : "./input";
    const char* optFile  = (argc > 2) ? argv[2] : "./solver_options";
    const char* resFile  = (argc > 3) ? argv[3] : "./result";

    //import data
    MatrixXd restV;
    Matrix2Xd initV;
    Matrix3Xi F;
    VectorXi handles;

    bool succeed = importData(dataFile,restV,initV,F,handles);
    if (!succeed) {
        std::cout << "usage: findInjective [inputFile] [solverOptionsFile] [resultFile]" << std::endl;
        return -1;
    }

    //import options
    NloptOptionManager options(optFile);
    //std::cout << "--- options ---" << std::endl;
    //options.printOptions();


    //init
    Optimization_Data data(restV, initV, F, handles, options.form, options.alphaRatio, options.alpha, options.theta);

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
    opt.set_maxeval(options.maxeval);

    //pass relevant options to LiftedData
    data.stopCode = options.stopCode;
    data.set_record_flags(options.record);

    //
    opt.set_min_objective(ojbective_func, &data);

    std::vector<double> x(data.x0.size());
    for (int i = 0; i < data.x0.size(); ++i) {
        x[i] = data.x0(i);
    }
    double minf;

    //optimize
    try{
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        nlopt::result result = opt.optimize(x, minf);
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Time difference: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " [microseconds]" << std::endl;
        std::cout << "result: ";
        switch(result) {
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
        if (data.solutionFound) std::cout << "yes" << std::endl;
        else std::cout << "no" << std::endl;
        //
        //std::cout << data.nb_feval << " function evalations, ";
        //std::cout << data.nb_geval << " gradient evalations." << std::endl;

    }
    catch(std::exception &e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }

    //export result
    data.exportResult(resFile);

    return 0;
}
