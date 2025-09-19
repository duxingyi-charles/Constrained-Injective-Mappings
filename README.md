# Constrained Injective Mappings

## Programs

This codebase provides various programs for computing injective mappings under different types of constraints. All programs follow a consistent `[energy]_[dimension]_[solver]` naming convention:

- **Energy formulations**:
  - **[TLC](https://duxingyi-charles.github.io/publication/lifting-simplices-to-find-injectivity/)** (Total Lifted Content): Energy for mapping optimization with prescribed boundaries
  - **[SEA](https://duxingyi-charles.github.io/publication/optimizing-global-injectivity-for-constrained-parameterization/)** (Smooth Excess Area): Energy for mapping optimization with general positional constraints
  - **sTGC**: (subtracted Total Generalized Content) Generalized variant of TLC for optimizing map distortions. 
  The 's' prefix means subtracting total signed area/volume, enabling optimization for free boundary mappings.
  **[sTLC_Iso](https://duxingyi-charles.github.io/publication/isometric-energies-for-recovering-injectivity-in-constrained-mapping/)** is a specialized variant of sTGC that encourages isometric mappings.
  - **SEA_Generalized**: Generalized variant of SEA for computing injective mappings while optimizing map distortions under general positional constraints. 
  **[SEA_Iso](https://duxingyi-charles.github.io/publication/isometric-energies-for-recovering-injectivity-in-constrained-mapping/)** is a specialized variant that encourages isometric mappings.
- **Dimension**: 2D for triangle meshes, 3D for tetrahedral meshes
- **Solver**: QN (quasi-Newton), PN (projected Newton)

### For Mappings with Prescribed Boundaries

Use these programs for mapping a mesh into a target domain with a prescribed boundary shape and want to optimize the interior mapping while keeping the boundary unchanged.

#### 2D Triangle Meshes
- **sTGC_2D_QN**: Recommended for most 2D problems. Supports both standard TLC and isometric TLC variants
- **sTGC_2D_PN**: Improved convergence rate. Recommended for hard problems where QN has difficulty converging.

#### 3D Tetrahedral Meshes
- **TLC_QN**: Standard TLC optimization for 3D meshes
- **sTGC_3D_QN**: Generalized variant that supports isometric TLC
- **sTGC_3D_PN**: Improved convergence rate for isometric TLC

### For Mappings with General Positional Constraints

Use this program when you have arbitrary positional constraints scattered throughout your mesh.

- **SEA_Generalized_2D_QN**: Handles both standard SEA and isometric SEA formulations for 2D triangle meshes

## Build


### Mac

Install [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse) (needed by PN solvers) by homebrew:

    brew install suite-sparse

Change the include/lib path of SuiteSparse in CMakeLists.txt if needed. Then, build with the following commands: 

    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make -j

### Windows

Install SuiteSparse using vcpkg or follow the instructions on the [SuiteSparse GitHub Repo](https://github.com/DrTimothyAldenDavis/SuiteSparse).

Change the include/lib path of SuiteSparse in CMakeLists.txt if needed. Then, use CMake to generate a Visual Studio project and build it.

### Add support for MaxIteration termination criterion

:bell:  **Important**:
We use [NLopt](https://nlopt.readthedocs.io/en/latest/) L-BFGS quasi-Newton solver for optimization. This implementation doesn't track L-BFGS iterations, so it doesn't natively support the maxIteration termination criterion. To support this termination criterion, you need to replace the NLopt source file `_deps/nlopt-src/src/algs/luksan/plis.c` by the file [`plis.c`](LBFGS_iteration_count/plis.c) we provide in the folder `LBFGS_iteration_count`. You need to rebuild the program for this change to take effect.

## How to use

**All mapping programs in this codebase follow the same command-line interface structure.** Here we demonstrate using `sTGC_2D_QN` as an example:

The executable takes 3 arguments: path to an input data file, path to a solver options file, and path to the file to store the result.

    ./sTGC_2D_QN [input_file] [solver_options_file] [result_file]

An example is provided in the `example` subdirectory. Test it by:

    ./sTGC_2D_QN example/input example/solver_options example/my_result

The result will be written to `example/my_result`.

In the 3 arguments, `input_file` is mandatory, while the rest two are optional. If `solver_options_file` is not specified, the program will look for a file named `solver_options` in the same directory as the binary. If that file is not found, the program will fall back to default options. If `result_file` is not given, results will be written to a file named `result` in the directory of the binary.


## File format

### input_file

_Input file_ contains vertices and simplices (triangles or tetrahedra) information about the source mesh and initial embedding, as well as the indices of constrained vertices (called handles or positional constraints). Vertices are indexed from 0.


    [num_sourceVert] [dimension_sourceVert]
    ... (num_sourceVert * dimension_sourceVert) Matrix ...
    [num_initVert]   [dimension_initVert]
    ... (num_initVert * dimension_initVert) Matrix ...
    [num_simplex]    [simplex_size]
    ... (num_simplex * simplex_size) Matrix ...
    [num_handles]
    ... (num_handles * 1) Matrix ...
 
 See input files under [`example`](example/) folder for a concrete example.


:tada: **It's possible to use your own mesh formats.** We provide two python scripts in directory [`IO`](IO/) to convert common mesh formats to our `input_file` format.

To use the two scripts, install [meshio](https://github.com/nschloe/meshio) with

     pip install meshio

To convert triangle meshes to our input format, run

    ./convert_input_2D.py [inputObjFile] [handleFile] [outFile]

Currently, we only support OBJ file with initial mesh as uv coordinates. Check out our [dataset](#dataset) for some concrete OBJ and handle files.
The generated `outFile` will have the format of our `input_file`.

We also provide a script in directory `IO` to generate a `handleFile` containing all the boundary vertex indices for a given input mesh. The script works for both triangle/tetrahedron mesh.

     ./extract_boundary_vert.py [inputMeshFile] [outputHandleFile] 

To convert tetrahedron rest(source) and initial meshes to our input format, run

    ./convert_input_3D.py [restFile] [initFile] [handleFile] [outFile]

All tet-mesh formats supported by `meshio` should be handled by this script. We have tested the VTK format. For more examples in VTK format, please check out our [dataset](#dataset).

 
### solver_options_file

_Solver options file_ contains parameters for the energy, options for QN/PN solver, and a list of intermediate status to record during optimization. They follow a common format but different programs may have different parameters.

See [`example`](example/) folder for concrete examples.


**sTGC_2D_QN**

    form
    [harmonic OR tutte-uniform]
    alpha
    [val]
    lambda1
    [val]
    lambda2
    [val]
    k
    [val]
    scale_rest_mesh
    [0 OR 1]
    subtract_total_signed_area
    [0 OR 1]
    aspect_ratio_threshold
    [val]
    ftol_abs
    [val]
    ftol_rel
    [val]
    xtol_abs
    [val]
    xtol_rel
    [val]
    algorithm
    [LBFGS]
    maxeval
    [val]
    stopCode
    [none OR no_flip_degenerate OR locally_injective]
    record
    vert    [0 OR 1]
    energy  [0 OR 1]
    minArea [0 OR 1]
    nbWindVert [0 OR 1]
    grad [0 OR 1]
    gradNorm [0 OR 1]
    initSingularValues [0 OR 1]
    resultSingularValues [0 OR 1]

Notes:
- 'form' determines the rest mesh used in the energy formulation. 'harmonic' means using the input rest mesh as-is. 'tutte-uniform' uses a rest mesh that consists of equilateral triangle/tetrahedron whose total area or volume is the same as the input rest mesh.
- 'lambad1', 'lambda2', 'k' determine the kind of map distortion to be minimized. They must satisfy $$lambda1 + lambda2 \times k = 1/2$$
- For standard TLC, lambda1 = 1/2, lambda2 = 0, k = 1. 
- For isometric TLC, lambda1 = 1/4, lambda2 = 1/4, k = 1.
- 'scale_rest_mesh': whether to scale the rest mesh to have the same total area/volume as the initial mesh.
- 'subtract_total_signed_area': whether to subtract the total signed area/volume from the energy.
- 'aspect_ratio_threshold': replace rest triangles whose aspect ratio is larger than this threshold with equilateral triangles with the same area. If set smaller than 1, no replacement will be performed.
- 'ftol_abs', 'ftol_rel', 'xtol_abs', 'xtol_rel': stop criteria for NLopt solver. Negative values mean disabled.
- 'stopCode' allows the solver to stop when the map becomes flip-free or locally injective.
- record options allow the solver to record intermediate status and save them to the result file.
  - 'vert': record target mesh vertices at each iteration.
  - 'energy': record energy value at each iteration.
  - 'minArea': record smallest simplex signed content (area or volume) at each iteration.
  - 'nbWindVert': record number of overwound vertices at each iteration.
  - 'grad': record gradient at each iteration.
  - 'gradNorm': record the norm of gradient at each iteration.
  - 'initSingularValues': record singular values of the initial map.
  - 'resultSingularValues': record singular values of the result map.

**sTGC_2D_PN**

    form
    [harmonic OR tutte-uniform]
    alpha
    [val]
    lambda1
    [val]
    lambda2
    [val]
    k
    [val]
    scale_rest_mesh
    [0 OR 1]
    subtract_total_signed_area
    [0 OR 1]
    aspect_ratio_threshold
    [val]
    ftol_abs
    [val]
    ftol_rel
    [val]
    xtol_abs
    [val]
    xtol_rel
    [val]
    gtol_abs
    [val]
    algorithm
    [Projected_Newton]
    maxeval
    [val]
    lineSearchGamma
    [val]
    stopCode
    [none OR no_flip_degenerate OR locally_injective]
    record
    vert    [0 OR 1]
    energy  [0 OR 1]
    minArea [0 OR 1]
    nbWindVert [0 OR 1]
    grad [0 OR 1]
    gradNorm [0 OR 1]
    searchDirection [0 OR 1]
    searchNorm [0 OR 1]
    stepSize [0 OR 1]
    stepNorm [0 OR 1]
    initSingularValues [0 OR 1]
    resultSingularValues [0 OR 1]
    lastNonFlip [0 OR 1]

Most options are the same as QN solvers, with the following additional options:
- 'gtol_abs': gradient norm stop criteria. Negative values mean disabled.
- 'lineSearchGamma': parameter for backtracking line search, see [here](https://en.wikipedia.org/wiki/Backtracking_line_search). Default value 1e-4 is recommended.
- additional record options:
  - 'searchDirection': record search direction at each iteration
  - 'searchNorm': record norm of search direction at each iteration
  - 'stepSize': record step size at each iteration
  - 'stepNorm': record norm of step size at each iteration
  - 'lastNonFlip': record the last non-flip iteration

**TLC_QN**

See [TLC repo](https://github.com/duxingyi-charles/lifting_simplices_to_find_injectivity?tab=readme-ov-file#solver_options_file) for solver opiton details.

**sTGC_3D_QN**

    form
    [harmonic OR tutte-uniform]
    alpha
    [val]
    lambda1
    [val]
    lambda2
    [val]
    k
    [val]
    scale_rest_mesh
    [0 OR 1]
    subtract_total_signed_area
    [0 OR 1]
    ftol_abs
    [val]
    ftol_rel
    [val]
    xtol_abs
    [val]
    xtol_rel
    [val]
    algorithm
    [LBFGS]
    maxeval
    [val]
    stopCode
    [none OR no_flip_degenerate]
    record
    vert    [0 OR 1]
    energy  [0 OR 1]
    minArea [0 OR 1]
    grad [0 OR 1]
    gradNorm [0 OR 1]
    initSingularValues [0 OR 1]
    resultSingularValues [0 OR 1]

Most options are the same as sTGC_2D_QN solver. 
- 'lambad1', 'lambda2', 'k' determine the kind of map distortion to be minimized. They must satisfy $$6 k^{1/3} \times lambda1 + 2k  \times lambda2 = 1$$

- For isometric TLC in 3D, set lambda1 = 1/6, lambda2 = 0, k = 1.

**sTGC_3D_PN**

    form
    [harmonic OR tutte-uniform]
    alpha
    [val]
    lambda1
    [val]
    lambda2
    [val]
    k
    [val]
    scale_rest_mesh
    [0 OR 1]
    subtract_total_signed_area
    [0 OR 1]
    ftol_abs
    [val]
    ftol_rel
    [val]
    xtol_abs
    [val]
    xtol_rel
    [val]
    gtol_abs
    [val]
    algorithm
    [Projected_Newton]
    maxeval
    [val]
    lineSearchGamma
    [val]
    stopCode
    [none OR no_flip_degenerate]
    record
    vert    [0 OR 1]
    energy  [0 OR 1]
    minArea [0 OR 1]
    grad [0 OR 1]
    gradNorm [0 OR 1]
    searchDirection [0 OR 1]
    searchNorm [0 OR 1]
    stepSize [0 OR 1]
    stepNorm [0 OR 1]
    initSingularValues [0 OR 1]
    resultSingularValues [0 OR 1]
    lastNonFlip [0 OR 1]

Most options are the same as sTGC_2D_PN solver.
- For isometric TLC in 3D, set lambda1 = 1/6, lambda2 = 0, k = 1.

**SEA_Generalized_2D_QN**

    form
    [harmonic OR tutte-uniform]
    alpha
    [val]
    lambda1
    [val]
    lambda2
    [val]
    k
    [val]
    theta
    [val]
    scale_rest_mesh
    [0 OR 1]
    aspect_ratio_threshold
    [val]
    ftol_abs
    [val]
    ftol_rel
    [val]
    xtol_abs
    [val]
    xtol_rel
    [val]
    algorithm
    [LBFGS]
    maxeval
    [val]
    stopCode
    [none OR no_flip_degenerate OR locally_injective OR globally_injective]
    record
    vert    [0 OR 1]
    energy  [0 OR 1]
    minArea [0 OR 1]
    nbWindVert [0 OR 1]
    grad [0 OR 1]
    gradNorm [0 OR 1]
    initSingularValues [0 OR 1]
    resultSingularValues [0 OR 1]
    lastInjective [0 OR 1]

Most options are the same as QN solvers:
- 'lambad1', 'lambda2', 'k' determine the kind of map distortion to be minimized. They must satisfy $$lambda1 + lambda2 \times k = 1/2$$
  - For standard SEA, lambda1 = 1/2, lambda2 = 0, k = 1.
  - For isometric SEA, lambda1 = 1/4, lambda2 = 1/4, k = 1.
- 'theta': the center angle of circular arcs replacing straight edges of the mesh boundary. See [SEA repo](https://github.com/duxingyi-charles/Smooth-Excess-Area).
- 'stopCode' allows the solver to stop when the map becomes flip-free, locally injective or globally injective.

### result_file

_Result file_ stores the vertices of result mesh, and also intermediate records as specified in solver options file.

    
    name dims
    data
    ...

See result files under[`example`](example/) folder for concrete examples.

We provide a script to convert a `result_file` to a mesh file in directory [`IO`](IO/).

Usage

    ./get_result_mesh.py [inputFile] [resultFile] [outFile]

For example, 

    ./get_result_mesh.py example/input example/result result.obj
    ./get_result_mesh.py example/input example/result result.vtk


## Dataset

For more test data, check out the two datasets we released with our papers:

- [Locally Injective Mappings Benchmark](https://github.com/duxingyi-charles/Locally-Injective-Mappings-Benchmark): Over 10000 triangle and tetrahedral meshes for mapping into domains with prescribed boundaries.
  **Download**: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3827969.svg)](https://doi.org/10.5281/zenodo.3827969)
- [Constrained Injective Mapping Benchmark](https://github.com/duxingyi-charles/Smooth-Excess-Area?tab=readme-ov-file#dataset): Over 1500 triangle meshes for mapping with arbitrary positional constraints.
  **Download**: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5547887.svg)](https://doi.org/10.5281/zenodo.5547887)

## TODO

- [ ] Add projected Newton solver for 3D TLC 