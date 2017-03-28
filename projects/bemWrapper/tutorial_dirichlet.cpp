// Copyright (C) 2011 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include "assembly/assembly_options.hpp"
#include "assembly/boundary_operator.hpp"
#include "assembly/context.hpp"
#include "assembly/evaluation_options.hpp"
#include "assembly/grid_function.hpp"
#include "assembly/l2_norm.hpp"
#include "assembly/numerical_quadrature_strategy.hpp"
#include "assembly/surface_normal_independent_function.hpp"

#include "assembly/identity_operator.hpp"
#include "assembly/laplace_3d_single_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_double_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_single_layer_potential_operator.hpp"
#include "assembly/laplace_3d_double_layer_potential_operator.hpp"

#include "common/boost_make_shared_fwd.hpp"

#include "grid/grid.hpp"
#include "grid/grid_factory.hpp"

#include "linalg/default_iterative_solver.hpp"

#include "space/piecewise_linear_continuous_scalar_space.hpp"
#include "space/piecewise_constant_scalar_space.hpp"

#include <iostream>
#include <fstream>

using namespace std;

string inputFilename;
string outputFilename;
int outputRes;

typedef double BFT; // basis function type
typedef double RT; // result type (type used to represent discrete operators)
typedef double CT; // coordinate type

int xRes, yRes, zRes;
double xCenter, yCenter, zCenter;
double xLength, yLength, zLength;

// read in a FIELD_3D file just to strip out the dimensions
void readFieldDimensions(const string& filename)
{
  FILE* file;
  file = fopen(filename.c_str(), "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " FIELD_3D read failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }

  fread((void*)&xRes, sizeof(int), 1, file);
  fread((void*)&yRes, sizeof(int), 1, file);
  fread((void*)&zRes, sizeof(int), 1, file);
  
  fread((void*)&xCenter, sizeof(double), 1, file);
  fread((void*)&yCenter, sizeof(double), 1, file);
  fread((void*)&zCenter, sizeof(double), 1, file);

  fread((void*)&xLength, sizeof(double), 1, file);
  fread((void*)&yLength, sizeof(double), 1, file);
  fread((void*)&zLength, sizeof(double), 1, file);

  fclose(file);

  cout << " Field res:     " << xRes << " " << yRes << " " << zRes << endl;
  cout << " Field center:  " << xCenter << " " << yCenter << " " << zCenter << endl;
  cout << " Field lengths: " << xLength << " " << yLength << " " << zLength << endl;
}

/*
def evalDirichletTrace(point):
    return 1
    */
class DirichletDataTrace
{
public:
    // Type of the function's values (e.g. float or std::complex<double>)
    typedef RT ValueType;
    // Type of coordinates (must be the "real part" of ValueType)
    typedef CT CoordinateType;

    // Number of components of the function's argument
    int argumentDimension() const { return 3; }
    // Number of components of the function's value
    int resultDimension() const { return 1; }

    // Evaluate the function at the point "point" and store result in
    // the array "result"
    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
        result(0) = 1;
    }
};
/*
def evalDirichletData(point):
    x, y, z = point
    r = np.sqrt(x**2 + y**2 + z**2)
    return 2 * x * z / r**5 - y / r**3
    */
class DirichletData
{
public:
    // Type of the function's values (e.g. float or std::complex<double>)
    typedef RT ValueType;
    // Type of coordinates (must be the "real part" of ValueType)
    typedef CT CoordinateType;

    // Number of components of the function's argument
    int argumentDimension() const { return 3; }
    // Number of components of the function's value
    int resultDimension() const { return 1; }

    // Evaluate the function at the point "point" and store result in
    // the array "result"
    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
        CoordinateType x = point(0), y = point(1), z = point(2);
        CoordinateType r = sqrt(point(0) * point(0) +
                point(1) * point(1) +
                point(2) * point(2));
        result(0) = 2 * x * z / (r * r * r * r * r) - y / (r * r * r);
        //result(0) = 1;
    }
};

/*
def evalExactNeumannData(point):
    x, y, z = point
    r = np.sqrt(x**2 + y**2 + z**2)
    return -6 * x * z / r**6 + 2 * y / r**4
    */
class ExactNeumannData
{
public:
    // Type of the function's values (e.g. float or std::complex<double>)
    typedef RT ValueType;
    // Type of coordinates (must be the "real part" of ValueType)
    typedef CT CoordinateType;

    // Number of components of the function's argument
    int argumentDimension() const { return 3; }
    // Number of components of the function's value
    int resultDimension() const { return 1; }

    // Evaluate the function at the point "point" and store result in
    // the array "result"
    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
        CoordinateType x = point(0), y = point(1), z = point(2);
        CoordinateType r = sqrt(point(0) * point(0) +
                point(1) * point(1) +
                point(2) * point(2));
        result(0) = -6 * x * z / (r * r * r * r * r * r) + 2 * y / (r * r * r * r);
        //result(0) = 0;
    }
};

int dirichletExample()
{
    // Import symbols from namespace Bempp to the global namespace

    using namespace Bempp;

    // Load mesh

    /*
    // PYTHON
    grid = createGridFactory().importGmshGrid(
    "triangular", "../../examples/meshes/sphere-h-0.2.msh")
    */
    const char* meshFile = "meshes/sphere-h-0.2.msh";
    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(params, meshFile);

    // Define the quadrature strategy
    /*
    // PYTHON
    accuracyOptions = createAccuracyOptions()
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2)
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(2)
    quadStrategy = createNumericalQuadratureStrategy( "float64", "float64", accuracyOptions)
    */
    AccuracyOptions accuracyOptions;
    // Increase by 2 the order of quadrature rule used to approximate
    // integrals of regular functions on pairs on elements
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2);
    // Increase by 2 the order of quadrature rule used to approximate
    // integrals of regular functions on single elements
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(2);
    NumericalQuadratureStrategy<BFT, RT> quadStrategy(accuracyOptions);

    // Specify the assembly method. We want to use ACA
    /*
    // PYTHON
    assemblyOptions = createAssemblyOptions()
    assemblyOptions.switchToAcaMode(createAcaOptions())
    context = createContext(quadStrategy, assemblyOptions)
    */
    AssemblyOptions assemblyOptions;
    AcaOptions acaOptions; // Default parameters for ACA
    assemblyOptions.switchToAcaMode(acaOptions);

    // Initialize the spaces
    /*
    // PYTHON
    pwiseConstants = createPiecewiseConstantScalarSpace(context, grid)
    pwiseLinears = createPiecewiseLinearContinuousScalarSpace(context, grid)
    */
    PiecewiseLinearContinuousScalarSpace<BFT> pwiseLinears(grid);
    PiecewiseConstantScalarSpace<BFT> pwiseConstants(grid);


    // Create the assembly context

    Context<BFT, RT> context(make_shared_from_ref(quadStrategy), assemblyOptions);

    // Construct elementary operators
    /*
    // PYTHON
    slpOp = createLaplace3dSingleLayerBoundaryOperator(
        context, pwiseConstants, pwiseLinears, pwiseConstants)
    dlpOp = createLaplace3dDoubleLayerBoundaryOperator(
        context, pwiseLinears, pwiseLinears, pwiseConstants)
    idOp = createIdentityOperator(
        context, pwiseLinears, pwiseLinears, pwiseConstants)
    */
    BoundaryOperator<BFT, RT> slpOp =
            laplace3dSingleLayerBoundaryOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseConstants),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseConstants));
    BoundaryOperator<BFT, RT> dlpOp =
            laplace3dDoubleLayerBoundaryOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseConstants));
    BoundaryOperator<BFT, RT> idOp =
            identityOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseConstants));

    // Form the right-hand side sum
    /*
    // PYTHON - where is lhs in the C++?
    lhsOp = slpOp
    rhsOp = -0.5 * idOp + dlpOp
    */
    BoundaryOperator<BFT, RT> rhsOp = -0.5 * idOp + dlpOp;

    // Construct the grid function representing the (input) Dirichlet data
    /*
    // PYTHON
    dirichletData = createGridFunction( context, pwiseLinears, pwiseLinears, evalDirichletData)
    */
    GridFunction<BFT, RT> dirichletData(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseLinears),
                surfaceNormalIndependentFunction(DirichletData()));

    // Construct the right-hand-side grid function
    /*
    // PYTHON
    rhs = rhsOp * dirichletData
    */
    GridFunction<BFT, RT> rhs = rhsOp * dirichletData;

    // Initialize the solver
    /*
    // PYTHON
    solver = createDefaultIterativeSolver(lhsOp)
    solver.initializeSolver(defaultGmresParameterList(1e-5))
    */
    DefaultIterativeSolver<BFT, RT> solver(slpOp);
    solver.initializeSolver(defaultGmresParameterList(1e-5));

    // Solve the equation
    /*
    // PYTHON
    solution = solver.solve(rhs)
    print solution.solverMessage()
    */
    Solution<BFT, RT> solution = solver.solve(rhs);
    std::cout << solution.solverMessage() << std::endl;

    // Extract the solution in the form of a grid function
    // and export it in VTK format
    /*
    // PYTHON
    solFun = solution.gridFunction()
    solFun.exportToVtk("cell_data", "neumann_data", "solution")
    */
    const GridFunction<BFT, RT>& solFun = solution.gridFunction();
    exportToVtk(solFun, VtkWriter::CELL_DATA, "Neumann_data", "solution");

/*
    // Compare the numerical and analytical solution on the grid

    // GridFunction<BFT, RT> exactSolFun(
    //             make_shared_from_ref(context),
    //             make_shared_from_ref(pwiseConstants),
    //             make_shared_from_ref(pwiseConstants),
    //             surfaceNormalIndependentFunction(ExactNeumannData()));
    CT absoluteError, relativeError;
    estimateL2Error(
                solFun, surfaceNormalIndependentFunction(ExactNeumannData()),
                quadStrategy, absoluteError, relativeError);
    std::cout << "Relative L^2 error: " << relativeError << std::endl;

    // GridFunction<BFT, RT> diff = solFun - exactSolFun;
    // double relativeError = diff.L2Norm() / exactSolFun.L2Norm();
    // std::cout << "Relative L^2 error: " << relativeError << std::endl;
    */

    // Prepare to evaluate the solution on an annulus outside the sphere

    // Create potential operators

    /*
    // PYTHON
    slPotOp = createLaplace3dSingleLayerPotentialOperator(context)
    dlPotOp = createLaplace3dDoubleLayerPotentialOperator(context)
    */
    Laplace3dSingleLayerPotentialOperator<BFT, RT> slPotOp;
    Laplace3dDoubleLayerPotentialOperator<BFT, RT> dlPotOp;

    // Construct the array 'evaluationPoints' containing the coordinates
    // of points where the solution should be evaluated
    /*
    // PYTHON
    rCount = 51;
    thetaCount = 361;
    r, theta, z = np.mgrid[1:2:rCount*1j, 0:2*np.pi:thetaCount*1j, 0:0:1j]
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    evaluationPoints = np.vstack((x.ravel(), y.ravel(), z.ravel()))
    */
    const int rCount = 51;
    const int thetaCount = 361;
    const CT minTheta = 0., maxTheta = 2. * M_PI;
    const CT minR = 1., maxR = 2.;
    const int dimWorld = 3;
    arma::Mat<CT> evaluationPoints(dimWorld, rCount * thetaCount);
    for (int iTheta = 0; iTheta < thetaCount; ++iTheta) {
        CT theta = minTheta + (maxTheta - minTheta) *
            iTheta / (thetaCount - 1);
        for (int iR = 0; iR < rCount; ++iR) {
            CT r = minR + (maxR - minR) * iR / (rCount - 1);
            evaluationPoints(0, iR + iTheta * rCount) = r * cos(theta); // x
            evaluationPoints(1, iR + iTheta * rCount) = r * sin(theta); // y
            evaluationPoints(2, iR + iTheta * rCount) = 0.;             // z
        }
    }

    // Use the Green's representation formula to evaluate the solution
    /*
    // PYTHON
    evaluationOptions = createEvaluationOptions()
    field = (-slPotOp.evaluateAtPoints(solFun, evaluationPoints,
                                       evaluationOptions) +
              dlPotOp.evaluateAtPoints(dirichletData, evaluationPoints,
                                       evaluationOptions))
    */
    EvaluationOptions evaluationOptions;
    arma::Mat<RT> field =
        -slPotOp.evaluateAtPoints(solFun, evaluationPoints,
                                  quadStrategy, evaluationOptions) +
         dlPotOp.evaluateAtPoints(dirichletData, evaluationPoints,
                                  quadStrategy, evaluationOptions);

    /*
    // can we retrieve the points here?
    typedef Fiber::EvaluatorForIntegralOperators<RT> Evaluator;
    std::auto_ptr<Evaluator> evaluator = slPotOp.makeEvaluator(solFun, quadStrategy, evaluationOptions);
    typedef Fiber::EvaluatorForIntegralOperators<RT> Base;
    //typedef typename ScalarTraits<RT>::RealType CoordinateType;
    typedef ScalarTraits<RT>::RealType CoordinateType;
    //typedef Base::CoordinateType CoordinateType;

    //Fiber::GeometricalData<CoordinateType>& trialGeomData = evaluator->m_farFieldTrialGeomData;
    evaluator->m_farFieldTrialGeomData;
    //const Fiber::GeometricalData<RT>& trialGeomData = evaluator->m_farFieldTrialGeomData;
    */


#if 0
    // Export the solution into text file
    std::ofstream out("solution.txt");
    out << "# x y z u\n";
    for (int i = 0; i < rCount * thetaCount; ++i)
        out << evaluationPoints(0, i) << ' '
            << evaluationPoints(1, i) << ' '
            << evaluationPoints(2, i) << ' '
            << field(0, i) << '\n';
    out.close();

    // Export the solution into Matlab file
    std::ofstream mout("solution.m");
    mout << " x = [ " << endl;
    for (int i = 0; i < rCount * thetaCount; ++i)
      mout << evaluationPoints(0, i) << ' ';
    mout << "];" << endl;
    mout << " y = [ " << endl;
    for (int i = 0; i < rCount * thetaCount; ++i)
      mout << evaluationPoints(1, i) << ' ';
    mout << "];" << endl;
    mout << " z = [ " << endl;
    for (int i = 0; i < rCount * thetaCount; ++i)
      mout << evaluationPoints(2, i) << ' ';
    mout << "];" << endl;
    mout << " field = [ " << endl;
    for (int i = 0; i < rCount * thetaCount; ++i)
      mout << field(0, i) << ' ';
    mout << "];" << endl;
    mout.close();

    // compute the solution over a regular grid
    int xRes = 100;
    int yRes = 100;
    double length = 4;
    arma::Mat<CT> evaluationPointsGrid(dimWorld, xRes * yRes);
    int index = 0;
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++, index++)
      {
        double xReal = ((double)x / xRes) * length - length / 2;
        double yReal = ((double)y / yRes) * length - length / 2;
        evaluationPointsGrid(0, index) = xReal;
        evaluationPointsGrid(1, index) = yReal;
        evaluationPointsGrid(2, index) = 0;
      }
    arma::Mat<RT> fieldGrid =
        -slPotOp.evaluateAtPoints(solFun, evaluationPointsGrid,
                                  quadStrategy, evaluationOptions) +
         dlPotOp.evaluateAtPoints(dirichletData, evaluationPointsGrid,
                                  quadStrategy, evaluationOptions);

    std::ofstream gout("solutionGrid.m");
    gout << " field = [ " << endl;
    index = 0;
    for (int y = 0; y < yRes; y++)
    {
      for (int x = 0; x < xRes; x++, index++)
        gout << fieldGrid(0, index) << " ";
      gout << ";" << endl;
    }
    gout << "];" << endl;
    gout.close();

    // output to a FIELD_3D
    {
      cout << " Generating 3D data ... " << flush;
      int res = 300;
      int xRes = res;
      int yRes = res;
      int zRes = res;
      double length = 4;
      arma::Mat<CT> evaluationPoints3D(dimWorld, xRes * yRes * zRes);
      int index = 0;
      for (int z = 0; z < zRes; z++)
        for (int y = 0; y < yRes; y++)
          for (int x = 0; x < xRes; x++, index++)
          {
            double xReal = ((double)x / xRes) * length - length / 2;
            double yReal = ((double)y / yRes) * length - length / 2;
            double zReal = ((double)z / zRes) * length - length / 2;
            evaluationPoints3D(0, index) = xReal;
            evaluationPoints3D(1, index) = yReal;
            evaluationPoints3D(2, index) = zReal;
          }
      arma::Mat<RT> field3D=
          -slPotOp.evaluateAtPoints(solFun, evaluationPoints3D,
                                    quadStrategy, evaluationOptions) +
           dlPotOp.evaluateAtPoints(dirichletData, evaluationPoints3D,
                                    quadStrategy, evaluationOptions);
      cout << " done. " << endl;

      string filename("bempp.field3d");
      cout << " Writing file " << filename.c_str() << " .... " << flush;
      FILE* file;
      file = fopen(filename.c_str(), "wb");
      if (file == NULL)
      {
        cout << " Couldn't open file " << filename.c_str() << endl;
        exit(0);
      }
      fwrite((void*)&xRes, sizeof(int), 1, file);
      fwrite((void*)&yRes, sizeof(int), 1, file);
      fwrite((void*)&zRes, sizeof(int), 1, file);

      double center[] = {0,0,0};
      double lengths[] = {4,4,4};
      fwrite((void*)&center[0], sizeof(double), 1, file);
      fwrite((void*)&center[1], sizeof(double), 1, file);
      fwrite((void*)&center[2], sizeof(double), 1, file);
      fwrite((void*)&lengths[0], sizeof(double), 1, file);
      fwrite((void*)&lengths[1], sizeof(double), 1, file);
      fwrite((void*)&lengths[2], sizeof(double), 1, file);

      int totalCells = xRes * yRes * zRes;
      double* data = new double[totalCells];
      for (int x = 0; x < totalCells; x++)
        data[x] = 1.0 - field3D(0, x);
      fwrite((void*)data, sizeof(double), totalCells, file);

      delete[] data;
      cout << " done. " << endl;
      fclose(file);
    }
#endif
}

int boltExample()
{
    // Import symbols from namespace Bempp to the global namespace

    using namespace Bempp;

    // Load mesh
    cout << " Loading mesh file ... " << inputFilename.c_str() << flush;
    const char* meshFile = inputFilename.c_str();
    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(params, meshFile);
    cout << " done. " << endl;

    // Define the quadrature strategy
    AccuracyOptions accuracyOptions;
    // Increase by 2 the order of quadrature rule used to approximate
    // integrals of regular functions on pairs on elements
    //accuracyOptions.doubleRegular.setRelativeQuadratureOrder(4);
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(8);
    // Increase by 2 the order of quadrature rule used to approximate
    // integrals of regular functions on single elements
    //accuracyOptions.singleRegular.setRelativeQuadratureOrder(2);
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(8);
    NumericalQuadratureStrategy<BFT, RT> quadStrategy(accuracyOptions);

    // Specify the assembly method. We want to use ACA
    AssemblyOptions assemblyOptions;
    AcaOptions acaOptions; // Default parameters for ACA
    acaOptions.eps = 1e-6;
    acaOptions.maximumBlockSize = (int)(grid->leafView()->entityCount(0) / 8);
    assemblyOptions.switchToAcaMode(acaOptions);

    // Create the assembly context
    // PYTHON
    // context = createContext(quadStrategy, assemblyOptions)
    Context<BFT, RT> context(make_shared_from_ref(quadStrategy), assemblyOptions);

    // Initialize the spaces
    PiecewiseLinearContinuousScalarSpace<BFT> pwiseLinears(grid);
    PiecewiseConstantScalarSpace<BFT> pwiseConstants(grid);

    // Construct elementary operators
    BoundaryOperator<BFT, RT> slpOp =
            laplace3dSingleLayerBoundaryOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseConstants),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseConstants), "SLP");
    BoundaryOperator<BFT, RT> dlpOp =
            laplace3dDoubleLayerBoundaryOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseConstants), "DLP");
    BoundaryOperator<BFT, RT> idOp =
            identityOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseConstants), "Id");

    // Form the right-hand side sum
    BoundaryOperator<BFT, RT> rhsOp = -0.5 * idOp + dlpOp;

    // Construct the grid function representing the (input) Dirichlet data
    GridFunction<BFT, RT> dirichletTrace(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseLinears),
                surfaceNormalIndependentFunction(DirichletDataTrace()));

    // Construct the right-hand-side grid function
    GridFunction<BFT, RT> rhs = rhsOp * dirichletTrace;

    // Initialize the solver
    DefaultIterativeSolver<BFT, RT> solver(slpOp);
    solver.initializeSolver(defaultGmresParameterList(1e-8));

    // Solve the equation
    Solution<BFT, RT> solution = solver.solve(rhs);
    std::cout << solution.solverMessage() << std::endl;

    // Extract the solution in the form of a grid function
    // and export it in VTK format
    const GridFunction<BFT, RT>& solFun = solution.gridFunction();
    exportToVtk(solFun, VtkWriter::CELL_DATA, "Neumann_data", "solution");

    // Create potential operators
    Laplace3dSingleLayerPotentialOperator<BFT, RT> slPotOp;
    Laplace3dDoubleLayerPotentialOperator<BFT, RT> dlPotOp;

    // output to a FIELD_3D
    {
      const int dimWorld = 3;
      cout << " Generating 3D data ... " << flush;
      double length = 1;
      double xHalfLength = 0.5 * xLength;
      double yHalfLength = 0.5 * yLength;
      double zHalfLength = 0.5 * zLength;

      double dx = xLength / xRes;
      double dy = yLength / yRes;
      double dz = zLength / zRes;

      arma::Mat<CT> evaluationPoints3D(dimWorld, xRes * yRes * zRes);
      for (int z = 0; z < zRes; z++)
        for (int y = 0; y < yRes; y++)
          for (int x = 0; x < xRes; x++)
          {
            int index = x + y * xRes + z * xRes * yRes;
            double xReal = xCenter - xHalfLength + x * dx + 0.5 * dx;
            double yReal = yCenter - yHalfLength + y * dy + 0.5 * dy;
            double zReal = zCenter - zHalfLength + z * dz + 0.5 * dz;

            evaluationPoints3D(0, index) = xReal;
            evaluationPoints3D(1, index) = yReal;
            evaluationPoints3D(2, index) = zReal;
          }
      EvaluationOptions evaluationOptions;
      cout << " Evaluating points ... " << flush;
      arma::Mat<RT> field3D=
          -slPotOp.evaluateAtPoints(solFun, evaluationPoints3D,
                                    quadStrategy, evaluationOptions) +
           dlPotOp.evaluateAtPoints(dirichletTrace, evaluationPoints3D,
                                    quadStrategy, evaluationOptions);
      cout << " done. " << endl;

      //string filename("bempp.field3d");
      string filename = outputFilename;
      cout << " Writing file " << filename.c_str() << " .... " << flush;
      FILE* file;
      file = fopen(filename.c_str(), "wb");
      if (file == NULL)
      {
        cout << " Couldn't open file " << filename.c_str() << endl;
        exit(0);
      }
      fwrite((void*)&xRes, sizeof(int), 1, file);
      fwrite((void*)&yRes, sizeof(int), 1, file);
      fwrite((void*)&zRes, sizeof(int), 1, file);

      fwrite((void*)&xCenter, sizeof(double), 1, file);
      fwrite((void*)&yCenter, sizeof(double), 1, file);
      fwrite((void*)&zCenter, sizeof(double), 1, file);
      fwrite((void*)&xLength, sizeof(double), 1, file);
      fwrite((void*)&yLength, sizeof(double), 1, file);
      fwrite((void*)&zLength, sizeof(double), 1, file);

      int totalCells = xRes * yRes * zRes;
      double* data = new double[totalCells];
      for (int x = 0; x < totalCells; x++)
        data[x] = 1.0 - field3D(0, x);
      fwrite((void*)data, sizeof(double), totalCells, file);

      delete[] data;
      cout << " done. " << endl;
      fclose(file);
    }
}

int main(int argc, char** argv)
{
  if (argc != 4)
  {
    cout << " USAGE: " << argv[0] << " <input Gmsh filename> <input distance field (just for the dimensions)> <output field3D filename> " << endl;
    exit(0);
  }

  inputFilename = string(argv[1]);
  string distanceFilename = string(argv[2]);
  outputFilename = string(argv[3]);

  cout << " Using the Gmsh file:                       " << inputFilename.c_str() << endl;
  cout << " Using the dimensions in the FIELD_3D file: " << distanceFilename.c_str() << endl;
  cout << " Outputting to:                             " << outputFilename.c_str() << endl;

  readFieldDimensions(distanceFilename);

  //dirichletExample();
  boltExample();
}
