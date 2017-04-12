#include <iostream>
#include <cstdio>
#include <ctime>
#include <math.h> 
#include <iomanip>
#include <cstdlib>
#include <stdlib.h>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCholesky>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/filesystem.hpp>
#include "imageLoadSave.h"
#include "editImages.h"
#include "chambollePockAlgorithm.h"
#include "csvLoadSafe.h"
#include "triangMesh.h"
#include "triangle.h"
#include "quadrature.h"
#include "configurators.h"
#include "FEOps.h"
#include "FETools.h"
#include "adaptiveTriangMeshWithTangentSpace.h"
#include "uniformSimplexMesh.h"

typedef double RealType;
typedef Eigen::Array < RealType, Eigen::Dynamic, Eigen::Dynamic > ImageType;
typedef Eigen::Matrix < RealType, Eigen::Dynamic, 1 > VectorNd;
typedef Eigen::Matrix < RealType, Eigen::Dynamic, Eigen::Dynamic > MatrixN;
typedef Eigen::SparseMatrix < RealType > SparseMatrixN;
typedef convexOpt::ForwardFD < RealType, VectorNd > OperatorType;
typedef convexOpt::ROFDataResolvent < RealType, VectorNd > ROFResolventDataTerm;
typedef convexOpt::KProjector < RealType, VectorNd > ProjectorOntoK;
typedef convexOpt::CProjector < RealType, VectorNd > ProjectorOntoC;

struct DataTypeContainerShellFE {
public:
  typedef double RealType;
  typedef Eigen::Vector2d DomVecType;
  typedef Eigen::Vector3d TangentVecType;
  typedef Eigen::Vector3d Point3DType;
  typedef Eigen::Vector2i Indices2DType;
  typedef Eigen::Vector3i Indices3DType;
  typedef Eigen::Matrix < RealType, 2, 2 > Matrix22;
  typedef Eigen::Matrix < RealType, 3, 2 > Matrix32;
  typedef Eigen::Matrix < RealType, 3, 3 > Matrix33;
  typedef Eigen::VectorXd VectorType;
  typedef std::vector < bool > MaskType;
  typedef Eigen::MatrixXd FullMatrixType;
  typedef Eigen::SparseMatrix < RealType, 0, int > SparseMatrixType;
  typedef Eigen::Triplet < RealType > TripletType;
};

typedef shellFE::ShellElementWithTangentSpaceAtVertex < DataTypeContainerShellFE > TriangleType;
typedef co::UniformSimplexMesh < DataTypeContainerShellFE > MeshType;
typedef shellFE::DuodecQuadrature < RealType, typename DataTypeContainerShellFE::DomVecType > QuadType;
typedef shellFE::UnitTriangMeshConfiguratorP1 < DataTypeContainerShellFE, MeshType, QuadType > ConfiguratorType;

namespace Eigen{

  template<class Matrix>
  void write_binary(const char* filename, const Matrix& matrix){
      std::ofstream out(filename,ios::out | ios::binary | ios::trunc);
      typename Matrix::Index rows=matrix.rows(), cols=matrix.cols();
      out.write((char*) (&rows), sizeof(typename Matrix::Index));
      out.write((char*) (&cols), sizeof(typename Matrix::Index));
      out.write((char*) matrix.data(), rows*cols*sizeof(typename Matrix::Scalar) );
      out.close();
  }

  template<class Matrix>
  void read_binary(const char* filename, Matrix& matrix){
      std::ifstream in(filename,ios::in | std::ios::binary);
      typename Matrix::Index rows=0, cols=0;
      in.read((char*) (&rows),sizeof(typename Matrix::Index));
      in.read((char*) (&cols),sizeof(typename Matrix::Index));
      matrix.resize(rows, cols);
      in.read( (char *) matrix.data() , rows*cols*sizeof(typename Matrix::Scalar) );
      in.close();
  }
}

template < typename RealType, typename VectorNd >
void getShiftedUR ( VectorNd& uRShifted, const VectorNd& uR, const VectorNd& thresholdedSolution, const int N ) {
  for ( int i = 0; i < N; ++i ) {
    RealType tmp = (thresholdedSolution ( i ) * (RealType) N);
    int j = static_cast < int > ( tmp );
    j += i;
    if ( j >= N )
      uRShifted ( i ) = uR ( N - 1 );
    else
      uRShifted ( i ) = uR ( j );
  }
}

template < typename RealType, typename VectorNd >
RealType getDataEnergy ( const VectorNd& DataFunction ) {
  RealType val = static_cast < RealType > ( 0 );
  const int length = DataFunction.size ();
  for ( int i = 0; i < length - 1; ++i ) {
    val += abs ( 0.5 * (std::abs ( DataFunction ( i + 1 ) ) + std::abs ( DataFunction ( i ) )) * static_cast < RealType > ( 1 ) / static_cast < RealType > ( length - 1 ) );
  }
  return val;
}

template < typename RealType, typename VectorNd >
RealType getBVEnergy ( const VectorNd& ThresholdedSolution ) {
  RealType val = static_cast < RealType > ( 0 );
  const int length = ThresholdedSolution.size ();
  for ( int i = 0; i < length - 1; ++i ) {
    val += abs ( ThresholdedSolution ( i + 1 ) - ThresholdedSolution ( i ) );
  }
  return val;
}

template < typename RealType, typename VectorNd >
class Saver {
  const boost::property_tree::ptree& _pt;
  std::string _saveDirectory;

public:
  Saver ( const boost::property_tree::ptree& Pt )
: _pt ( Pt ) {
    _saveDirectory = _pt.get < std::string > ( "Directory.saveDirectory" );
    _saveDirectory += "-" + std::to_string ( co::addCounterToSaveDirectory ( "counter.txt" ) ) + "/";
    boost::filesystem::create_directory ( _saveDirectory );
    std::cout << "save directory: " << _saveDirectory << std::endl;
  }

  std::string getSaveDirectory () const {
    return _saveDirectory;
  }

  void saveImage ( const VectorNd& Image, const int NumRows, const int NumColumns, const std::string FilenameBase, const co::IMAGE_SCALE_MODE ScaleMode ) const {
    const RealType minValue = Image.minCoeff ();
    const RealType maxValue = Image.maxCoeff ();

    if ( minValue < static_cast < RealType > ( 0 ) )
      co::displayWarning ( "Smallest entry in image below 0! Smallest entry: " + std::to_string ( minValue ) + " in " + FilenameBase );
    if ( maxValue > static_cast < RealType > ( 1 ) )
      co::displayWarning ( "Greatest entry in image greater than 1! Greatest entry: " + std::to_string ( maxValue ) + " in " + FilenameBase );
    ImageType outputImage ( NumRows, NumColumns );
    co::vectorToImage < RealType > ( Image, outputImage );
    switch ( ScaleMode ) {
    case co::SCALE_TO_FULL:
      co::scaleToFull < RealType, ImageType > ( outputImage );
      break;

    case co::MULTIPLY_TO_FULL:
      outputImage *= static_cast < RealType > ( 255 );
      break;

    default:
      throw std::runtime_error ( "Wrong ScaleMode!" );
      break;
    }
    const std::string saveName = _saveDirectory + FilenameBase;
    co::saveBitmap ( saveName, outputImage );
  }

  void saveFDImage ( const int N, const VectorNd& Image, const std::string FilenameBase ) const {
    ImageType outputImage ( N, N );
    co::vectorToImage < RealType > ( Image, outputImage );
    co::scaleToFull < RealType, ImageType > ( outputImage );
    const std::string saveName = _saveDirectory + FilenameBase;
    co::saveBitmap ( saveName, outputImage );
  }

  void saveConfigFile ( const std::string FilenameBase = "parameter.ini" ) const {
    const std::string iniSave = _saveDirectory + FilenameBase;
    boost::property_tree::write_ini ( iniSave, _pt );
  }

  void writeCSVData ( const VectorNd& PrimalSolution, const VectorNd& uL, const VectorNd& uR ) const {
    const int N = uL.size ();
    VectorNd thresholdedSolution ( N );
    const RealType threshold = _pt.get < RealType > ( "General.threshold" );
    co::thresholding < RealType, VectorNd > ( PrimalSolution, thresholdedSolution, threshold );

    VectorNd uRShifted ( N );
    getShiftedUR < RealType, VectorNd > ( uRShifted, uR, thresholdedSolution, N );
    const std::string uRshiftedSave = _saveDirectory + "uRshifted.csv";
    co::safeSignal < RealType, VectorNd > ( uRShifted, uRshiftedSave );

    VectorNd dataFunction = uL - uRShifted;
    const std::string dataFunctionSave = _saveDirectory + "dataFunction.csv";
    co::safeSignal < RealType, VectorNd > ( dataFunction, dataFunctionSave );

    const RealType dataEnergy = getDataEnergy < RealType > ( dataFunction );
    const RealType BVEnergy = getBVEnergy < RealType > ( thresholdedSolution );

    const std::string generalInfoSave = _saveDirectory + "generalInfo.txt";
    std::ofstream file ( generalInfoSave );
    if ( file.is_open () ) {
      file << "DataEnergy: " << dataEnergy << std::endl;
      file << "BVEnergy: " << BVEnergy << std::endl;
    }
    else
      throw std::runtime_error ( "Cannot open file for writing!" );
    file.close ();
  }

  template < typename ConfiguratorType, typename VectorType >
  void plotFEFunction ( const typename ConfiguratorType::InitType& Mesh, const VectorType& Dofs, const std::string FilenameBase ) const {
    typename ConfiguratorType::InitType plotMesh ( Mesh );
    for ( int i = 0; i < plotMesh.getNumVertices (); ++i ) {
      const typename Eigen::Vector3d& currentVertex = plotMesh.getVertex ( i );
      typename Eigen::Vector3d coords;
      coords << currentVertex[0], currentVertex[1], Dofs[i];
      plotMesh.setVertex ( i, coords );
    }
    const std::string filename = _saveDirectory + FilenameBase + ".vtk";
    shellFE::MeshWithData < typename ConfiguratorType::InitType > ( plotMesh ).saveAsLegacyVTK ( filename );
  }

  void createPlot () {
    const std::string filename = _saveDirectory + "plot.gp";
    std::ofstream out ( filename );
    out << "set terminal eps;" << std::endl;
    out << "set output '" << _saveDirectory << "results.eps';" << std::endl;
    out << "unset key;" << std::endl;
    out << "set key off;" << std::endl;
    out << "set loadpath '" << _saveDirectory << "';" << std::endl;
    out << "set yrange [-0.2:1];" << std::endl;
    out << "set size 1,1;" << std::endl;
    out << "set origin 0,0;" << std::endl;
    out << "set grid ytics lc rgb '#bbbbbb' lw 1 lt 0;" << std::endl;
    out << "set grid xtics lc rgb '#bbbbbb' lw 1 lt 0;" << std::endl;
    out << "set multiplot layout 2,2 columnsfirst scale 1.1,0.9;" << std::endl;
    out << "plot '" << _pt.get < std::string > ( "Images.uL" ) << "' with lines;" << std::endl;
    out << "plot '" << _saveDirectory << "uRshifted.csv' with lines, '" << _saveDirectory << "dataFunction.csv' with lines;" << std::endl;
    out << "plot '" << _pt.get < std::string > ( "Images.uR" ) << "' with lines;" << std::endl;
    out << "plot 'thresholdedSolution.csv' with lines;" << std::endl;
    out.close ();
    const std::string syscommand = "gnuplot " + filename;
    if ( system ( syscommand.c_str () ) ) {
      std::cerr << "system command: " << syscommand << std::endl;
      throw std::runtime_error ( "Gnuplot command not working!" );
    }
  }
};

enum SMOOTHING_TYPE {
  NO_SMOOTHING = 0, GAUSSIAN_BLUR = 1, PERONA_MALIK_DIFFUSION = 2
};

int main () {
  try {
    //get Parameters_Chambolle
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini ( "../../../Input/config.ini", pt );

#ifdef _OPENMP
    omp_set_num_threads ( 4 );
#endif

    //load signals and create image function g
    VectorNd uL;
    VectorNd uR;
    VectorNd disparity;
    const std::string fileName = pt.get < std::string > ( "Images.uL" );
    co::loadSignal < RealType, VectorNd > ( uL, fileName );
    co::loadSignal < RealType, VectorNd > ( uR, pt.get < std::string > ( "Images.uR" ) );
    const RealType lambda = pt.get < RealType > ( "General.lambda" );
    const BOUNDARY_TYPE_FEOPS boundaryType = static_cast < BOUNDARY_TYPE_FEOPS > ( pt.get < int > ( "General.boundaryType" ) );

    int N = uL.size ();
    convexOpt::createDisparity ( uL, uR, disparity, lambda );
    std::cout << "N: " << N << std::endl;
/*
    //create discrete operators for gradient and the projections
    OperatorType Fd ( N );
    VectorNd primalSolution ( disparity.size () );
    VectorNd dualSolution ( 2 * primalSolution.size () );
    std::vector < Eigen::Ref < VectorNd > > dualSolutionComponents;
    for ( short i = 0; i < 2; ++i )
      dualSolutionComponents.push_back ( dualSolution.segment ( i * co::sqr ( N ), co::sqr ( N ) ) );

    auto startTime = std::clock ();
    convexOpt::applyConvexOptimizationAlgorithm < RealType, VectorNd, convexOpt::ForwardFD < RealType, VectorNd > > ( primalSolution, dualSolution, disparity, N, Fd, pt, convexOpt::CHAMBOLLE_POCK_ALGORITHM1 );
    std::cout << "Duration: " << (std::clock () - startTime) / static_cast < double > ( CLOCKS_PER_SEC ) << std::endl;

    std::cout << "Primal solution min/max: " << primalSolution.minCoeff () << "\t" << primalSolution.maxCoeff () << std::endl;
    std::cout << "Dual solution 0 min/max: " << dualSolutionComponents[0].minCoeff () << "\t" << dualSolutionComponents[0].maxCoeff () << std::endl;
    std::cout << "Dual solution 1 min/max: " << dualSolutionComponents[1].minCoeff () << "\t" << dualSolutionComponents[1].maxCoeff () << std::endl;
*/
    // save everything
    Saver < RealType, VectorNd > saver ( pt );
    saver.saveConfigFile ();
/*
    saver.saveFDImage ( N, primalSolution, "primalSolution.bmp" );
    saver.saveFDImage ( N, disparity, "disparity.bmp" );
    for ( short i = 0; i < 2; ++i ) {
      const std::string dualSave = "dualSolution" + std::to_string ( i ) + ".bmp";
      saver.saveFDImage ( N, dualSolutionComponents[i], dualSave );
    }
    const RealType threshold = pt.get < RealType > ( "General.threshold" );
    VectorNd thresholdedSolution ( N );
    co::thresholding < RealType, VectorNd > ( primalSolution, thresholdedSolution, threshold );
    const std::string thresholdedSolutionSave = saver.getSaveDirectory () + "thresholdedSolution.csv";
    co::safeSignal < RealType, VectorNd > ( thresholdedSolution, thresholdedSolutionSave );
    saver.writeCSVData ( primalSolution, uL, uR );
        
    const std::string str1 = saver.getSaveDirectory () + "primal.txt";
    const char* PrimalSolutionSave = str1.c_str();
    write_binary(PrimalSolutionSave, primalSolution);
    const std::string str2 = saver.getSaveDirectory () + "dual1.txt";
    const char* DualSolution1Save = str2.c_str();
    write_binary(DualSolution1Save, dualSolutionComponents[0]);
    const std::string str3 = saver.getSaveDirectory () + "dual2.txt";
    const char* DualSolution2Save = str3.c_str();
    write_binary(DualSolution2Save, dualSolutionComponents[1]);

    //test
    VectorNd readerPrimal(primalSolution.size());
    read_binary(PrimalSolutionSave, readerPrimal);
    saver.saveFDImage(N, readerPrimal, "primalSolutionAfterReading.bmp");
    VectorNd readerDual1(primalSolution.size());
    read_binary(DualSolution1Save, readerDual1);
    saver.saveFDImage(N, readerDual1, "dualSolution1AfterReading.bmp");
    VectorNd readerDual2(primalSolution.size());
    read_binary(DualSolution2Save, readerDual2);
    saver.saveFDImage(N, readerDual2, "dualSolution2AfterReading.bmp");

   */

    VectorNd readerPrimal(disparity.size());
    const std::string fileName1 = pt.get < std::string > ( "Solutions.primal" );
    const char* read1 = fileName1.c_str();
    read_binary(read1, readerPrimal);
    VectorNd readerDual0(disparity.size());
    const std::string fileName2 = pt.get < std::string > ( "Solutions.dual1" );
    const char* read2 = fileName2.c_str();
    read_binary(read2, readerDual0);
    VectorNd readerDual1(disparity.size());
    const std::string fileName3 = pt.get < std::string > ( "Solutions.dual2" );
    const char* read3 = fileName3.c_str();
    read_binary(read3, readerDual1);

    saver.saveFDImage(N, readerPrimal, "primalSolution.bmp");
    saver.saveFDImage(N, readerDual0, "dual1Solution.bmp");
    saver.saveFDImage(N, readerDual1, "dual2Solution.bmp");

    std::ofstream fileError ( saver.getSaveDirectory () + "errorInfo.txt" );
    fileError << "iteration / refined elements / number of elements / localErrorFStarSum / localErrorGSum / localErrorGStarSum / totalErrorSum" << std::endl;
    for ( int iteration = pt.get < int > ( "General.startGridLevel" ); iteration < pt.get < int > ( "General.stopGridLevel" ); ++iteration ) {
      const int numDofsInDirection = std::pow ( 2, iteration ) + 1;
      MeshType grid ( numDofsInDirection, numDofsInDirection );

      ConfiguratorType config ( grid );
      std::vector < bool > dirichletMask;
      co::createBoundaryMaskDirichletOnUnitSquare < RealType, MeshType, ConfiguratorType::Point3DType > ( grid, dirichletMask );

      // compute FE functions
      ProjectFiniteDifferenceOntoFE < ConfiguratorType > projector ( config, N, boundaryType, &dirichletMask );

      const RealType hAverage = std::sqrt ( static_cast < RealType > ( 1 ) / static_cast < RealType > ( grid.getNumTriangs () ) );
      // caution: relative smoothing parameter
      const RealType smoothingParameterProjection = pt.get < RealType > ( "Adaptive_Mesh.smoothingParameterProjection" ) * hAverage;


      VectorNd primalSolutionFE ( grid.getNumVertices () );
      projector.apply ( readerPrimal, primalSolutionFE, smoothingParameterProjection );

      VectorNd dualSolutionFE0 ( grid.getNumVertices () );
      VectorNd dualSolutionFE1 ( grid.getNumVertices () );
      projector.apply  ( readerDual0, dualSolutionFE0 , smoothingParameterProjection );
      projector.apply  ( readerDual1, dualSolutionFE1 , smoothingParameterProjection );

      // smoothing
      // caution: relative smoothing parameter
      const RealType smoothingParameterPrimalSmoothing = pt.get < RealType > ( "Adaptive_Mesh.smoothingParameterPrimal" ) * hAverage;
      const RealType smoothingParameterDualSmoothing = pt.get < RealType > ( "Adaptive_Mesh.smoothingParameterDual" ) * hAverage;
      const RealType weightParameterPeronaMalik = pt.get < RealType > ( "Adaptive_Mesh.weightParameterPeronaMalik" );

      const SMOOTHING_TYPE smoothingType = static_cast < SMOOTHING_TYPE > ( pt.get < int > ( "Adaptive_Mesh.smoothingType" ) );
      switch ( smoothingType ) {
      case GAUSSIAN_BLUR: {
        LinearSmoothOp < ConfiguratorType > smoothOp ( config, boundaryType, &dirichletMask );
        smoothOp.apply ( primalSolutionFE, smoothingParameterPrimalSmoothing );
        smoothOp.apply ( dualSolutionFE0, smoothingParameterDualSmoothing );
        smoothOp.apply ( dualSolutionFE1, smoothingParameterDualSmoothing );
      }
      break;

      case PERONA_MALIK_DIFFUSION: {
        PeronaMalikDiffusionFiltering < ConfiguratorType > smoothOp ( config, boundaryType, &dirichletMask );
        smoothOp.apply ( primalSolutionFE, weightParameterPeronaMalik, smoothingParameterPrimalSmoothing );
        smoothOp.apply ( dualSolutionFE0, weightParameterPeronaMalik, smoothingParameterDualSmoothing );
        smoothOp.apply ( dualSolutionFE1, weightParameterPeronaMalik, smoothingParameterDualSmoothing );
      }

      break;

      default:
        throw std::runtime_error ( "Bad smoothingType!" );
        break;
      }

      // compute a posteriori error estimator
      VectorNd localErrorFStar ( grid.getNumTriangs () );
      VectorNd localErrorG ( grid.getNumTriangs () );
      VectorNd localErrorGStar ( grid.getNumTriangs () );

      EnergyFStar < ConfiguratorType > ( config, primalSolutionFE, dualSolutionFE0, dualSolutionFE1 ).apply ( localErrorFStar );
      EnergyG < ConfiguratorType > ( config, dualSolutionFE0, dualSolutionFE1 ).apply ( localErrorG );
      EnergyGStar < ConfiguratorType > ( config, primalSolutionFE ).apply ( localErrorGStar );
      VectorNd totalError = localErrorFStar + localErrorG + localErrorGStar;

      VectorNd localErrorFStarNode, localErrorGNode, localErrorGStarNode, totalErrorNode;
      co::distributeElementwiseErrorToNodesOnUnitSquare < RealType > ( localErrorFStar, localErrorFStarNode, grid );
      co::distributeElementwiseErrorToNodesOnUnitSquare < RealType > ( localErrorG, localErrorGNode, grid );
      co::distributeElementwiseErrorToNodesOnUnitSquare < RealType > ( localErrorGStar, localErrorGStarNode, grid );
      co::distributeElementwiseErrorToNodesOnUnitSquare < RealType > ( totalError, totalErrorNode, grid );

      //       // refine grid
      //       int numRefinedElements = 0;
      //       //const RealType thresholdRefinement = pt.get < RealType > ( "General.alpha" ) * totalError.maxCoeff ();
      //       for ( int i = 0; i < grid.getNumTriangs (); ++i ) {
      //         // TODO: change for adaptive grid
      //           grid.mark ( i );
      //           ++numRefinedElements;
      //       }
      //       grid.refineMarkedTriangles ();

      // plot routines
      saver.saveImage ( primalSolutionFE, grid.getNumX (), grid.getNumY (), "primalFE_iteration" + std::to_string ( iteration ), co::SCALE_TO_FULL );
      const std::string dualFilename0 = "dualFE" + std::to_string ( 0 ) + "_iteration" + std::to_string ( iteration );
      saver.saveImage ( dualSolutionFE0, grid.getNumX (), grid.getNumY (), dualFilename0, co::SCALE_TO_FULL );
      const std::string dualFilename1 = "dualFE" + std::to_string ( 1 ) + "_iteration" + std::to_string ( iteration );
      saver.saveImage ( dualSolutionFE1, grid.getNumX (), grid.getNumY (), dualFilename1, co::SCALE_TO_FULL );

      saver.saveImage ( localErrorFStarNode, grid.getNumX (), grid.getNumY (), "localErrorFStar_iteration" + std::to_string ( iteration ), co::SCALE_TO_FULL );
      saver.saveImage ( localErrorGNode, grid.getNumX (), grid.getNumY (), "localErrorG_iteration" + std::to_string ( iteration ), co::SCALE_TO_FULL );
      saver.saveImage ( localErrorGStarNode, grid.getNumX (), grid.getNumY (), "localErrorGStar_iteration" + std::to_string ( iteration ), co::SCALE_TO_FULL );
      saver.saveImage ( totalErrorNode, grid.getNumX (), grid.getNumY (), "localErrorTotal_iteration" + std::to_string ( iteration ), co::SCALE_TO_FULL );

      const RealType localErrorFStarSum = localErrorFStar.sum ();
      const RealType localErrorGSum = localErrorG.sum ();
      const RealType localErrorGStarSum = localErrorGStar.sum ();
      const RealType totalErrorSum = totalError.sum ();

      std::cout << "iteration         : " << iteration << std::endl;
      //       std::cout << "refined elements  : " << numRefinedElements << std::endl;
      std::cout << "number of elements: " << grid.getNumTriangs () << std::endl;
      std::cout << "F^*   error       : " << localErrorFStarSum << std::endl;
      std::cout << "G     error       : " << localErrorGSum << std::endl;
      std::cout << "G^*   error       : " << localErrorGStarSum << std::endl;
      std::cout << "Total error       : " << totalErrorSum << std::endl << std::endl << std::endl;

      fileError << iteration << "\t" /*<< numRefinedElements << "\t"*/<< grid.getNumTriangs () << "\t" << localErrorFStarSum << "\t" << localErrorGSum << "\t" << localErrorGStarSum << "\t" << totalErrorSum << std::endl;

    }
    fileError.close ();
    saver.createPlot ();
  }
  catch ( std::exception& error ) {
    std::cerr << "exception caught: " << error.what () << '\n';
  }
  return 0;
}
