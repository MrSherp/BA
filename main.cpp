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
  typedef Eigen::Vector3i Indices3DType;
  typedef Eigen::Matrix < RealType, 2, 2 > Matrix22;
  typedef Eigen::Matrix < RealType, 3, 2 > Matrix32;
  typedef Eigen::Matrix < RealType, 3, 3 > Matrix33;
  typedef Eigen::VectorXd VectorType;
  typedef std::vector < RealType > MaskType;
  typedef Eigen::MatrixXd FullMatrixType;
  typedef Eigen::SparseMatrix < RealType, 0, int > SparseMatrixType;
  typedef Eigen::Triplet < RealType > TripletType;
};

typedef shellFE::ShellElementWithTangentSpaceAtVertex < DataTypeContainerShellFE > TriangleType;
typedef shellFE::AdaptiveTriangMesh < DataTypeContainerShellFE, TriangleType > MeshType;
typedef shellFE::DuodecQuadrature < RealType, typename DataTypeContainerShellFE::DomVecType > QuadType;
typedef shellFE::UnitTriangMeshConfiguratorP1 < DataTypeContainerShellFE, MeshType, QuadType > ConfiguratorType;

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
  const int _N;
  const boost::property_tree::ptree& _pt;
  std::string _saveDirectory;

public:
  Saver ( const int N, const boost::property_tree::ptree& Pt )
: _N ( N ), _pt ( Pt ) {
    _saveDirectory = _pt.get < std::string > ( "Directory.saveDirectory" );
    _saveDirectory += "-" + std::to_string ( addCounterToSaveDirectory ( "counter.txt" ) ) + "/";
    boost::filesystem::create_directory ( _saveDirectory );
    std::cout << "save directory: " << _saveDirectory << std::endl;
  }

  std::string getSaveDirectory () const {
    return _saveDirectory;
  }

  void saveFDImages ( const VectorNd& Image, const std::string FilenameBase ) const {
    ImageType outputImage ( _N, _N );
    vectorToImage < VectorNd, ImageType, RealType > ( Image, outputImage );
    scaleToFull < RealType, ImageType > ( outputImage );
    const std::string saveName = _saveDirectory + FilenameBase;
    saveBitmap ( saveName, outputImage );
  }

  void saveConfigFile ( const std::string FilenameBase = "parameter.ini" ) const {
    const std::string iniSave = _saveDirectory + FilenameBase;
    boost::property_tree::write_ini ( iniSave, _pt );
  }

  void writeCSVData ( const VectorNd& PrimalSolution, const VectorNd& uL, const VectorNd& uR ) const {
    VectorNd thresholdedSolution ( _N );
    const RealType threshold = _pt.get < RealType > ( "General.threshold" );
    thresholding < RealType, VectorNd > ( PrimalSolution, thresholdedSolution, threshold );

    VectorNd uRShifted ( _N );
    getShiftedUR < RealType, VectorNd > ( uRShifted, uR, thresholdedSolution, _N );
    const std::string uRshiftedSave = _saveDirectory + "uRshifted.csv";
    safeSignal < RealType, VectorNd > ( uRShifted, uRshiftedSave );

    VectorNd dataFunction = uL - uRShifted;
    const std::string dataFunctionSave = _saveDirectory + "dataFunction.csv";
    safeSignal < RealType, VectorNd > ( dataFunction, dataFunctionSave );

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

  void createPlot ( ) {
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
    out.close ( );
    const std::string syscommand = "gnuplot " + filename;
    if ( system ( syscommand.c_str () ) ) {
      std::cerr << "system command: " << syscommand << std::endl;
      throw std::runtime_error ( "Gnuplot command not working!" );
    }
  }
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
    const std::string fileName = pt.get < std::string > ( "Images.uL" ).c_str ();
    loadSignal < RealType, VectorNd > ( uL, fileName );
    loadSignal < RealType, VectorNd > ( uR, pt.get < std::string > ( "Images.uR" ) );
    const RealType lambda = pt.get < RealType > ( "General.lambda" );

    int N = uL.size ();
    convexOpt::createDisparity ( uL, uR, disparity, lambda );
    std::cout << "N: " << N << std::endl;

    //create discrete operators for gradient and the projections
    OperatorType Fd ( N );
    VectorNd primalSolution ( disparity.size () );
    VectorNd dualSolution ( 2 * primalSolution.size () );
    std::vector < Eigen::Ref < VectorNd > > dualSolutionComponents;
    for ( short i = 0; i < 2; ++i )
      dualSolutionComponents.push_back ( dualSolution.segment ( i * co::sqr ( N ), co::sqr ( N ) ) );


    auto startTime = std::clock ();
    convexOpt::applyChambollePockAccelerated < RealType, VectorNd, convexOpt::ForwardFD < RealType, VectorNd > > ( primalSolution, dualSolution, disparity, N, Fd, pt );
    std::cout << "Duration: " << (std::clock () - startTime) / static_cast < double > ( CLOCKS_PER_SEC ) << std::endl;
    Saver < RealType, VectorNd > saver ( N, pt );

    // save everything
    saver.saveConfigFile ();
    saver.saveFDImages ( primalSolution, "primalSolution.bmp" );
    saver.saveFDImages ( disparity, "disparity.bmp" );
    for ( short i = 0; i < 2; ++i ) {
      const std::string dualSave = "dualSolution" + std::to_string ( i ) + ".bmp";
      saver.saveFDImages ( dualSolutionComponents[i], dualSave );
    }
    const RealType threshold = pt.get < RealType > ( "General.threshold" );
    VectorNd thresholdedSolution ( N );
    thresholding < RealType, VectorNd > ( primalSolution, thresholdedSolution, threshold );
    const std::string thresholdedSolutionSave = saver.getSaveDirectory () + "thresholdedSolution.csv";
    safeSignal < RealType, VectorNd > ( thresholdedSolution, thresholdedSolutionSave );
    saver.writeCSVData ( primalSolution, uL, uR );

    MeshType grid ( pt.get < std::string > ( "Adaptive_Mesh.meshInput" ) );
    std::ofstream fileError ( saver.getSaveDirectory () + "errorInfo.txt");
    fileError << "iteration / localErrorFStarSum / localErrorGSum / localErrorGStarSum / totalErrorSum" << std::endl;
    for ( int iteration = 0; iteration < std::stoi ( pt.get < std::string > ( "General.numberRefinementSteps" ) ); ++iteration ) {
      ConfiguratorType config ( grid );

      // compute FE functions
      ProjectFiniteDifferenceOntoFE < ConfiguratorType > projector ( config, N );
      const RealType smoothingParameterProjection = std::stod ( pt.get < std::string > ( "Adaptive_Mesh.smoothingParameterProjection" ) );
      VectorNd primalSolutionFE ( grid.getNumVertices () );
      projector.apply ( primalSolution, primalSolutionFE, smoothingParameterProjection );
      VectorNd dualSolutionFE ( 2 * grid.getNumVertices () );
      std::vector < Eigen::Ref < VectorNd > > dualSolutionFEComponents;
      for ( short i = 0; i < 2; ++i )
        dualSolutionFEComponents.push_back ( dualSolutionFE.segment ( i * grid.getNumVertices (), grid.getNumVertices () ) );
      for ( short i = 0; i < 2; ++i )
        projector.apply < Eigen::Ref < VectorNd > > ( dualSolutionComponents[i], dualSolutionFEComponents[i], smoothingParameterProjection );

      // compute a posteriori error estimator
      VectorNd localErrorFStar ( grid.getNumTriangs () );
      VectorNd localErrorG ( grid.getNumTriangs () );
      VectorNd localErrorGStar ( grid.getNumTriangs () );
      EnergyFStar < ConfiguratorType > (config, primalSolutionFE, dualSolutionFEComponents).apply ( localErrorFStar );
      EnergyG < ConfiguratorType > (config, dualSolutionFEComponents).apply ( localErrorG );
      EnergyGStar < ConfiguratorType > (config, primalSolutionFE).apply ( localErrorGStar );
      VectorNd totalError = localErrorFStar + localErrorG + localErrorGStar;

      // plot routines
      saver.plotFEFunction < ConfiguratorType, VectorNd > ( grid, primalSolutionFE, "primalFE_iteration" + std::to_string ( iteration ) );
      for ( short i = 0; i < 2; ++i ) {
        const std::string dualFilename = "dualFE" + std::to_string ( i ) + "_iteration" + std::to_string ( iteration );
        saver.plotFEFunction < ConfiguratorType, VectorNd > ( grid, dualSolutionFEComponents[i], dualFilename );
      }
      saver.plotFEFunction < ConfiguratorType > ( grid, localErrorFStar, "localErrorFStar_iteration" + std::to_string ( iteration ) );
      saver.plotFEFunction < ConfiguratorType > ( grid, localErrorG, "localErrorG_iteration" + std::to_string ( iteration ) );
      saver.plotFEFunction < ConfiguratorType > ( grid, localErrorGStar, "localErrorGStar_iteration" + std::to_string ( iteration ) );
      saver.plotFEFunction < ConfiguratorType > ( grid, totalError, "localErrorTotal_iteration" + std::to_string ( iteration ) );

      const RealType localErrorFStarSum = localErrorFStar.sum ();
      const RealType localErrorGSum = localErrorG.sum ();
      const RealType localErrorGStarSum = localErrorGStar.sum ();
      const RealType totalErrorSum = totalError.sum ();

      std::cout << "iteration   : " << iteration << std::endl;
      std::cout << "F^*   error : " << localErrorFStarSum << std::endl;
      std::cout << "G     error : " << localErrorGSum << std::endl;
      std::cout << "G^*   error : " << localErrorGStarSum << std::endl;
      std::cout << "Total error : " << totalErrorSum << std::endl << std::endl << std::endl;

      fileError << iteration << "\t" << localErrorFStarSum << "\t" << localErrorGSum << "\t" << localErrorGStarSum << "\t" << totalErrorSum << std::endl;

      // refine grid
      const RealType thresholdRefinement = pt.get < RealType > ( "General.alpha" ) * totalError.maxCoeff ();
      for ( int i = 0; i < grid.getNumTriangs (); ++i ) {
        if ( totalError[i] > thresholdRefinement )
          grid.mark ( i );
      }
      grid.refineMarkedTriangles ();
    }
    fileError.close ( );
    saver.createPlot ( );
  }
  catch ( std::exception& error ) {
    std::cerr << "exception caught: " << error.what () << '\n';
  }
  return 0;
}
