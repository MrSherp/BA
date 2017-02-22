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
#include "triangleShellFE.h"
#include "quadratureShellFE.h"
#include "configuratorsShellFE.h"
#include "FEOps.h"
#include "FETools.h"

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
typedef shellFE::TriangMesh < DataTypeContainerShellFE, TriangleType > MeshType;
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
void getDataEnergy ( RealType& dataEnergy, const VectorNd& dataFunction, const int& N ) {
  for ( int i = 0; i < N - 2; ++i ) {
    dataEnergy += abs ( 0.5 * (abs ( dataFunction ( i + 1 ) ) + abs ( dataFunction ( i ) )) * static_cast < RealType > ( 1 ) / static_cast < RealType > ( N - 1 ) );
  }
}

template < typename RealType, typename VectorNd >
void getBVEnergy ( RealType& BVEnergy, const VectorNd& thresholdedSolution, const int& N ) {
  for ( int i = 0; i < N + 2; ++i ) {
    BVEnergy += abs ( thresholdedSolution ( i + 1 ) - thresholdedSolution ( i ) );
  }
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

  void saveFDImages ( const VectorNd& PrimalSolution, const VectorNd& DualSolution, const VectorNd& Disparity ) const {
    ImageType outputImagePrimalSolution ( _N, _N );
    vectorToImage < VectorNd, ImageType, RealType > ( PrimalSolution, outputImagePrimalSolution, _N );
    scaleToFull < RealType, ImageType > ( outputImagePrimalSolution );
    const std::string primalSave = _saveDirectory + "primalSolution.bmp";
    saveBitmap ( primalSave, outputImagePrimalSolution );

    ImageType outputImageDisparity ( _N, _N );
    vectorToImage < VectorNd, ImageType, RealType > ( Disparity, outputImageDisparity, _N );
    scaleToFull < RealType, ImageType > ( outputImageDisparity );
    const std::string disparitySave = _saveDirectory + "disparity.bmp";
    saveBitmap ( disparitySave, outputImageDisparity );

    const RealType threshold = std::stod ( _pt.get < std::string > ( "General.threshold" ) );
    VectorNd thresholdedSolution ( _N );
    thresholding < RealType, VectorNd > ( PrimalSolution, thresholdedSolution, threshold );
    const std::string thresholdedSolutionSave = _saveDirectory + "thresholdedSolution.csv";
    safeSignal < RealType, VectorNd > ( thresholdedSolution, thresholdedSolutionSave );

    const std::string iniSave = _saveDirectory + "parameter.ini";
    boost::property_tree::write_ini ( iniSave, _pt );
  }

  void writeCSVData ( const VectorNd& PrimalSolution, const VectorNd& uL, const VectorNd& uR ) const {
    VectorNd thresholdedSolution ( _N );
    const RealType threshold = std::stod ( _pt.get < std::string > ( "General.threshold" ) );
    thresholding < RealType, VectorNd > ( PrimalSolution, thresholdedSolution, threshold );

    VectorNd uRShifted ( _N );
    getShiftedUR < RealType, VectorNd > ( uRShifted, uR, thresholdedSolution, _N );
    const std::string uRshiftedSave = _saveDirectory + "uRshifted.csv";
    safeSignal < RealType, VectorNd > ( uRShifted, uRshiftedSave );

    VectorNd dataFunction = uL - uRShifted;
    const std::string dataFunctionSave = _saveDirectory + "dataFunction.csv";
    safeSignal < RealType, VectorNd > ( dataFunction, dataFunctionSave );

    RealType dataEnergy = static_cast < RealType > ( 0 );
    getDataEnergy ( dataEnergy, dataFunction, _N );

    RealType BVEnergy = static_cast < RealType > ( 0 );
    getBVEnergy ( BVEnergy, thresholdedSolution, _N );

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
};

int main () {
  try {
    //get Parameters_Chambolle
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini ( "../../../Input/config.ini", pt );

#ifdef _OPENMP
    omp_set_num_threads ( 8 );
#endif

    //load signals and create image function g
    VectorNd uL;
    VectorNd uR;
    VectorNd disparity;
    const std::string fileName = pt.get < std::string > ( "Images.uL" ).c_str ();
    loadSignal < RealType, VectorNd > ( uL, fileName );
    loadSignal < RealType, VectorNd > ( uR, pt.get < std::string > ( "Images.uR" ) );
    const RealType lambda = std::stod ( pt.get < std::string > ( "General.lambda" ) );

    int N = uL.size ();
    convexOpt::createDisparity ( uL, uR, disparity, lambda );
    std::cout << "N: " << N << std::endl;

    //create discrete operators for gradient and the projections
    OperatorType Fd ( N );
    VectorNd primalSolution ( disparity.size () );
    VectorNd dualSolution ( 2 * primalSolution.size () );

    auto startTime = std::clock ();
    convexOpt::applyChambollePockAccelerated < RealType, VectorNd, convexOpt::ForwardFD < RealType, VectorNd > > ( primalSolution, dualSolution, disparity, N, Fd, pt );
    std::cout << "Duration: " << (std::clock () - startTime) / static_cast < double > ( CLOCKS_PER_SEC ) << std::endl;

    Saver < RealType, VectorNd > saver ( N, pt );
    saver.saveFDImages ( primalSolution, dualSolution, disparity );
    saver.writeCSVData ( primalSolution, uL, uR );

    MeshType grid ( pt.get < std::string > ( "Adaptive_Mesh.meshInput" ) );

    std::cerr << "0" << std::endl;
    for ( int iteration = 0; iteration < std::stoi ( pt.get < std::string > ( "General.numberRefinementSteps" ) ); ++iteration ) {
      ConfiguratorType config ( grid );
      VectorNd primalSolutionFE ( grid.getNumVertices () );
std::cerr << "1" << std::endl;
      ProjectFiniteDifferenceOntoFE < ConfiguratorType > projector ( config, N );
      projector.apply ( primalSolution, primalSolutionFE, std::stod ( pt.get < std::string > ( "Adaptive_Mesh.smoothingParameterProjection" ) ) );
std::cerr << "2" << std::endl;
      saver.plotFEFunction < ConfiguratorType, VectorNd > ( grid, primalSolutionFE, "primalFE" );
      std::cerr << "3" << std::endl;
    }




    /*
     MeshType grid ( pt.get < std::string > ( "Adaptive_Mesh.meshInput" ) );
     ConfiguratorType config ( grid );
     const ConfiguratorType::ElementType& El ( config.getInitializer ().getTriang ( 1 ) );
     //    ConfiguratorType::DomVecType dest;
     //    config.getGlobalCoords ( El, 0, dest );
     //    std::cout << dest << std::endl;

     const int NN = 3;

     const int numDofs = config.getNumLocalDofs ( );
     ConfiguratorType::DomVecType globalRef;


     //for ( int elementIdx = 0; elementIdx < config.getInitializer ().getNumTriangs (); ++elementIdx ) {

     for ( unsigned int i = 0; i < 3; ++i ) {
     //std::cout << El.getNode ( i );
     std::cout << El.getNode ( i ).segment ( 0, 2 ) << std::endl << std::endl;
     }
     ConfiguratorType::DomVecType refCoords;
     refCoords << 0,1;

     ConfiguratorType::Point3DType baryCentric;
     ConfiguratorType::DomVecType Dest;

     config.convertRefCoordsToBaryCentric ( refCoords, baryCentric);
     std::cout <<  std::endl << "baryCentric" << baryCentric << std::endl;

     config.getGlobalCoordsFromBarycentric ( El, baryCentric, Dest );
     std::cout <<  std::endl << "Dest" << Dest << std::endl;

     std::cout << std::endl << getGlobalIndex ( refCoords, NN ) << std::endl;
     */

  }
  catch ( std::exception& error ) {
    std::cerr << "exception caught: " << error.what () << '\n';
  }
  return 0;

}

