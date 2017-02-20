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

#include <configuratorsShellFE.h>
#include <quadratureShellFE.h>
#include <legacyVtkWriter.h>
#include <shellHandler.h>
#include <triangMesh.h>
#include <triangleShellFE.h>
#include <basefunctionSetShellFE.h>
#include <EigenDataContainer.h>
#include <aol.h>

#include "imageLoadSave.h"
#include "editImages.h"
#include "pockChambolleAlgorithm.h"
#include "csvLoadSafe.h"
#include "stiffnessMatrixIntegrator.h"
#include "massMatrixIntegrator.h"
#include "characteristicToFEMatrixIntegrator.h"


typedef double RealType;
typedef Eigen::Array <RealType, Dynamic, Dynamic> ImageType;
typedef Matrix<RealType, Dynamic, 1> VectorNd;
typedef Eigen::Matrix <RealType, Dynamic, Dynamic> MatrixN;
typedef ForwardFD < RealType, VectorNd > OperatorType;
typedef ROFDataResolvent < RealType, VectorNd > ROFResolventDataTerm;
typedef KProjector < RealType, VectorNd > ProjectorOntoK;
typedef CProjector2 < RealType, VectorNd > ProjectorOntoC;

typedef shellFE::ShellElementWithTangentSpaceAtVertex<DataTypeContainerShellFE> TriangleType;
typedef shellFE::TriangMesh<DataTypeContainerShellFE, TriangleType >            MeshType;
typedef shellFE::CenterOfMassQuadrature<RealType, typename DataTypeContainerShellFE::DomVecType> QuadType;
typedef shellFE::UnitTriangMeshConfiguratorP1<DataTypeContainerShellFE, MeshType, QuadType> ConfiguratorType;
typedef typename ConfiguratorType::VectorType VectorType;
typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;

    
int main()
{
    
    //get Parameters_Chambolle
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini("/home/staff/scherping/Input/config.ini", pt);    
    RealType lambda = std::stod( pt.get<std::string>("Parameters_Chambolle.lambda") );
    const RealType gamma = std::stod( pt.get<std::string>("Parameters_Chambolle.gamma") );
    RealType tau = std::stod( pt.get<std::string>("Parameters_Chambolle.tau") );
    RealType sigma = std::stod( pt.get<std::string>("Parameters_Chambolle.sigma") );
    const int maxIter = std::stoi( pt.get<std::string>("Parameters_Chambolle.maxIter") );
    const RealType stopEpsilon = std::stod( pt.get<std::string>("Parameters_Chambolle.stopEpsilon") );
    int endIter = std::stoi( pt.get<std::string>("Parameters_Chambolle.endIter") );
    RealType endEpsilon = std::stod( pt.get<std::string>("Parameters_Chambolle.endEpsilon") );
    RealType threshold = std::stod( pt.get<std::string>("Parameters_Others.threshold") );
    const char *plot1 = pt.get<std::string>("Images.uL").c_str();
    const char *plot2 = pt.get<std::string>("Images.uR").c_str();    

    
    
 /*   //create directory for saving 
    string saveDirectory = pt.get<std::string>("Directory.saveDirectory");
    boost::filesystem::path directory (saveDirectory);
    boost::filesystem::create_directory(directory);
    boost::filesystem::copy_file("/home/staff/scherping/Build/" ,saveDirectory);
    
    */
 

    //load signals and create imagefunction g
    VectorNd uL;
    VectorNd uR;
    VectorNd G; 
    loadSignal < RealType, VectorNd >(uL, plot1);
    loadSignal < RealType, VectorNd >(uR, plot2);
    int N = uL.size();
    createG ( uL, uR, G, lambda );    
    std::cout << "N: " << N <<std::endl;
    

    //create discrete operators for gradient and the projections
    ForwardFD < RealType, VectorNd > Fd ( N );
    KProjector < RealType, VectorNd > K ( G );
    CProjector2 < RealType, VectorNd > C ( N );   
    
    
    //create meshes and FE projectors
    MeshType mesh("/home/staff/scherping/Input/square.vtk");std::cout << "Ok" << std::endl;
    ConfiguratorType conf(mesh);std::cout << "Ok" << std::endl;
    shellFE::ShellHandler<ConfiguratorType> shellHandler(conf);std::cout << "Ok" << std::endl;
    SparseMatrixType systemMat(conf.getNumGlobalDofs(), conf.getNumGlobalDofs());std::cout << "Ok" << std::endl;
    MassMatrixIntegrator<ConfiguratorType> integrand(conf);std::cout << "Ok" << std::endl;
    integrand.assembleDirichlet<SparseMatrixType> (systemMat, shellHandler.getDirichletMask());std::cout << "Ok" << std::endl;
    
    std::cout << systemMat << std::endl;
    std::cout << "Ok" << std::endl;
    
    
    
    
    LegacyVtkWriter <MeshType> VWriter(mesh);
    
    
    //
    ImageType outimageSol1;
    ImageType outimageSol2;   
    ImageType outimageG;
    ImageType outimageUL;    
    VectorNd v ( G.size() );
    VectorNd phi ( 2 * v.size() );    
    startVectors < RealType, VectorNd >( v, phi, G, N );
    VectorNd w ( v );
    VectorNd psi ( phi );
    /*
    
    //start clock for and Algo1
    std::clock_t start1;
    double duration1;
    start1 = std::clock();   
    ChambollePockAlgorithm1 ( v, phi, N, Fd, C, K, tau, sigma, maxIter, stopEpsilon, endIter, endEpsilon );    
    std::cout << std::endl << "endIter: " << endIter << std::endl << "endEpsilon: " << endEpsilon << std::endl;

    
    //save result and stop clock
    outimageSol1.resize ( N, N);
    outimageSol1 = vectorToImage < VectorNd, ImageType, RealType >( v , N );
    scaleToFull < RealType, ImageType > ( outimageSol1 );
    string name1 = "solution1.bmp";
    saveBitmap(name1, outimageSol1);       
    duration1 = ( std::clock() - start1 ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Duration1: "<< duration1 << std::endl;
    
    *//*
    //start clock for and Algo2
    std::clock_t start2;
    double duration2;
    start2 = std::clock();        
    ChambollePockAlgorithm2 ( w, psi, N, Fd, C, K, gamma, tau, sigma, maxIter, stopEpsilon, endIter, endEpsilon );
    std::cout << std::endl << "endIter: " << endIter << std::endl << "endEpsilon: " << endEpsilon << std::endl;
        
    
    //save result and stop clock
    char s[20];
    outimageSol2.resize ( N, N);
    outimageSol2 = vectorToImage < VectorNd, ImageType, RealType >( w , N );
    scaleToFull < RealType, ImageType > ( outimageSol2 );  
    string solutionName = pt.get<std::string>("Parameters_Chambolle.lambda") + "_" + pt.get<std::string>("Parameters_Chambolle.gamma") + "_" + pt.get<std::string>("Parameters_Chambolle.tau") + "_" +
    pt.get<std::string>("Parameters_Chambolle.sigma") + "_" + pt.get<std::string>("Parameters_Chambolle.maxIter") + "_" + std::to_string(N);
    string solutionNameBmp = solutionName + ".bmp";
    saveBitmap(solutionNameBmp, outimageSol2);     
    duration2 = ( std::clock() - start2 ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Duration2: "<< duration2 << std::endl;
    //sprintf ( s, "TVout%05d.png", maxIter );
    //saveBitmap(s, outimageSol2);
    
     
 
    // save g image
    outimageG.resize ( N, N);
    outimageG = vectorToImage < VectorNd, ImageType, RealType >( G , N );
    scaleToFull < RealType, ImageType > ( outimageG );  
    string name4 = "G.bmp"; 
    saveBitmap(name4, outimageG); 
    
    
    //save result as csv
    VectorNd result (N);
    const string solutionNameCsv = solutionName + "Result.csv";
    const char *filename3 = solutionNameCsv.c_str();
    thresholding < RealType, VectorNd > ( w, result, threshold );
    safeSignal < RealType, VectorNd > ( result, filename3 );
    /*
     */
    
}




