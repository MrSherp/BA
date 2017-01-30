#include <iostream>
#include <cstdio>
#include <ctime>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <imageLoadSave.h>
#include <editImages.h>
#include <pockChambolleAlgorithm.h>
#include <csvLoadSafe.h>
#include <math.h> 
#include <iomanip>
#include <cstdlib>
#include <stdlib.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/filesystem.hpp>


typedef double RealType;
typedef Eigen::Array <RealType, Dynamic, Dynamic> ImageType;
typedef Matrix<RealType, Dynamic, 1> VectorNd;
typedef Eigen::Matrix <RealType, Dynamic, Dynamic> MatrixN;
typedef ForwardFD < RealType, VectorNd > OperatorType;
typedef ROFDataResolvent < RealType, VectorNd > ROFResolventDataTerm;
typedef KProjector < RealType, VectorNd > ProjectorOntoK;
typedef CProjector2 < RealType, VectorNd > ProjectorOntoC;

    
int main()
{
    //get Parameters_Chambolle
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini("config.ini", pt);    
    RealType lambda = std::stod( pt.get<std::string>("Parameters_Chambolle.lambda") );
    const RealType gamma = std::stod( pt.get<std::string>("Parameters_Chambolle.gamma") );
    RealType tau = std::stod( pt.get<std::string>("Parameters_Chambolle.tau") );
    RealType sigma = std::stod( pt.get<std::string>("Parameters_Chambolle.sigma") );
    const int maxIter = std::stoi( pt.get<std::string>("Parameters_Chambolle.maxIter") );
    const RealType stopEpsilon = std::stod( pt.get<std::string>("Parameters_Chambolle.stopEpsilon") );
    int endIter = std::stoi( pt.get<std::string>("Parameters_Chambolle.endIter") );
    RealType endEpsilon = std::stod( pt.get<std::string>("Parameters_Chambolle.endEpsilon") );
    RealType threshold = std::stod( pt.get<std::string>("Parameters_Others.threshold") );
    
    
    
    //create directory for saving 
    string saveDirectory = pt.get<std::string>("Directory.safeDirectory");
    boost::filesystem::path directory (saveDirectory);
    boost::filesystem::create_directory(directory);
    boost::filesystem::copy_file("/home/staff/scherping/Build/" ,saveDirectory);
    
    

    //load signals and create imagefunction g
    VectorNd uL;
    VectorNd uR;
    VectorNd G; 
    const char *filename1 = "plot1.csv";
    const char *filename2 = "plot2.csv";
    loadSignal < RealType, VectorNd >(uL, filename1);
    loadSignal < RealType, VectorNd >(uR, filename2);
    int N = uL.size();
    createG ( uL, uR, G, lambda );    
    std::cout << "N: " << N <<std::endl;
    

    //create discrete operators for gradient and the projections
    ForwardFD < RealType, VectorNd > Fd ( N );
    KProjector < RealType, VectorNd > K ( G, lambda );
    CProjector2 < RealType, VectorNd > C ( N );   
    //ROFDataResolvent < RealType, VectorNd > Resolvent ( G, lambda );    
    
    
    //test mit N = 3
    
    /*VectorNd testv (9);
    VectorNd testimage (9);
    testimage << 0,0,0,0.2,0.8,0,0,0.6,0;
    VectorNd testphi (18);
    testphi.setZero();
    testv << -1,-2,-3,-4,-6,-8,-9,-12,-15;
    ForwardFD < RealType, VectorNd > testfd ( 3 );
    KProjector < RealType, VectorNd > testk ( testimage );
    CProjector1 < RealType, VectorNd > testc ( 3 );   
    int j = 0;
   
    VectorNd xBar ( testv ) ;    
    VectorNd oldPrimalSolution ( testv );                                  
    VectorNd adjointOfDual ( testv );
    VectorNd gradientOfXBar ( testphi );
    
    testfd.apply ( xBar, gradientOfXBar );
    std::cout << "gradientOfXBar: " << gradientOfXBar << std::endl;
    gradientOfXBar *= sigma;
    testphi += gradientOfXBar;
    std::cout << "testphi: " << testphi << std::endl;
    testk.projectOntoK ( testphi, j);
    std::cout << "testphi: " << testphi << std::endl;    
    xBar = testv;
    testfd.applyAdjoint ( testphi, adjointOfDual );
    std::cout << "adjointOfDual: " << adjointOfDual << std::endl;
    testv -= tau * adjointOfDual; 
    std::cout << "testv: " << testv << std::endl;
    testc.projectOntoC ( testv );
    std::cout << "testv: " << testv << std::endl;
    xBar *= -1.;                                                               
    xBar += 2. * testv;    
    std::cout << "xBar: " << xBar << std::endl;
    */ 
    
    
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
    
    */
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
    string name2 = "solution2.bmp"; 
    saveBitmap(name2, outimageSol2);     
    duration2 = ( std::clock() - start2 ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Duration2: "<< duration2 << std::endl;
    sprintf ( s, "TVout%05d.png", maxIter );
    saveBitmap(s, outimageSol2);
    /*
     */
 
    // save g image
    outimageG.resize ( N, N);
    outimageG = vectorToImage < VectorNd, ImageType, RealType >( G , N );
    scaleToFull < RealType, ImageType > ( outimageG );  
    string name4 = "G.bmp"; 
    saveBitmap(name4, outimageG); 
    
    
    //save result as csv
    VectorNd result (N);
    const char *filename3 = "result.csv"; 
    thresholding < RealType, VectorNd > ( w, result, threshold );
    safeSignal < RealType, VectorNd > ( result, filename3 );
    /*
     */
}




