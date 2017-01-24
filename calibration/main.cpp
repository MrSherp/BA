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


typedef double RealType;
typedef Eigen::Array <RealType, Dynamic, Dynamic> ImageType;
typedef Matrix<RealType, Dynamic, 1> VectorNd;
typedef Eigen::Matrix <RealType, Dynamic, Dynamic> MatrixN;
typedef ForwardFD < RealType, VectorNd > OperatorType;
typedef ROFDataResolvent < RealType, VectorNd > ROFResolventDataTerm;
typedef KProjector < RealType, VectorNd > ProjectorOntoK;
typedef CProjector1 < RealType, VectorNd > ProjectorOntoC;

    
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
    RealType threshhold = std::stod( pt.get<std::string>("Parameters_Others.threshhold") );
    

    //load signals and create imagefunction g
    VectorNd uL;
    VectorNd uR;
    VectorNd G; 
    VectorNd Test;
    const char *filename1 = "plot1.csv";
    const char *filename2 = "plot2.csv";
    loadSignal < RealType, VectorNd >(uL, filename1);
    loadSignal < RealType, VectorNd >(uR, filename2);
    int N = uL.size();
    createG ( uL, uR, G, lambda );    
    std::cout << "N: " << N <<std::endl;
    

    //create discrete operators for gradient and the projections
    ForwardFD < RealType, VectorNd > Fd ( N );
    KProjector < RealType, VectorNd > K ( G );
    CProjector1 < RealType, VectorNd > C ( N );   
    ROFDataResolvent < RealType, VectorNd > Resolvent ( G, lambda );    
    
    
    //
    ImageType outimageSol1;
    ImageType outimageSol2;   
    ImageType outimageG;
    ImageType outimageUL;    
    VectorNd v ( G.size() );
    v.setConstant ( 0.); 
    VectorNd w ( v );
    VectorNd phi ( 2 * v.size() );
    phi.setConstant(1); 
    VectorNd psi ( phi );
    
    
    //start clock for and Algo1
    std::clock_t start1;
    double duration1;
    start1 = std::clock();   
    ChambollePockAlgorithm1 ( v, phi, N, Fd, C, K, tau, sigma, maxIter, stopEpsilon, endIter, endEpsilon );    
    std::cout << std::endl << "endIter: " << endIter << "\nendEpsilon: " << endEpsilon << std::endl;

    
    //save result and stop clock
    outimageSol1.resize ( N, N);
    outimageSol1 = vectorToImage < VectorNd, ImageType, RealType >( v , N );
    scaleToFull < RealType, ImageType > ( outimageSol1 );
    string name1 = "solution1.bmp";
    saveBitmap(name1, outimageSol1);        
    duration1 = ( std::clock() - start1 ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Duration1: "<< duration1 <<'\n';
    
    
    //start clock for and Algo2
    std::clock_t start2;
    double duration2;
    start2 = std::clock();        
    ChambollePockAlgorithm2 ( w, psi, N, Fd, C, K, gamma, tau, sigma, maxIter, stopEpsilon, endIter, endEpsilon );
    std::cout << "\nendIter: " << endIter << "\nendEpsilon: " << endEpsilon << std::endl;
        
    
    //save result and stop clock
    outimageSol2.resize ( N, N);
    outimageSol2 = vectorToImage < VectorNd, ImageType, RealType >( w , N );
    scaleToFull < RealType, ImageType > ( outimageSol2 );  
    string name2 = "solution2.bmp"; 
    saveBitmap(name2, outimageSol2);     
    duration2 = ( std::clock() - start2 ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Duration2: "<< duration2 <<'\n';
    
 
    // save g image
    outimageG.resize ( N, N);
    outimageG = vectorToImage < VectorNd, ImageType, RealType >( G , N );
    scaleToFull < RealType, ImageType > ( outimageG );  
    string name4 = "G.bmp"; 
    saveBitmap(name4, outimageG); 
    
    
    //save result as csv
    VectorNd result (N);
    const char *filename3 = "result.csv"; 
    threshholding < RealType, VectorNd > ( v, result, threshhold );
    safeSignal < RealType, VectorNd > ( result, filename3 );
    /*
     */
}




