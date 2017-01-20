#include <iostream>
#include <cstdio>
#include <ctime>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <imageLoadSave.h>
#include <editImages.h>
#include <pockChambolleAlgorithm.h>
#include <csvLoad.h>
#include <math.h> 
#include <iomanip>
#include <cstdlib>
#include <stdlib.h>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/config.hpp>
#include <boost/program_options/environment_iterator.hpp>
#include <boost/program_options/eof_iterator.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/option.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/version.hpp>



typedef double RealType;
typedef Eigen::Array <RealType, Dynamic, Dynamic> ImageType;
typedef Matrix<RealType, Dynamic, 1> VectorNd;
typedef Eigen::Matrix <RealType, Dynamic, Dynamic> MatrixN;
typedef ForwardFD < RealType, VectorNd > OperatorType;
typedef ROFDataResolvent < RealType, VectorNd > ROFResolventDataTerm;
typedef KProjector < RealType, VectorNd > ProjectorOntoK;

    
    
int main(int argc, char **argv)
{
    
    namespace po = boost::program_options;

// Declare the supported options.
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
    ("help", "produce help message")
    ("compression", po::value<int>(), "set compression level");

    boost::program_options::variables_map vm;

    VectorNd uL;
    VectorNd uR;
    const char *filename1 = "plot1.csv";
    const char *filename2 = "plot2.csv";
    loadSignal < RealType, VectorNd >(uL, filename1);
    loadSignal < RealType, VectorNd >(uR, filename2);
    
    int N = uL.size();
    std::cout << "N: " << N <<std::endl;
    signalToIndicatorFunction < VectorNd >( uL ); 
    signalToIndicatorFunction < VectorNd >( uR ); 
    VectorNd G (uL - uR);    

    ImageType outimage1;
    ImageType outimage2;   
    ImageType outimage3;
    
    VectorNd v ( G.size() );
    v.setConstant ( 1.); 
    VectorNd w ( v );
    
    ForwardFD < RealType, VectorNd > Fd ( N );
    KProjector < RealType, VectorNd > K ( G );
    
    VectorNd phi ( 2 * v.size() );
    phi.setZero(); 
    VectorNd psi ( phi );
    
    RealType lambda = 0.3;
    const RealType gamma = 0.001;
    RealType tau = 0.01;
    RealType sigma = 0.03;//1. / ( 16. * tau * static_cast < RealType > ( co::sqr ( N ) ) );
    const int maxIter = 10000;
    const RealType stopEpsilon = 0.0000001;
    int endIter = 0;
    RealType endEpsilon = 0;
    ROFDataResolvent < RealType, VectorNd > Resolvent ( G, lambda );    
    
   
    std::clock_t start1;
    double duration1;
    start1 = std::clock();
   
    ChambollePockAlgorithm1 ( v, phi, Fd, Resolvent, tau, sigma, maxIter, stopEpsilon, endIter, endEpsilon );    
    std::cout << "\nendIter: " << endIter << "\nendEpsilon: " << endEpsilon << std::endl;
    
    outimage1.resize ( N, N);
    outimage1 = vectorToImage < VectorNd, ImageType, RealType >( v , N );
    scaleToFull ( outimage1 );
    string name1 = "test1.bmp";
    saveBitmap(name1, outimage1);    
    
    duration1 = ( std::clock() - start1 ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Duration1: "<< duration1 <<'\n';
    
    
    std::clock_t start2;
    double duration2;
    start2 = std::clock();    
    
    ChambollePockAlgorithm2 ( w, psi, Fd, Resolvent, gamma, tau, sigma, maxIter, stopEpsilon, endIter, endEpsilon );
    std::cout << "\nendIter: " << endIter << "\nendEpsilon: " << endEpsilon << std::endl;
        
    outimage2.resize ( N, N);
    outimage2 = vectorToImage < VectorNd, ImageType, RealType >( w , N );
    scaleToFull ( outimage2 );  
    string name2 = "test2.bmp"; 
    saveBitmap(name2, outimage2); 
    
    duration2 = ( std::clock() - start2 ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Duration2: "<< duration2 <<'\n';
    /*
     */
    
}




