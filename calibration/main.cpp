#include <iostream>
#include <cstdio>
#include <ctime>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <imageLoadSave.h>
#include <editImages.h>
#include <pockChambolleAlgorithm.h>
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


    
    
int main(int argc, char **argv)
{
    
 namespace po = boost::program_options;

// Declare the supported options.
boost::program_options::options_description desc("Allowed options");
desc.add_options()
    ("help", "produce help message")
    ("compression", po::value<int>(), "set compression level");

boost::program_options::variables_map vm;


    
    
   
  
   /* VectorNd m(2);
    m << 0,1;
    VectorNd n(6);
    n << 1,2,3,4,5,6;*/
   // std::cout << "v:\n" << partition(v, 1, 3) << std::endl;
    
  
    /*int xMax =1000;
    int yMax =1000;
    
    ImageType inputImage (xMax,yMax);
    ImageType outimage1 (xMax,yMax);
    ImageType outimage2 (xMax,yMax);
    ImageType outimage3 (xMax,yMax);   
    ImageType outimage4 (xMax,yMax);  
    
    loadBitmap("g1.bmp", inputImage);

    mirrorImage(inputImage, outimage1);    
    fourInOneSameSize(inputImage, outimage2);       
    copy(inputImage, outimage3);
    changeSize(outimage3, 255, 255);
    copy(inputImage, outimage4);
    invert(outimage4);
    
    string name1 = "test1.bmp";
    string name2 = "test2.bmp";
    string name3 = "test3.bmp";
    string name4 = "test4.tmp";
    //saveBitmap(name1, outimage1); 
    //saveBitmap(name2, outimage2);
    //saveBitmap(name3, outimage3);
    //saveBitmap(name4, outimage4); 
*/

    
 /*   ImageType image(5,5);                                                                                    
    image << 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 8, 4, 4, 4, 4, 1, 1, 1, 1, 1;
    std::cout << "image:\n" << image << std::endl;
    
    imageToVector < VectorNd, ImageType, RealType >( image );
    VectorNd v (image);
    std::cout << "v:\n" << v << std::endl;
    
    VectorNd phi(50);
    phi.setZero();
    ForwardFD < RealType, VectorNd > K ( 5 );
    VectorNd w (v);

    K.apply ( v, phi );
    VectorNd psi (phi); 
    std::cout << "phi:\n" << phi << std::endl;
    
    divideByMaximumOfOneOrNorm <RealType> ( phi );
    std::cout << "phi:\n" << phi << std::endl;      
    */
 
    ImageType inputImage;
    ImageType outimage1;
    ImageType outimage2;   
    ImageType outimage3;
    
    loadBitmap ( "a.bmp" , inputImage );
    scaleToZeroOne ( inputImage );
    int N = inputImage.rows();
    imageToVector < VectorNd, ImageType, RealType >( inputImage );
    
    VectorNd u ( inputImage );
    VectorNd v ( u.size() );
    v.setConstant ( 1.); 
    VectorNd w ( v );
    
    ForwardFD < RealType, VectorNd > K ( N );
    
    VectorNd phi ( 2 * v.size() );
    phi.setZero(); 
    VectorNd psi ( phi );
    
    RealType lambda = 0.3;
    const RealType gamma = 0.001;
    RealType tau = 0.01;
    RealType sigma = 0.03;//1. / ( 16. * tau * static_cast < RealType > ( co::sqr ( N ) ) );
    const int maxIter = 10000;
    const RealType stopEpsilon = 0.000001;
    int endIter = 0;
    RealType endEpsilon = 0;
    ROFDataResolvent < RealType, VectorNd > Resolvent ( u, lambda );    
    
    
   /* ChambollePockAlgorithm1 ( v, phi, K, Resolvent, tau, sigma, maxIter, stopEpsilon, endIter, endEpsilon );    
    std::cout << "v: " << v << "\nphi: \n" << phi << "\nendIter: " << endIter << "\nendEpsilon: " << endEpsilon << std::endl;
    
    ChambollePockAlgorithm2 ( w, psi, K, Resolvent, gamma, tau, sigma, maxIter, stopEpsilon, endIter, endEpsilon );
    std::cout << "w: " << w << "\npsi: \n" << psi << "\nendIter: " << endIter << "\nendEpsilon: " << endEpsilon << std::endl;
   
    */
   
    std::clock_t start1;
    double duration1;
    start1 = std::clock();
   
    ChambollePockAlgorithm1 ( v, phi, K, Resolvent, tau, sigma, maxIter, stopEpsilon, endIter, endEpsilon );    
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
    
    ChambollePockAlgorithm2 ( w, psi, K, Resolvent, gamma, tau, sigma, maxIter, stopEpsilon, endIter, endEpsilon );
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




