#include <iostream>
#include <cstdio>
#include <ctime>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <imageLoadSave.h>
#include <editImages.h>
#include <pockChambolleAlgorithm.h>
#include <math.h> 

typedef double RealType;
typedef Eigen::Array <RealType, Dynamic, Dynamic> ImageType;
typedef Matrix<RealType, Dynamic, 1> VectorNd;
typedef Eigen::Matrix <RealType, Dynamic, Dynamic> MatrixN;
typedef ForwardFD < RealType, VectorNd > OperatorType;
typedef ROFDataResolvent < RealType, VectorNd > ResolventDataTerm;

int main()
{
    std::clock_t start;
    double duration;
    start = std::clock();
  
    /*int xMax =1000;
    int yMax =1000;
    
    ImageType inimage (xMax,yMax);
    ImageType outimage1 (xMax,yMax);
    ImageType outimage2 (xMax,yMax);
    ImageType outimage3 (xMax,yMax);   
    ImageType outimage4 (xMax,yMax);  
    
    loadBitmap("g1.bmp", inimage);

    mirrorImage(inimage, outimage1);    
    fourInOneSameSize(inimage, outimage2);       
    copy(inimage, outimage3);
    changeSize(outimage3, 255, 255);
    copy(inimage, outimage4);
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
 
    ImageType inimage ( 1000, 1000 );
    ImageType outimage1 ( 1000, 1000 );
    ImageType outimage2 ( 1000, 1000 );   

    
    loadBitmap ( "g1.bmp" , inimage );
    changeSize ( inimage, 50, 50 );
    ImageType outimage3 ( 1000, 1000 );
    copy ( inimage, outimage3);
    int N = inimage.rows();
    imageToVector < VectorNd, ImageType, RealType >( inimage );
    VectorNd u ( inimage );
    VectorNd v ( u.size() );
    v.setZero();
    VectorNd w ( v );
    ForwardFD < RealType, VectorNd > K ( N );
    VectorNd phi ( 2 * v.size() );
    phi.setZero(); 
    VectorNd psi ( phi );
    
    RealType lambda = 0.5;
    const RealType gamma = 0.1;
    RealType tau = 0.04;
    RealType sigma = 0.06;
    const int maxIter = 100000;
    const RealType stopEpsilon = 0.00001;
    int endIter = 0;
    RealType endEpsilon = 0;
    ROFDataResolvent < RealType, VectorNd > Resolvent ( u, lambda );    
    
    
   /* ChambollePockAlgorithm1 ( v, phi, K, Resolvent, tau, sigma, maxIter, stopEpsilon, endIter, endEpsilon );    
    std::cout << "v: " << v << "\nphi: \n" << phi << "\nendIter: " << endIter << "\nendEpsilon: " << endEpsilon << std::endl;
    
    ChambollePockAlgorithm2 ( w, psi, K, Resolvent, gamma, tau, sigma, maxIter, stopEpsilon, endIter, endEpsilon );
    std::cout << "w: " << w << "\npsi: \n" << psi << "\nendIter: " << endIter << "\nendEpsilon: " << endEpsilon << std::endl;
   
    */
   
    ChambollePockAlgorithm1 ( v, phi, K, Resolvent, tau, sigma, maxIter, stopEpsilon, endIter, endEpsilon );    
    std::cout << "\nendIter: " << endIter << "\nendEpsilon: " << endEpsilon << std::endl;
    
    ChambollePockAlgorithm2 ( w, psi, K, Resolvent, gamma, tau, sigma, maxIter, stopEpsilon, endIter, endEpsilon );
    std::cout << "\nendIter: " << endIter << "\nendEpsilon: " << endEpsilon << std::endl;
    
    outimage1.resize ( N, N);
    outimage2.resize ( N, N);
    outimage1 = vectorToImage < VectorNd, ImageType, RealType >( v , N );
    outimage2 = vectorToImage < VectorNd, ImageType, RealType >( w , N );   
    
    changeSize ( outimage1, 1000, 1000);
    changeSize ( outimage2, 1000, 1000);
    changeSize ( outimage3, 1000, 1000);
    
    string name1 = "test1.bmp";
    string name2 = "test2.bmp";
    string name3 = "test3.bmp";   
    
    saveBitmap(name1, outimage1); 
    saveBitmap(name2, outimage2);
    saveBitmap(name3, outimage3);   
    
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"printf: "<< duration <<'\n';
    
}





