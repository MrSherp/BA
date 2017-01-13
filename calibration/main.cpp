#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <imageLoadSave.h>
#include <editImages.h>
#include <pockAlgorithm.h>

using namespace Eigen;
using namespace std;

typedef double RealType;
typedef Eigen::Array <RealType, Dynamic, Dynamic> ImageType;
typedef Matrix<RealType, Dynamic, 1> VectorNd;
typedef Eigen::Matrix <RealType, Dynamic, Dynamic> MatrixN;


int main()
{
  //MatrixXd m(2,2);m(0,0) = 3;m(1,0) = 2.5;m(0,1) = -1;m(1,1) = m(1,0) + m(0,1);std::cout << "Here is the matrix m:\n" << m << std::endl;VectorXd v(2);v(0) = 4;v(1) = v(0) - 1;std::cout << "Here is the vector v:\n" << v << std::endl;
  
    int N =200;
    RealType h =1./N;
    cout << h<<std::endl;
  
    int xMax =1000;
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
    
    VectorNd v(9);
    v << 0, 0, 0, 5, 4, 3, 1, 1, 1;
    VectorNd phi(18);
    phi.setZero();
    
    std::cout << "v:\n" << v <<  std::endl;
  
  
    ForwardFD < RealType, VectorNd > K (  5  );
    
    K.apply(v, phi);
    
    std::cout << "phi:\n" << phi <<std::endl;
}

