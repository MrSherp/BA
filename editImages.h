#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>

#include "imageLoadSave.h"

using namespace Eigen;
using namespace std;




template < typename ImageType>
void loadBitmap(const char* Filename, ImageType& Image){    
    loadColorBitmap(Filename, Image, Image, Image);
}



template < typename ImageType>
void copy(const ImageType& input, ImageType& output){    
    const unsigned int xMax = input.rows();
    const unsigned int yMax = input.cols();

    output.resize(xMax, yMax);
    
    for(int y=0; y < yMax; ++y) {
        
        for(int x=0; x < xMax; ++x) {
            
            output(x,y) = input(x,y);   
        }
    }
}



template < typename ImageType >
void scaleToZeroOne ( ImageType& Arg ){
    Arg /= 255; 
}



template < typename RealType, typename ImageType >
void scaleToFull ( ImageType& Arg ){
    RealType min = Arg.minCoeff();
    if( min < 0 )
        Arg -= min;
    RealType max = Arg.maxCoeff();
    RealType scale = 1./max;
    if( max > 1 )
        Arg *= scale;
    Arg *= 255;
}



template < typename ImageType>
void mirrorImage(const ImageType& input, ImageType& output) {    
    const unsigned int xMax = input.rows();
    const unsigned int yMax = input.cols();
    output.resize(xMax, yMax);
    
    for(int y=0; y < yMax; ++y) {        
        for(int x=0; x < xMax; ++x) {            
            output(x,y) = input(xMax-1-x,y);    //static_cast < int > ( sqrt ( co::sqr ( input( xMax-x-1, y ) )));
        }
    }
} 



template <typename ImageType>
void changeSize(ImageType& Image, int X, int Y){    
    const unsigned int xMax = Image.rows();
    const unsigned int yMax = Image.cols();
    
    if(X != xMax || Y != yMax){        
        ImageType tmp(X,Y);        
        for(double y=0; y < Y; y+=1.){            
            for(double x=0; x < X; x+=1.){                
                tmp(x,y)=Image(x*xMax/X,y*yMax/Y);
            }
        }        
        Image.resize(X,Y);        
        copy(tmp, Image);
    }
}
    

    
template <typename ImageType>
void fourInOneSameSize(const ImageType& input, ImageType& output) { 
    const unsigned int xMax = input.rows();
    const unsigned int yMax = input.cols();    
    output.resize( 2 * xMax, 2 * yMax);
    
    for(int y=0; y < yMax; ++y) {        
        for(int x=0; x < xMax; ++x) {            
            output(x,y) = input(x,y);    
            output(x+xMax,y) = input(x,y);    
            output(x,y+yMax) = input(x,y);    
            output(x+xMax,y+yMax) = input(x,y);                
        }
    }
    changeSize(output,xMax,yMax);
}



template <typename ImageType>
void invert(ImageType& input){    
    const unsigned int xMax = input.rows();
    const unsigned int yMax = input.cols();    
    ImageType tmp(xMax,yMax);
    
    for(int y=0; y < yMax; ++y) {        
        for(int x=0; x < xMax; ++x) {            
            tmp(x,y) = 255.-input(x,y);   
        }
    }
    
    copy(tmp, input);
}



template < typename VectorNd, typename ImageType, typename RealType >
void imageToVector ( ImageType& Image){
    VectorNd _vector( Image.cols() * Image.rows(), 1 );
    for ( int i = 0; i < Image.rows(); ++i ){
        _vector.segment( i * Image.cols(), Image.cols() ) = Image.row( i ); 
    }
    Image.resize( _vector.size(), 1);
    Image = _vector;
}



template < typename VectorNd, typename ImageType, typename RealType >
ImageType vectorToImage ( VectorNd& Vector, const int Cols){
    ImageType _image( Cols, Vector.size() / Cols);
    for ( int i = 0; i < Cols; ++i){
        _image.col(i) = Vector.segment( i * _image.rows(), _image.rows() );
    }
    return _image;
}


