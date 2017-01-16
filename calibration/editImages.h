#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <imageLoadSave.h>

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
    
    output.resize(2*xMax,2*yMax);
    
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

