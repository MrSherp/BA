#ifndef GABOR_H
#define GABOR_H

#include <core.h>

template < typename RealType >
class GaborFnc {
  RealType _theta;
  RealType _lambda;
  RealType _gamma;
  RealType _sigma;

public:
  GaborFnc ( const RealType Theta = 0.0, const RealType Lambda = 0.1, const RealType Gamma = 1.0, const RealType Sigma = 1.0 )
  : _theta ( Theta ), _lambda ( Lambda ), _gamma ( Gamma ), _sigma ( Sigma ) { }

  RealType operator ( ) ( const RealType x, const RealType y ) const {
    const RealType xx = x * cos ( _theta ) +  y * sin ( _theta );
    const RealType yy = - x * sin( _theta ) + y * cos ( _theta );
    return exp ( - ( co::sqr ( xx ) + co::sqr ( _gamma *yy ) ) / ( static_cast < RealType > ( 2 ) * co::sqr ( _sigma ) ) ) * cos ( static_cast < RealType > ( 2 ) * co::Core::Pi * xx / _lambda );      
  }
};


#endif