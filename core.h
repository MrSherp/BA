#ifndef CORE_H
#define CORE_H

#include <sys/stat.h>
#include <math.h> 

namespace co {

// mathematical helper functions
template < typename RealType >
RealType sqr ( const RealType Arg ) {
  return Arg * Arg;
}

template < typename RealType >
RealType min ( const RealType A, const RealType B ) {
  if ( A > B )
    return B;
  else
    return A;
}

template < typename RealType >
RealType max ( const RealType A, const RealType B ) {
  if ( A < B )
    return B;
  else
    return A;
}


// core struct
struct Core {
  static constexpr long double Pi = 3.14159265358979323846;
    
  static void displayWarning ( const std::string& TextMassage ) {
    std::cout << "WARNING: " << TextMassage << std::endl << std::endl;   
  }
};

}

#endif
