#ifndef CORE_H
#define CORE_H

#include <sys/stat.h>
#include <math.h> 

namespace co {

enum IMAGE_SCALE_MODE {
  SCALE_TO_FULL = 0, MULTIPLY_TO_FULL = 1   
};
    
namespace color {
const std::string reset       = "\033[0;0m";
const std::string invert      = "\033[0;7m";
const std::string black       = "\033[0;30m";
const std::string red         = "\033[0;31m";
const std::string green       = "\033[0;32m";
const std::string brown       = "\033[0;33m";
const std::string blue        = "\033[0;34m";
const std::string purple      = "\033[0;35m";
const std::string cyan        = "\033[0;36m";
const std::string light_grey  = "\033[0;37m";
const std::string dark_grey   = "\033[1;30m";
const std::string light_red   = "\033[1;31m";
const std::string light_green = "\033[1;32m";
const std::string yellow      = "\033[1;33m";
const std::string light_blue  = "\033[1;34m";
const std::string pink        = "\033[1;35m";
const std::string light_cyan  = "\033[1;36m";
const std::string white       = "\033[1;37m";
const std::string beep        = "\007";
const std::string error       = beep + red;
const std::string ok          = green;
const std::string residuum    = blue;
}

// mathematical helper functions
template < typename RealType >
RealType sqr ( const RealType Arg ) {
  return Arg * Arg;
}

template < typename RealType >
RealType min ( const RealType Arg1, const RealType Arg2 ) {
  if ( Arg1 > Arg2 )
    return Arg2;
  else
    return Arg1;
}

template < typename RealType >
RealType max ( const RealType Arg1, const RealType Arg2 ) {
  if ( Arg1 < Arg2 )
    return Arg2;
  else
    return Arg1;
}

template < typename RealType >
bool approxEqual ( const RealType Arg1, const RealType Arg2, const RealType Tolerance = 1e-12 ) {
  return ( std::abs ( Arg1 - Arg2 ) < Tolerance ) ? true : false;   
}

static constexpr long double Pi = 3.14159265358979323846;

void displayWarning ( const std::string& TextMassage ) {
  std::cout << color::red << "WARNING: " << color::reset << TextMassage << std::endl << std::endl;
}


}

#endif
