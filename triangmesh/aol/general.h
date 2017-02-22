#ifndef __GENERALLIB_H
#define __GENERALLIB_H

// C standard libraries
#include <cmath>
#include <cstdio>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <memory>

// STL
#include <complex>
#include <limits>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <algorithm>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <vector>

#include <typeinfo>

#include <chrono>
#include <unistd.h> //for sleep function

using namespace std;

#include <stdexcept>


//! platform dependent
#include <sys/stat.h>

#ifdef __GNUC__
// With this define, GCC_VERSION is 40301 for GCC 4.3.1 and can be conveniently checked (e.g. __WARNING_ON)
#define GCC_VERSION ( __GNUC__ * 10000 \
                      + __GNUC_MINOR__ * 100 \
                      + __GNUC_PATCHLEVEL__ )
#if (__GNUC__ >= 3) && !defined(_LIBCPP_VERSION)
#define HAVE_STDIO_FILEBUF 1
#endif
#endif

#ifdef HAVE_STDIO_FILEBUF
#include <ext/stdio_filebuf.h>
using  namespace __gnu_cxx;
#endif

#include <stdint.h>

// Mechanism to disable and enable warnings
// Usage example: WARNING_OFF ( old-style-cast ) and WARNING_ON ( old-style-cast )
// Note: Support for diagnostic pragmas was added in GCC 4.2.0.
#if ( defined ( __GNUC__ ) ) && ( GCC_VERSION >= 40200 ) && !( defined ( __INTEL_COMPILER ) ) && !( defined ( __NVCC__ ) )
#define __PRAGMA(P) _Pragma(#P)
#define __WARNING_OFF(WARN) __PRAGMA(GCC diagnostic ignored #WARN)
#define __WARNING_ON(WARN) __PRAGMA(GCC diagnostic warning #WARN)
#define WARNING_OFF(WARN) __WARNING_OFF(-W ## WARN)
#define WARNING_ON(WARN) __WARNING_ON(-W ## WARN)
#else
#define WARNING_OFF(WARN)
#define WARNING_ON(WARN)
#endif


namespace aol {
    
template <typename _DomainType, typename _RangeType = _DomainType>
class Op {

public:
  typedef _DomainType DomainType;
  typedef _RangeType  RangeType;

  Op() { }

  // Destroy polymorphic Ops correctly, important!
  virtual ~Op () {}

  virtual void apply ( const DomainType &Arg, RangeType &Dest ) const = 0;
};
    


template <typename _DomainType, typename _RangeType = _DomainType>
class PreparedOp {

public:
  typedef _DomainType DomainType;
  typedef _RangeType  RangeType;

  PreparedOp() { }

  // Destroy polymorphic Ops correctly, important!
  virtual ~PreparedOp () {}

  virtual void apply ( const DomainType &Arg, RangeType &Dest ) const = 0;
  virtual void prepare ( const DomainType &Arg ) const = 0;
  virtual void applyIntoDirection ( const DomainType &Arg, DomainType &Dest ) const = 0;
};

//! table of ansi color codes
namespace color {
const string reset       = "\033[0;0m";
const string invert      = "\033[0;7m";
const string black       = "\033[0;30m";
const string red         = "\033[0;31m";
const string green       = "\033[0;32m";
const string brown       = "\033[0;33m";
const string blue        = "\033[0;34m";
const string purple      = "\033[0;35m";
const string cyan        = "\033[0;36m";
const string light_grey  = "\033[0;37m";
const string dark_grey   = "\033[1;30m";
const string light_red   = "\033[1;31m";
const string light_green = "\033[1;32m";
const string yellow      = "\033[1;33m";
const string light_blue  = "\033[1;34m";
const string pink        = "\033[1;35m";
const string light_cyan  = "\033[1;36m";
const string white       = "\033[1;37m";
const string beep        = "\007";
const string error       = beep + red;
const string ok          = green;
const string residuum    = blue;
}

// Returns minimum/maximium
template<class T> inline T Min ( const T a, const T b ) { return ( ( a < b ) ? a : b ); }
template<class T> inline T Max ( const T a, const T b ) { return ( ( a < b ) ? b : a ); }
// Returns Value clamped into [Min,Max].
template<class T> inline T Clamp ( const T Value, const T Min, const T Max ) { return ( aol::Max ( aol::Min ( Value, Max ), Min ) ); }

template<class T> inline T Sqr (const T a) { return a * a; }
template<class T> inline T Cub (const T a) { return a * a * a; }

//! Signum function template, signum ( 0 ) = 0
template<class T> inline T signum ( const T x ) {
  if (x > 0.0) return 1;
  else if (x == 0.0) return 0;
  else if (x < 0.0) return -1;
  return x; // NaN
}
// int countDigitsOfNumber ( const int N ) {  return ( N != 0 ) ? static_cast<int>(floor(log10(static_cast<double>(std::abs(N)))))+1 : 1; }

//! Give back formatted string, analogously to sprintf, but save the long way 'round with char arrays.
string strprintf(const char * format, ...);

void makeDirectory ( const char *DirectoryName, bool verbose = true );
bool fileExists ( std::string filename );



} // namespace aol

#endif

