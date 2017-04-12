#include <iostream>
#include <cstdio>
#include <ctime>
#include <math.h> 
#include <iomanip>
#include <cstdlib>
#include <stdlib.h>


#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/filesystem.hpp>
#include "core.h"

typedef double RealType;

int main () {
  try {
    //get Parameters_Chambolle
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini ( "../../../Input/config.ini", pt );
    
    const std::string filename1 = pt.get < std::string > ( "CSV.saveName1" );
    const std::string filename2 = pt.get < std::string > ( "CSV.saveName2" );

    if ( boost::filesystem::exists ( filename1 ) || boost::filesystem::exists ( filename2 ) )
      throw std::runtime_error ( "File already exists!" );
    std::ofstream out1 ( filename1 );
    std::ofstream out2 ( filename2 );

    const unsigned int N = pt.get < unsigned int > ( "CSV.N" );
    if ( N < 2 )
      throw std::runtime_error ( "CSV.N must be larger than 1!" );
    
    for ( unsigned int i = 0; i < N; ++i ) {
      const RealType h = static_cast < RealType > ( i ) / static_cast < RealType > ( N - 1 );
      
      // change here
      const RealType func1 = co::max ( static_cast < RealType > ( 0 ), - static_cast < RealType > ( 10 ) * co::sqr ( h - static_cast < RealType > ( 0.6 ) ) + static_cast < RealType > ( 0.5 ) );
      const RealType func2 = co::max ( static_cast < RealType > ( 0 ), - static_cast < RealType > ( 10 ) * co::sqr ( h - static_cast < RealType > ( 0.4 ) ) + static_cast < RealType > ( 0.5 ) ); 

      out1 << func1 << std::endl;
      out2 << func2 << std::endl;
    }    
    out1.close ( );
    out2.close ( );
}
  catch ( std::exception& error ) {
    std::cerr << "exception caught: " << error.what () << '\n';
  }
  return 0;
}
