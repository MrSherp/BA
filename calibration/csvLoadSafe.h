#ifndef CSVLOADSAVE_H
#define CSVLOADSAVE_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>

namespace co {

template < typename RealType, typename VectorNd >
  void loadSignal ( VectorNd& Destination, const std::string Filename ) {
    std::string line;
    std::ifstream file;
    file.open ( Filename );
    int N = 0;
    if ( file.is_open () ) {
      while ( getline ( file, line ) ) {
        ++N;
      }
    }
    else {
      throw std::runtime_error ( "unable to open csv file " + Filename );
    }
    file.close ();
    Destination.resize ( N );
    int i = 0;
    file.open ( Filename );
    while ( getline ( file, line ) ) {
      Destination ( i ) = std::stod ( line );
      ++i;
    }
    file.close ();

    RealType max = Destination.maxCoeff ();
    RealType min = Destination.minCoeff ();
    RealType scale = 1. / (max - min);
    for ( int i = 0; i < N; ++i ) {
      Destination ( i ) -= min;
      Destination ( i ) *= scale;
    }
  }

template < typename RealType, typename VectorNd >
  void safeSignal ( const VectorNd& Signal, const std::string Filename ) {

    std::ofstream myfile;
    myfile.open ( Filename.c_str () );
    int N = Signal.size ();
    for ( int i = 0; i < N; ++i ) {
      myfile << Signal ( i ) << std::endl;
    }
    myfile.close ();

  }

template < typename RealType, typename VectorNd >
  void thresholding ( const VectorNd& Image, VectorNd& Signal, const RealType threshhold ) {
    const int N = Signal.size ();
    for ( int i = 0; i < N; ++i ) {
#ifdef DEBUG
      if( threshhold > 1 || threshhold < 0 ) {
        std::cout << "Threshold not between 0 and 1!" << std::endl;
        break;
      }
#endif        
      int j = 0;
      while ( Image ( i + (N - 1 - j) * N ) >= threshhold ) {
        ++j;
      }
      Signal ( i ) = static_cast < RealType > ( j - 1 ) / static_cast < RealType > ( N - 1 );
    }
  }

}
  
#endif
