#ifndef EDITIMAGES_H
#define EDITIMAGES_H

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>

#include "imageLoadSave.h"

template < typename RealType, typename ImageType >
void scaleToFull ( ImageType& Arg ) {
  const RealType min = Arg.minCoeff ();
  const RealType max = Arg.maxCoeff ();
  Arg *= static_cast < RealType > ( 255 ) / ( max - min );
  Arg -= ( static_cast < RealType > ( 255 ) / ( max - min ) ) * min;
}

template < typename VectorNd, typename ImageType, typename RealType >
void vectorToImage ( const VectorNd& Arg, ImageType& Dest ) {
  const auto rows = Dest.rows ();
  const auto columns = Dest.cols ();
  if ( rows * columns != Arg.size () )
    throw std::runtime_error ( "Size mismatch!" );
  for ( int i = 0; i < columns; ++i ) {
    Dest.col ( i ) = Arg.segment ( i * rows, rows );
  }
}

#endif
