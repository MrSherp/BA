#ifndef EDITIMAGES_H
#define EDITIMAGES_H

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>

#include "imageLoadSave.h"

template < typename ImageType >
void scaleToZeroOne ( ImageType& Arg ) {
  Arg /= 255;
}

template < typename RealType, typename ImageType >
void scaleToFull ( ImageType& Arg ) {
  const RealType min = Arg.minCoeff ();
  const RealType max = Arg.maxCoeff ();
  Arg *= static_cast < RealType > ( 255 ) / ( max - min );
  Arg -= ( static_cast < RealType > ( 255 ) / ( max - min ) ) * min;
}


template < typename VectorNd, typename ImageType, typename RealType >
void imageToVector ( ImageType& Image ) {
  VectorNd _vector ( Image.cols () * Image.rows (), 1 );
  for ( int i = 0; i < Image.rows (); ++i ) {
    _vector.segment ( i * Image.cols (), Image.cols () ) = Image.row ( i );
  }
  Image.resize ( _vector.size (), 1 );
  Image = _vector;
}

template < typename VectorNd, typename ImageType, typename RealType >
void vectorToImage ( const VectorNd& Arg, ImageType& dest, const int Cols ) {
  ImageType _image ( Cols, Arg.size () / Cols );
  for ( int i = 0; i < Cols; ++i ) {
    _image.col ( i ) = Arg.segment ( i * _image.rows (), _image.rows () );
  }
}

#endif
