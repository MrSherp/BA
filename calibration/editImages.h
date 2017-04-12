#ifndef EDITIMAGES_H
#define EDITIMAGES_H

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>

#include "imageLoadSave.h"

namespace co {

template < typename RealType, typename ImageType >
void scaleToFull ( ImageType& Arg ) {
  const RealType min = Arg.minCoeff ();
  const RealType max = Arg.maxCoeff ();
  Arg *= static_cast < RealType > ( 255 ) / ( max - min );
  Arg -= ( static_cast < RealType > ( 255 ) / ( max - min ) ) * min;
}

template < typename RealType, typename VectorType, typename ImageType >
void vectorToImage ( const VectorType& Arg, ImageType& Dest ) {
  const auto rows = Dest.rows ();
  const auto columns = Dest.cols ();
  if ( rows * columns != Arg.size () )
    throw std::runtime_error ( "Size mismatch in vectorToImage!" );
  for ( int i = 0; i < rows; ++i ) {
    Dest.row ( i ) = Arg.segment ( i * columns, columns );
  }
}

template < typename RealType, typename ImageType, typename VectorType >
void imageToVector ( const ImageType& Arg, VectorType& Dest ) {
  const auto rows = Arg.rows ();
  const auto columns = Arg.cols ();
  if ( rows * columns != Dest.size () )
    throw std::runtime_error ( "Size mismatch in imageToVector!" );
  for ( int i = 0; i < rows; ++i ) {
    Dest.segment ( i * columns, columns ) = Arg.row ( i ) ;
  }
}

template < typename RealType, typename MeshType, typename VectorType >
void distributeElementwiseErrorToNodesOnUnitSquare ( const VectorType& Arg, VectorType& Dest, const MeshType& Mesh ) {
  if ( Arg.size () != Mesh.getNumTriangs () )
    throw std::runtime_error ( "Size mismatch in distributeElementwiseErrorToNodesOnUnitSquare!" );
  Dest.resize ( Mesh.getNumVertices () );
  VectorType multiplicityNode ( Dest.size ( ) );
  for ( int elementIdx = 0; elementIdx < Mesh.getNumTriangs (); ++elementIdx ) {
    const typename MeshType::ElementType& El ( Mesh.getTriang ( elementIdx ) );
    const auto globalElementIndex = El.getGlobalElementIdx ( );
    for ( unsigned short int i = 0; i < 3; ++i ) {
      const auto globalNodeIndex = El.getGlobalNodeIdx ( i );
      Dest[globalNodeIndex] = Arg[globalElementIndex];
      multiplicityNode[globalNodeIndex] =+ static_cast < RealType > ( 1 );
    }
  }
  for ( int i = 0; i < Dest.size ( ); ++i ) {
    Dest[i] /= static_cast < RealType > ( multiplicityNode[i] );
  }
}

}

#endif
