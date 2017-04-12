#ifndef __FETOOLS_H_
#define __FETOOLS_H_

#include "meshWithData.h"
#include "core.h"

namespace co {


template < typename RealType, typename VectorType >
bool checkIfBoundaryNodeOnUnitSquare ( const VectorType& Arg ) {

    if ( co::approxEqual ( static_cast < RealType > ( 0 ), Arg[0] ) || co::approxEqual ( static_cast < RealType > ( 1 ), Arg[0] )
    || co::approxEqual ( static_cast < RealType > ( 0 ), Arg[1] ) || co::approxEqual ( static_cast < RealType > ( 1 ), Arg[1] ) )
    return true;
else
    return false;
}
    
template < typename RealType, typename MeshType, typename VectorType >
void createBoundaryMaskDirichletOnUnitSquare ( const MeshType& Mesh, std::vector < bool >& DirichletMask ) {
    const auto numVertices = Mesh.getNumVertices ( );
    DirichletMask.resize ( numVertices );
    for ( int i = 0; i < numVertices; ++i ) {
       DirichletMask[i] = checkIfBoundaryNodeOnUnitSquare < RealType > ( Mesh.getVertex ( i ) ); 
    }
}


}

#endif
