#ifndef __MESHWITHDATA_H
#define __MESHWITHDATA_H

#include <bzipiostream.h>

namespace shellFE {

//! data can either belong to vertices or to faces
enum DataSupp { VERTEX_DATA, FACE_DATA };
//! vector-valued data can be saved (in VTK legacy format) as 3-vectors, normals or texture coordinates (the file format also supports color scalars and lookup tables, which we will not use).
enum VectorSpec { VECTORS, NORMALS };

//! Use this class like:  MeshWithData<> ( mesh ) -> .addData ( result, "color", VERTEX_DATA ) ->  .saveAsLegacyVTK ( "result.vtk" );
//! This class is a container for a mesh plus data vectors on vertices and faces.  We will not store any of the data vectors here, but only keep pointers to them.
template <class MeshType>
class MeshWithData {
public:
  typedef typename MeshType::RealType RealType;
  typedef typename MeshType::Point3DType Point3DType;
  typedef typename MeshType::Indices3DType Indices3DType;
  typedef typename MeshType::VectorType VectorType;

  MeshWithData ( const MeshType & mesh ) : _mesh ( mesh ), _precision( 8 ) {}
  MeshWithData ( const MeshType & mesh, int precision ) : _mesh ( mesh ), _precision( precision ) {}

  MeshWithData & addScalarData ( const VectorType & data, string dataDescr, DataSupp supp ) {
    
    ScalarData entry = { dataDescr, &data };

    switch ( supp ) {
    case VERTEX_DATA:
      if ( data.size() != static_cast<unsigned>( _mesh.getNumVertices () ) ) throw std::invalid_argument( aol::strprintf ( "Wrong size. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
      _scalarVertexData.push_back ( entry );
      break;

    case FACE_DATA:
      if ( data.size() != static_cast<unsigned>( _mesh.getNumTriangs () ) ) throw std::invalid_argument( aol::strprintf ( "Wrong size. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
      _scalarFaceData.push_back ( entry );
      break;
      
    default:
      throw std::invalid_argument( aol::strprintf ( "Unknown DataSupp. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
    return *this;
  }

  MeshWithData& addVectorData ( const VectorType & data, const int numComponents, string dataDescr, DataSupp supp, VectorSpec vSpec = VECTORS ) {
    VectorData entry = { dataDescr, vSpec, &data, numComponents };
    switch ( supp ) {
    case VERTEX_DATA:
      _vectorVertexData.push_back ( entry );
      break;
    case FACE_DATA:
      _vectorFaceData.push_back ( entry );
      break;
    default:
      throw std::invalid_argument( aol::strprintf ( "Unknown DataSupp. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
    return *this;
  }

  //! set precision in saving methods
  void setPrecisionTo( int prec ){ _precision = prec; }

  void saveAsLegacyVTK ( string filename ) const {

    const int numVertices = _mesh.getNumVertices();
      
    aol::Bzipofstream out ( filename.c_str() );                         // ctor checks if file could be opened
    //TODO same as std::ofstream out ( filename.c_str() );    

    out << "# vtk DataFile Version 2.0" << endl
        << "written by method MeshWithData::saveAsLegacyVTK" << endl
        << "ASCII" << endl
        << "DATASET POLYDATA" << endl
        << "POINTS " << numVertices << " float" << endl;

    // vertex coordinates
    for ( int nodeIter = 0; nodeIter < _mesh.getNumVertices(); ++nodeIter ) {
      const Point3DType& coords ( _mesh.getVertex(nodeIter) );
      for ( short i = 0; i < 3; ++i )
        out << ( i == 0 ? "" : " " ) << coords[i];
      out << endl;
    }

    out << "POLYGONS " << _mesh.getNumTriangs() << " " << 4 * _mesh.getNumTriangs() << endl;
    // triangles' vertex indices
    for ( int elementIter = 0; elementIter < _mesh.getNumTriangs(); ++elementIter ) {
      out << "3 ";
      for ( short i = 0; i < 3; ++i )
        out << ( i == 0 ? "" : " " ) << _mesh.getTriangNodeIdx( elementIter, i );
      out << endl;
    }

    if ( _scalarVertexData.size() > 0 || _vectorVertexData.size() > 0 )
      out << "POINT_DATA " << numVertices << endl;
    
    // scalar data on vertices
    for (int i = 0; i < _scalarVertexData.size(); ++i) {
      out << "SCALARS " <<  _scalarVertexData[i]._descr << " float" << endl;
      out << "LOOKUP_TABLE default" << endl;
      for ( unsigned vx = 0; vx < (*_scalarVertexData[i]._data).size(); ++vx )
        out << (*_scalarVertexData[i]._data)[vx] << endl;
    }
    
    // vector data on vertices
    for (int i = 0; i < _vectorVertexData.size(); ++i) {
      string spec;
      if ( _vectorVertexData[i]._spec == VECTORS ) spec = "VECTORS";
      if ( _vectorVertexData[i]._spec == NORMALS ) spec = "NORMALS";
      out << spec << " " << _vectorVertexData[i]._descr << " float" << endl;
      for ( int vx = 0; vx < numVertices; ++vx ) {
        for ( int comp = 0; comp < _vectorVertexData[i]._numComponents; ++comp )
          out << ( comp == 0 ? "" : " " ) << (*_vectorVertexData[i]._data)(vx + comp * numVertices); //TODO index mapper
        out << endl;
      }
    }

    // scalar data on faces
    for( int i=0; i < _scalarFaceData.size(); ++i ) {
        out << "CELL_DATA " << _mesh.getNumTriangs() << endl;
        out << "SCALARS " << _scalarFaceData[i]._descr << " float" << endl;
        out << "LOOKUP_TABLE default" << endl;
        for ( int vx = 0; vx < (*_scalarFaceData[i]._data).size(); ++vx )
            out << (*_scalarFaceData[i]._data)[vx] << endl;
    }
    
    // vector data on faces
    for (int i = 0; i < _vectorFaceData.size(); ++i) {
      string spec;
      if ( _vectorFaceData[i]._spec == VECTORS ) spec = "VECTORS";
      if ( _vectorFaceData[i]._spec == NORMALS ) spec = "NORMALS";
      out << spec << " " << _vectorFaceData[i]._descr << " float" << endl;
      
      for ( int vx = 0; vx < numVertices; ++vx ) {
        for ( int comp = 0; comp < _vectorFaceData[i]._numComponents; ++comp )
          out << ( comp == 0 ? "" : " " ) << (*_vectorFaceData[i]._data)(vx + comp * numVertices); //TODO index mapper
        out << endl;
      }
    }
  }

protected:

  struct ScalarData {
    string             _descr;
    const VectorType * _data;
  };

  struct VectorData {
    string             _descr;
    VectorSpec         _spec;
    const VectorType * _data;
    const int          _numComponents;
  };

  std::vector<ScalarData> _scalarVertexData;
  std::vector<VectorData> _vectorVertexData;
  std::vector<ScalarData> _scalarFaceData;
  std::vector<VectorData> _vectorFaceData;

  const MeshType _mesh;
  int _precision;      
};

} // end of namespace

#endif
