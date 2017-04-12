#ifndef __UNIFORMSIMPLEXMESH_H
#define __UNIFORMSIMPLEXMESH_H

namespace co {

template < typename DataTypeContainer >
class UniformSimplexElement {

  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::DomVecType DomVecType;
  typedef typename DataTypeContainer::Indices3DType Indices3DType;
  typedef std::vector < DomVecType > VertexIterator;

  const int _globalIndex;
  const Indices3DType _globalNodeIndex;
  const VertexIterator& _nodes;

public:
  UniformSimplexElement ( const int GlobalIndex, const Indices3DType& GlobalNodeIndex, const VertexIterator& Nodes )
: _globalIndex ( GlobalIndex ), _globalNodeIndex ( GlobalNodeIndex ), _nodes ( Nodes ) { }

  void printNodes () const {
    std::cout << "node0 = " << getNode ( 0 ).transpose () << std::endl << "node1 = " << getNode ( 1 ).transpose () << std::endl << "node2 = " << getNode ( 2 ).transpose () << std::endl;
  }

  void print () const {
    std::cout << std::endl << "Element " << _globalIndex << std::endl;
    printNodes ();
  }

  int getGlobalElementIdx () const {
    return _globalIndex;
  }

  const Indices3DType & getGlobalNodeIdx () const {
    return _globalNodeIndex;
  }

  int getGlobalNodeIdx ( int localIndex ) const {
    return _globalNodeIndex[localIndex];
  }

  const DomVecType& getNode ( int i ) const {
    return _nodes[_globalNodeIndex[i]];
  }
};

template < typename DataTypeContainer >
class UniformSimplexMesh {
public:
  typedef UniformSimplexElement < DataTypeContainer > ElementType;
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::DomVecType DomVecType;
  typedef typename DataTypeContainer::Indices3DType Indices3DType;

protected:
  const int _numX;
  const int _numY;
  std::vector < DomVecType > _vertexIterator;
  std::vector < ElementType > _triangIterator;

  void init () {
    _vertexIterator.reserve ( _numX * _numY );
    for ( int y = 0; y < _numY; ++y ) {
      for ( int x = 0; x < _numX; ++x ) {
        _vertexIterator.emplace_back ( static_cast < RealType > ( x ) / static_cast < RealType > ( _numX - 1 ), static_cast < RealType > ( y ) / static_cast < RealType > ( _numY - 1 ) );
      }
    }
    int currentIndex = 0;
    const int numTriangs = 2 * (_numX - 1) * (_numY - 1);
    _triangIterator.reserve ( numTriangs );
    for ( int y = 0; y < _numY - 1; ++y ) {
      for ( int x = 0; x < _numX - 1; ++x ) {
        //   UniformSimplexElement( const int GlobalIndex, const Indices3DType& GlobalNodeIndex, const VertexIterator& Nodes  )
        // lower left
        Indices3DType lowerLeftIndex, upperRightIndex;
        lowerLeftIndex << (y * _numX + x), (y * _numX + x + 1), ((y + 1) * _numX + x);
        upperRightIndex << (y * _numX + x + 1), ((y + 1) * _numX + x), ((y + 1) * _numX + x + 1);
        _triangIterator.emplace_back ( currentIndex, lowerLeftIndex, _vertexIterator );
        ++currentIndex;
        _triangIterator.emplace_back ( currentIndex, upperRightIndex, _vertexIterator );
        ++currentIndex;
      }
    }
  }

public:
  UniformSimplexMesh ( const int NumX, const int NumY )
: _numX ( NumX ), _numY ( NumY ), _vertexIterator (), _triangIterator () {
    init ();
  }

  int getNumX () const {
    return _numX;
  }

  int getNumY () const {
    return _numY;
  }

  int getNumVertices () const {
    return (_vertexIterator.size ());
  }
  int getNumTriangs () const {
    return (static_cast < int > ( _triangIterator.size () ));
  }

  const DomVecType& getVertex ( const int num ) const {
    return _vertexIterator[num];
  }

  const ElementType& getTriang ( const int num ) const {
    return _triangIterator[num];
  }

};

}         //end namespace

#endif
