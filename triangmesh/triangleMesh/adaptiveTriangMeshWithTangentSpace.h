#ifndef __ADAPTIVETRIANGMESHWITHTANGENTSPACE_H
#define __ADAPTIVETRIANGMESHWITHTANGENTSPACE_H

#include <general.h>
#include <triangMesh.h>

//! defines for local node and edge numbering
typedef short LocalIndex;
typedef int GlobalIndex;
const int IndexNotSet = -1;

template <typename DomVecType>
struct ParentInformation{ 
      int globalIndices[2]; 
      int ElementIndex;
      int EdgeNumber;
      DomVecType RefCoords;
};

namespace shellFE {

template< typename DataTypeContainer, typename TriangleType >
class AdaptiveTriangMesh : public TriangMesh<DataTypeContainer, TriangleType> {
  
public:   
  typedef TriangleType ElementType;
  typedef typename DataTypeContainer::RealType      RealType;
  typedef typename DataTypeContainer::DomVecType    DomVecType;
  typedef typename DataTypeContainer::Point3DType   Point3DType;
  typedef typename DataTypeContainer::Indices3DType Indices3DType;
  typedef typename DataTypeContainer::VectorType    VectorType;
  typedef typename DataTypeContainer::MaskType      MaskType;
  typedef std::vector< Indices3DType >              VertexIndicesType;
  typedef std::vector< Point3DType >                VertexCoordsType;
  
protected:  
// half-edge iterator (only for internal refinement!)
class DartIterator{
  
  typedef AdaptiveTriangMesh<DataTypeContainer, TriangleType> MeshType;
    
  protected:
    const MeshType & _mesh;
    GlobalIndex _triangle;
    LocalIndex  _node;
    LocalIndex  _edge;

  public:
    DartIterator( const MeshType & Mesh, GlobalIndex triangleIndex, LocalIndex localEdgeIndex, LocalIndex localNodeIndex ) : 
    _mesh( Mesh ), _triangle (triangleIndex), _node(localNodeIndex), _edge(localEdgeIndex){}    
    DartIterator( const MeshType & Mesh, GlobalIndex triangleIndex, LocalIndex localEdgeIndex ) : 
    _mesh( Mesh ), _triangle (triangleIndex), _node((localEdgeIndex + 1) % 3), _edge(localEdgeIndex){}

    void set ( const GlobalIndex triangle, const LocalIndex edge, const LocalIndex node )  { _triangle = triangle; _node = node; _edge = edge;}
    GlobalIndex getGlobalTriangleIndex() const { return _triangle; }
    GlobalIndex getGlobalNodeIndex()     const { return _mesh.getTriangNodeIdx( _triangle, _node );}
    LocalIndex  getLocalNodeIndex()      const { return _node; }
    LocalIndex  getLocalEdgeIndex()      const { return _edge;}
    LocalIndex  getNextNodeLocalIndex()  const { return 3 - (_node + _edge);}
    LocalIndex  getNextEdgeLocalIndex()  const { return 3 - (_node + _edge); }  
    GlobalIndex getNextNodeGlobalIndex() const { return _mesh.getTriangNodeIdx( _triangle, getNextNodeLocalIndex() );}
    GlobalIndex getNextTriangleIndex()   const { return _mesh.getNeighbour( _triangle, _edge ); }
 
    // returns local index of common node (first checking if "d" also refers to _triangle, but to a different node)
    template<typename DartIteratorType>
    LocalIndex getCommonNodeLocalIndex( const DartIteratorType & d ) const {
      if(_triangle != d.getGlobalTriangleIndex()) throw std::invalid_argument( aol::strprintf("T=%d, but this->tringle = %d. In File %s at line %d.", d.getGlobalTriangleIndex(), _triangle, __FILE__, __LINE__).c_str() );
      if ( _edge == d.getLocalEdgeIndex() ) throw std::invalid_argument( aol::strprintf("points to same edge=. In File %s at line %d.", __FILE__, __LINE__).c_str() );
      return 3 - ( getLocalEdgeIndex() + d.getLocalEdgeIndex() );
    }

    // returns global index of common node (also if "d" refers to a different triangle than _triangle )
    template<typename DartIteratorType>
    GlobalIndex getCommonNodeGlobalIndex( const DartIteratorType & d ) const {
      if ( getGlobalNodeIndex() == d.getGlobalNodeIndex() || getGlobalNodeIndex() == d.getNextNodeGlobalIndex() ) return getGlobalNodeIndex();
      if ( getNextNodeGlobalIndex() == d.getGlobalNodeIndex() || getNextNodeGlobalIndex() == d.getNextNodeGlobalIndex() ) return getNextNodeGlobalIndex();
      throw std::invalid_argument( aol::strprintf("Dart d and *this have no node in common. In File %s at line %d.", __FILE__, __LINE__).c_str() );
      return -1;
    }
    
    void print() const{ cerr << "triangle = " << _triangle << ", node = " << _node << ", edge = " << _edge << endl; }

    // does _triangle have a neighbour across _edge?
    bool canFlipTriangle() const { return this->getNextTriangleIndex() != IndexNotSet ; }

    // moves to neighbouring triangle across _edge
    void flipTriangle() {
          if( !canFlipTriangle() ) throw std::invalid_argument( aol::strprintf("Cannot flip triangle. In File %s at line %d.", __FILE__, __LINE__).c_str() );
        
          GlobalIndex newTriangleIndex = getNextTriangleIndex();
          GlobalIndex globalNodeIndex = _mesh.getTriangNodeIdx( _triangle, _node );    
          GlobalIndex globalOtherNodesIndex = _mesh.getTriangNodeIdx( _triangle, 3 - (_node + _edge) );

          // find out which node of the new triangle is ours
          LocalIndex newNodeIndex = IndexNotSet;
          for (LocalIndex i = 0; i < 3; ++i){
            if ( _mesh.getTriangNodeIdx( newTriangleIndex, i) == globalNodeIndex ){ newNodeIndex = i; break;}
          }
            
          // find out which node of the new triangle lies on our edge
          LocalIndex newOtherNodesIndex = IndexNotSet;
          for (LocalIndex i = 1; i < 3; ++i){
            if ( _mesh.getTriangNodeIdx( newTriangleIndex,(newNodeIndex + i) % 3 ) == globalOtherNodesIndex ){ newOtherNodesIndex = (newNodeIndex + i) % 3; break; }
          }

          if ( min( newNodeIndex, newOtherNodesIndex ) < 0 ) 
              throw std::invalid_argument( aol::strprintf( "newTriang should neighbour of _triangle but has no common node .In File %s at line %d.", __FILE__, __LINE__ ).c_str() );

          LocalIndex newEdgeIndex = (-(newNodeIndex + newOtherNodesIndex)) % 3;
          if (newEdgeIndex < 0) newEdgeIndex += 3;

          // all information found, now set:
          _triangle = newTriangleIndex; _node = newNodeIndex; _edge = newEdgeIndex; 
    }
    
    // moves to other node along _edge (inside current triangle)
    void flipNode() { _node = getNextNodeLocalIndex();}
    // moves to other edge with _node (inside current triangle)
    void flipEdge() {_edge = getNextEdgeLocalIndex();}

};
  

protected:  
  std::vector<bool> _markedForRefinement;
  std::map< int, ParentInformation<DomVecType> > _interpolationMap;
  
public:
  
  AdaptiveTriangMesh ( const string& fileName ) : 
      TriangMesh<DataTypeContainer, TriangleType >( fileName ), 
      _markedForRefinement( this->getNumTriangs(), false ) {
    this->makeNeighbour();
  }   

  const std::map< int, ParentInformation<DomVecType> > & getInterpolationMap ( ) const { return _interpolationMap;}
  
  // mark
  void mark( int element ){ _markedForRefinement[ element ] = true;}
  void markAll(){ for( unsigned i = 0; i < _markedForRefinement.size(); i++ ) _markedForRefinement[i] = true; }
  void unmark( int element ){ _markedForRefinement[ element ] = false;}
  void unmarkAll(){ for( unsigned i = 0; i < _markedForRefinement.size(); i++ ) _markedForRefinement[i] = false; }
  bool isMarkedForRefinement( int element ) const {return _markedForRefinement[element];}
  
  // when adding a new face we have to expand the array _markedForRefinement correspondingly
  int pushBackTriang ( const Indices3DType newTriang ) {
    _markedForRefinement.push_back( false );
    return TriangMesh<DataTypeContainer, TriangleType >::pushBackTriang( newTriang );
  }
  
  void refineMarkedTriangles( ) {
      const int numElements = this->getNumTriangs();
      for ( int elementIndex = 0; elementIndex < numElements; ++elementIndex  )
        if ( isMarkedForRefinement( elementIndex ) )  refine( elementIndex );
      this->makeOrientationConsistent();
      // generate tangent space and update projection coefficients TODO should be done only for new elements!
      //this->generateApproximativeTangentSpaceAtNodes();
      //for ( int elementIndex = 0; elementIndex < this->getNumTriangs(); ++elementIndex  )
      //this->_triangIterator[elementIndex].updateProjectionCoefficients( this->_tangentSpaceVec1, this->_tangentSpaceVec2 );
      const int numVertices = this->getNumVertices();
      this->_tangentSpaceVec1.resize ( numVertices ); this->_tangentSpaceVec2.resize ( numVertices ); this->_normalSpaceVec.resize ( numVertices );
  }
  
  // prolongate by linear interpolation
  void prolongateLinearly( VectorType& function ) const {
      VectorType oldFunction = function;
      int oldSize = function.size();
      function.resize( this->getNumVertices() );
      for( int i = 0; i < oldSize; ++i ) function[i] = oldFunction[i];
      typename std::map< int, ParentInformation<DomVecType> >::const_iterator iter;
      for( int i = oldSize; i < function.size(); ++i ){
        iter = _interpolationMap.find( i );
        if( iter == _interpolationMap.end() ) throw std::invalid_argument ( aol::strprintf ( "unknown vertex in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
        function[i] = 0.5 * ( function[iter->second.globalIndices[0]] + function[iter->second.globalIndices[1]] );
      }  
  }
  
public:
  // internal refinement functions
  GlobalIndex refine( GlobalIndex triangleToRefine ) {
    DartIterator d( *this, triangleToRefine, getLongestEdgeIndex( triangleToRefine ) );
    return refine( d );
  }
  
  GlobalIndex refine( const DartIterator& d) {
    // first cut the neighbour on edge l in two triangles:
    GlobalIndex P_new = addEdgeMidpoint( d ); 
    DartIterator d_prime = d;
    bool hasNeighbourOn_d_prime = d_prime.canFlipTriangle();
    GlobalIndex T_l_new = IndexNotSet;

    if ( hasNeighbourOn_d_prime ){
        d_prime.flipTriangle();
        T_l_new = refineOnlyThis(d_prime, P_new);
    }
    
    // create new triangle T_new and connect with neighbours
    GlobalIndex T_new = refineOnlyThis(d, P_new);

    // take care of neighbours
    if (hasNeighbourOn_d_prime){
        // set d_prime to edge which lies in T_l_new, opposite to edge d2.edge of T_new
        d_prime.flipNode();
        d_prime.flipEdge();
        // flipTriangle() here will always work, as we have just created this triangle. Thus, no canFlipTriangle() call is needed.
        d_prime.flipTriangle();
        d_prime.flipEdge();
        // set neighbour GlobalIndex in T_l_new
        this->setNeighbour ( d_prime.getGlobalTriangleIndex(), d_prime.getLocalEdgeIndex(), T_new );
        // move d1 to T_new
        DartIterator d1 = d; d1.flipNode(); d1.flipEdge();
        d1.flipTriangle();
        d1.flipEdge();
        // set neighbour in T_new
        this->setNeighbour ( d1.getGlobalTriangleIndex(), d1.getLocalEdgeIndex(), T_l_new );
    }

    return T_new;
  }
    
  GlobalIndex refineOnlyThis( const DartIterator& d, GlobalIndex midpoint) {
    // in comments for this function, we write T for *this.
    LocalIndex l = getLongestEdgeIndex( d.getGlobalTriangleIndex() );
    if (l != d.getLocalEdgeIndex()) {
        DartIterator d_l(*this, d.getGlobalTriangleIndex(), l);
        LocalIndex commonNode = d.getCommonNodeLocalIndex(d_l);
        DartIterator d_prime( *this, d.getGlobalTriangleIndex(), l, commonNode);
        // bisect T along the longest edge
        refine(d_prime);
        // because d_prime.node also lies on d.edge, T is now modifies in such a way that d.edge has remained unchanged. Proceed by recursion until d.edge is the longest edge in T.
        return refineOnlyThis(d, midpoint);
    }
    // create DartIterators to all edges, starting from l
    DartIterator d1 = d;  d1.flipNode(); d1.flipEdge();
    DartIterator d2 = d1; d2.flipNode(); d2.flipEdge();
    // create new triangle T_new and connect with neighbours
    GlobalIndex T_newIndex = this->pushBackTriang ( Indices3DType( midpoint, d1.getGlobalNodeIndex(), d2.getGlobalNodeIndex() ) );
    this->_neighbour_.push_back( Indices3DType ( d1.getNextTriangleIndex(), d.getGlobalTriangleIndex(), d.getNextTriangleIndex() ) );
    // set neighbour GlobalIndex in triangle lying opposite to d1.edge to T_new
    DartIterator d1_prime = d1;
    if (d1_prime.canFlipTriangle()){
        d1_prime.flipTriangle();
        this->setNeighbour ( d1_prime.getGlobalTriangleIndex(), d1_prime.getLocalEdgeIndex(), T_newIndex );
    }
    // connect ourselves (i. e., connect T) to T_new:
    this->setNeighbour( d1.getGlobalTriangleIndex(), d1.getLocalEdgeIndex(), T_newIndex );
    // do not have to refine this triangle again
    unmark( d1.getGlobalTriangleIndex() );
    // old triangle now has a new vertex, ie the midpoint
    this->getTriang( d1.getGlobalTriangleIndex() ).setGlobalNodeIdx( d1.getLocalNodeIndex(), midpoint);

    return T_newIndex;
  }
  
  // get longest edge index (starting search possibly with a preferred edge)
  LocalIndex getLongestEdgeIndex( GlobalIndex triangle ) const {
      // get global node indices of triangle
      const Indices3DType &triangIdx ( this->getTriang( triangle ).getGlobalNodeIdx() );
      // assume that edge 0 is longest edge
      LocalIndex localIdx = 0;
      Point3DType edge = this->getVertex( triangIdx[1] ); edge -= this->getVertex( triangIdx[2] ); RealType maxLength = edge.squaredNorm();  
      // now check if edge 1 is longer
      edge = this->getVertex( triangIdx[2] ); edge -= this->getVertex( triangIdx[0] ); RealType candidate = edge.squaredNorm();
      if( candidate > maxLength ){ localIdx = 1; maxLength = candidate; }
      // now check if edge 2 is longer
      edge = this->getVertex( triangIdx[0] ); edge -= this->getVertex( triangIdx[1] ); candidate = edge.squaredNorm();
      if( candidate > maxLength ) localIdx = 2;
      
      return localIdx;  
   }
  
  // add edge midpoint on d.edge, update _interpolationMap and return global index of node
  GlobalIndex addEdgeMidpoint( const DartIterator& d ) {    
    DartIterator d1 = d; d1.flipNode();
    
    ParentInformation<DomVecType> parents;
    parents.globalIndices[0] = d.getGlobalNodeIndex(); parents.globalIndices[1] = d1.getGlobalNodeIndex();
    parents.ElementIndex = d.getGlobalTriangleIndex(); parents.EdgeNumber = getLocalEdgeNumber( parents );
    DomVecType RefCoords = getRefCoordsOfEdgeMidpoint( parents );
    parents.RefCoords = RefCoords;
    
    // get edge midpoint between nodes n1 and n2
    Point3DType newVertex; newVertex.setZero();
    newVertex += this->getVertex ( d.getGlobalNodeIndex() );
    newVertex += this->getVertex ( d1.getGlobalNodeIndex() );
    newVertex /= 2.;
    
    // add new vertex
    int newNodeIdx = this->pushBackVertex( newVertex ); 
    _interpolationMap.insert( std::pair< int, ParentInformation<DomVecType> >( newNodeIdx, parents ) );
    return newNodeIdx;
  }
  
  int getLocalEdgeNumber ( const ParentInformation<DomVecType> & parents ) const {
      const Indices3DType &GlobalNodeIndizesOfEl ( this->getTriang( parents.ElementIndex ).getGlobalNodeIdx() );
      for( int i=0; i<3; ++i ){
        if( (GlobalNodeIndizesOfEl(i) != parents.globalIndices[0]) && (GlobalNodeIndizesOfEl(i) != parents.globalIndices[1]) ) return i;
      }
      throw std::invalid_argument ( aol::strprintf ( "something wrong with global indices in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
      return 23;
  }
  
  DomVecType getRefCoordsOfEdgeMidpoint ( const ParentInformation<DomVecType> & parents ) const {
      DomVecType RefCoords;
      const Indices3DType &GlobalNodeIndizesOfEl ( this->getTriang( parents.ElementIndex ).getGlobalNodeIdx() );
      if( (GlobalNodeIndizesOfEl(0) != parents.globalIndices[0]) && (GlobalNodeIndizesOfEl(0) != parents.globalIndices[1]) ){
        RefCoords << 0.5, 0.5;
        return RefCoords;
      }
      if( (GlobalNodeIndizesOfEl(1) != parents.globalIndices[0]) && (GlobalNodeIndizesOfEl(1) != parents.globalIndices[1]) ){
        RefCoords << 0.0, 0.5;
        return RefCoords;
      }
      if( (GlobalNodeIndizesOfEl(2) != parents.globalIndices[0]) && (GlobalNodeIndizesOfEl(2) != parents.globalIndices[1]) ){
        RefCoords << 0.5, 0.0;
        return RefCoords;
      }
      throw std::invalid_argument ( aol::strprintf ( "something wrong with global indices in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
      return RefCoords;
  }

};


} // end namespace 

#endif
