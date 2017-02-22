#ifndef __TRIANGMESH_H
#define __TRIANGMESH_H

#include <queue>

namespace shellFE {

template< typename DataTypeContainer, typename TriangleType >
class TriangMesh {
public:
  typedef TriangleType ElementType;
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::Point3DType Point3DType;
  typedef typename DataTypeContainer::TangentVecType    TangentVecType;
  typedef typename DataTypeContainer::Indices3DType Indices3DType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::MaskType MaskType;

protected:
  std::vector< Point3DType > _vertexIterator;
  std::vector< TriangleType > _triangIterator;
  mutable std::vector<TangentVecType> _tangentSpaceVec1;
  mutable std::vector<TangentVecType> _tangentSpaceVec2;
  mutable std::vector<TangentVecType> _normalSpaceVec;
  mutable std::vector< Indices3DType > _neighbour_;

public:
  //! Create empty TriangMesh
  TriangMesh ( ) : _vertexIterator(), _triangIterator()  { }

  TriangMesh ( const string& fileName ){
    this->loadFromLegacyVTK ( fileName );
    generateApproximativeTangentSpaceAtNodes();
    updateAllProjectionCoefficients();
  }

  virtual ~TriangMesh ( ) {} 

public:
  int getNumVertices ( ) const { return ( _vertexIterator.size() );}
  int getNumTriangs ( ) const { return ( static_cast<int> ( _triangIterator.size() ) );}

  //! insert new vertex and return global index
  int pushBackVertex ( const Point3DType newVertex ) {
    _vertexIterator.push_back ( newVertex );
    return getNumVertices() - 1;
  }
  //! insert new triangle and return global index
  int pushBackTriang ( const Indices3DType nodeIdx ) {
    int globalIdx = getNumTriangs();
    _triangIterator.push_back ( TriangleType( globalIdx, nodeIdx, _vertexIterator ) );
    return globalIdx;
  }

  const Point3DType& getVertex ( const int num ) const { return _vertexIterator[num];}
  void setVertex ( const int num, const Point3DType Arg ) { _vertexIterator[num] = Arg;}

  int getNeighbour ( const int elementID, const int acrossLocalNode ) const { return _neighbour_[elementID][acrossLocalNode];}
  void setNeighbour ( const int elementID, const int acrossLocalNode, const int value ) const {  _neighbour_[elementID][acrossLocalNode] = value;}

  const TriangleType& getTriang ( const int num ) const {return _triangIterator[num];}
  TriangleType & getTriang ( const int num ) {return _triangIterator[num];}
  void setTriang ( const int num, const TriangleType Arg ) { _triangIterator[num] = Arg;}

  int getTriangNodeIdx  ( const int num, const int localNode ) const { return _triangIterator[num].getGlobalNodeIdx(localNode);}
  void setTriangNodeIdx ( const int num, const int localNode, const int newIdx ) {_triangIterator[num].setGlobalNodeIdx(localNode, newIdx );}

  RealType getTangentVec1 ( int numComp, int VertexIndex ) const { return _tangentSpaceVec1[VertexIndex][numComp]; }
  RealType getTangentVec2 ( int numComp, int VertexIndex ) const { return static_cast<RealType>( _tangentSpaceVec2[VertexIndex][numComp] );  }
  const TangentVecType & getTangentVec1 (const int& vertexIndex ) const{ return _tangentSpaceVec1[vertexIndex];}
  const TangentVecType & getTangentVec2 (const int& vertexIndex ) const { return _tangentSpaceVec2[vertexIndex];}
  const TangentVecType & getNormalVec (const int& vertexIndex ) const{ return _normalSpaceVec[vertexIndex];}
  void setTangentVec1 (const int& vertexIndex, const TangentVecType& vec ) { _tangentSpaceVec1[vertexIndex] = vec;}
  void setTangentVec2 (const int& vertexIndex, const TangentVecType& vec ) { _tangentSpaceVec2[vertexIndex] = vec;}
  void setNormalVec (const int& vertexIndex, const TangentVecType& vec ) { _normalSpaceVec[vertexIndex] = vec;}

  // refinement operations only update NodeIndex informations, hence have to update coords, edges and projection coefficients of elements
  void updateAllTriangles ( ) {
    for ( int elementIndex = 0; elementIndex < getNumTriangs(); ++elementIndex  ) _triangIterator[elementIndex].updateNodesAndEdges( _vertexIterator );
    for ( int elementIndex = 0; elementIndex < getNumTriangs(); ++elementIndex  ) _triangIterator[elementIndex].updateProjectionCoefficients( _tangentSpaceVec1, _tangentSpaceVec2 );   
  }

  void updateAllProjectionCoefficients ( ) {
    for ( int elementIndex = 0; elementIndex < getNumTriangs(); ++elementIndex  ) _triangIterator[elementIndex].updateProjectionCoefficients( _tangentSpaceVec1, _tangentSpaceVec2 );   
  }

  void print ( ) {
    for ( int nodeIndex = 0; nodeIndex < getNumVertices(); ++nodeIndex ){
      cout << endl << "node = " << nodeIndex << endl;
      cout << "tv1   = " << _tangentSpaceVec1[nodeIndex].transpose() << endl;
      cout << "tv2   = " << _tangentSpaceVec2[nodeIndex].transpose() << endl;
      cout << "n     = " << _normalSpaceVec[nodeIndex].transpose() << endl;
      cout << "v1xv2 = " << (_tangentSpaceVec1[nodeIndex].cross( _tangentSpaceVec2[nodeIndex] ) ).transpose() << endl;
    }
    for ( int elementIndex = 0; elementIndex < getNumTriangs(); ++elementIndex  ) _triangIterator[elementIndex].print();
  }

public:

  //! load from file in the .vtk file format. Currently only loads geometric information.
  void loadFromLegacyVTK( const string& filename ) {
    /*std::ifstream vtkFile ( filename.c_str() ); */
    std::ifstream vtkFile;
    vtkFile.open ( filename.c_str() );
    if ( vtkFile.fail() ) throw std::invalid_argument( "Cannot open file" );
    bool readVertices = false; bool readFaces = false;
    while ( !vtkFile.eof() ) {
      char line[256];
      vtkFile.getline ( line, 256 );
      // first: expect vertices
      if( !strncmp ( line, "POINTS ", strlen ( "POINTS " ) ) ){
        readVertices = true;
        continue;
      }
      // second: expect triangles (i.e. starting with "3 " )
      if( !strncmp ( line, "POLYGONS ", strlen ( "POLYGONS " ) ) ){
        readVertices = false;
        readFaces = true;
        continue;
      }
      // geometric information ends with the first line that does not start with "3 " anymore
      if( readFaces && strncmp ( line, "3 ", strlen ( "3 " ) ) ) break;
      // read in the vertex coordinates and add it to the mesh
      if( readVertices ){
        Point3DType vertex;
        char * substring = strtok ( line, " " );
        for ( int i = 0; i < 3; i++ ){
          vertex[i] = atof( substring );
          substring = strtok (NULL, " ");
        }
        pushBackVertex ( vertex );
      }
      // read in the face and add it to the mesh
      if( readFaces ){
        Indices3DType triangle;
        char * substring = strtok ( line, " " );
        for ( int i = 0; i < 3; i++ ){
          substring = strtok (NULL, " ");
          triangle[i] = atof( substring );
        }
        pushBackTriang( triangle );
      }
    }
    vtkFile.close();
  }


  void generateApproximativeTangentSpaceAtNodes (){
    const int numVertices = this->getNumVertices ();
    _tangentSpaceVec1.resize ( numVertices ); _tangentSpaceVec2.resize ( numVertices ); _normalSpaceVec.resize( numVertices );
    std::fill( _normalSpaceVec.begin(), _normalSpaceVec.end(), TangentVecType(0.,0.,0.) );
    for ( int elementIndex = 0; elementIndex < this->getNumTriangs(); ++elementIndex  ) {
      TangentVecType normalizedNormal;
      _triangIterator[elementIndex].getNormalizedNormalForFlattenedTriangle( normalizedNormal );
      RealType areaOfElement = _triangIterator[elementIndex].getAreaOfFlattenedTriangle ();
      for (int idx = 0; idx < 3; ++idx)
        _normalSpaceVec[_triangIterator[elementIndex].getGlobalNodeIdx(idx)] += areaOfElement * normalizedNormal;
    }
    for (int vertexIndex = 0; vertexIndex < numVertices; ++vertexIndex) _normalSpaceVec[vertexIndex].normalize();
    for (int vertexIndex = 0; vertexIndex < numVertices; ++vertexIndex) getTangentSpaceFromNormal ( _tangentSpaceVec1[vertexIndex], _tangentSpaceVec2[vertexIndex], _normalSpaceVec[vertexIndex] );
  }

  // for boundary node only averge both neighboring triangle normals s.t. the corresponding triangle has another boundary node
  void generateApproximativeTangentSpaceAtNodes ( const MaskType &boundaryMask ){
    const int numVertices = this->getNumVertices ();
    _tangentSpaceVec1.resize ( numVertices ); _tangentSpaceVec2.resize ( numVertices ); _normalSpaceVec.resize( numVertices );
    std::fill( _normalSpaceVec.begin(), _normalSpaceVec.end(), TangentVecType(0.,0.,0.) );
    for ( int elementIndex = 0; elementIndex < this->getNumTriangs(); ++elementIndex  ) {
      TangentVecType normalizedNormal; _triangIterator[elementIndex].getNormalizedNormalForFlattenedTriangle( normalizedNormal );
      RealType areaOfElement = _triangIterator[elementIndex].getAreaOfFlattenedTriangle ();
      for (int idx = 0; idx < 3; ++idx){
        int globIdx = _triangIterator[elementIndex].getGlobalNodeIdx(idx);
        if( boundaryMask[globIdx] == false ){ _normalSpaceVec[ globIdx ] += areaOfElement * normalizedNormal;
        }else {
          // node is on boundary, so ask if element has another boundary node
          if( boundaryMask[ _triangIterator[elementIndex].getGlobalNodeIdx( (idx + 1) % 3) ] || boundaryMask[ _triangIterator[elementIndex].getGlobalNodeIdx( (idx + 2) % 3) ] )
            _normalSpaceVec[ globIdx ] += areaOfElement * normalizedNormal;
        }
      }
    }
    for (int vertexIndex = 0; vertexIndex < numVertices; ++vertexIndex) _normalSpaceVec[vertexIndex].normalize();
    for (int vertexIndex = 0; vertexIndex < numVertices; ++vertexIndex) getTangentSpaceFromNormal ( _tangentSpaceVec1[vertexIndex], _tangentSpaceVec2[vertexIndex], _normalSpaceVec[vertexIndex] );
  }


  //Compute orthonormal basis (v1, v2, n ) of (tangent space + normal space);   //since |n|^2 = nx^2 + ny^2 + nz^2 = 1  => at least one |ni| >= 1/sqrt(3) > 0.5
  void getTangentSpaceFromNormal ( TangentVecType &tv1, TangentVecType &tv2, const TangentVecType& normal ) const {
    if ( (std::abs ( normal(2) ) > 0.5) ) tv1 = TangentVecType( (-1.) * normal(2), 0., normal(0) );
    else                                  tv1 = TangentVecType( (-1.) * normal(1), normal(0), 0. );
    tv1.normalize();
    tv2 = normal.cross( tv1 );
    tv2.normalize();
  }


  void makeNeighbour() const {

    int  *n_to_t, *num_of_n, *start_of_n;
    const int numNodes = getNumVertices(); const int numElements = getNumTriangs();

    if ( int(_neighbour_.size()) != numElements ) _neighbour_.resize ( numElements );

    n_to_t     = new int[numElements*3];
    num_of_n   = new int[numNodes];
    start_of_n = new int[numNodes+1];

    // iterate over all vertices and set num_of_n to zero, i.e. initialize
    for ( int nodeIndex = 0; nodeIndex < numNodes; nodeIndex++ ) num_of_n[nodeIndex] = 0;

    // iterate over all triangles and all neighbours, i.e. 3, since every triangle has 3 neighbours for closed meshes (TODO: what happens for non closed? )
    for ( int elementIndex = 0; elementIndex < numElements; elementIndex++ )
      for ( int j = 0; j < 3; j++ ){
        num_of_n[getTriangNodeIdx ( elementIndex,j ) ]++;
        setNeighbour ( elementIndex, j, -1 );
      }

    // get Startindex for node
    int nodeIdx = 0; start_of_n[nodeIdx++] = 0;
    while ( nodeIdx < numNodes ) {
      start_of_n[nodeIdx] = start_of_n[nodeIdx-1] + num_of_n[nodeIdx-1];
      nodeIdx++;
    }
    start_of_n[numNodes] = 3 * numElements;

    // initialize reference
    for ( int elementIndex = 0; elementIndex < numElements*3; ++elementIndex ) n_to_t[elementIndex] = -1;

    // build reference
    for ( int elementIndex = 0; elementIndex < numElements; elementIndex++ )
      for ( int j = 0; j < 3; j++ ) {
        int k = start_of_n[getTriangNodeIdx ( elementIndex,j ) ];
        while ( n_to_t[k] > -1 ) k++;
        n_to_t[k] = elementIndex;
      }
    // find neighbour
    for ( int elementIndex = 0; elementIndex < numElements; elementIndex++ )
      for ( int j = 0; j < 3; j++ ) {
        int node1 = getTriangNodeIdx ( elementIndex, j        );
        int node2 = getTriangNodeIdx ( elementIndex, ( j + 1 ) % 3 );
        int node3 = getTriangNodeIdx ( elementIndex, ( j + 2 ) % 3 );

        for ( int k = start_of_n[node1]; k < start_of_n[node1+1]; k++ ) {
          int n = n_to_t[k]; // Tetraeder
          if ( elementIndex < n ) // set neighborhood only once
            for ( int l = 0; l < 3; l++ ) {
              if ( node3 == getTriangNodeIdx ( n, l ) ) {
                setNeighbour ( elementIndex, ( j + 1 ) % 3, n );
                for ( int v = 0; v < 3; v++ )
                  if ( v != l && getTriangNodeIdx ( n, v ) != node1 ) setNeighbour ( n, v, elementIndex );
              }
              else {
                if ( node2 == getTriangNodeIdx ( n, l ) ) {
                  setNeighbour ( elementIndex, ( j + 2 ) % 3, n );
                  for ( int v = 0; v < 3; v++ )
                    if ( v != l && getTriangNodeIdx ( n, v ) != node1 ) setNeighbour ( n, v, elementIndex );
                }
              }
            }
        }
      }

    delete[] n_to_t; delete[] num_of_n; delete[] start_of_n;
  }



  void makeOrientationConsistent() {
    if ( int(_neighbour_.size()) != this->getNumTriangs() ) makeNeighbour();
    // true for all triangles T, whose neighboring triangles have already been oriented consistent with T
    MaskType alreadyHandled( this->getNumTriangs() );
    // true for all triangles who have already been dealt with or who are already waiting in queue
    MaskType inQueue( this->getNumTriangs() );
    // contains all triangles which are already oriented and whose neighbors will be dealt with next (may contain triangles twice for simplicity)
    std::queue<int> activeTriangles;
    activeTriangles.push( 0 );
    inQueue[0] = true;
    // the triangle whose neighbors are currently handled
    int currentTriangle;
    // while there are triangles left whose neighbors are not consistently oriented...
    while( !activeTriangles.empty() ){
      currentTriangle = activeTriangles.front();
      activeTriangles.pop();
      // deal with all three neighbors of currentTriangle, i.e. orient them and add them to the list to deal with their neighbors
      for ( int i = 0; i < 3; i++ ){
        int neighbor = getNeighbour( currentTriangle, i );
        if ( neighbor >= 0 && neighbor < this->getNumTriangs() && !alreadyHandled[neighbor] ){
          // compute the nodes "currentTriangle" and "neighbor" have in common
          int node1 = getTriangNodeIdx ( currentTriangle, ( i + 1 ) % 3 );
          int node2 = getTriangNodeIdx ( currentTriangle, ( i + 2 ) % 3 );
          // check whether common nodes occur in reversed order in "neighbor", if not, change order
          int j = 0;
          while ( getTriangNodeIdx ( neighbor, j ) != node2 )
            j++;
          if ( getTriangNodeIdx ( neighbor, ( j + 1 ) % 3 ) != node1 ){
            // change order of nodes
            int exchangeCache = getTriangNodeIdx ( neighbor, ( j + 1 ) % 3 );
            setTriangNodeIdx( neighbor, ( j + 1 ) % 3, getTriangNodeIdx ( neighbor, ( j + 2 ) % 3 ) );
            setTriangNodeIdx( neighbor, ( j + 2 ) % 3, exchangeCache );
            // change order of corresponding neighbours
            exchangeCache = getNeighbour( neighbor, ( j + 1 ) % 3 );
            setNeighbour( neighbor, ( j + 1 ) % 3, getNeighbour( neighbor, ( j + 2 ) % 3 ) );
            setNeighbour( neighbor, ( j + 2 ) % 3, exchangeCache );
          }
          if ( !inQueue[neighbor] ){
            activeTriangles.push( neighbor );
            inQueue[neighbor] = true;
          }
        }
      }
      alreadyHandled[currentTriangle] = true;
    }
  }

};


}//end namespace

#endif
