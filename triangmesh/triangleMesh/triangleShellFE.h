#ifndef __TRIANGLESHELLFE_H
#define __TRIANGELSHELLFE_H

namespace shellFE {

//! Triangle which has a tangent space at each node 
template<typename DataTypeContainer>
class ShellElementWithTangentSpaceAtVertex {

  typedef typename DataTypeContainer::RealType        RealType;
  typedef typename DataTypeContainer::DomVecType      DomVecType;
  typedef typename DataTypeContainer::Point3DType     Point3DType;
  typedef typename DataTypeContainer::TangentVecType  TangentVecType;
  typedef typename DataTypeContainer::Indices3DType   Indices3DType;
  typedef std::vector<Point3DType>                    VertexIterator;
  typedef std::vector<TangentVecType>                 TangentVecIterator;

  public :   
    // global indices of element an nodes
    int _globIdx;
    Indices3DType _globNodeIdx;
    Point3DType _nodes[3];
    TangentVecType _edges[3]; // Notation:   e_i = N_{i-1} - N_{i+1} = N_{i+2} - N_{i+1}, i = 0,1,2

    // corresponding coefficients: projVec = projCoeff1 v_1 + projCoeff2 v_2. Save as projcoeff(localNodeIndex)(Direction)
    RealType _projCoeff1[3][2];
    RealType _projCoeff2[3][2];
    
public:
  
  ShellElementWithTangentSpaceAtVertex () : _globIdx(-1) {}
  
  ShellElementWithTangentSpaceAtVertex( const int globalIdx, const Indices3DType globalNodeIndex, const VertexIterator &nodes  ) : 
     _globIdx(globalIdx),
     _globNodeIdx ( globalNodeIndex )
    {
         // fill nodes and edges
        for ( int i = 0; i < 3; ++i ) _nodes[i] = nodes[ _globNodeIdx[i] ];
        for ( int i = 0; i < 3; ++i ) _edges[i] = getNode((i+2)%3) - getNode((i+1)%3);
    }
  
  ~ShellElementWithTangentSpaceAtVertex(){}
  
  //Let v1, v2 basis of tangent space, write w as linear combination of v1, v2 and (v1 x v2): w = a1 * v1 + a2 * v2 + a3 * (v1 x v2), return coefficient a1 and a2
  void computeCoefficientsOfProjection (const TangentVecType& v1, const TangentVecType& v2, const TangentVecType& w, 
                                        RealType& a1, RealType& a2 ){
    RealType normV1Squared = v1.dot(v1); RealType normV2Squared = v2.dot(v2);
    RealType scalarProductV1V2 = v1.dot(v2); RealType scalarProductV1W = v1.dot(w); RealType scalarProductV2W = v2.dot(w);
    RealType scalarProductV1V2Squared = scalarProductV1V2 * scalarProductV1V2;
    
    a1 = normV2Squared * scalarProductV1W - scalarProductV1V2 * scalarProductV2W;
    a1 /= (normV1Squared * normV2Squared - scalarProductV1V2Squared);
    
    a2 = scalarProductV2W * normV1Squared - scalarProductV1W * scalarProductV1V2;
    a2 /= (normV1Squared * normV2Squared - scalarProductV1V2Squared);
  }

  
  void updateNodesAndEdges( const VertexIterator &nodes ){
      for ( int i = 0; i < 3; ++i ) _nodes[i] = nodes[ _globNodeIdx[i] ];
      for ( int i = 0; i < 3; ++i ) _edges[i] = getNode((i+2)%3) - getNode((i+1)%3);
  }
  
  void updateProjectionCoefficients( const TangentVecIterator &tangentVec1, const TangentVecIterator &tangentVec2 ){ 
    
    // project edges to tangent space
    RealType tmpCoeff1, tmpCoeff2;
    TangentVecType vecToProject;
    
    // at node N_0 = (0,0)
    // vec corresponding to (1,0)
    computeCoefficientsOfProjection ( tangentVec1 [ _globNodeIdx[0] ], tangentVec2 [ _globNodeIdx[0] ], _edges [2], _projCoeff1[0][0], _projCoeff2 [0][0] );
    // vec corresponding to (0,1)
    vecToProject = -1.0 * _edges[1];
    computeCoefficientsOfProjection ( tangentVec1 [ _globNodeIdx[0] ], tangentVec2 [ _globNodeIdx[0] ], vecToProject, _projCoeff1 [0][1], _projCoeff2 [0][1] );
    
    // at node N_1 = (1,0)
    // vec corresponding to (1,0)
    vecToProject = -1.0 * _edges[2];
    computeCoefficientsOfProjection ( tangentVec1 [ _globNodeIdx[1] ], tangentVec2 [ _globNodeIdx[1] ], vecToProject, tmpCoeff1 , tmpCoeff2  );
    _projCoeff1 [1][0] = - 1.0 * tmpCoeff1;
    _projCoeff2 [1][0] = - 1.0 * tmpCoeff2;
    // vec corresponding to (0,1)
    _projCoeff1 [1][1] = - 1.0 * tmpCoeff1;
    _projCoeff2 [1][1] = - 1.0 * tmpCoeff2;
    computeCoefficientsOfProjection ( tangentVec1 [ _globNodeIdx[1] ], tangentVec2 [ _globNodeIdx[1] ], _edges[0], tmpCoeff1 , tmpCoeff2 );
    _projCoeff1 [1][1] +=  tmpCoeff1;
    _projCoeff2 [1][1] +=  tmpCoeff2;

    // at node N_2 = (0,1)
    // vec corresponding to (0,1)
    computeCoefficientsOfProjection ( tangentVec1 [ _globNodeIdx[2] ], tangentVec2 [ _globNodeIdx[2] ], _edges[1], tmpCoeff1 , tmpCoeff2  );
    _projCoeff1 [2][1] = - 1.0 * tmpCoeff1;
    _projCoeff2 [2][1] = - 1.0 * tmpCoeff2;
    // vec corresponding to (1,0)
    _projCoeff1 [2][0] = -1.0 * tmpCoeff1;
    _projCoeff2 [2][0] = -1.0 * tmpCoeff2;
    vecToProject = -1.0 * _edges[0];
    computeCoefficientsOfProjection ( tangentVec1 [ _globNodeIdx[2] ], tangentVec2 [ _globNodeIdx[2] ], vecToProject, tmpCoeff1 , tmpCoeff2 );
    _projCoeff1 [2][0] +=  tmpCoeff1;
    _projCoeff2 [2][0] +=  tmpCoeff2;
    
  }

  void printNodes() const { cout << "node0 = " << getNode(0).transpose() << endl << "node1 = " << getNode(1).transpose() << endl << "node2 = " << getNode(2).transpose() << endl;}

  void print() const {
      cout << endl << "Element " << _globIdx << endl;
      printNodes();
      cout << "projcoeff = " << endl;
      for( int i=0; i<3; ++i )
          for( int j=0; j<2; ++j ){
           cout << _projCoeff1[i][j] << "  ,    " << _projCoeff2[i][j] << endl;   
          }
  }
  
  //The following methods compute the weighted/normalized normal of the triangle elements and its area under the assumption that the triangle is flat
  void getWeightedNormalForFlattenedTriangle ( TangentVecType &normal ) const {
    TangentVecType dx1 = getNode(1), dx2 = getNode(2);
    dx1 -= getNode(0);
    dx2 -= getNode(0);
    normal = TangentVecType( dx1.cross(dx2) );
  }

  void getNormalizedNormalForFlattenedTriangle ( TangentVecType &normal ) const {
    getWeightedNormalForFlattenedTriangle ( normal );
    normal.normalize();
  }

  inline RealType getAreaOfFlattenedTriangle ( ) const {
    TangentVecType n; getWeightedNormalForFlattenedTriangle ( n );
    return n.norm() / 2.;
  }
    
  // get and set functions
  int getGlobalElementIdx(  ) const { return _globIdx;} //TODO old: globIdx, getIndex
  const Indices3DType & getGlobalNodeIdx( ) const{ return _globNodeIdx;}
  int getGlobalNodeIdx(int localIndex) const{ return _globNodeIdx[localIndex];}
  void setGlobalNodeIdx(int localIndex, int globalIndex) { _globNodeIdx[localIndex] = globalIndex;}
  
  const Point3DType& getNode ( int i ) const { return _nodes[i];}
  Point3DType& getNode ( int i ) { return _nodes[i];}
  void setNode ( int i, const Point3DType& node ) {_nodes[i] = node;}
  
  const TangentVecType& operator[] ( int i ) const { return getNode(i);}
  TangentVecType& operator[] ( int i ) {return getNode(i);}
  
  
  RealType getProjCoeff1 ( const int BaseFuncNum ) const {
    int localNodeIndex = std::floor(BaseFuncNum/3);
    int direction = BaseFuncNum%3;
    direction -= 1;
    return _projCoeff1[localNodeIndex][direction];
  }
  
  RealType getProjCoeff2 ( const int BaseFuncNum ) const {
    int localNodeIndex = std::floor(BaseFuncNum/3);
    int direction = BaseFuncNum%3;
    direction -= 1;
    return _projCoeff2[localNodeIndex][direction];
  }

};


} //end namespace


#endif
