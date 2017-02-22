#ifndef __CONFIGURATORSSHELLFE_H
#define __CONFIGURATORSSHELLFE_H

#include <basefunctionSetShellFE.h>

namespace shellFE {

enum ShellFEType {
  NodalValuedDofs,
  C1Dofs,
};

//! each triangle of the mesh is the unit triangle embedded in R^2, given by the three positions (0,0), (1,0) and (0,1) with area 0.5
template <typename DataTypeContainer, typename MeshType, typename QuadType >
class UnitTriangMeshConfiguratorP1 {
protected:
  const MeshType &_mesh;

public:

  static const int maxNumLocalDofs = 3;
  static const ShellFEType _ShellFEType = NodalValuedDofs;

  typedef DataTypeContainer                                             DTContainer;
  typedef MeshType                                                      InitType;               //!< that's the type that is needed by the constructor of the configurator
  typedef typename MeshType::ElementType                                ElementType;
  typedef typename DataTypeContainer::RealType                          RealType;
  typedef typename DataTypeContainer::DomVecType                        DomVecType;
  typedef typename DataTypeContainer::TangentVecType                    TangentVecType;
  typedef typename DataTypeContainer::Point3DType                       Point3DType;
  typedef typename DataTypeContainer::Matrix22                          Matrix22;
  typedef typename DataTypeContainer::Matrix32                          Matrix32;
  typedef typename DataTypeContainer::Matrix33                          Matrix33;
  typedef typename DataTypeContainer::VectorType                        VectorType;
  typedef typename DataTypeContainer::FullMatrixType                    FullMatrixType;
  typedef typename DataTypeContainer::TripletType                       TripletType;
  typedef typename DataTypeContainer::SparseMatrixType                  SparseMatrixType;
  typedef typename DataTypeContainer::MaskType                          MaskType;
  typedef Eigen::Matrix<RealType, maxNumLocalDofs, 1 >                  LocalVectorType;
  typedef Eigen::Matrix<RealType, maxNumLocalDofs, maxNumLocalDofs >    LocalMatrixType;

  typedef UnitTriangMeshBaseFunctionSetP1<DataTypeContainer, QuadType, ElementType> BaseFuncSetType;

  UnitTriangMeshConfiguratorP1 ( const InitType &Mesh ) : _mesh ( Mesh ) {}

  const InitType& getInitializer( ) const { return this->_mesh; }

  mutable BaseFuncSetType _baseFuncSet;

  inline int getNumLocalDofs ( const ElementType & ) const { return 3;}
  inline int getNumLocalDofs (  ) const { return 3;}
  const MeshType& getMesh ( ) const { return _mesh; }
  int getNumGlobalDofs( ) const {return this->_mesh.getNumVertices();}
  int maxNumQuadPoints( ) const { return QuadType::numQuadPoints;}
  const BaseFuncSetType& getBaseFunctionSet ( const ElementType &/*T*/ ) const {return _baseFuncSet;}
  const BaseFuncSetType& getBaseFunctionSet ( ) const {return _baseFuncSet;}

  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const ElementType &T, int localIndex ) const {
    return T.getGlobalNodeIdx( localIndex );
  }

  void convertRefCoordsToBaryCentric ( const DomVecType& RefCoords, Point3DType& BaryCentricCoords ) const {
    BaryCentricCoords[1] = RefCoords[0];
    BaryCentricCoords[2] = RefCoords[1];
    BaryCentricCoords[0] = 1. - BaryCentricCoords[1] - BaryCentricCoords[2];
  }

  void getGlobalCoordsFromBarycentric ( const ElementType &El, const Point3DType& BaryCentric, DomVecType& Dest ) const {
    Dest.setZero ();
    for ( unsigned int i = 0; i < 3; ++i ) {
      Dest += BaryCentric[i] * El.getNode ( i ).segment ( 0, 2 );
    }
  }

  void getGlobalCoords ( const ElementType &El, const int QuadPoint, DomVecType& Dest ) const {
    DomVecType refCoords ( _baseFuncSet.getRefCoord ( QuadPoint ) );
    Point3DType baryCentric;
    convertRefCoordsToBaryCentric ( refCoords, baryCentric);
    getGlobalCoordsFromBarycentric ( El, baryCentric, Dest );
  }

};

} // end namespace

#endif
