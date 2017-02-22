#ifndef __SHELLHANDLER_H
#define __SHELLHANDLER_H

#include <PParser.h>
#include <discreteFunctionShellFE.h>
#include <configuratorsShellFE.h> // for ShellFEType
#include <meshWithData.h>

namespace shellFE {

enum ShellBoundaryType {
      ALLBOUNDARY = 100,
      NOBOUNDARY = 0,
  //Plate
      PlateLeft = 1,
      PlateLeftTop = 2,
      PlateAll = 3,
  //Cylinder
      CylinderTopBottom = 11,
      CylinderBottom = 12,
  //Sphere
      SphereTopBottom = 21,
      SphereEquator = 22,
  //PipeConnection
      PipeConnectionAll = 51
};


template< typename ConfiguratorType >
void generateDirichletBoundaryMaskUponShellBoundaryType ( const ShellBoundaryType ShellType, const typename ConfiguratorType::InitType &mesh,
                                     typename ConfiguratorType::MaskType & mask, int & numBoundaryNodes,
                                     const bool clampedBoundaryCondition = false ) {
      
  typedef typename ConfiguratorType::Point3DType    Point3DType;
  
  numBoundaryNodes = 0;
  const int numVertices = mesh.getNumVertices();
  
  switch ( ShellType ){
      
      case NOBOUNDARY : {
      } break;
      
      case ALLBOUNDARY : {
            mesh.makeNeighbour();
            for( int ElementIndex = 0; ElementIndex < mesh.getNumTriangs(); ++ElementIndex ){
                for( int i = 0; i < 3; i++ ){
                    if ( mesh.getNeighbour( ElementIndex , i ) == -1 ){
                        mask[ mesh.getTriangNodeIdx( ElementIndex, (i+1)%3) ] = true;
                        mask[ mesh.getTriangNodeIdx( ElementIndex, (i+2)%3) ] = true;
                    }
                }
            }
      } break;
      
      case PlateLeft : {
            for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
              if (coords [0] == 0. ){
                ++numBoundaryNodes;
                mask[nodeIdx] = true;
              }
            }
      } break;
      
      case PlateLeftTop : {
            for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
              if (coords [0] == 0. || coords[1] == 1. ){
                ++numBoundaryNodes;
                mask[nodeIdx] = true;
              }
            }
      } break;
      
     case PlateAll : {
            for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
              if (coords [0] == 0. || coords [0] == 1. || coords [1] == 0. || coords[1] == 1. ){
                ++numBoundaryNodes;
                mask[nodeIdx] = true;
              }
            }
      } break;
      
      case CylinderTopBottom : {
            for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
              if (coords [2] == 0. || coords [2] == 1.){
                ++numBoundaryNodes;
                mask[nodeIdx] = true;
              }
            } 
      } break;     
      
      case CylinderBottom : {
            for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
              if (coords [2] == 0. ){
                ++numBoundaryNodes;
                mask[nodeIdx] = true;
              }
            } 
      } break;   
      
      case SphereTopBottom : {
            for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
              if (coords [2] >= 0.9 || coords [2] <= -0.9){
                ++numBoundaryNodes;
                mask[nodeIdx] = true;
              }
            } 
      } break;
      
      case SphereEquator : {
            for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
              if (coords [2] <= 0.2 && coords [2] >= -0.2){
                ++numBoundaryNodes;
                mask[nodeIdx] = true;
              }
            } 
      } break;
      
      
     case PipeConnectionAll : {
         for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
             const Point3DType& coords ( mesh.getVertex(nodeIdx) );
             if (coords [2] == 0.5 || coords [2] == -0.5){
                 ++numBoundaryNodes;
                 mask[nodeIdx] = true;
             }
             if (coords [1] == 0.5 || coords [1] == -0.5){
                 ++numBoundaryNodes;
                 mask[nodeIdx] = true;
             }
         }
      } break;
      
      
      default :
        throw std::invalid_argument( aol::strprintf ( "Wrong boundary condition. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
        break;
    }
    
    //clamped boundary condition
    if( (ConfiguratorType::_ShellFEType == C1Dofs) && ( clampedBoundaryCondition ) ){
        for( int i=0; i<numVertices; ++i ){
          if( mask[i] ){
            mask[ i + numVertices] = true;
            mask[ i + 2. * numVertices] = true;
          }
        }
    }
    
}






template< typename ConfiguratorType >
class ShellHandler{
  
public:
  
  typedef typename ConfiguratorType::RealType       RealType;
  typedef typename ConfiguratorType::InitType       MeshType;
  typedef typename ConfiguratorType::MaskType       MaskType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Point3DType    Point3DType;
  typedef typename ConfiguratorType::VectorType     VectorType;
  
  const aol::PParser &_parser;
  const ConfiguratorType &_conf;
  const MeshType &_mesh;
  const int _numVertices, _numGlobalDofs;
  MaskType _DirichletMask;
  int _numBoundaryNodes;
  mutable VectorType _xA;
  const DiscreteVectorFunctionStorage <ConfiguratorType> *_xACachePtr;
  
public:
  
  ShellHandler( const aol::PParser &Parser, const ConfiguratorType &conf ) : 
  _parser ( Parser),
  _conf ( conf ),
  _mesh( conf.getInitializer() ),
  _numVertices ( _mesh.getNumVertices() ),
  _numGlobalDofs ( conf.getNumGlobalDofs() ),
  _xA ( 3 * conf.getNumGlobalDofs() )
  {
    generateChart_xA ( );
    _xACachePtr = new DiscreteVectorFunctionStorage<ConfiguratorType> ( _conf, _xA, 3 );
    generateDirichletBoundaryMask( _DirichletMask, _numBoundaryNodes );
  }

  void generateChart_xA( ) const{
      for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
        const Point3DType& coords ( _mesh.getVertex(nodeIdx) );
        
        for( int comp=0; comp<3; ++comp )
          _xA[ nodeIdx + _numGlobalDofs * comp ] = coords[comp];
               
        if( ConfiguratorType::_ShellFEType == C1Dofs ){
            TangentVecType firstTangentVecAtNode ( _mesh.getTangentVec1 (nodeIdx) );
            TangentVecType secondTangentVecAtNode ( _mesh.getTangentVec2 (nodeIdx) );
            
            for( int comp=0; comp<3; ++comp ){
              _xA[ nodeIdx +     _numVertices + _numGlobalDofs * comp ] = firstTangentVecAtNode  [comp];
              _xA[ nodeIdx + 2 * _numVertices + _numGlobalDofs * comp ] = secondTangentVecAtNode [comp];
            }
        }
      }
  }

  VectorType &getChartToUndeformedShell (  ) const {  return _xA;}
  const DiscreteVectorFunctionStorage <ConfiguratorType> &getChartToUndeformedShell_Cache () const { return *_xACachePtr;}
  
  void generateDirichletBoundaryMask ( MaskType & mask, int & numBoundaryNodes ) const{
    mask.resize( _numGlobalDofs, false );
    if( _parser.get<int> ( "ClampedBoundaryCondition" ) == 1 )
        generateDirichletBoundaryMaskUponShellBoundaryType<ConfiguratorType> ( static_cast<ShellBoundaryType>( _parser.get<int>( "ShellType" ) ), _mesh, mask, numBoundaryNodes, true  );
    else
        generateDirichletBoundaryMaskUponShellBoundaryType<ConfiguratorType> ( static_cast<ShellBoundaryType>( _parser.get<int>( "ShellType" ) ), _mesh, mask, numBoundaryNodes );
  }
  
  const MaskType & getDirichletMask ( ) const { return _DirichletMask;}
  const int & getNumBoundaryNodes() const { return _numBoundaryNodes; }
  
 //===============================================================================================================================================================================
  void constructConstantMaterial( VectorType &material, RealType materialConstant ){
    for( int i=0; i<_numVertices; ++i )
      material[i] = materialConstant;
  }
  
  void constructHardSoftMaterial( VectorType &material,  const RealType sizeFirstMaterial, 
                                  const int direction = 0, const bool firstHard = true, 
                                  const RealType mhard = 1.0, const RealType msoft = -1.0 ){
    for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
        const Point3DType& coords ( _mesh.getVertex(nodeIdx) );
        if( coords[direction] < sizeFirstMaterial )  material[nodeIdx] = mhard;
        else material[nodeIdx] = msoft;
    }
    if( firstHard == false ) material *= -1.0;        
  }
  
  void switchMaterialType( VectorType &material ){
//       material.reallocate( _numVertices );
      switch( _parser.get<int>( "initMaterialType" ) ){
        case 1:
          constructConstantMaterial( material, _parser.get<double>( "materialConstant" ) );
          break;
        case 2:
          constructHardSoftMaterial( material, _parser.get<double>( "sizeFirstMaterial" ), 0, true );
          break;
        case 3:
          constructHardSoftMaterial( material, _parser.get<double>( "sizeFirstMaterial" ), 0, false );
          break;
        default:
          throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
          break;
      }
  }
  
 //===============================================================================================================================================================================
  //==================================  Plot     ============================================================================
 //===============================================================================================================================================================================
 void plotShell ( const string dispOrDeform,
                  const VectorType & dispOrDeformVec,
                  const string outfile_base_name, const int iter = 0 ) const{
    MeshType mesh ( _mesh );
    
    if( dispOrDeform == "deform" ){
        for( int i = 0; i < mesh.getNumVertices(); i++ ){
            TangentVecType coords ( dispOrDeformVec[i], dispOrDeformVec[i+_numGlobalDofs], dispOrDeformVec[i+2*_numGlobalDofs] );
            mesh.setVertex( i, coords );
        }    
    }
    
    if( dispOrDeform == "deform" ){    
        for( int i = 0; i < mesh.getNumVertices(); ++i ){
            TangentVecType coords = mesh.getVertex( i );
            for( int comp = 0; comp < 3; ++comp ) coords[comp] += dispOrDeformVec[i + comp * _numGlobalDofs];
            mesh.setVertex( i, coords );
        }
    }

    MeshWithData<MeshType> meshSaver ( mesh );
    meshSaver.saveAsLegacyVTK( aol::strprintf ( "%s/%s_%d.vtk", _parser.get<std::string> ( "saveDirectory" ).c_str (), outfile_base_name.c_str(), iter ) );
  }

  void plotShellWithMaterial ( const string dispOrDeform,
                               const VectorType & dispOrDeformVec,
                               const VectorType & material, 
                               const string outfile_base_name, const int iter = 0 ) const{
    MeshType mesh ( _mesh );
    
    if( dispOrDeform == "deform" ){
        for( int i = 0; i < mesh.getNumVertices(); ++i ){
            TangentVecType coords;
            for( int comp = 0; comp < 3; ++comp ) coords[comp] = dispOrDeformVec[i + comp * _numGlobalDofs];
            mesh.setVertex( i, coords );
        }
    }
    
    if( dispOrDeform == "disp" ){
        for( int i = 0; i < mesh.getNumVertices(); ++i ){
            TangentVecType coords = mesh.getVertex( i );
            for( int comp = 0; comp < 3; ++comp ) coords[comp] += dispOrDeformVec[i + comp * _numGlobalDofs];
            mesh.setVertex( i, coords );
        }
    }
    
    VectorType materialAtNodes ( _numVertices );
    for( int i=0; i<_numVertices; ++i )
      materialAtNodes[i] = material[i];
    
    MeshWithData<MeshType> meshSaver ( mesh );
    meshSaver.addScalarData ( materialAtNodes, "material", VERTEX_DATA );
    meshSaver.saveAsLegacyVTK( aol::strprintf ( "%s/%s_%d.vtk", _parser.get<std::string> ( "saveDirectory" ).c_str (), outfile_base_name.c_str(), iter ) );
  }
  
  void plotUndeformedShellWithMaterial ( const VectorType & material,  const string outfile_base_name, const int iter = 0 ) const{
    
    VectorType materialAtNodes ( _numVertices );
    for( int i=0; i<_numVertices; ++i ) materialAtNodes[i] = material[i];
    
    MeshWithData<MeshType> meshSaver ( _mesh );
    meshSaver.addScalarData ( materialAtNodes, "material", VERTEX_DATA );
    meshSaver.saveAsLegacyVTK( aol::strprintf ( "%s/%s_%d.vtk", _parser.get<std::string> ( "saveDirectory" ).c_str (), outfile_base_name.c_str(), iter ) );
  }
  
  void plotNormalField ( const int iter = 0 ) const{   
    MeshType mesh (_mesh);
    int numVertices = _mesh.getNumVertices ();
    VectorType normals ( 3 * numVertices );;
    for( int nodeIter=0; nodeIter < numVertices; ++nodeIter){
        const TangentVecType &firstTangentVec (_mesh.getTangentVec1(nodeIter) ); const TangentVecType &secondTangentVec (_mesh.getTangentVec2(nodeIter) );
        TangentVecType normal = firstTangentVec.cross(secondTangentVec);
        for( int comp = 0; comp<3; ++comp ) normals[nodeIter + comp * numVertices] = normal[comp];
    }
    MeshWithData<MeshType> meshSaver ( mesh );
    meshSaver.addVectorData ( normals, 3, "normalField", VERTEX_DATA, VECTORS );
    meshSaver.saveAsLegacyVTK( aol::strprintf ( "%s/normalField_%d.vtk", _parser.get<std::string> ( "saveDirectory" ).c_str (), iter ) );
  }
  
  
  ~ShellHandler() {
     delete _xACachePtr;
  };
  
};

} //end namespace


#endif //__SHANDLER_H
