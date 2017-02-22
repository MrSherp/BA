#ifndef __DISCRETEFUNCTIONSHELLFE_H
#define __DISCRETEFUNCTIONSHELLFE_H

#include <initializer_list>
#include <EigenIncludes.h>
#include <tensor.h>

namespace shellFE {
  
template <typename ConfiguratorType, typename VectorType, typename NLTYPE >
struct DiscreteFunctionLookup {
  
  void getLocalDof ( const ConfiguratorType & conf, const VectorType &dofs, const typename ConfiguratorType::ElementType &El,  
                     const int bfNum,  const NLTYPE &bfValue, NLTYPE &aux ) {
      aux = bfValue;
      aux *= dofs[ conf.localToGlobal ( El, bfNum ) ];
  }

};
  

template <typename ConfiguratorType, typename VectorType = typename ConfiguratorType::VectorType >
class DiscreteFunctionDefaultShellFE {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType  DomVecType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Matrix22  Matrix22;
  typedef typename ConfiguratorType::ElementType ElementType;
//   typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::BaseFuncSetType BaseFuncSetType;
  //typedef typename ConfiguratorType::ApproxGradientBaseFuncSetType ApproxGradientBaseFuncSetType; //TODO: only for DKT

  const ConfiguratorType & _conf;
  const VectorType & _dofs;

  DiscreteFunctionDefaultShellFE ( const ConfiguratorType & config, const VectorType & Dofs )
    : _conf ( config ),
      _dofs ( Dofs ) { }
    
  RealType evaluate ( const ElementType &El, const DomVecType& RefCoord ) const {
    const BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
    RealType w = 0.;
    RealType v, aux;
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      v = bfs.evaluate ( b, RefCoord );
      DiscreteFunctionLookup<ConfiguratorType, VectorType, RealType>().getLocalDof( _conf, _dofs, El, b, v, aux );
      w += aux;
    }
    return w;
  }

  RealType evaluateAtQuadPoint ( const ElementType &El, int QuadPoint ) const {
    const BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
    RealType w = 0.;
    RealType v, aux;
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      v = bfs.evaluate ( b, QuadPoint );
      DiscreteFunctionLookup<ConfiguratorType, VectorType, RealType>().getLocalDof ( _conf, _dofs, El, b, v, aux );
      w += aux;
    }
    return w;
  }
  
  void evaluateGradient ( const ElementType &El, const DomVecType& RefCoord, DomVecType& Grad ) const {
    const BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
    DomVecType v, aux;
    Grad.setZero();
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      bfs.evaluateGradient ( b, RefCoord, v );
      DiscreteFunctionLookup<ConfiguratorType, VectorType, DomVecType>().getLocalDof ( _conf, _dofs, El, b, v, aux );
      Grad += aux;
    }
  }
  
  void evaluateGradientAtQuadPoint ( const ElementType &El, int QuadPoint, DomVecType& Grad ) const {
    const BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
    DomVecType v, aux;
    Grad.setZero();
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      v = bfs.evaluateGradient ( b, QuadPoint );
      DiscreteFunctionLookup<ConfiguratorType, VectorType, DomVecType >().getLocalDof ( _conf, _dofs, El, b, v, aux );
      Grad += aux;
    }
  }
  
  //! evaluation of approximative gradient and its derivative
  
  void evaluateApproxGradient ( const ElementType &El, const DomVecType& RefCoord, DomVecType& Grad ) const {
    const typename ConfiguratorType::ApproxGradientBaseFuncSetType &bfs = _conf.getApproxGradientBaseFunctionSet ( El );
    DomVecType v, aux;
    Grad.setZero();
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      bfs.evaluateApproxGradient ( b, RefCoord, v );
      DiscreteFunctionLookup<ConfiguratorType, VectorType, DomVecType>().getLocalDof ( _conf, _dofs, El, b, v, aux );
      Grad += aux;
    }
  }
  
  void evaluateApproxGradientAtQuadPoint ( const ElementType &El, int QuadPoint, DomVecType& Grad ) const {
    const typename ConfiguratorType::ApproxGradientBaseFuncSetType &bfs = _conf.getApproxGradientBaseFunctionSet ( El );
    DomVecType v, aux;
    Grad.setZero();
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      v = bfs.evaluateApproxGradient ( b, QuadPoint );
      DiscreteFunctionLookup<ConfiguratorType, VectorType, DomVecType >().getLocalDof ( _conf, _dofs, El, b, v, aux );
      Grad += aux;
    }
  }
  
  
  //returns HessianVec = ( d_x theta_1 , d_y theta_2, d_y theta_1 + d_x theta_2  )  
  void evaluateApproxHessianAsVecAtQuadPoint ( const ElementType &El, int QuadPoint, TangentVecType & Hessian ) const {
    const typename ConfiguratorType::ApproxGradientBaseFuncSetType &bfs = _conf.getApproxGradientBaseFunctionSet ( El );
    TangentVecType v, aux;
    Hessian.setZero();
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      v = bfs.evaluateApproxHessianAsVec ( b, QuadPoint );
      DiscreteFunctionLookup<ConfiguratorType, VectorType, TangentVecType >().getLocalDof ( _conf, _dofs, El, b, v, aux );
      Hessian += aux;
    }
  }
  
  //returns Hessian = ( d_x theta_1 ,  d_y theta_1 )
  //                    d_x theta_2 ,  d_y theta_2 )  
  void evaluateApproxHessianAtQuadPoint ( const ElementType &El, int QuadPoint, Matrix22 & Hessian ) const {
    const typename ConfiguratorType::ApproxGradientBaseFuncSetType &bfs = _conf.getApproxGradientBaseFunctionSet ( El );
    Matrix22 v, aux;
    Hessian.setZero();
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      v = bfs.evaluateApproxHessian ( b, QuadPoint );
      DiscreteFunctionLookup<ConfiguratorType, VectorType, Matrix22 >().getLocalDof ( _conf, _dofs, El, b, v, aux );
      Hessian += aux;
    }
  }
  
  //returns Hessian = ( d_x theta_1,                       1/2 (d_y theta_1 + d_x theta_2)
  //                    1/2( d_y theta_1 + d_x theta_2) ,  d_y theta_2,  )  
  void evaluateApproxHessianSymAtQuadPoint ( const ElementType &El, int QuadPoint, Matrix22 & Hessian ) const {
    const typename ConfiguratorType::ApproxGradientBaseFuncSetType &bfs = _conf.getApproxGradientBaseFunctionSet ( El );
    Matrix22 v, aux;
    Hessian.setZero();
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      v = bfs.evaluateApproxHessianSym ( b, QuadPoint );
      DiscreteFunctionLookup<ConfiguratorType, VectorType, Matrix22 >().getLocalDof ( _conf, _dofs, El, b, v, aux );
      Hessian += aux;
    }
  }
  
  const VectorType& getDofs (  ) const {return _dofs;}
};

template <typename ConfiguratorType>
class DiscreteVectorFunctionDefaultShellFE {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType  DomVecType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::Matrix22  Matrix22;
  typedef typename ConfiguratorType::Matrix32  Matrix32;
  typedef typename ConfiguratorType::Matrix33  Matrix33;
  
  typedef typename ConfiguratorType::Tensor222Type  Tensor222Type;
  typedef typename ConfiguratorType::Tensor322Type  Tensor322Type;
  
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef DiscreteFunctionDefaultShellFE<ConfiguratorType, Eigen::Ref<const VectorType> > DiscFuncType;

  const int _numComponents;
  mutable std::vector<Eigen::Ref<const VectorType> > _refs;
  mutable std::vector<DiscFuncType> _discrFuncs;

public:
  
  DiscreteVectorFunctionDefaultShellFE ( const ConfiguratorType &Configurator,  const VectorType &Dofs, const int numComponents ) :
      _numComponents ( numComponents )
  {
     const int numGlobalDofs = Configurator.getNumGlobalDofs();
     _refs.reserve( _numComponents );
     _discrFuncs.reserve( _numComponents );

     _refs.push_back ( Dofs.segment( 0, numGlobalDofs) );
     _refs.push_back ( Dofs.segment( numGlobalDofs, numGlobalDofs) );
     _refs.push_back ( Dofs.segment( 2 * numGlobalDofs, numGlobalDofs) );
     
     _discrFuncs.push_back ( DiscFuncType ( Configurator, _refs[0] ) );
     _discrFuncs.push_back ( DiscFuncType ( Configurator, _refs[1] ) );
     _discrFuncs.push_back ( DiscFuncType ( Configurator, _refs[2] ) );
 
  }

  void evaluate( const ElementType &El, const DomVecType& RefCoord, TangentVecType &Value ) const {
    for ( int c = 0; c < 3; ++c )
      Value[c] = _discrFuncs[c].evaluate ( El, RefCoord );
  }
  
  void evaluateAtQuadPoint ( const ElementType &El, int QuadPoint, TangentVecType &Value ) const {
    for ( int c = 0; c < 3; ++c )
      Value[c] = _discrFuncs[c].evaluateAtQuadPoint ( El, QuadPoint );
  }
  
  // e.g. for Lagrange-Multiplier
//   void evaluateAsSymMatAtQuadPoint ( const ElementType &El, int QuadPoint, Matrix22 &SymMat ) const {
//     for ( int c = 0; c < 3; c++ ) {
//       SymMat[0][0] = _discrFuncs[0].evaluateAtQuadPoint ( El, QuadPoint );
//       SymMat[0][1] = SymMat[1][0] = _discrFuncs[0].evaluateAtQuadPoint ( El, QuadPoint );
//       SymMat[1][1] = _discrFuncs[2].evaluateAtQuadPoint ( El, QuadPoint );
//     }
//   }

  void evaluateGradient ( const ElementType &El, const DomVecType& RefCoord, Matrix32 &Dx ) const {
    DomVecType v;
    for ( int c = 0; c < 3; c++ ) {
      _discrFuncs[c].evaluateGradient ( El, RefCoord, v );
      Dx( c, 0 ) = v[0];
      Dx( c, 1 ) = v[1];
    }
  }

  void evaluateGradientAtQuadPoint ( const ElementType &El, int QuadPoint, Matrix32 &Dx ) const {
    DomVecType v;
    for ( int c = 0; c < 3; ++c ) {
      _discrFuncs[c].evaluateGradientAtQuadPoint ( El, QuadPoint, v );
      Dx( c, 0 ) = v[0];
      Dx( c, 1 ) = v[1];
    }
  }
  
  void evaluateFirstFundamentalForm ( const ElementType &El, const DomVecType& RefCoord, Matrix22 &g ) const {
    Matrix32 Dx;
    evaluateGradient( El, RefCoord, Dx );
    g = Dx.transpose() * Dx;
  }
  
  
  void evaluateFirstFundamentalFormAtQuadPoint ( const ElementType &El, int QuadPoint, Matrix22 &g ) const {
    Matrix32 Dx;
    evaluateGradientAtQuadPoint( El, QuadPoint, Dx );
    g = Dx.transpose() * Dx;
  }
  
  void evaluateFirstFundamentalForm ( const Matrix32 &Dx, Matrix22 &g ) const {
    g = Dx.transpose() * Dx;
  }
  
  // given chart X -> compute n \circ X  = d_1 X \cross d_2 X / |.|
  void evaluateNormalAtQuadPointOnShell ( const ElementType &El, int QuadPoint, TangentVecType &normal ) const {
    Matrix32 Dx;
    evaluateGradientAtQuadPoint ( El, QuadPoint, Dx );
    TangentVecType unNormalizedNormal ( (Dx.col(0)).cross(Dx.col(1)) );
    normal = unNormalizedNormal.normalized();
  }
  
  void evaluateNormal ( const Matrix32 &Dx, TangentVecType &normal ) const {
    TangentVecType unNormalizedNormal ( (Dx.col(0)).cross(Dx.col(1)) );
    normal = unNormalizedNormal.normalized();
  }
  

  //! evaluation of approximative gradient and its derivative
  
  void evaluateApproxGradient ( const ElementType &El, const DomVecType& RefCoord, Matrix32 &Dx ) const {
    DomVecType v;
    for ( int c = 0; c < 3; c++ ) {
      _discrFuncs[c].evaluateApproxGradient ( El, RefCoord, v );
      Dx( c, 0 ) = v[0];
      Dx( c, 1 ) = v[1];
    }
  }
  
  void evaluateApproxGradientAtQuadPoint ( const ElementType &El, int QuadPoint, Matrix32 &Dx ) const {
    DomVecType v;
    for ( int c = 0; c < 3; c++ ) {
      _discrFuncs[c].evaluateApproxGradientAtQuadPoint ( El, QuadPoint, v );
      Dx( c, 0 ) = v[0];
      Dx( c, 1 ) = v[1];
    }
  }
  
  void evaluateApproxFirstFundamentalFormAtQuadPoint ( const ElementType &El, int QuadPoint, Matrix22 &g ) const {
    Matrix32 Dx;
    evaluateApproxGradientAtQuadPoint( El, QuadPoint, Dx );
    g = Dx.transpose() * Dx;
  }

  
  // ddX = ( HVec_comp0      i.e. component corresponds to column
  //         HVec_comp1
  //         HVec_comp2 )
  void evaluateApproxHessianAsVecAtQuadPoint ( const ElementType &El, int QuadPoint, Matrix33 & ddX ) const {
    TangentVecType hessian_c;
    for ( int c = 0; c < 3; ++c ) {
      _discrFuncs[c].evaluateApproxHessianAsVecAtQuadPoint ( El , QuadPoint, hessian_c );
      for( int j=0; j<3; ++j )
          ddX( c, j ) = hessian_c[j];
    }
  }

  void evaluateApproxHessianAtQuadPoint ( const ElementType &El, int QuadPoint, Tensor322Type & ddX ) const {
    Matrix22 hessian_c;
    for ( int c = 0; c < 3; ++c ) {
      _discrFuncs[c].evaluateApproxHessianAtQuadPoint ( El , QuadPoint, hessian_c );
      for( int j=0; j<2; ++j )
        for( int k=0; k<2; ++k )
          ddX.set ( c, j, k, hessian_c( j, k ) );
    }
  }
  
  void evaluateApproxHessianSymAtQuadPoint ( const ElementType &El, int QuadPoint, Tensor322Type & ddX ) const {
    Matrix22 hessian_c;
    for ( int c = 0; c < 3; ++c ) {
      _discrFuncs[c].evaluateApproxHessianSymAtQuadPoint ( El , QuadPoint, hessian_c );
      for( int j=0; j<2; ++j )
        for( int k=0; k<2; ++k )
          ddX.set ( c, j, k, hessian_c( j, k ) );
    }
  }

  //! compute the shape tensor (second fundamental form), the product of the second derivatives with the unit normal
  //TODO choise of Hessian? AsVec? Sym?
  void evaluateSecondFundamentalFormAtQuadPoint ( const ElementType &El, int QuadPoint, Matrix22& shapeTensor ) const {
      Tensor322Type ddX;
      evaluateApproxHessianSymAtQuadPoint ( El, QuadPoint, ddX );
      TangentVecType normal;
      evaluateNormalAtQuadPointOnShell( El, QuadPoint, normal );
      for ( int i = 0; i < 2; ++i )
        for ( int j = 0; j < 2; ++j ) {
          shapeTensor( i , j ) = normal.dot( ddX.getVector ( i, j ) );
        }
  }
  
  void evaluateSecondFundamentalForm ( const Tensor322Type &ddX, const TangentVecType &normal, Matrix22& shapeTensor ) const {
      for ( int i = 0; i < 2; ++i )
        for ( int j = 0; j < 2; ++j ) {
          TangentVecType ddX_ij;
          ddX.getVector( ddX_ij, i, j );
          shapeTensor( i , j ) = normal.dot( ddX_ij );
        }
  }
  
  
  //! computes derivative of g: Deriv_ijk = \partial_k g_ij = \sum_comp  \partial_ik x^comp \partial_j x^comp + \partial_i x^comp \partial_jk x^comp
  // TODO choise of Hessian? approximation of first derivative?
  void evaluateApproxDerivativeOfFirstFundamentalFormAtQuadPoint ( const ElementType &El, int QuadPoint, Tensor222Type& Deriv ) const {
      Matrix32 dX;
      evaluateGradientAtQuadPoint( El, QuadPoint, dX );
      Tensor322Type ddX;
      evaluateApproxHessianSymAtQuadPoint ( El, QuadPoint, ddX );
      for ( int i = 0; i < 2; ++i )
        for ( int j = 0; j < 2; ++j )
          for( int k=0; k<2; ++k ){
            RealType tmp = 0.0;
            for( int comp=0; comp<3; ++comp ){
              tmp += ddX( comp, i, k ) * dX( comp, j ) + dX( comp, i ) * ddX( comp, j, k );  
            }
            Deriv.set ( i , j , k, tmp );
        }
  }
  
  void evaluateApproxDerivativeOfFirstFundamentalForm ( const Matrix32 &dX, const Tensor322Type &ddX, Tensor222Type& Deriv ) const {
      for ( int i = 0; i < 2; ++i )
        for ( int j = 0; j < 2; ++j )
          for( int k=0; k<2; ++k ){
            RealType tmp = 0.0;
            for( int comp=0; comp<3; ++comp ){
              tmp += ddX( comp, i, k ) * dX( comp, j ) + dX( comp, i ) * ddX( comp, j, k );  
            }
            Deriv.set ( i , j , k, tmp );
        }
  }
  
  //! Gamma_ijk = 1/2 ( \partial_j g_ki + \partial_i g_kj - \partial_k g_ij )
  void evaluateApproxChristoffelSymbolsOfFirstKindAtQuadPoint ( const ElementType &El, int QuadPoint, Tensor222Type& ChristoffelSym ) const {
      Tensor222Type Dg;
      evaluateApproxDerivativeOfFirstFundamentalFormAtQuadPoint ( El, QuadPoint, Dg );
      for ( int i = 0; i < 2; ++i )
        for ( int j = 0; j < 2; ++j )
          for( int k=0; k<2; ++k ){
            ChristoffelSym.set ( i , j , k, 0.5 * ( Dg( k, i, j) + Dg(k, j, i) - Dg(i,j,k) )   );
        }
  }
  
  void evaluateApproxChristoffelSymbolsOfFirstKind ( const Matrix32 &dX, const Tensor322Type &ddX, Tensor222Type& ChristoffelSym ) const {
      Tensor222Type Dg;
      evaluateApproxDerivativeOfFirstFundamentalForm ( dX, ddX, Dg );
      for ( int i = 0; i < 2; ++i )
        for ( int j = 0; j < 2; ++j )
          for( int k=0; k<2; ++k ){
            ChristoffelSym.set ( i , j , k, 0.5 * ( Dg( k, i, j) + Dg(k, j, i) - Dg(i,j,k) )   );
        }
  }
  
  //! Gamma_ij^m = \sum_k g^m,k Gamma_i,j,k = \sum_k 1/2 g^m,k ( \partial_j g_ki + \partial_i g_kj - \partial_k g_ij )
  void evaluateApproxChristoffelSymbolsOfSecondKindAtQuadPoint ( const ElementType &El, int QuadPoint, Tensor222Type& ChristoffelSym2 ) const {
    Matrix22 g, ginv;
    evaluateFirstFundamentalFormAtQuadPoint (El, QuadPoint, g );
    ginv = g.inverse();
    Tensor222Type ChristoffelSym1;
    evaluateApproxChristoffelSymbolsOfFirstKindAtQuadPoint ( El, QuadPoint, ChristoffelSym1 );
    for ( int i = 0; i < 2; ++i )
        for ( int j = 0; j < 2; ++j )
          for( int m=0; m<2; ++m ){
            RealType tmp = 0.0;
            for( int k=0; k<2; ++k ){
              tmp += 0.5 * ginv( m, k ) * ChristoffelSym1( i, j, k );  
            }
            ChristoffelSym2.set ( i , j , m, tmp  );
        }
  }
  
  void evaluateApproxChristoffelSymbolsOfSecondKind ( const Matrix22 &ginv, const Matrix32 &dX, const Tensor322Type &ddX, Tensor222Type &ChristoffelSym2 ) const {
    Tensor222Type ChristoffelSym1;
    evaluateApproxChristoffelSymbolsOfFirstKind ( dX, ddX, ChristoffelSym1 );
    for ( int i = 0; i < 2; ++i )
        for ( int j = 0; j < 2; ++j )
          for( int m=0; m<2; ++m ){
            RealType tmp = 0.0;
            for( int k=0; k<2; ++k ){
              tmp += 0.5 * ginv( m, k ) * ChristoffelSym1( i, j, k );  
            }
            ChristoffelSym2.set ( i , j , m, tmp  );
        }
  }
  
  //! g^-1 Gamma
  void evaluateVectorForLaplacianAtQuadPoint ( const ElementType &El, int QuadPoint, DomVecType& vec ) const {
    vec.setZero();
    Matrix22 g, ginv;
    evaluateFirstFundamentalFormAtQuadPoint (El, QuadPoint, g );
    ginv = g.inverse();
    Tensor222Type ChristoffelSym2;
    evaluateApproxChristoffelSymbolsOfSecondKindAtQuadPoint ( El, QuadPoint, ChristoffelSym2 );
    for( int l=0; l<2; ++l )
      for ( int j=0; j<2; ++j )
        for( int k=0; k<2; ++k )
          vec[l] += ginv( j, k ) * ChristoffelSym2( j, k, l );  

  }
  
  void evaluateVectorForLaplacian ( const Matrix22 &ginv, const Matrix32 &dX, const Tensor322Type &ddX, DomVecType& vec ) const {
    vec.setZero();
    Tensor222Type ChristoffelSym2;
    evaluateApproxChristoffelSymbolsOfSecondKind ( ginv, dX, ddX, ChristoffelSym2 );
    for( int l=0; l<2; ++l )
      for ( int j=0; j<2; ++j )
        for( int k=0; k<2; ++k )
          vec[l] += ginv( j, k ) * ChristoffelSym2( j, k, l );  

  }
  
  const DiscFuncType& operator[] ( int i ) const { return _discrFuncs[i];}
  DiscFuncType& operator[] ( int i ) { return _discrFuncs[i];}

};


} // end namespace

#endif
