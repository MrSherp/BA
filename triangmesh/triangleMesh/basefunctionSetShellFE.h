#ifndef __BASEFUNCTIONSETSHELLFE_H
#define __BASEFUNCTIONSETSHELLFE_H

namespace shellFE {

//! Inteface
template <typename DataTypeContainer, class QuadRuleType >
class UnitTriangleBaseFunctionSetInterface  {
    
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::DomVecType DomVecType;
    
public:
  UnitTriangleBaseFunctionSetInterface( ) {}

  int numQuadPoints( ) const {
//     return _quadRule.getNumQuadPoints();
    return QuadRuleType::numQuadPoints;
  }

  inline RealType getWeight ( int QuadPoint ) const {
    return _quadRule.getWeight ( QuadPoint );
  }

  inline const DomVecType& getRefCoord ( int QuadPoint ) const {
    return _quadRule.getRefCoord ( QuadPoint );
  }

protected:
  QuadRuleType _quadRule;
};
  
  
  
//! Base function set for unit triangle. 
//! Unit triangle embedded in R^2 is given by the three positions (0,0), (1,0) and (0,1).
template <typename DataTypeContainer, typename QuadType, typename TriangleType>
class UnitTriangMeshBaseFunctionSetP1 : public UnitTriangleBaseFunctionSetInterface< DataTypeContainer, QuadType >  {

  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::DomVecType DomVecType;
    
  static RealType _b1   ( const DomVecType &c ) { return 1. - c[0] - c[1]; }
  static RealType _b2   ( const DomVecType &c ) { return c[0]; }
  static RealType _b3   ( const DomVecType &c ) { return c[1]; }

  typedef RealType ( *BASIS_FUNC_TYPE ) ( const DomVecType &RefCoord );
  BASIS_FUNC_TYPE _basis[3];
  
  const typename DataTypeContainer::TangentVecType _d1b, _d2b;

public:
  UnitTriangMeshBaseFunctionSetP1(  ) : _d1b(-1.,1.,0.), _d2b(-1.,0.,1.) {
    _basis[0] = _b1;
    _basis[1] = _b2;
    _basis[2] = _b3;
  }
  
  enum { numBaseFuncs = 3 };

  void setTriangle ( const TriangleType &/*T*/ ) { }
  
  void evaluateGradient ( int BaseFuncNum, const DomVecType &/*RefCoord*/, DomVecType &Gradient ) const {
    Gradient[0] = _d1b[BaseFuncNum];
    Gradient[1] = _d2b[BaseFuncNum];
  }

  inline DomVecType evaluateGradient ( int BaseFuncNum, int /*QuadPoint*/ ) const {
    return DomVecType( _d1b[BaseFuncNum], _d2b[BaseFuncNum] );
  }

  RealType evaluate ( int BaseFuncNum, const DomVecType &RefCoord ) const {
    return _basis[BaseFuncNum] ( RefCoord );
  }

  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return evaluate ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ) );
  }
};



} // end namespace




#endif //__BASEFUNCTIONSETSHELLFE_H
