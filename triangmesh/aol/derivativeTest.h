#ifndef __DERIVATIVETESTER_H
#define __DERIVATIVETESTER_H

#include <general.h>
#include <EigenIncludes.h>

// Compare DE(testPoint)(testDirection) with diffQuotient 1/h ( E(testPoint + h testDirection ) - E(testPoint) )
template<typename RealType, typename VectorType >
class FirstDerivativeTester {

protected:

  const aol::Op<VectorType, RealType > &_E;
  const aol::Op<VectorType, VectorType> &_DE;
  const RealType _stepSize;
  const VectorType _testPoint;
  const int _numSteps;

public:
  FirstDerivativeTester ( const aol::Op<VectorType, RealType > & E, const aol::Op<VectorType, VectorType> & DE,
                          const RealType stepSize, const int numSteps,
                          const VectorType & testPoint ) 
  : _E ( E ),
    _DE ( DE ),
    _stepSize ( stepSize ),
    _testPoint ( testPoint ),
    _numSteps ( numSteps ) {}        
   
  RealType testSingleDirection ( const VectorType & testDirection ) const {
    
    // evalualte J(m)
    RealType energy;
    _E.apply( _testPoint, energy );
    
    //evaluate J(m + h testDirection )
    RealType energyShifted;
    VectorType shiftedPoint = _testPoint + _stepSize * testDirection;    
    _E.apply( shiftedPoint, energyShifted );
    RealType diffQuotient = ( energyShifted - energy ) / (_stepSize ) ;
    
    return diffQuotient;

  }
  
  
  void testAllDirections( ){
      unsigned numDirections = _testPoint.size();
      
      //evaluate DJ(m)
      VectorType derivative ( numDirections );
      _DE.apply( _testPoint, derivative );
      
      // diffQuotient
      VectorType approxDerivative ( numDirections );
      
      VectorType testDirection( numDirections );
      for( unsigned i=0; i<numDirections; ++i ){
        testDirection.setZero(); testDirection[i] = 1.;
        approxDerivative[i] = testSingleDirection ( testDirection );
      }
      
//       cout << "gradient = " << endl << derivative.transpose() << endl << endl;
      printVector<VectorType> ( derivative, 9, "gradient" );
//       cout << "diffQuotient = " << endl << approxDerivative.transpose() << endl;
      printVector<VectorType> ( approxDerivative, 9, "diffQuotient" );
      VectorType diff = derivative - approxDerivative;
      printVector<VectorType> ( diff, 9, "grad - diffQuotient" );
//       cout << endl << endl << "grad - diffQuotient = " << endl << diff.transpose() << endl << endl;
  }
                      
};



template<typename RealType, typename VectorType, typename MatrixType >
class SecondDerivativeTester {

protected:

  const aol::Op<VectorType, VectorType> &_DE;
  //const aol::Op<VectorType, MatrixType > &_D2E;
  const aol::PreparedOp<VectorType, MatrixType > &_D2E;
  const RealType _stepSize;
  const VectorType _testPoint;
  const int _numSteps;

public:
  SecondDerivativeTester ( const aol::Op<VectorType, VectorType> & DE, const aol::PreparedOp<VectorType, MatrixType> & D2E,
                           const RealType stepSize, const int numSteps,
                           const VectorType & testPoint ) 
  : _DE ( DE ),
    _D2E ( D2E ),
    _stepSize ( stepSize ),
    _testPoint ( testPoint ),
    _numSteps ( numSteps ) {}        
   
  RealType testSingleDirection ( const VectorType & testDirection1, const VectorType & testDirection2 ) const {
      
    // evalualte DJ(m)
    VectorType derivative ( testDirection1.size() );
    _DE.apply( _testPoint, derivative );
    
    //evaluate DJ(m + h testDirection )
    VectorType derivativeShifted ( testDirection1.size() );
    VectorType shiftedPoint = _testPoint + _stepSize * testDirection2;    
    _DE.apply( shiftedPoint, derivativeShifted );
    RealType diffQuotient = ( derivativeShifted.dot(testDirection1) - derivative.dot(testDirection1) ) / (_stepSize ) ;
    
    return diffQuotient;
  }
  
  
  void testAllDirections( ) const {
      unsigned numDirections = _testPoint.size();
      //evaluate D2J(m)
      MatrixType hessian ( numDirections, numDirections ); 
      Eigen::MatrixXd approxHessian (numDirections, numDirections );
      _D2E.apply( _testPoint, hessian );
      VectorType testDirection1( numDirections ), testDirection2( numDirections );
      for( unsigned i=0; i<numDirections; ++i ){
          testDirection1.setZero();
          testDirection1[i] = 1.;
          for( unsigned j=0; j<numDirections; ++j ){
              testDirection2.setZero();
              testDirection2[j] = 1.;
              approxHessian(i,j) = testSingleDirection ( testDirection1, testDirection2 );
          }
      }
      cout << "hessian = " << endl << Eigen::MatrixXd( hessian ) << endl << endl;
      cout << "diffQuotient = " << endl << approxHessian << endl;
      cout << endl << endl << endl << "hessian - diffQuotient = " << endl << Eigen::MatrixXd( hessian ) - approxHessian << endl;
  }
              
              
  void testAllDirectionsWithPreparedSolver( ) const {
      unsigned numDirections = _testPoint.size();
      //evaluate D2J(m)
      Eigen::MatrixXd approxHessian (numDirections, numDirections );
      VectorType testDirection1( numDirections ), testDirection2( numDirections );
      
      _D2E.prepare( _testPoint );
      cout << endl << endl << "compute with prepared solver " << endl;
      for( unsigned i=0; i<numDirections; ++i ){
          testDirection1.setZero();
          testDirection1[i] = 1.;
           _D2E.applyIntoDirection( testDirection1, testDirection2  );
           cout << testDirection2.transpose() << endl;
      }
      
      for( unsigned i=0; i<numDirections; ++i ){
          testDirection1.setZero();
          testDirection1[i] = 1.;
          for( unsigned j=0; j<numDirections; ++j ){
              testDirection2.setZero();
              testDirection2[j] = 1.;
              approxHessian(i,j) = testSingleDirection ( testDirection1, testDirection2 );
          }
      }
      cout << "diffQuotient = " << endl << approxHessian << endl;
  }
};






template<typename RealType, typename VectorType, typename MatrixType >
class SecondDerivativeTesterOld {

protected:

  const aol::Op<VectorType, VectorType> &_DE;
  const aol::Op<VectorType, MatrixType > &_D2E;
  //const aol::PreparedOp<VectorType, MatrixType > &_D2E;
  const RealType _stepSize;
  const VectorType _testPoint;
   
public:
  SecondDerivativeTesterOld ( const aol::Op<VectorType, VectorType> & DE, const aol::Op<VectorType, MatrixType> & D2E,
                              const RealType stepSize, const VectorType & testPoint ) 
  : _DE ( DE ),
    _D2E ( D2E ),
    _stepSize ( stepSize ),
    _testPoint ( testPoint ) {}        
   
  RealType testSingleDirection ( const VectorType & testDirection1, const VectorType & testDirection2 ) const {
      
    // evalualte DJ(m)
    VectorType derivative ( testDirection1.size() );
    _DE.apply( _testPoint, derivative );
    
    //evaluate DJ(m + h testDirection )
    VectorType derivativeShifted ( testDirection1.size() );
    VectorType shiftedPoint = _testPoint + _stepSize * testDirection2;    
    _DE.apply( shiftedPoint, derivativeShifted );
    RealType diffQuotient = ( derivativeShifted.dot(testDirection1) - derivative.dot(testDirection1) ) / (_stepSize ) ;
    
    return diffQuotient;
  }
  
  
  void testAllDirections( ) const {
      unsigned numDirections = _testPoint.size();
      //evaluate D2J(m)
      MatrixType hessian ( numDirections, numDirections ); 
      Eigen::MatrixXd approxHessian ( numDirections, numDirections );
      _D2E.apply( _testPoint, hessian );
      VectorType testDirection1( numDirections ), testDirection2( numDirections );
      for( unsigned i=0; i<numDirections; ++i ){
          testDirection1.setZero();
          testDirection1[i] = 1.;
          for( unsigned j=0; j<numDirections; ++j ){
              testDirection2.setZero();
              testDirection2[j] = 1.;
              approxHessian(i,j) = testSingleDirection ( testDirection1, testDirection2 );
          }
      }
      
      //output
//       for( int argComp = 0; argComp < 3; ++argComp )
//           for( int destComp = 0; destComp < 3; ++destComp ){
//             cout << endl << endl << "===============================================================================================================" << endl
//                 << "argComp = " << argComp << " , destComp = " << destComp << endl << endl;;
//             cout << "hessian = " << endl << Eigen::MatrixXd( hessian ).block( argComp * _numGlobalDofs, destComp * _numGlobalDofs, _numGlobalDofs, _numGlobalDofs ) << endl << endl;
//             cout << "diffQuotient = " << endl << approxHessian.block( argComp * _numGlobalDofs, destComp * _numGlobalDofs, _numGlobalDofs, _numGlobalDofs ) << endl;
//             cout << endl << endl << endl << "hessian - diffQuotient = " << endl << ( Eigen::MatrixXd( hessian ) - approxHessian ).block( argComp * _numGlobalDofs, destComp * _numGlobalDofs, _numGlobalDofs, _numGlobalDofs ) << endl;
//           }

      cout << "hessian = " << endl << Eigen::MatrixXd( hessian ) << endl << endl;
      cout << "diffQuotient = " << endl << approxHessian << endl;
      cout << endl << endl << endl << "hessian - diffQuotient = " << endl << Eigen::MatrixXd( hessian ) - approxHessian << endl;
      
  }

};









// template <typename DataTypeContainer >
// void DerivativeTest ( const int DerivativeTestNumber ) {
// 
//     VectorType testMaterial = VectorType::Random( initMaterial.size() );
//         
//     // First Derivative Test
//     if( DerivativeTestNumber == 1 ){
//         FirstDerivativeTester<RealType, VectorType> ( energyOp, gradientOp, _parser.get<double>( "StepSizeDerivativeTest" ), _parser.get<int>( "NumStepsDerivativeTest" ), testMaterial ).testAllDirections();
//         return;
//      }
//      
//      if( DerivativeTestNumber == 2 ){
//          SecondDerivativeTester<RealType, VectorType, SparseMatrixType> ( gradientOp, hessianOp, _parser.get<double>( "StepSizeDerivativeTest" ), _parser.get<int>( "NumStepsDerivativeTest" ), testMaterial ).testAllDirections();
//           return;
//      }
//         if( DerivativeTestNumber == 3 ){
//           cout << endl << endl << endl << "==================================================================================" << endl
//                                  << "=====  First Derivative Test ================== "        << endl
//                                  << "==================================================================================" << endl << endl;
//           FirstDerivativeTester<RealType, VectorType> ( energyOp, gradientOp, _parser.get<double>( "StepSizeDerivativeTest" ), _parser.get<int>( "NumStepsDerivativeTest" ), testMaterial ).testAllDirections();
//          
//           cout << endl << endl << endl << "==================================================================================" << endl
//                                  << "=====  Second Derivative Test ================== "        << endl
//                                  << "==================================================================================" << endl << endl;
//           SecondDerivativeTester<RealType, VectorType, SparseMatrixType> ( gradientOp, hessianOp, _parser.get<double>( "StepSizeDerivativeTest" ), _parser.get<int>( "NumStepsDerivativeTest" ), testMaterial ).testAllDirections();
//           
//           cout << endl << endl << endl << "==================================================================================" << endl
//                                  << "=====  Second Derivative Test with prepared solver ================== "        << endl
//                                  << "==================================================================================" << endl << endl;
//           SecondDerivativeTester<RealType, VectorType, SparseMatrixType> ( gradientOp, hessianOp, _parser.get<double>( "StepSizeDerivativeTest" ), _parser.get<int>( "NumStepsDerivativeTest" ), testMaterial ).testAllDirectionsWithPreparedSolver();
//           
//           return;
//         }
// 
// }


#endif //__DERIVATIVETESTER_H