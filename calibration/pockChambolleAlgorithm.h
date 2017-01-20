#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <math.h> 

template < typename RealType, typename VectorNd >
void signalToIndicatorFunction( VectorNd& Arg){
    VectorNd tmp (Arg);
    const unsigned int N = Arg.size();
    Arg.resize( co::sqr( N ));
    for ( int j = N-2; j>0 ; --j){
        double value =   (double)  (j);
        value /= (double) (N-1);
        for ( int i = 0; i < N; ++i){
            if ( tmp( i ) >= value ){
                Arg( i + (N-j-1) * N ) = 1.;
            }
            else{
                Arg( i + (N-j-1) * N ) = 0.;
            }
        }
    }
    for ( int i = 0; i < N; ++i){
        Arg( i ) = 0;
        Arg( i + ( N - 1 ) * N) = 1;        
    }
}


template < typename RealType, typename VectorNd >
void divideByMaximumOfOneOrNorm( VectorNd& Arg ){
    int numX = Arg.rows() / 2;
    for ( int i = 0; i < numX; ++i){
        RealType scale = 1.;   
        VectorNd tmp(2);
        tmp ( 0 ) = Arg ( i );
        tmp ( 1 ) = Arg ( i + numX ); 
        if( tmp.norm() > 1 ){
            scale /= tmp.norm();
            Arg( i ) *= scale;
            Arg( i + numX ) *= scale;
        }    
    }
}



template < typename RealType, typename VectorNd >
class KProjector {
    VectorNd& _image;
  
public:
    KProjector ( VectorNd& Image )
:   _image ( Image ){ }

void projectOntoK( VectorNd& Arg ){
    int numX = Arg.rows() / 2;
        for ( int i = 0; i < numX; ++i){
            RealType scale = 1.;   
            VectorNd tmp(2);
            tmp ( 0 ) = Arg ( i );
            tmp ( 1 ) = Arg ( i + numX ); 
            if( tmp(0) > 1 ){
                scale /= tmp(0);
                Arg( i ) *= scale;
            }   
            if( tmp(1) < -_image( i ) ){
                Arg( i + numX ) = -_image( i );
            }
        }
    }
};



template < typename RealType, typename VectorNd >
class ROFDataResolvent {
    VectorNd& _image;
    const RealType _lambda;
    const RealType _hscale;
  
public:
    ROFDataResolvent ( VectorNd& Image, RealType Lambda )
:   _image ( Image ), _lambda ( Lambda ), _hscale ( static_cast < RealType > (1) ) { }
  
    void apply ( VectorNd& Arg, const RealType Tau ) const { 
        Arg += Tau * _lambda * _hscale * _image;
        Arg /= (1. + Tau * _lambda * _hscale);
    }
};



template < typename RealType, typename VectorNd >
void projectOntoC( VectorNd& Arg ){
    int N = (int) sqrt( (double) Arg.size());
    if( sqrt ( (double) Arg.size()) > (double) N){
        std::cout << "projectOntoC doesn't work - Arg not quadratic" << std::endl;
    }
    else{
        for( int i = 0; i < N; ++i){
            Arg(i) = 0;
            Arg(i + (N-1) * N) = 1;
        }
    }        
}



template < typename RealType, typename VectorNd >
VectorNd partition(VectorNd& V, int I, int Dim){
    if ( V.size() % Dim != 0 ){
        std::cout << "Dimension does not fit partition!\n" << std::endl;
        return V;
    }
    else{
        int j = V.size() / Dim;
        return V.segment( I  * j, j );
    }

}
    


template < typename RealType, typename OperatorType, typename ResolventDataTerm, typename VectorNd >
void ChambollePockAlgorithm1 ( VectorNd& PrimalSolution,
    VectorNd& DualSolution, const OperatorType& K, ResolventDataTerm& ResolventOfG , const RealType tau, const RealType sigma, const int maxIter, const RealType StopEpsilon, int& enditer, RealType& endepsilon ) {
  VectorNd xBar ( PrimalSolution ) ;                                               //( PrimalSolution, aol::DEEP_COPY );
  VectorNd oldPrimalSolution ( PrimalSolution );                                  //( PrimalSolution, aol::STRUCT_COPY );
  VectorNd adjointOfDual ( PrimalSolution );
  VectorNd gradientOfXBar ( DualSolution );

  int iter;
  for ( iter = 0; iter < maxIter; ++iter ) { 
    
    oldPrimalSolution = PrimalSolution;
    K.apply ( xBar, gradientOfXBar );
    DualSolution += sigma * gradientOfXBar;
    divideByMaximumOfOneOrNorm < RealType > ( DualSolution );
    
    //xBar not needed beyond this point, so an update is fine
    xBar = PrimalSolution;

    K.applyAdjoint ( DualSolution, adjointOfDual );
    PrimalSolution -= tau * adjointOfDual; 
    ResolventOfG.apply ( PrimalSolution, tau );
    xBar *= -1.;                                                                // xBar.scaleAndAddMultiple( - aol::ZOTrait<RealType>::one, PrimalSolution, aol::ZOTrait<RealType>::one + aol::ZOTrait<RealType>::one ); 
    xBar += 2. * PrimalSolution; 
    
    if ( iter % 2500 == 0 && iter != 0 ) {
        VectorNd epsilon = oldPrimalSolution - PrimalSolution;
        cout << "\nIterations: " << iter << "\nEpsilon: " << epsilon.squaredNorm() << std::endl;
    }
    
    if ( iter % 10 == 0 && iter != 0 ) {
      oldPrimalSolution -= PrimalSolution;
     if (oldPrimalSolution.squaredNorm ()  < StopEpsilon ) { //(l_2^2-Norm)
	endepsilon = oldPrimalSolution.squaredNorm ();
        enditer = iter;
        cout << "Number of required iterations: " << iter << endl;
        break;
      }
      oldPrimalSolution = PrimalSolution;
    }
  }
  if ( iter == maxIter){ 
    enditer = maxIter;
    oldPrimalSolution -= PrimalSolution;
    endepsilon = oldPrimalSolution.squaredNorm ();
  }
  cout << "\ntau  : " << tau << "\nsigma: " << sigma << endl;
}



template < typename RealType, typename OperatorType, typename ResolventDataTerm, typename VectorNd >
void ChambollePockAlgorithm2 ( VectorNd& PrimalSolution,
    VectorNd& DualSolution, const OperatorType& K, ResolventDataTerm& ResolventOfG, const RealType gamma, RealType& tau, RealType& sigma, const int maxIter, const RealType StopEpsilon, int& enditer, RealType& endepsilon ) {
  VectorNd xBar ( PrimalSolution );
  VectorNd oldPrimalSolution ( PrimalSolution );
  VectorNd adjointOfDual ( PrimalSolution );
  VectorNd gradientOfXBar ( DualSolution );
  RealType theta = 1.;                                                          // aol::ZOTrait<RealType>::one;

  int iter;
  for ( iter = 0; iter < maxIter; ++iter ) {
      
    oldPrimalSolution = PrimalSolution;
    K.apply ( xBar, gradientOfXBar );
    DualSolution += sigma * gradientOfXBar;
    divideByMaximumOfOneOrNorm < RealType > ( DualSolution );

    //xBar not needed beyond this point, so an update is fine
    xBar = PrimalSolution;

    K.applyAdjoint ( DualSolution, adjointOfDual );
    PrimalSolution -= tau * adjointOfDual;
    ResolventOfG.apply ( PrimalSolution, tau );
    
    theta = 1. / std::sqrt ( 1. + static_cast < RealType > ( 2 ) * gamma * tau );
    tau *= theta;
    sigma /= theta;
    xBar *= -theta;
    xBar += ( 1. + theta ) * PrimalSolution;                                    //( - theta, PrimalSolution, aol::ZOTrait<RealType>::one + theta );
    
    if ( iter % 10 == 0 && iter != 0 ) {
      oldPrimalSolution -= PrimalSolution;
     if ( oldPrimalSolution.squaredNorm () < StopEpsilon ) { //(l_2^2-Norm)
	endepsilon = oldPrimalSolution.squaredNorm ();
	enditer = iter;
        cout << "Number of required iterations: " << iter << endl;
        break;
      }
      oldPrimalSolution = PrimalSolution;
    }
  }
  if ( iter == maxIter){ 
    enditer = maxIter;
    oldPrimalSolution -= PrimalSolution;
    endepsilon = oldPrimalSolution.squaredNorm ();
  }
  cout << "\ntheta: " << theta << "\ntau  : " << tau << "\nsigma: " << sigma << "\ngamma: " << gamma << endl;
  cout << "\ntau*gamma:" << tau*gamma << endl;
}



template < typename RealType, typename VectorNd > 
class ForwardFD { 
  typedef Triplet<RealType> TripletType;

  const int _numX;
  Eigen::SparseMatrix < RealType>  _matX ;
  Eigen::SparseMatrix < RealType>  _matY ;

public: 
  
  ForwardFD ( const int NumX ) : _numX ( NumX ){ 
    
    generateX ( _matX );                                                                                 
    generateY ( _matY );
  }
  
  void applyAdd ( const VectorNd& Arg, VectorNd& Dest ) const {  
    Dest.head( _numX * _numX) += _matX * Arg;                                                                                             //_matX->applyAdd ( Arg, Dest[0] );
    Dest.tail( _numX * _numX) += _matY * Arg;
  }
  
  void apply ( const VectorNd& Arg, VectorNd& Dest ) const {
    Dest.setZero();
    this->applyAdd ( Arg, Dest );
  }

  void applyAddAdjoint ( const VectorNd& Arg, VectorNd& Dest ) const {
    Dest += _matX.transpose() * Arg.head( _numX * _numX );
    Dest += _matY.transpose() * Arg.tail( _numX * _numX );
  }
  
  void applyAdjoint ( const VectorNd& Arg, VectorNd& Dest ) const {
    Dest.setZero ( );
    this->applyAddAdjoint ( Arg, Dest );
  }
  
  const RealType getX ( const int i, const int j) const {
    return _matX.coeff(i,j);  
  }
  
  const RealType getY ( const int i, const int j) const {
    return _matY.coeff(i,j);  
  }
  
private:
void generateX ( Eigen::SparseMatrix < RealType >& MatX  ) {
    std::vector < TripletType > tripletList;
    tripletList.reserve ( 2 * co::sqr( _numX ) );
    int nonZeros = 0;
    const RealType oneOverH = 1. / static_cast < RealType > ( _numX - 1 );            
    
    for ( int x = 0; x < _numX - 1; ++x ) {
      for ( int y = 0; y < _numX - 1; ++y ) {
        const int currentPosition = y * _numX + x;
        tripletList.push_back( TripletType ( currentPosition, currentPosition, -oneOverH ) ); 
        tripletList.push_back( TripletType ( currentPosition, currentPosition + 1, oneOverH ) );
        nonZeros += 2;
      }
    }
    for ( int y = 0; y < _numX -2; ++y ) {
      const int currentPosition = y * _numX + (_numX - 1);
      tripletList.push_back ( TripletType (currentPosition, currentPosition, -oneOverH ));
      tripletList.push_back ( TripletType (currentPosition, y * _numX, oneOverH ));
      nonZeros += 2;
    }
    MatX.resize(co::sqr(_numX),co::sqr(_numX));
    MatX.data().squeeze();
    MatX.reserve(Eigen::VectorXi::Constant( co::sqr(_numX), 6));
    MatX.setFromTriplets( tripletList.begin(), tripletList.end() );    
    MatX.makeCompressed(); 
}

  void generateY ( Eigen::SparseMatrix < RealType >& MatY ) {
    
    std::vector< TripletType > tripletList;
    tripletList.reserve(6*_numX*_numX);
    int nonZeros = 0; 
    const RealType oneOverH = 1. / static_cast < RealType > ( _numX - 1 );
    
    for ( int y = 0; y < _numX - 1; ++y ) {
      for ( int x = 0; x < _numX ; ++x ) {
        const int currentPosition = y * _numX + x;
        tripletList.push_back (TripletType ( currentPosition, currentPosition, -oneOverH ));
        tripletList.push_back (TripletType ( currentPosition, currentPosition + _numX, oneOverH ));
        nonZeros += 2;
      }
    }
    MatY.resize(co::sqr(_numX),co::sqr(_numX));
    MatY.data().squeeze();
    MatY.reserve(Eigen::VectorXi::Constant( co::sqr(_numX) , 6));
    MatY.setFromTriplets( tripletList.begin(), tripletList.end() );    
    MatY.makeCompressed();        
  }
};

/* template < typename VectorNd, typename ImageType, typename RealType >
class Converter {
    VectorNd _vector;
    ImageType _image;
    
public:
    
    Converter (){}
    
    void imageToVector ( ImageType& Image){
        _vector.resize( Image.cols() * Image.rows(), 1);
        for ( int i = 0; i < Image.rows(); ++i ){
            _vector.segment( i * Image.cols(), Image.cols() ) = Image.row( i ); 
        }
        Image.resize( _vector.size(), 1);
        Image = _vector;
    }
    
    void vectorToImage ( VectorNd& Vector, const int Rows){
        _image.resize( Rows, Vector.size() / Rows );
        for ( int i = 0; i < Rows; ++i){
            _image.row(i) = Vector.segment( i * _image.cols(), _image.cols() );
        }
        Vector.resize( Rows, _image.cols() );
        Vector = _image;
    }
    
};
*/

template < typename VectorNd, typename ImageType, typename RealType >
void imageToVector ( ImageType& Image){
    VectorNd _vector( Image.cols() * Image.rows(), 1 );
    for ( int i = 0; i < Image.rows(); ++i ){
        _vector.segment( i * Image.cols(), Image.cols() ) = Image.row( i ); 
    }
    Image.resize( _vector.size(), 1);
    Image = _vector;
}


template < typename VectorNd, typename ImageType, typename RealType >
ImageType vectorToImage ( VectorNd& Vector, const int Cols){
    ImageType _image( Cols, Vector.size() / Cols);
    for ( int i = 0; i < Cols; ++i){
        _image.col(i) = Vector.segment( i * _image.rows(), _image.rows() );
    }
    return _image;
}





