#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <math.h> 

/*template < typename RealType, typename VectorNd >
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
*/

template < typename RealType, typename VectorNd >
void divideByMaximumOfOneOrNorm( VectorNd& Arg ){
    int numX = Arg.rows() / 2;
    for ( int i = 0; i < numX; ++i){
        RealType scale = static_cast < RealType > ( 1 );   
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
void createG ( VectorNd UL, VectorNd UR, VectorNd& Dest, RealType Lambda ){
    const unsigned int N = UL.size();
    Dest.resize( co::sqr( N ) );
    for (int i = 0; i < N; ++i){
        for ( int j = 0; j < N; ++j){
            /*RealType scale = (N-1) / (N - 1 + i);
            int k = (int) ((RealType) (j + i) * scale);
            Dest ( j + (N - 1 - i) * N ) = abs ( Lambda * ( UL ( j ) - UR ( k ) ));
            */
            if ( j + i < N)    
                Dest ( j + (N - 1 - i) * N ) = abs (Lambda * ( UL ( j ) - UR ( j + i ) ));
            else
                Dest ( j + (N - 1 - i) * N ) = abs (Lambda * ( UL ( j ) - UR ( N-1 ) ));
            
        }
    }   
}



template < typename RealType, typename VectorNd >
void startVectors( VectorNd& V, VectorNd& Phi, VectorNd G, const int N ){
    const unsigned int Nsquare = V.size();
    
    
    
   /*

    VectorNd tmp (Nsquare);
    tmp.setZero();
    Phi << tmp, -G;
    V.setZero();
    */
    
    V.setZero();
    Phi.setZero();
   /* for (int i = 0; i < N; ++i){
        for ( int j = N/2; j < N; ++j){
            V( i + j * N ) = 1.;
        }
    }   
    */
}



template < typename RealType, typename VectorNd >
class KProjector {
    VectorNd _g;
    const RealType _regularizationWeight;
  
public:
    KProjector ( VectorNd G, RealType Regularizer )
:   _g ( G ), _regularizationWeight ( Regularizer ){ }

void projectOntoK( VectorNd& Arg, int& anzahl ){
    int numX = Arg.size() / 2;
        for ( int i = 0; i < numX; ++i){ 
            if( Arg( i ) > 1 ){
                Arg( i ) = static_cast < RealType > ( 1 );
            }   
            if( Arg( i ) < -1 ){
                Arg( i ) = static_cast < RealType > ( -1 );
            }   
            if( Arg( i + numX ) < -_g( i ) ){
                Arg( i + numX ) = -_g( i );
                anzahl++;
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

/*

template < typename RealType, typename VectorNd >
class CProjector1 {
    int _N;
public:
    CProjector1 ( int M )
:   _N ( M ) {}

void projectOntoC( VectorNd& Arg ){
    RealType min = Arg.minCoeff();
    if( min < 0 ){
        for ( int i = 0; i < _N; ++i){
            for ( int j = 0; j < _N; ++j){
                Arg ( i + j * _N ) -= min;
            }
        }        
    }
    RealType max = Arg.maxCoeff(); 
    RealType scale = 1. / max;  
    if( max > 1 ){
        Arg *= scale;
    }
    for( int i = 0; i < _N; ++i){
        Arg ( i ) = 0;
        Arg ( i + (_N-1) * _N ) = 1;
    }
}

};

*/

template < typename RealType, typename VectorNd >
class CProjector2 {
    int _N;
public:
    CProjector2 ( int M )
:   _N ( M ) {}

void projectOntoC( VectorNd& Arg ){
    for ( int i = 0; i < _N; ++i){
        for ( int j = 0; j < _N; ++j){
            if(Arg ( i + j * _N ) < 0)
                Arg ( i + j * _N ) = static_cast < RealType > ( 0 );
            if(Arg ( i + j * _N ) > 1)
                Arg ( i + j * _N ) = static_cast < RealType > ( 1 );          
        }
    }        
    for( int i = 0; i < _N; ++i){
        Arg ( i ) = static_cast < RealType > ( 0 );
        Arg ( i + (_N-1) * _N ) = static_cast < RealType > ( 1 );
    }
}

};



template < typename RealType, typename VectorNd >
void partitionInDimension(VectorNd Arg, VectorNd& Dest, int I, int Dim){
    if ( Arg.size() % Dim != 0 ){
        std::cout << "Dimension does not fit partition!\n" << std::endl;
    }
    else{
        int j = Arg.size() / Dim;
        Dest = Arg.segment( I  * j, j );
    }

}
    


template < typename RealType, typename VectorNd, typename OperatorType, typename ProjectorOntoC, typename ProjectorOntoK >
void ChambollePockAlgorithm1 ( VectorNd& PrimalSolution,
    VectorNd& DualSolution, const int N, const OperatorType& K, ProjectorOntoC& ResolventOfG, ProjectorOntoK& ResolventOfH,  const RealType tau, const RealType sigma, const int maxIter, const RealType StopEpsilon, int& enditer, RealType& endepsilon ) {
  VectorNd xBar ( PrimalSolution ) ;                                               
  VectorNd oldPrimalSolution ( PrimalSolution );                                  
  VectorNd adjointOfDual ( PrimalSolution );
  VectorNd gradientOfXBar ( DualSolution );

  int anzahl = 0;
  for ( int iter = 0; iter < maxIter; ++iter ) { 
    
    oldPrimalSolution = PrimalSolution;
    K.apply ( xBar, gradientOfXBar );
    DualSolution += sigma * gradientOfXBar;
    ResolventOfH.projectOntoK ( DualSolution, anzahl );
    
    //xBar not needed beyond this point, so an update is fine
    xBar = PrimalSolution;

    K.applyAdjoint ( DualSolution, adjointOfDual );
    PrimalSolution -= tau * adjointOfDual; 
    ResolventOfG.projectOntoC ( PrimalSolution );
    xBar *= -static_cast < RealType > ( 1 );                                                               
    xBar += static_cast < RealType > ( 2 ) * PrimalSolution; 
    
    if ( iter % 100 == 0 && iter != 0 ) {
        VectorNd epsilon = oldPrimalSolution - PrimalSolution;
        cout << endl << "Iterations: " << iter << endl << "Anzahl Resolvente Projektionen: " << anzahl << endl << "Epsilon: " << epsilon.squaredNorm() << std::endl;
    }
    
    if ( iter % 10 == 0 && iter != 0 ) {
      oldPrimalSolution -= PrimalSolution;
     if (oldPrimalSolution.squaredNorm ()  < StopEpsilon ) { //(l_2^2-Norm)
	endepsilon = oldPrimalSolution.squaredNorm ();
        enditer = iter;
        cout << "Number of required iterations: " << iter << endl;
        break;
      }
      oldPrimalSolution += PrimalSolution;
    }
  }
    enditer = maxIter;
    oldPrimalSolution -= PrimalSolution;
    endepsilon = oldPrimalSolution.squaredNorm ();
    cout << endl << "tau  : " << tau << endl << "sigma: " << sigma << endl << "Anzahl Projektionen: " << anzahl << endl;

}



template < typename RealType, typename VectorNd, typename OperatorType, typename ProjectorOntoC, typename ProjectorOntoK >
void ChambollePockAlgorithm2 ( VectorNd& PrimalSolution,
    VectorNd& DualSolution, const int N, const OperatorType& K, ProjectorOntoC& ResolventOfG, ProjectorOntoK& ResolventOfH, const RealType gamma, RealType& tau, RealType& sigma, const int maxIter, const RealType StopEpsilon, int& enditer, RealType& endepsilon ) {
  VectorNd xBar ( PrimalSolution );
  VectorNd oldPrimalSolution ( PrimalSolution );
  VectorNd adjointOfDual ( PrimalSolution );
  VectorNd gradientOfXBar ( DualSolution );
  RealType theta = static_cast < RealType > ( 1 );                                                          

  int anzahl = 0;
  int iter;
  for ( iter = 0; iter < maxIter; ++iter ) {
      
    oldPrimalSolution = PrimalSolution;
    K.apply ( xBar, gradientOfXBar );
    DualSolution += sigma * gradientOfXBar;
    ResolventOfH.projectOntoK ( DualSolution, anzahl );

    //xBar not needed beyond this point, so an update is fine
    xBar = PrimalSolution;

    K.applyAdjoint ( DualSolution, adjointOfDual );
    PrimalSolution -= tau * adjointOfDual;
    ResolventOfG.projectOntoC ( PrimalSolution );
    
    theta = static_cast < RealType > ( 1 ) / std::sqrt ( static_cast < RealType > ( 1 ) + static_cast < RealType > ( 2 ) * gamma * tau );
    tau *= theta;
    sigma /= theta;
    xBar *= -theta;
    xBar += ( static_cast < RealType > ( 1 ) + theta ) * PrimalSolution;                                    //( - theta, PrimalSolution, aol::ZOTrait<RealType>::one + theta );
    
    if ( iter % 250 == 0 && iter != 0 ) {
        VectorNd epsilon = oldPrimalSolution - PrimalSolution;
        cout << endl << "Iterations: " << iter << endl << "Anzahl Resolvente Projektionen: " << anzahl << endl << "Epsilon: " << epsilon.squaredNorm() << std::endl;
    }
    
    if ( iter % 10 == 0 && iter != 0 ) {
      oldPrimalSolution -= PrimalSolution;
     if ( oldPrimalSolution.squaredNorm () < StopEpsilon ) { //(l_2^2-Norm)
	endepsilon = oldPrimalSolution.squaredNorm ();
	enditer = iter;
        cout << "Number of required iterations: " << iter << endl;
        break;
      }
      oldPrimalSolution += PrimalSolution;
    }
  }
    if ( iter == maxIter )
    {
        enditer = maxIter;
        oldPrimalSolution -= PrimalSolution;
        endepsilon = oldPrimalSolution.squaredNorm ();
    }
    cout << endl << "theta: " << theta << endl << "tau  : " << tau << endl << "sigma: " << sigma << endl << "gamma: " << gamma << endl;
    cout << endl << "tau*gamma:" << tau*gamma << endl << "Anzahl Projektionen: " << anzahl << endl;
    
   // cout << endl << "PrimalSolution " << endl << PrimalSolution << endl; //<< "DualSolution: " << endl << DualSolution << endl;
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
      for ( int y = 0; y < _numX ; ++y ) {
        const int currentPosition = y * _numX + x;
        tripletList.push_back( TripletType ( currentPosition, currentPosition, -oneOverH ) ); 
        tripletList.push_back( TripletType ( currentPosition, currentPosition + 1, oneOverH ) );
        nonZeros += 2;
      }
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





