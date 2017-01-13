#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>



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
    //Dest += _matXtranspose * Arg.head( _numX * _numX );
    //Dest += _matYtranspose * Arg.tail( _numX * _numX );
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
    tripletList.reserve ( 6*_numX*_numX );
      
    const RealType oneOverH = 1. / static_cast < RealType > ( _numX - 1. );                                                                              //aol::ZOTrait<RealType>::one / static_cast < RealType > ( _numX - 1 );
    for ( int x = 0; x < _numX - 1; ++x ) {
      for ( int y = 0; y < _numX; ++y ) {
        const int currentPosition = y * _numX + x;
        tripletList.push_back( TripletType ( currentPosition, currentPosition, -oneOverH ) ); 
        tripletList.push_back( TripletType ( currentPosition, currentPosition + 1, oneOverH ) );
      }
    }
    
    for ( int y = 0; y < _numX; ++y ) {
      const int currentPosition = y * _numX + (_numX - 1);
      tripletList.push_back ( TripletType (currentPosition, currentPosition, -oneOverH ));
      tripletList.push_back ( TripletType (currentPosition, y * _numX, oneOverH ));
      
    }
    
    MatX.reserve(Eigen::VectorXi::Constant( MatX.cols(), 6));
    for ( int i = 0; i < MatX.nonZeros(); ++i ){ 
        if( tripletList[i].value() != 0 ){
            MatX.insert(tripletList[i].row(),tripletList[i].col()) = tripletList[i].value();
        }
        else{
           // i--;
        }
    }
    MatX.makeCompressed();
}

  void generateY ( Eigen::SparseMatrix < RealType >& MatY ) {
    
    std::vector< TripletType > tripletList;
    tripletList.reserve(6*_numX*_numX);
      
    const RealType oneOverH = 1. / static_cast < RealType > ( _numX - 1 );
    for ( int y = 0; y < _numX - 1; ++y ) {
      for ( int x = 0; x < _numX; ++x ) {
        const int currentPosition = y * _numX + x;
        tripletList.push_back (TripletType ( currentPosition, currentPosition, -oneOverH ));
        tripletList.push_back (TripletType ( currentPosition, currentPosition + _numX, oneOverH ));
      }
    }
    
    for ( int x = 0; x < _numX; ++x ) {
      const int currentPosition = (_numX - 1 ) * _numX + x;
      tripletList.push_back ( TripletType (currentPosition, currentPosition, -oneOverH ));
      tripletList.push_back ( TripletType (currentPosition, x, oneOverH));
    }
    
    MatY.reserve(Eigen::VectorXi::Constant( MatY.cols() , 6));
    for ( int i = 0; i < MatY.nonZeros() ; ++i ){                                                   //co::sqr(_numX)
        if( tripletList[i].value() != 0 ){
            MatY.insert(tripletList[i].row(),tripletList[i].col()) = tripletList[i].value();
        }
        else{
            //i--;
        }
    }
    
    MatY.makeCompressed();        
  }
};











template <typename VectorNd>
void algorithm1(int N, VectorNd uL, VectorNd uR, Eigen::MatrixXd v, double eta){
    
    
    
    
    
    
    
    
    

}




void algorithm2(){
    
    
}