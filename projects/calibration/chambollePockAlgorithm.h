#ifndef CHAMBOLLEPOCKALGORITHM_H
#define CHAMBOLLEPOCKALGORITHM_H

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <math.h> 

namespace convexOpt {

template < typename RealType, typename VectorNd >
void divideByMaximumOfOneOrNorm ( VectorNd& Arg ) {
  int numX = Arg.rows () / 2;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for ( int i = 0; i < numX; ++i ) {
    RealType scale = static_cast < RealType > ( 1 );
    VectorNd tmp ( 2 );
    tmp ( 0 ) = Arg ( i );
    tmp ( 1 ) = Arg ( i + numX );
    if ( tmp.norm () > 1 ) {
      scale /= tmp.norm ();
      Arg ( i ) *= scale;
      Arg ( i + numX ) *= scale;
    }
  }
}

template < typename RealType, typename VectorNd >
void createDisparity ( const VectorNd& UL, const VectorNd& UR, VectorNd& G, const RealType Lambda ) {
  const unsigned int N = UL.size ();
  G.resize ( co::sqr ( N ) );
  for ( int i = 0; i < N; ++i ) {
    for ( int j = 0; j < N; ++j ) {
      /*RealType scale = (N-1) / (N - 1 + i);
         int k = (int) ((RealType) (j + i) * scale);
         G ( j + (N - 1 - i) * N ) = abs ( Lambda * ( UL ( j ) - UR ( k ) ));
       */
      if ( j + i < N )
        G ( j + (N - 1 - i) * N ) = abs ( Lambda * (UL ( j ) - UR ( j + i )) );
      else
        G ( j + (N - 1 - i) * N ) = abs ( Lambda * (UL ( j ) - UR ( N - 1 )) );

    }
  }
}

template < typename RealType, typename VectorNd >
class KProjector {
  const VectorNd& _g;
  const RealType _regularizationWeight;

public:
  KProjector ( const VectorNd& G ) : _g ( G ), _regularizationWeight () { }

  void projectOntoK ( VectorNd& Arg, int& Iterations ) {
    Iterations = 0;
    int numX = Arg.size () / 2;
    for ( int i = 0; i < numX; ++i ) {
      if ( Arg ( i ) > 1 ) {
        Arg ( i ) = static_cast < RealType > ( 1 );
      }
      if ( Arg ( i ) < -1 ) {
        Arg ( i ) = static_cast < RealType > ( -1 );
      }
      if ( Arg ( i + numX ) > _g ( i ) ) {
        Arg ( i + numX ) = _g ( i );
        ++Iterations;
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
  ROFDataResolvent ( const VectorNd& Image, const RealType Lambda )
: _image ( Image ), _lambda ( Lambda ), _hscale ( static_cast < RealType > ( 1 ) ) {
  }

  void apply ( VectorNd& Arg, const RealType Tau ) const {
    Arg += Tau * _lambda * _hscale * _image;
    Arg /= (1. + Tau * _lambda * _hscale);
  }
};

template < typename RealType, typename VectorNd >
class CProjector {
  int _N;
public:
  CProjector ( int N ) : _N ( N ) { }

  void projectOntoC ( VectorNd& Arg ) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < _N; ++i ) {
      Arg ( i ) = static_cast < RealType > ( 0 );
      Arg ( i + (_N - 1) * _N ) = static_cast < RealType > ( 1 );
    }
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < co::sqr ( _N ); ++i )
      Arg ( i ) = co::min ( static_cast < RealType > ( 1 ), co::max ( static_cast < RealType > ( 0 ), Arg ( i )  ) );
  }

};

template < typename RealType, typename VectorNd >
void partitionInDimension ( VectorNd Arg, VectorNd& Dest, int I, int Dim ) {
  if ( Arg.size () % Dim != 0 ) {
    std::cout << "Dimension does not fit partition!\n" << std::endl;
  }
  else {
    int j = Arg.size () / Dim;
    Dest = Arg.segment ( I * j, j );
  }

}

template < typename RealType, typename VectorNd, typename OperatorType, typename ProjectorOntoC, typename ProjectorOntoK >
void ChambollePockAlgorithm1 ( VectorNd& PrimalSolution, VectorNd& DualSolution, const int N, const OperatorType& K, ProjectorOntoC& ResolventOfG, ProjectorOntoK& ResolventOfH, const RealType tau, const RealType sigma, const int maxIter, const RealType StopEpsilon, int& enditer, RealType& endepsilon ) {
  VectorNd xBar ( PrimalSolution );
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
      std::cout << std::endl << "Iterations: " << iter << std::endl << "Anzahl Resolvente Projektionen: " << anzahl << std::endl << "Epsilon: " << epsilon.squaredNorm () << std::endl;
    }

    if ( iter % 10 == 0 && iter != 0 ) {
      oldPrimalSolution -= PrimalSolution;
      if ( oldPrimalSolution.squaredNorm () < StopEpsilon ) {   //(l_2^2-Norm)
        endepsilon = oldPrimalSolution.squaredNorm ();
        enditer = iter;
        std::cout << "Number of required iterations: " << iter << std::endl;
        break;
      }
      oldPrimalSolution += PrimalSolution;
    }
  }
  enditer = maxIter;
  oldPrimalSolution -= PrimalSolution;
  endepsilon = oldPrimalSolution.squaredNorm ();
  std::cout << std::endl << "tau  : " << tau << std::endl << "sigma: " << sigma << std::endl << "Anzahl Projektionen: " << anzahl << std::endl;

}

template < typename RealType, typename VectorNd, typename OperatorType, typename ProjectorOntoC, typename ProjectorOntoK >
void ChambollePockAlgorithm2 ( VectorNd& PrimalSolution, VectorNd& DualSolution, const int N, const OperatorType& K, ProjectorOntoC& ResolventOfG, ProjectorOntoK& ResolventOfH, const RealType gamma, RealType& tau, RealType& sigma, const int maxIter, const RealType StopEpsilon, int& enditer, RealType& endepsilon ) {
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
    xBar += (static_cast < RealType > ( 1 ) + theta) * PrimalSolution;                                    //( - theta, PrimalSolution, aol::ZOTrait<RealType>::one + theta );

    if ( iter % 250 == 0 && iter != 0 ) {
      VectorNd epsilon = oldPrimalSolution - PrimalSolution;
      std::cout << std::endl << "Iterations: " << iter << std::endl << "Anzahl Resolvente Projektionen: " << anzahl << std::endl << "Epsilon: " << epsilon.squaredNorm () << std::endl;
    }

    if ( iter % 10 == 0 && iter != 0 ) {
      oldPrimalSolution -= PrimalSolution;
      if ( oldPrimalSolution.squaredNorm () < StopEpsilon ) {   //(l_2^2-Norm)
        endepsilon = oldPrimalSolution.squaredNorm ();
        enditer = iter;
        std::cout << "Number of required iterations: " << iter << std::endl;
        break;
      }
      oldPrimalSolution += PrimalSolution;
    }
  }
  if ( iter == maxIter ) {
    enditer = maxIter;
    oldPrimalSolution -= PrimalSolution;
    endepsilon = oldPrimalSolution.squaredNorm ();
  }
  std::cout << std::endl << "theta: " << theta << std::endl << "tau  : " << tau << std::endl << "sigma: " << sigma << std::endl << "gamma: " << gamma << std::endl;
  std::cout << std::endl << "tau*gamma:" << tau * gamma << std::endl << "Anzahl Projektionen: " << anzahl << std::endl;
}

template < typename RealType, typename VectorNd >
class ForwardFD {
  typedef Eigen::Triplet < RealType > TripletType;

  const int _numX;
  Eigen::SparseMatrix < RealType > _matX;
  Eigen::SparseMatrix < RealType > _matY;

public:

  ForwardFD ( const int NumX )
: _numX ( NumX ) {

    generateX ( _matX );
    generateY ( _matY );
  }

  void applyAdd ( const VectorNd& Arg, VectorNd& Dest ) const {
    Dest.head ( _numX * _numX ) += _matX * Arg;                                                                                             //_matX->applyAdd ( Arg, Dest[0] );
    Dest.tail ( _numX * _numX ) += _matY * Arg;
  }

  void apply ( const VectorNd& Arg, VectorNd& Dest ) const {
    Dest.setZero ();
    this->applyAdd ( Arg, Dest );
  }

  void applyAddAdjoint ( const VectorNd& Arg, VectorNd& Dest ) const {
    Dest += _matX.transpose () * Arg.head ( _numX * _numX );
    Dest += _matY.transpose () * Arg.tail ( _numX * _numX );
  }

  void applyAdjoint ( const VectorNd& Arg, VectorNd& Dest ) const {
    Dest.setZero ();
    this->applyAddAdjoint ( Arg, Dest );
  }

  const RealType getX ( const int i, const int j ) const {
    return _matX.coeff ( i, j );
  }

  const RealType getY ( const int i, const int j ) const {
    return _matY.coeff ( i, j );
  }

private:
  void generateX ( Eigen::SparseMatrix < RealType >& MatX ) {
    std::vector < TripletType > tripletList;
    tripletList.reserve ( 2 * co::sqr ( _numX ) );
    int nonZeros = 0;
    const RealType oneOverH = 1. / static_cast < RealType > ( _numX - 1 );

    for ( int x = 0; x < _numX - 1; ++x ) {
      for ( int y = 0; y < _numX; ++y ) {
        const int currentPosition = y * _numX + x;
        tripletList.push_back ( TripletType ( currentPosition, currentPosition, -oneOverH ) );
        tripletList.push_back ( TripletType ( currentPosition, currentPosition + 1, oneOverH ) );
        nonZeros += 2;
      }
    }

    MatX.resize ( co::sqr ( _numX ), co::sqr ( _numX ) );
    MatX.data ().squeeze ();
    MatX.reserve ( Eigen::VectorXi::Constant ( co::sqr ( _numX ), 6 ) );
    MatX.setFromTriplets ( tripletList.begin (), tripletList.end () );
    MatX.makeCompressed ();
  }

  void generateY ( Eigen::SparseMatrix < RealType >& MatY ) {

    std::vector < TripletType > tripletList;
    tripletList.reserve ( 6 * _numX * _numX );
    int nonZeros = 0;
    const RealType oneOverH = 1. / static_cast < RealType > ( _numX - 1 );

    for ( int y = 0; y < _numX - 1; ++y ) {
      for ( int x = 0; x < _numX; ++x ) {
        const int currentPosition = y * _numX + x;
        tripletList.push_back ( TripletType ( currentPosition, currentPosition, -oneOverH ) );
        tripletList.push_back ( TripletType ( currentPosition, currentPosition + _numX, oneOverH ) );
        nonZeros += 2;
      }
    }
    MatY.resize ( co::sqr ( _numX ), co::sqr ( _numX ) );
    MatY.data ().squeeze ();
    MatY.reserve ( Eigen::VectorXi::Constant ( co::sqr ( _numX ), 6 ) );
    MatY.setFromTriplets ( tripletList.begin (), tripletList.end () );
    MatY.makeCompressed ();
  }
};

template < typename RealType, typename VectorNd, typename OpType >
void applyChambollePockAccelerated ( VectorNd& PrimalSolution, VectorNd& DualSolution, VectorNd& Disparity, const int N, const OpType& Fd, const boost::property_tree::ptree& Pt ) {

  convexOpt::KProjector < RealType, VectorNd > resolventData ( Disparity );
  convexOpt::CProjector < RealType, VectorNd > resolventBV ( N );
  int endIter = 0;
  RealType endEpsilon = 0.;

  const RealType gamma = std::stod ( Pt.get < std::string > ( "Parameters_Chambolle.gamma" ) );
  RealType tau = std::stod ( Pt.get < std::string > ( "Parameters_Chambolle.tau" ) );
  RealType sigma = std::stod ( Pt.get < std::string > ( "Parameters_Chambolle.sigma" ) );
  const int maxIter = std::stoi ( Pt.get < std::string > ( "Parameters_Chambolle.maxIter" ) );
  const RealType stopEpsilon = std::stod ( Pt.get < std::string > ( "Parameters_Chambolle.stopEpsilon" ) );

  ChambollePockAlgorithm2 ( PrimalSolution, DualSolution, N, Fd, resolventBV, resolventData, gamma, tau, sigma, maxIter, stopEpsilon, endIter, endEpsilon );
  std::cout << std::endl << "endIter: " << endIter << std::endl << "endEpsilon: " << endEpsilon << std::endl;

}

}

#endif
