#ifndef CHAMBOLLEPOCKALGORITHM_H
#define CHAMBOLLEPOCKALGORITHM_H

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <math.h> 

namespace convexOpt {

enum ALGORITHM_CHOICE {
  CHAMBOLLE_POCK_ALGORITHM1 = 0, CHAMBOLLE_POCK_ALGORITHM2 = 1
};

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
  const int N = UL.size ();
  G.resize ( co::sqr ( N ) );
  for ( int i = 0; i < N; ++i ) {
    for ( int j = 0; j < N; ++j ) {
      if ( j + i < N )
        G ( j + (N - 1 - i) * N ) = co::sqr ( Lambda * (UL ( j ) - UR ( j + i )) );
      else
        G ( j + (N - 1 - i) * N ) = co::sqr ( Lambda * (UL ( j ) - UR ( N - 1 )) );

    }
  }
}

template < typename RealType, typename VectorNd >
class KProjector {
  const VectorNd& _g;
  const RealType _regularizationWeight;

public:
  KProjector ( const VectorNd& G )
: _g ( G ), _regularizationWeight () {
  }

  void projectOntoK ( VectorNd& Arg ) const {
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
  CProjector ( int N )
: _N ( N ) {
  }

  void projectOntoC ( VectorNd& Arg ) const {
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
  Arg ( i ) = co::min ( static_cast < RealType > ( 1 ), co::max ( static_cast < RealType > ( 0 ), Arg ( i ) ) );
  }

};

template < typename RealType, typename VectorType >
RealType squaredL2NormOfDifference ( const VectorType& Arg1, const VectorType& Arg2 ) {
  VectorType diff = Arg1;
  diff -= Arg2;
  return diff.squaredNorm ();
}

template < typename RealType, typename VectorNd, typename OperatorType, typename ProjectorOntoC, typename ProjectorOntoK >
void ChambollePockAlgorithm1 ( VectorNd& PrimalSolution, VectorNd& DualSolution, const OperatorType& K, const ProjectorOntoC& ResolventOfG, const ProjectorOntoK& ResolventOfH, const RealType Tau, const RealType Sigma, const int MaxIter, const RealType StopEpsilon, int& EndIter, RealType& EndEpsilon ) {
  VectorNd xBar ( PrimalSolution );
  VectorNd oldPrimalSolution ( PrimalSolution );
  VectorNd adjointOfDual ( PrimalSolution );
  VectorNd gradientOfXBar ( DualSolution );

  for ( EndIter = 0; EndIter < MaxIter; ++EndIter ) {

    oldPrimalSolution = PrimalSolution;
    K.apply ( xBar, gradientOfXBar );
    DualSolution += Sigma * gradientOfXBar;
    ResolventOfH.projectOntoK ( DualSolution );

    //xBar not needed beyond this point, so an update is fine
    xBar = PrimalSolution;

    K.applyAdjoint ( DualSolution, adjointOfDual );
    PrimalSolution -= Tau * adjointOfDual;
    ResolventOfG.projectOntoC ( PrimalSolution );
    xBar *= -static_cast < RealType > ( 1 );
    xBar += static_cast < RealType > ( 2 ) * PrimalSolution;

    if ( EndIter % 1000 == 0 && EndIter != 0 ) {
      EndEpsilon = squaredL2NormOfDifference < RealType > ( oldPrimalSolution, PrimalSolution );
      std::cout << std::endl << "Iterations: " << EndIter << std::endl << "Epsilon: " << EndEpsilon << std::endl;
      if ( EndEpsilon < StopEpsilon )
        break;
    }
  }
  std::cout << std::endl << "tau  : " << Tau << std::endl << "sigma: " << Sigma << std::endl;
}

template < typename RealType, typename VectorNd, typename OperatorType, typename ProjectorOntoC, typename ProjectorOntoK >
void ChambollePockAlgorithm2 ( VectorNd& PrimalSolution, VectorNd& DualSolution, const OperatorType& K, const ProjectorOntoC& ResolventOfG, const ProjectorOntoK& ResolventOfH, const RealType gamma, RealType& tau, RealType& sigma, const int MaxIter, const RealType StopEpsilon, int& EndIter, RealType& EndEpsilon ) {
  VectorNd xBar ( PrimalSolution );
  VectorNd oldPrimalSolution ( PrimalSolution );
  VectorNd adjointOfDual ( PrimalSolution );
  VectorNd gradientOfXBar ( DualSolution );
  RealType theta = static_cast < RealType > ( 1 );

  for ( EndIter = 0; EndIter < MaxIter; ++EndIter ) {

    oldPrimalSolution = PrimalSolution;
    K.apply ( xBar, gradientOfXBar );
    DualSolution += sigma * gradientOfXBar;
    ResolventOfH.projectOntoK ( DualSolution );

    //xBar not needed beyond this point, so an update is fine
    xBar = PrimalSolution;

    K.applyAdjoint ( DualSolution, adjointOfDual );
    PrimalSolution -= tau * adjointOfDual;
    ResolventOfG.projectOntoC ( PrimalSolution );

    theta = static_cast < RealType > ( 1 ) / std::sqrt ( static_cast < RealType > ( 1 ) + static_cast < RealType > ( 2 ) * gamma * tau );
    tau *= theta;
    sigma /= theta;
    xBar *= -theta;
    xBar += (static_cast < RealType > ( 1 ) + theta) * PrimalSolution;

    if ( EndIter % 1000 == 0 && EndIter != 0 ) {
      EndEpsilon = squaredL2NormOfDifference < RealType > ( oldPrimalSolution, PrimalSolution );
      std::cout << std::endl << "Iterations: " << EndIter << std::endl << "Epsilon: " << EndEpsilon << std::endl;
      if ( EndEpsilon < StopEpsilon )
        break;
    }
  }
  std::cout << std::endl << "theta: " << theta << std::endl << "tau  : " << tau << std::endl << "sigma: " << sigma << std::endl << "gamma: " << gamma << std::endl;
  std::cout << std::endl << "tau*gamma:" << tau * gamma << std::endl;
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
    Dest.head ( _numX * _numX ) += _matX * Arg;
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
    const RealType oneOverH = static_cast < RealType > ( 1 ) / static_cast < RealType > ( _numX - 1 );

    for ( int x = 0; x < _numX - 1; ++x ) {
      for ( int y = 0; y < _numX; ++y ) {
        const int currentPosition = y * _numX + x;
        tripletList.push_back ( TripletType ( currentPosition, currentPosition, -oneOverH ) );
        tripletList.push_back ( TripletType ( currentPosition, currentPosition + 1, oneOverH ) );
      }
    }

    MatX.resize ( co::sqr ( _numX ), co::sqr ( _numX ) );
    MatX.setFromTriplets ( tripletList.cbegin (), tripletList.cend () );
    MatX.makeCompressed ();
  }

  void generateY ( Eigen::SparseMatrix < RealType >& MatY ) {
    std::vector < TripletType > tripletList;
    tripletList.reserve ( 2 * co::sqr ( _numX ) );
    const RealType oneOverH = static_cast < RealType > ( 1 ) / static_cast < RealType > ( _numX - 1 );

    for ( int y = 0; y < _numX - 1; ++y ) {
      for ( int x = 0; x < _numX; ++x ) {
        const int currentPosition = y * _numX + x;
        tripletList.push_back ( TripletType ( currentPosition, currentPosition, -oneOverH ) );
        tripletList.push_back ( TripletType ( currentPosition, currentPosition + _numX, oneOverH ) );
      }
    }
    MatY.resize ( co::sqr ( _numX ), co::sqr ( _numX ) );
    MatY.setFromTriplets ( tripletList.cbegin (), tripletList.cend () );
    MatY.makeCompressed ();
  }
};

template < typename RealType, typename VectorNd, typename OpType >
void applyConvexOptimizationAlgorithm ( VectorNd& PrimalSolution, VectorNd& DualSolution, const VectorNd& Disparity, const int N, const OpType& Fd, const boost::property_tree::ptree& Pt, const ALGORITHM_CHOICE AlgorithmChoice, int& endIter ) {

  convexOpt::KProjector < RealType, VectorNd > resolventData ( Disparity );
  convexOpt::CProjector < RealType, VectorNd > resolventBV ( N );
  RealType endEpsilon = static_cast < RealType > ( 0 );

  const RealType gamma = Pt.get < RealType > ("Parameters_Chambolle.gamma");
  RealType tau = Pt.get < RealType > ("Parameters_Chambolle.tau");
  RealType sigma = Pt.get < RealType > ("Parameters_Chambolle.sigma");
  const int maxIter = Pt.get < int > ( "Parameters_Chambolle.maxIter" );
  const RealType stopEpsilon = Pt.get < RealType > ("Parameters_Chambolle.stopEpsilon");

  switch ( AlgorithmChoice ) {
  case CHAMBOLLE_POCK_ALGORITHM1: {
    std::cout << co::color::green << "Starting Chambolle-Pock Algorithm 1!" << co::color::reset << std::endl;
    ChambollePockAlgorithm1 ( PrimalSolution, DualSolution, Fd, resolventBV, resolventData, tau, sigma, maxIter, stopEpsilon, endIter, endEpsilon );
    std::cout << std::endl << "EndIter: " << endIter << std::endl << "EndEpsilon: " << endEpsilon << std::endl;
  }
  break;
  case CHAMBOLLE_POCK_ALGORITHM2: {
    std::cout << co::color::green << "Starting Chambolle-Pock Algorithm 2!" << co::color::reset << std::endl;
    ChambollePockAlgorithm2 ( PrimalSolution, DualSolution, Fd, resolventBV, resolventData, gamma, tau, sigma, maxIter, stopEpsilon, endIter, endEpsilon );
    std::cout << std::endl << "EndIter: " << endIter << std::endl << "EndEpsilon: " << endEpsilon << std::endl;
  }
  break;

  default:
    throw std::runtime_error ( "Wrong algorithm choice!" );
    break;
  }
}

}

#endif
