#ifndef __FEOPS_H
#define __FEOPS_H

//! General Interface for matrix valued integrators
template < typename ConfiguratorType, typename Imp >
class MatrixValuedIntegratorBase {
  const ConfiguratorType &_config;

public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::MaskType MaskType;

  explicit MatrixValuedIntegratorBase ( const ConfiguratorType &conf ) : _config ( conf ) { }

protected:
  void assembleTripletList ( std::vector < TripletType > & tripletList, const RealType Factor ) const {
    tripletList.reserve ( _config.getInitializer ().getNumTriangs () * co::sqr ( _config.getNumLocalDofs () ) );
    LocalMatrixType localMatrix;
    int globalDofs[ConfiguratorType::maxNumLocalDofs];
    for ( int elementIdx = 0; elementIdx < _config.getInitializer ().getNumTriangs (); ++elementIdx ) {
      const ElementType& El ( _config.getInitializer ().getTriang ( elementIdx ) );
      this->asImp ().prepareLocalMatrix ( El, localMatrix );
      const int numLocalDofs = _config.getNumLocalDofs ( El );

      for ( int i = 0; i < numLocalDofs; ++i )
        globalDofs[i] = _config.localToGlobal ( El, i );

      for ( int i = 0; i < numLocalDofs; ++i ) {
        int glob_i = globalDofs[i];
        for ( int j = 0; j < numLocalDofs; ++j ) {
          int glob_j = globalDofs[j];
          tripletList.push_back ( TripletType ( glob_i, glob_j, 0.5 * Factor * localMatrix ( i, j ) ) );
        }
      }
    }
  }

public:
  template < typename SparseMatrixType >
  void assemble ( SparseMatrixType &Dest, const RealType Factor = 1.0 ) const {
    std::vector < TripletType > tripletList;
    assembleTripletList ( tripletList, Factor );
    Dest.setFromTriplets ( tripletList.cbegin (), tripletList.cend () );
  }

  template < typename SparseMatrixType >
  void assembleDirichlet ( SparseMatrixType &Dest, const MaskType& boundaryMask, const RealType Factor = 1.0 ) const {
    std::vector < TripletType > tripletList;
    assembleTripletList ( tripletList, Factor );
    std::vector < TripletType > tripletListMasked;
    tripletListMasked.reserve ( _config.getInitializer ().getNumTriangs () * co::sqr ( _config.getNumLocalDofs () ) );
    for ( unsigned iter = 0; iter < tripletList.size (); ++iter ) {
      if ( (boundaryMask[tripletList[iter].row ()]) || (boundaryMask[tripletList[iter].col ()]) ) {
        //Boundary node!
      }
      else {
        tripletListMasked.push_back ( tripletList[iter] );
      }
    }

    for ( int i = 0; i < _config.getNumGlobalDofs (); ++i ) {
      if ( boundaryMask[i] )
        tripletListMasked.push_back ( TripletType ( i, i, 1.0 ) );
    }
    Dest.setFromTriplets ( tripletListMasked.begin (), tripletListMasked.end () );
  }

protected:
  // barton-nackman
  inline Imp& asImp () { return static_cast < Imp& > ( *this ); }
  inline const Imp& asImp () const { return static_cast < const Imp& > ( *this ); }
};

//! The corresponding matrix assembly yields \f$ \left(\int_\Omega  w(x) \phi_i \phi_j dx\right)_{ij} \f$ for FE basis functions \f$ \phi_i,\phi_j \f$.
template < typename ConfiguratorType, typename Imp > class UnitTriangleFELinWeightedMassIntegratorShellFE :
    public MatrixValuedIntegratorBase < ConfiguratorType, UnitTriangleFELinWeightedMassIntegratorShellFE < ConfiguratorType, Imp > > {
    protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  const ConfiguratorType &_config;

    public:
  UnitTriangleFELinWeightedMassIntegratorShellFE ( const ConfiguratorType& Config ) :
    MatrixValuedIntegratorBase < ConfiguratorType, UnitTriangleFELinWeightedMassIntegratorShellFE < ConfiguratorType, Imp > > ( Config ),
    _config ( Config ) {}

  //! this function has to be provided in the implementation (derived class) of the interface
  inline RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    return this->asImp().getNonlinearity ( El, QuadPoint );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, LocalMatrixType &LocalMatrix ) const {
    const int numDofs = _config.getNumLocalDofs ( El );
    for ( int i = 0; i < numDofs; ++i )
      for ( int j = 0; j < numDofs; ++j )
        LocalMatrix(i,j) = 0.;

    RealType nonlinearity;
    const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
    const int numQuadPoints = bfs.numQuadPoints( );
    for ( int q = 0; q < numQuadPoints; ++q ) {
      nonlinearity = getNonlinearity ( El, q );
      for ( int i = 0; i < numDofs; ++i ) {
        RealType b_i = bfs.evaluate( i, q );
        for ( int j = 0; j < numDofs; ++j ) {
          RealType b_j = bfs.evaluate( j, q );
          LocalMatrix(j,i) += nonlinearity * b_i * b_j * bfs.getWeight ( q );
        }
      }
    }
  }

    protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }
};

template < typename ConfiguratorType >
class UnitTriangleFEMassMatrix : public UnitTriangleFELinWeightedMassIntegratorShellFE < ConfiguratorType, UnitTriangleFEMassMatrix < ConfiguratorType > > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;

public:
  UnitTriangleFEMassMatrix ( const ConfiguratorType & Config ) : UnitTriangleFELinWeightedMassIntegratorShellFE < ConfiguratorType, UnitTriangleFEMassMatrix < ConfiguratorType > > ( Config ) { }

  inline RealType getNonlinearity ( const typename ConfiguratorType::ElementType&, int ) const {
    return 1.;
  }
};

//////////////////////////////////////////////


//! The corresponding matrix assembly yields \f$ \left(\int_\Omega  w(x) nabla\phi_i nabla\phi_j dx\right)_{ij} \f$ for FE basis functions \f$ \phi_i,\phi_j \f$.
template < typename ConfiguratorType, typename Imp >
class UnitTriangleFELinScalarWeightedStiffIntegratorShellFE : public MatrixValuedIntegratorBase < ConfiguratorType, UnitTriangleFELinScalarWeightedStiffIntegratorShellFE < ConfiguratorType, Imp > > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  const ConfiguratorType &_config;

public:
  UnitTriangleFELinScalarWeightedStiffIntegratorShellFE ( const ConfiguratorType & Config )
: MatrixValuedIntegratorBase < ConfiguratorType, UnitTriangleFELinScalarWeightedStiffIntegratorShellFE < ConfiguratorType, Imp > > ( Config ), _config ( Config ) {
  }

  //! this function has to be provided in the implementation (derived class) of the interface
  inline RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    return this->asImp ().getNonlinearity ( El, QuadPoint );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, LocalMatrixType &LocalMatrix ) const {

    const int numDofs = _config.getNumLocalDofs ( El );

    for ( int i = 0; i < numDofs; ++i )
      for ( int j = 0; j < numDofs; ++j )
        LocalMatrix ( i, j ) = 0.;

    RealType nonlinearity;

    const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
    const int numQuadPoints = bfs.numQuadPoints ();

    //TODO cache b_i, b_j? use symmetry?
    for ( int q = 0; q < numQuadPoints; ++q ) {
      nonlinearity = getNonlinearity ( El, q );
      for ( int i = 0; i < numDofs; ++i ) {
        DomVecType grad_i = bfs.evaluateGradient ( i, q );
        for ( int j = 0; j < numDofs; ++j ) {
          DomVecType grad_j = bfs.evaluateGradient ( j, q );
          LocalMatrix ( j, i ) += nonlinearity * grad_i.dot ( grad_j ) * bfs.getWeight ( q );
        }
      }
    }
  }

protected:
  inline Imp &asImp () { return static_cast < Imp& > ( *this ); }
  inline const Imp &asImp () const { return static_cast < const Imp& > ( *this ); }
};

template < typename ConfiguratorType >
class UnitTriangleFEStiffMatrix : public UnitTriangleFELinScalarWeightedStiffIntegratorShellFE < ConfiguratorType, UnitTriangleFEStiffMatrix < ConfiguratorType > > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;

public:
  UnitTriangleFEStiffMatrix ( const ConfiguratorType & Config )
: UnitTriangleFELinScalarWeightedStiffIntegratorShellFE < ConfiguratorType, UnitTriangleFEStiffMatrix < ConfiguratorType > > ( Config ) {
  }

  inline RealType getNonlinearity ( const typename ConfiguratorType::ElementType &, int ) const {
    return 1.;
  }
};

//////////////////////////////////////////////

//! General Interface for mixed matrix valued integrators
template < typename ConfiguratorType >
class MixedMatrixValuedIntegratorBase {
  const ConfiguratorType &_config;
  const int _N;

public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::ElementType ElementType;

  explicit MixedMatrixValuedIntegratorBase ( const ConfiguratorType& Config, const int N ) : _config ( Config ), _N ( N ) { }

protected:

  const int getGlobalIndex ( const DomVecType& GlobalRef ) const {
    const int x = co::min ( static_cast < int > ( _N * GlobalRef[0] ), _N - 1 ) % _N;
    const int y = co::min ( static_cast < int > ( _N * GlobalRef[1] ), _N - 1 ) % _N;
    return x + y * _N;
  }

  void assembleTripletList ( std::vector < TripletType > & TripletList, const RealType Factor ) const {
    const int numDofs = _config.getNumLocalDofs ( );
    DomVecType globalRef;

    const int numQuadPoints = _config.getBaseFunctionSet ( _config.getInitializer ().getTriang ( 0 ) ).numQuadPoints( );

    TripletList.reserve ( _config.getInitializer ().getNumTriangs () * numDofs * numQuadPoints );
    int globalDofs[ConfiguratorType::maxNumLocalDofs];
    for ( int elementIdx = 0; elementIdx < _config.getInitializer ().getNumTriangs (); ++elementIdx ) {
      const ElementType& El ( _config.getInitializer ().getTriang ( elementIdx ) );
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
      for ( int i = 0; i < numDofs; ++i ) {
        const int glob_i = _config.localToGlobal ( El, i );
        for ( int q = 0; q < numQuadPoints; ++q ) {
          _config.getGlobalCoords ( El, q, globalRef );
          const int glob_j = getGlobalIndex ( globalRef );
          TripletList.push_back ( TripletType ( glob_i, glob_j, 0.5 * Factor * bfs.evaluate( i, q ) * bfs.getWeight ( q ) ) );
        }
      }
    }
  }

public:
  template < typename SparseMatrixType >
  void assemble ( SparseMatrixType &Dest, const RealType Factor = 1.0 ) const {
    std::vector < TripletType > tripletList;
    assembleTripletList ( tripletList, Factor );
    Dest.setFromTriplets ( tripletList.cbegin (), tripletList.cend () );
  }
};

template < typename ConfiguratorType >
class ProjectFiniteDifferenceOntoFE {
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  const ConfiguratorType& _config;
  const typename ConfiguratorType::InitType& _mesh;
  SparseMatrixType _massMatrix;
  SparseMatrixType _stiffMatrix;
  SparseMatrixType _mixedMatrix;

public:
  ProjectFiniteDifferenceOntoFE ( const ConfiguratorType& Config, const int N ) : _config ( Config ), _mesh ( _config.getMesh ( ) ),
  _massMatrix ( _mesh.getNumVertices ( ), _mesh.getNumVertices ( ) ), _stiffMatrix ( _mesh.getNumVertices ( ), _mesh.getNumVertices ( ) ),
  _mixedMatrix ( _mesh.getNumVertices ( ), co::sqr ( N ) ) {
    UnitTriangleFEMassMatrix < ConfiguratorType > massOp ( _config );
    UnitTriangleFEStiffMatrix < ConfiguratorType > stiffOp ( _config );
    MixedMatrixValuedIntegratorBase < ConfiguratorType > mixedFEOp ( _config, N );
    massOp.assemble ( _massMatrix );
    stiffOp.assemble ( _stiffMatrix );
    mixedFEOp.assemble( _mixedMatrix );
  }

  template < typename VectorType >
  void apply ( const VectorType& FDVector, VectorType& FEVector, const RealType SmoothingFactor ) const {
    VectorType temp ( _mesh.getNumVertices ( ) );
    temp = _mixedMatrix * FDVector;
    SparseMatrixType systemMatrix = _massMatrix;
    if ( SmoothingFactor > 0 ) {
      systemMatrix += SmoothingFactor * _stiffMatrix;
    }
    Eigen::ConjugateGradient < SparseMatrixType, Eigen::Lower| Eigen::Upper > solver;
    solver.compute ( systemMatrix );
    FEVector = solver.solve ( temp );
  }
};


#endif
