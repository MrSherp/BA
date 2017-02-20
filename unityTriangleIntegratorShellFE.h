/*! \file unityTriangleIntegratorShellFE.h
 *  \brief Base classes for integration type operators defined using the unit triangle
 */

#ifndef __UNITTRIANGLEINTEGRATORSHELLFE_H
#define __UNITTRIANGLEINTEGRATORSHELLFE_H

namespace shellFE {

//! Integrator to compute \f$\int_\Omega f(\phi,x) dx\f$, where \f$\phi\f$ is the argument of the operator.
template <typename ConfiguratorType, typename Imp>
class UnitTriangleFENonlinIntegrationScalarIntegratorShellFE{
public:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::ElementType ElementType;
protected:
    const ConfiguratorType &_config;
public:

    UnitTriangleFENonlinIntegrationScalarIntegratorShellFE ( const ConfiguratorType & Config ) : _config( Config ) {}

    virtual ~UnitTriangleFENonlinIntegrationScalarIntegratorShellFE( ) {}

    void assembleAdd ( RealType &Dest ) const {

        RealType res = 0.;

        for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumTriangs(); ++elementIdx){
            const ElementType& El ( _config.getInitializer().getTriang( elementIdx ) );
            const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
            const int numQuadPoints = bfs.numQuadPoints( );

            RealType a = 0.;
            for ( int q = 0; q < numQuadPoints; ++q ) {
                a += this->asImp().evaluateIntegrand ( El, q ) * bfs.getWeight ( q );
            }
            res += a;
        }
        Dest += 0.5 * res;
    }

    //! interface function, has to be provided in derived classes.
    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
        throw std::invalid_argument( aol::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
        return this->asImp().evaluateIntegrand ( El, QuadPoint );
    }

protected:
    // barton-nackman
    inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
    inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};

//=============================================================================
// Vector-valued Intefaces
//=============================================================================

//! Integrator for \f$ (\int_\Omega s(x)  w_i(x) da )_{i} \f$, of some scalar valued function \f$ s\f$.
template <typename ConfiguratorType, typename Imp>
class UnitTriangleFENonlinOpIntegratorShellFE  {

protected:
	typedef typename ConfiguratorType::RealType RealType;
	typedef typename ConfiguratorType::VectorType VectorType;
	typedef typename ConfiguratorType::ElementType ElementType;
	const ConfiguratorType & _config;

public:
	UnitTriangleFENonlinOpIntegratorShellFE( const ConfiguratorType & Conf ) : _config( Conf){ }

	virtual ~UnitTriangleFENonlinOpIntegratorShellFE( ) {}

	void assembleAdd ( VectorType &Dest ) const {
		RealType *nl_cache = new RealType[ _config.maxNumQuadPoints() ];

		for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumTriangs(); ++elementIdx){
			const ElementType& El ( _config.getInitializer().getTriang( elementIdx ) );
			const int numLocalDofs = _config.getNumLocalDofs ( El );

			const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
			const int numQuadPoints = bfs.numQuadPoints( );

			for ( int q = 0; q < numQuadPoints; ++q )
				nl_cache[q] = this->asImp().getNonlinearity ( El, q );

			RealType aux;

			for ( int dof = 0; dof < numLocalDofs; ++dof ) {
				aux = 0.;

				for ( int q = 0; q < numQuadPoints; ++q )
					aux += nl_cache[q] * bfs.evaluate ( dof, q ) * bfs.getWeight ( q );

				Dest[ _config.localToGlobal ( El, dof ) ] += 0.5 * aux;
			}
		}
		delete[] nl_cache;
	}

	//! interface function, has to be provided in derived classes.
	RealType getNonlinearity ( const typename ConfiguratorType::ElementType & El, int QuadPoint ) const {
		throw std::invalid_argument( aol::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
		return this->asImp().getNonlinearity ( El, QuadPoint );
	}

protected:
	// barton-nackman
	inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
	inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};

//=============================================================================
// Matrix-valued Intefaces
//=============================================================================


//! \brief General interface for matrix valued integrators
template < typename ConfiguratorType, typename Imp >
class MatrixValuedIntegratorBase {
public:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::MaskType MaskType;

  explicit MatrixValuedIntegratorBase(const ConfiguratorType &conf, const RealType factor = 1.0)
  : _config ( conf ), _factor(factor), _tripletListAssembled(false)
  {}

protected:
  void assembleTripletList ( std::vector<TripletType> & tripletList, const RealType Factor ) const {
    tripletList.reserve(_config.getInitializer().getNumTriangs() * aol::Sqr( _config.getNumLocalDofs() ) );
    LocalMatrixType localMatrix;
    int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
    for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumTriangs(); ++elementIdx){
      const ElementType& El ( _config.getInitializer().getTriang( elementIdx ) );
      this->asImp().prepareLocalMatrix ( El, localMatrix );
      const int numLocalDofs = _config.getNumLocalDofs ( El );

      for ( int i = 0; i < numLocalDofs; ++i )
        globalDofs[ i ] = _config.localToGlobal ( El, i );

      for ( int i = 0; i < numLocalDofs; ++i ) {
        int glob_i = globalDofs[ i ];
        for ( int j = 0; j < numLocalDofs; ++j ) {
          int glob_j = globalDofs[ j ];
          tripletList.push_back( TripletType( glob_i, glob_j, Factor * localMatrix(i,j) ) );
        }
      }
    }
  }
  
  void assembleTripletListForN ( std::vector<TripletType> & tripletList, const RealType Factor ) const {
    tripletList.reserve(_config.getInitializer().getNumTriangs() * aol::Sqr( _config.getNumLocalDofs() ) );
    LocalMatrixType localMatrix;
    int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
    for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumTriangs(); ++elementIdx){
      const ElementType& El ( _config.getInitializer().getTriang( elementIdx ) );
      this->asImp().prepareLocalMatrix ( El, localMatrix );
      const int numLocalDofs = _config.getNumLocalDofs ( El );

      for ( int i = 0; i < numLocalDofs; ++i )
        globalDofs[ i ] = _config.localToGlobal ( El, i );

      for ( int i = 0; i < numLocalDofs; ++i ) {
        int glob_i = globalDofs[ i ];
        int glob_j = elementIdx;
        tripletList.push_back( TripletType( glob_i, glob_j, Factor * localMatrix(0,i) ) );        
      }
    }
  }
  
  
public:
  void assembleTripletListCached() const {
      if(!_tripletListAssembled) {
          assembleTripletList(_tripletList, _factor);
          _tripletListAssembled = true;
      }
  }

  void assembleTripletListCachedForN() const {
      if(!_tripletListAssembled) {
          assembleTripletListForN(_tripletList, _factor);
          _tripletListAssembled = true;
      }
  }
  
  template <typename SparseMatrixType>
  void assemble ( SparseMatrixType &Dest ) const {
      assembleTripletListCached();
      Dest.setFromTriplets( _tripletList.cbegin(), _tripletList.cend() );
  }
  
  template <typename SparseMatrixType>
  void assembleDirichlet ( SparseMatrixType &Dest, const MaskType& boundaryMask ) const {
    assembleTripletListCached();

    std::vector<TripletType> tripletListMasked;
    tripletListMasked.reserve(_config.getInitializer().getNumTriangs() * aol::Sqr( _config.getNumLocalDofs() ) );

    for( unsigned iter=0; iter < _tripletList.size(); ++iter ){
      if( (boundaryMask[_tripletList[iter].row()]) || (boundaryMask[_tripletList[iter].col()]) ){
       //Boundary node!        
      } else {
        tripletListMasked.push_back( _tripletList[iter] );
      }
    }
    
    for ( int i = 0; i < _config.getNumGlobalDofs(); ++i ){
      if ( boundaryMask[i] )
         tripletListMasked.push_back( TripletType( i, i, 1.0 ) );
    }
    
    Dest.setFromTriplets( tripletListMasked.begin(), tripletListMasked.end() );
  }
  
  template <typename SparseMatrixType>
  void assembleDirichletForN ( SparseMatrixType &Dest, const MaskType& boundaryMask ) const {
    assembleTripletListCachedForN();

    std::vector<TripletType> tripletListMasked;
    tripletListMasked.reserve(_config.getInitializer().getNumTriangs() * aol::Sqr( _config.getNumLocalDofs() ) );

    for( unsigned iter=0; iter < _tripletList.size(); ++iter ){
      if( (boundaryMask[_tripletList[iter].row()]) || (boundaryMask[_tripletList[iter].col()]) ){
       //Boundary node!        
      } else {
        tripletListMasked.push_back( _tripletList[iter] );
      }
    }
    
    for ( int i = 0; i < _config.getNumGlobalDofs(); ++i ){
      if ( boundaryMask[i] )
         tripletListMasked.push_back( TripletType( i, i, 1.0 ) );
    }
    
    Dest.setFromTriplets( tripletListMasked.begin(), tripletListMasked.end() );
  }
  
protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
  
  const ConfiguratorType &_config;
  mutable std::vector<TripletType> _tripletList;
  RealType _factor;
  mutable bool _tripletListAssembled;
};




//! \brief Provides an easy interface to Finite Element operators of the form \f$ \mbox{div}(A(x)\nabla u)\f$,
//!         where \f$A\f$ is an ASYMMETRIC coefficient matrix.
//! The corresponding matrix assembly yields \f$ \left(\int_\Omega \nabla\phi_i\cdot A(x)\nabla\phi_j dx\right)_{ij} \f$
//! for FE basis functions \f$ \phi_i,\phi_j \f$.
template <typename ConfiguratorType, typename Imp >
class UnitTriangleFELinWeightedStiffIntegrator :
      public MatrixValuedIntegratorBase< ConfiguratorType, UnitTriangleFELinWeightedStiffIntegrator<ConfiguratorType, Imp> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  const ConfiguratorType &_config;
  
public:
  UnitTriangleFELinWeightedStiffIntegrator ( const ConfiguratorType & Config ) :
   MatrixValuedIntegratorBase<  ConfiguratorType, UnitTriangleFELinWeightedStiffIntegrator<ConfiguratorType, Imp> > ( Config ),
  _config ( Config ) {}


  //! \brief This function has to be provided in the implementation (derived class) of the interface
  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint,
                               Matrix22 &Matrix ) const {
    this->asImp().getCoeffMatrix ( El, QuadPoint, Matrix );
  }

  //! \brief Computes the numerical quadrature of the bilinear form and saves the values locally
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, 
                                  LocalMatrixType &LocalMatrix ) const {
    
    
    const int numDofs = _config.getNumLocalDofs ( El );	
    
    for ( int i = 0; i < numDofs; ++i )
      for ( int j = 0; j < numDofs; ++j )
        LocalMatrix(i,j) = 0.;

    Matrix22 mat;
    DomVecType matgrad1;

    const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    for ( int q = 0; q < numQuadPoints; ++q ) {
      getCoeffMatrix ( El, q, mat );
      for ( int i = 0; i < numDofs; ++i ) {
          DomVecType grad_i = bfs.evaluateGradient( i, q );
          matgrad1 = mat * grad_i;
        for ( int j = 0; j < numDofs; ++j ) {
          DomVecType grad_j = bfs.evaluateGradient( j, q );
          LocalMatrix(j,i) += matgrad1.dot( grad_j ) * bfs.getWeight ( q );
        }
      }
    }
    for ( int i = 0; i < numDofs; ++i )
      for ( int j = 0; j < numDofs; ++j ) 
        LocalMatrix(i,j) *= 0.5;  // 0.5 is the volume of the unit triangle
  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

};







} // end namespace

#endif
