#ifndef __MASSMATRIXINTEGRATOR_H
#define __MASSMATRIXINTEGRATOR_H

#include <unityTriangleIntegratorShellFE.h>

using namespace shellFE;



template <typename ConfiguratorType, typename Imp >
class UnitTriangleFELinWeightedMassIntegrator :
      public MatrixValuedIntegratorBase< ConfiguratorType, UnitTriangleFELinWeightedMassIntegrator<ConfiguratorType, Imp> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  const ConfiguratorType &_config;
  
public:
  UnitTriangleFELinWeightedMassIntegrator ( const ConfiguratorType & Config ) :
   MatrixValuedIntegratorBase<  ConfiguratorType, UnitTriangleFELinWeightedMassIntegrator<ConfiguratorType, Imp> > ( Config ),
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

    const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
    const int numQuadPoints = bfs.numQuadPoints( );
    
    for ( int q = 0; q < numQuadPoints; ++q ) {
      for ( int i = 0; i < numDofs; ++i ) {
          RealType value_i = bfs.evaluate( i, q );
        for ( int j = 0; j < numDofs; ++j ) {
          RealType value_j = bfs.evaluate( j, q );
          LocalMatrix(j,i) += value_i * value_j * bfs.getWeight ( q );
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



//! \brief Assembles \f$ \sum_{EL} \int_T \sqrt{\det g} g^{-1} D v_i \cdot D v_j \f$
template<typename ConfiguratorType>
class MassMatrixIntegrator : public
    UnitTriangleFELinWeightedMassIntegrator
    <ConfiguratorType, MassMatrixIntegrator <ConfiguratorType> >
{
protected:
    typedef typename ConfiguratorType::RealType   RealType;
    typedef typename ConfiguratorType::Matrix22   Matrix22;

public:
    //! \brief Configurator
    //! \param[in] conf Reference to the configurator used with the mesh
    MassMatrixIntegrator(const ConfiguratorType &conf)
        : UnitTriangleFELinWeightedMassIntegrator<
          ConfiguratorType, MassMatrixIntegrator<ConfiguratorType>>(conf)
    {}

    //! \brief Computes the weight matrix
    //! \param[in] el The element for which matrix shall be computed
    //! \param[out] matrix The weight matrix
    inline void getCoeffMatrix(const typename ConfiguratorType::ElementType &el,
                               const int &/*quadrature point*/, Matrix22 &matrix) const
    {
        auto &edge1 = -1.0 * el.getEdge(1);
        auto &edge2 = el.getEdge(2);
        Eigen::Matrix<RealType, 3, 2> dx;
        // Flat triangles ~> Last entry 0.0
        dx << edge2[0], edge1[0], edge2[1], edge1[1], edge2[2], edge1[2];

        matrix = el.getAreaOfFlattenedTriangle() * (dx.transpose() * dx).inverse();
    }
};





#endif
