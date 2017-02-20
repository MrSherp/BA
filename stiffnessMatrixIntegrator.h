#ifndef __STIFFNESSMATRIXINTEGRATOR_H
#define __STIFFNESSMATRIXINTEGRATOR_H

#include <unityTriangleIntegratorShellFE.h>

using namespace shellFE;

//! \brief Assembles \f$ \sum_{EL} \int_T \sqrt{\det g} g^{-1} D v_i \cdot D v_j \f$
template<typename ConfiguratorType>
class StiffnessMatrixIntegrator : public
    UnitTriangleFELinWeightedStiffIntegrator
    <ConfiguratorType, StiffnessMatrixIntegrator <ConfiguratorType> >
{
protected:
    typedef typename ConfiguratorType::RealType   RealType;
    typedef typename ConfiguratorType::Matrix22   Matrix22;

public:
    //! \brief Configurator
    //! \param[in] conf Reference to the configurator used with the mesh
    StiffnessMatrixIntegrator(const ConfiguratorType &conf)
        : UnitTriangleFELinWeightedStiffIntegrator<
          ConfiguratorType, StiffnessMatrixIntegrator<ConfiguratorType>>(conf)
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
