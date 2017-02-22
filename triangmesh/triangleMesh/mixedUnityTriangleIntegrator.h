#ifndef __MIXEDUNITTRIANGLEINTEGRATORSHELLFE_H
#define __MIXEDUNITTRIANGLEINTEGRATORSHELLFE_H

//! General Interface for matrix valued integrators: assembles R^{u,m}-Matrix: m=NumVecCompsArg (ConfigPf), u=NumVecCompsDest (Config)
template < typename ConfiguratorType, typename ConfiguratorTypePf, typename Imp> //int NumVecCompsArg = 1, int NumVecCompsDest = 3>
class BlockMatrixValuedIntegratorMixedConfig {
public:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::MaskType MaskType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::LocalVectorType LocalVectorType;

  explicit BlockMatrixValuedIntegratorMixedConfig ( const ConfiguratorType &conf, const ConfiguratorTypePf& confpf ) : _config ( conf ), _configPf ( confpf ) {}

  
  //! this function has to be provided in the implementation (derived class) of the interface
  inline void getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint, LocalVectorType (&localVector)[3] ) const {
    this->asImp().getNonlinearity ( El, QuadPoint, localVector );
  }
  
protected:
    void assembleTripletList ( std::vector<TripletType> & tripletList, const RealType Factor ) const {
        tripletList.reserve( 3 * _config.getNumLocalDofs() *_config.getInitializer().getNumTriangs () * _configPf.getNumLocalDofs() *_configPf.getInitializer().getNumTriangs() );
       
        const int numGlobalDofs = _config.getNumGlobalDofs();
        //const int numGlobalDofsPf = _configPf.getNumGlobalDofs();
        int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
        int globalDofsPf[ ConfiguratorTypePf::maxNumLocalDofs ];
        for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumTriangs(); ++elementIdx){
            const ElementType& El ( _config.getInitializer().getTriang( elementIdx ) );
            const int numLocalDofs = _config.getNumLocalDofs ( El );
            const int numLocalDofsPf = _configPf.getNumLocalDofs ( El );
            for ( int i = 0; i < numLocalDofs;   ++i ) globalDofs[ i ]   = _config.localToGlobal ( El, i );
            for ( int i = 0; i < numLocalDofsPf; ++i ) globalDofsPf[ i ] = _configPf.localToGlobal ( El, i );
            const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet(El);
            const typename ConfiguratorTypePf::BaseFuncSetType &bfsPf = _configPf.getBaseFunctionSet(El);
            
            for (int quadPoint = 0; quadPoint < _config.maxNumQuadPoints(); ++quadPoint) {
                LocalVectorType localVector[3];
                this->asImp().getNonlinearity ( El, quadPoint, localVector );
                for ( int destComp = 0; destComp < 3; ++destComp ){
                  for ( int i = 0; i < numLocalDofs; ++i ) {
                    int glob_i = globalDofs[ i ];
                    RealType b_i = bfs.evaluate( i, quadPoint );
                    for ( int j = 0; j < numLocalDofsPf; ++j ) {
                      int glob_j = globalDofsPf[ j ];
                      RealType b_j = bfsPf.evaluate( j, quadPoint );
                      tripletList.push_back( TripletType( glob_i + destComp * numGlobalDofs, glob_j, 0.5 * Factor * localVector[destComp](i) * b_j * bfs.getWeight(quadPoint) ) );
                    }
                  }
                }
                
            }//end quadPoint-loop
        }//end element-loop
    }
  
public:

  template <typename BlockMatrixType>
  void assemble ( BlockMatrixType &Dest, const RealType Factor = 1.0 ) const {
    std::vector<TripletType> tripletList;
    assembleTripletList ( tripletList, Factor );
    Dest.setFromTriplets( tripletList.begin(), tripletList.end() ); 
  }
  
  template <typename BlockMatrixType>
  void assembleDirichlet ( BlockMatrixType &Dest, const MaskType& boundaryMask, const RealType Factor = 1.0 ) const {
    
    std::vector<TripletType> tripletList;
    assembleTripletList ( tripletList, Factor );

    // Boundary Mask
    const int numGlobalDofs = _config.getNumGlobalDofs();
    std::vector<TripletType> tripletListMasked;
    tripletListMasked.reserve( 3 * _config.getNumLocalDofs() *_config.getInitializer().getNumTriangs () * _configPf.getNumLocalDofs() *_configPf.getInitializer().getNumTriangs() );

    for( unsigned int iter=0; iter < tripletList.size(); ++iter ){
      if( boundaryMask[tripletList[iter].row() % numGlobalDofs] ){
       //Boundary node!        
      } else {
        tripletListMasked.push_back( tripletList[iter] );
      }
    }
    
    Dest.setFromTriplets( tripletListMasked.begin(), tripletListMasked.end() );
    
  }
  

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
  
  const ConfiguratorType &_config;
  const ConfiguratorTypePf &_configPf;
};



//! General Interface for matrix valued integrators: assembles R^{u,m}-Matrix: m=NumVecCompsArg (ConfigPf), u=NumVecCompsDest (Config)
// template < typename ConfiguratorType, typename ConfiguratorTypePf, typename Imp> //int NumVecCompsArg = 1, int NumVecCompsDest = 3>
// class VectorWeightedBlockMatrixValuedIntegratorMixedConfig {
// public:
//   typedef typename ConfiguratorType::RealType   RealType;
//   typedef typename ConfiguratorType::DomVecType DomVecType;
//   typedef typename ConfiguratorType::TangentVecType TangentVecType;
//   typedef typename ConfiguratorType::Matrix22 Matrix22;
//   typedef typename ConfiguratorType::ElementType ElementType;
//   typedef typename ConfiguratorType::MaskType MaskType;
//   typedef typename ConfiguratorType::TripletType TripletType;
// 
//   explicit VectorWeightedBlockMatrixValuedIntegratorMixedConfig ( const ConfiguratorType &conf, const ConfiguratorTypePf& confpf ) : _config ( conf ), _configPf ( confpf ) {}
// 
//   
//   //! this function has to be provided in the implementation (derived class) of the interface
//   inline void getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint, TangentVecType &NL ) const {
//     this->asImp().getNonlinearity ( El, QuadPoint, NL );
//   }
//   
// protected:
//     void assembleTripletList ( std::vector<TripletType> & tripletList, const RealType Factor ) const {
//         tripletList.reserve( 3 * _config.getNumLocalDofs() *_config.getInitializer().getNumTriangs () * _configPf.getNumLocalDofs() *_configPf.getInitializer().getNumTriangs() );
//        
//         const int numGlobalDofs = _config.getNumGlobalDofs();
//         //const int numGlobalDofsPf = _configPf.getNumGlobalDofs();
//         int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
//         int globalDofsPf[ ConfiguratorTypePf::maxNumLocalDofs ];
//         for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumTriangs(); ++elementIdx){
//             const ElementType& El ( _config.getInitializer().getTriang( elementIdx ) );
//             const int numLocalDofs = _config.getNumLocalDofs ( El );
//             const int numLocalDofsPf = _configPf.getNumLocalDofs ( El );
//             for ( int i = 0; i < numLocalDofs;   ++i ) globalDofs[ i ]   = _config.localToGlobal ( El, i );
//             for ( int i = 0; i < numLocalDofsPf; ++i ) globalDofsPf[ i ] = _configPf.localToGlobal ( El, i );
//             const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet(El);
//             const typename ConfiguratorTypePf::BaseFuncSetType &bfsPf = _configPf.getBaseFunctionSet(El);
//             
//             for (int quadPoint = 0; quadPoint < _config.maxNumQuadPoints(); ++quadPoint) {
//                 TangentVecType nl;
//                 this->asImp().getNonlinearity ( El, quadPoint, nl );
//                 for ( int destComp = 0; destComp < 3; ++destComp ){
//                   for ( int i = 0; i < numLocalDofs; ++i ) {
//                     int glob_i = globalDofs[ i ];
//                     RealType b_i = bfs.evaluate( i, quadPoint );
//                     for ( int j = 0; j < numLocalDofsPf; ++j ) {
//                       int glob_j = globalDofsPf[ j ];
//                       RealType b_j = bfsPf.evaluate( j, quadPoint );
//                       tripletList.push_back( TripletType( glob_i + destComp * numGlobalDofs, glob_j, 0.5 * Factor * nl[destComp] * b_i * b_j * bfs.getWeight(quadPoint) ) ); //
//                     }
//                   }
//                 }
//                 
//             }//end quadPoint-loop
//         }//end element-loop
//     }
//   
// public:
// 
//   template <typename BlockMatrixType>
//   void assemble ( BlockMatrixType &Dest, const RealType Factor = 1.0 ) const {
//     std::vector<TripletType> tripletList;
//     assembleTripletList ( tripletList, Factor );
//     Dest.setFromTriplets( tripletList.begin(), tripletList.end() ); 
//   }
//   
// //   template <typename BlockMatrixType>
// //   void assembleDirichlet ( BlockMatrixType &Dest, const MaskType& boundaryMask, const RealType Factor = 1.0 ) const {
// //     
// //     std::vector<TripletType> tripletList;
// //     assembleTripletList ( tripletList, Factor );
// // 
// //     // Boundary Mask
// //     const int numGlobalDofs = _config.getNumGlobalDofs();
// //     std::vector<TripletType> tripletListMasked;
// //     tripletListMasked.reserve( NumVecCompsArg * NumVecCompsDest * aol::Sqr( _config.getNumLocalDofs() ) *_config.getInitializer ().getNumTriangs ());
// // 
// //     for( unsigned int iter=0; iter < tripletList.size(); ++iter ){
// //       if( (boundaryMask[tripletList[iter].row() % numGlobalDofs] ) || (boundaryMask[tripletList[iter].col() % numGlobalDofs] ) ){
// //        //Boundary node!        
// //       } else {
// //         tripletListMasked.push_back( tripletList[iter] );
// //       }
// //     }
// //     
// //     for ( int i = 0; i < _config.getNumGlobalDofs(); ++i ){
// //       if ( boundaryMask[i] ){
// //         for ( int Comp = 0; Comp < 3; ++Comp )
// //             tripletListMasked.push_back( TripletType( i + Comp * numGlobalDofs, i + Comp * numGlobalDofs, 1.0 ) );
// //       }
// //     }
// //     
// //     Dest.setFromTriplets( tripletListMasked.begin(), tripletListMasked.end() );
// //     
// //   }
// 
// protected:
//   // barton-nackman
//   inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
//   inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
//   
//   const ConfiguratorType &_config;
//   const ConfiguratorTypePf &_configPf;
// };


#endif