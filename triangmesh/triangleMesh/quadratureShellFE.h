#ifndef __QUADRATURESHELLFE_H
#define __QUADRATURESHELLFE_H

namespace shellFE {


//! Different quadrature types for triangular meshes
    
    
// CenterQuadrature
template <typename RealType, typename DomVecType>
class CenterQuadrature {
public:
  enum { numQuadPoints = 1 };

  inline const DomVecType& getRefCoord ( int /*quadpoint*/ ) const  {
    static const DomVecType c ( 1. / 3., 1. / 3. );
    return c;
  }

  inline RealType getWeight ( int /*quadpoint*/ ) const { return 1.;}
};

// // TriQuadrature
template <typename RealType, typename DomVecType>
class TriQuadrature {
public:
  enum { numQuadPoints = 3 };

  inline const DomVecType& getRefCoord ( int quadpoint ) const  {
    static const DomVecType c[3] = { DomVecType ( 1. / 6., 1. / 6. ),
                                     DomVecType ( 1. / 6., 2. / 3. ),
                                     DomVecType ( 2. / 3., 1. / 6. )
                                    };
    return c[quadpoint];
  }

  inline RealType getWeight ( int /*quadpoint*/ ) const { return 1. / 3.;}
};


// EdgeQuadrature
template <typename RealType, typename DomVecType>
class EdgeQuadrature {
public:
  enum { numQuadPoints = 3 };

  inline const DomVecType& getRefCoord ( int quadpoint ) const  {
    static const DomVecType c[3] = {    DomVecType ( 1./2., 1./2. ),
                                        DomVecType ( 0., 1./2. ),
                                        DomVecType ( 1./2., 0. )
                                   };
    return c[quadpoint];
  }

  inline RealType getWeight ( int /*quadpoint*/ ) const { return 1./3.;}
};

// // EdgeQuadrature
// template <typename RealType>
// class CornerQuadrature {
// public:
//   enum { numQuadPoints = 3 };
// 
//   inline const DomVecType& getRefCoord ( int quadpoint ) const  {
//     static const DomVecType c[3] = {
//                                               DomVecType ( 0., 0. ),
//                                               DomVecType ( 0., 1. ),
//                                               DomVecType ( 1., 0. )
//                                             };
//     return c[quadpoint];
//   }
// 
//   inline RealType getWeight ( int ) const {
//     return 1./3.;
//   }
// };
// 

// // QuadriQuadrature, exact for polynomials of degree 3
// template <typename RealType>
// class QuadriQuadrature {
// public:
//   enum { numQuadPoints = 4 };
// 
//   inline const DomVecType& getRefCoord ( int quadpoint ) const  {
//     static const DomVecType c[4] = {
//                                               DomVecType ( 1. / 3., 1. / 3. ),
//                                               DomVecType ( 3. / 5., 1. / 5. ),
//                                               DomVecType ( 1. / 5., 1. / 5. ),
//                                               DomVecType ( 1. / 5., 3. / 5. )
//                                             };
//     return c[quadpoint];
//   }
// 
//   inline RealType getWeight ( int quadpoint) const {
//     static const RealType w[4] = { 
//                     RealType (-27./48. ),
//                     RealType ( 25./48. ),
//                     RealType ( 25./48. ),
//                     RealType ( 25./48. )
//                   };
//     return w[quadpoint];
//   }
// };
// 
// // SexaQuadrature, exact for polynomials of degree 4
// template <typename RealType>
// class SexaQuadrature {
// public:
//   enum { numQuadPoints = 6 };
// 
//   inline const DomVecType& getRefCoord ( int quadpoint ) const  {
//     static const DomVecType c[6] = {
//                                               DomVecType ( 0.44594849091597 , 0.44594849091597 ),
//                                               DomVecType ( 0.44594849091597 , 0.10810301816807 ),
//                                               DomVecType ( 0.10810301816807 , 0.44594849091597 ),
//                                               DomVecType ( 0.09157621350977 , 0.09157621350977 ),
//                                               DomVecType ( 0.09157621350977 , 0.81684757298046 ),
//                                               DomVecType ( 0.81684757298046 , 0.09157621350977 )
//                                             };
//     return c[quadpoint];
//   }
// 
//   inline RealType getWeight ( int quadpoint) const {
//     static const RealType w[6] = { 
//                     RealType ( 0.22338158967801 ),
//                     RealType ( 0.22338158967801 ),
//                     RealType ( 0.22338158967801 ),
//                     RealType ( 0.10995174365532 ),
//                     RealType ( 0.10995174365532 ),
//                     RealType ( 0.10995174365532 )
//                   };
//     return w[quadpoint];
//   }
// };
// 

// // SeptiQuadrature, exact for polynomials of degree 5
// template <typename RealType>
// class SeptiQuadrature {
// public:
//   enum { numQuadPoints = 7 };
// 
//   inline const DomVecType& getRefCoord ( int quadpoint ) const  {
//     static const DomVecType c[7] = {
//                                               DomVecType ( 0.33333333333333 , 0.33333333333333 ),
//                                               DomVecType ( 0.47014206410511 , 0.47014206410511 ),
//                                               DomVecType ( 0.47014206410511 , 0.05971587178977 ),
//                                               DomVecType ( 0.05971587178977 , 0.4701420641051 ),
//                                               DomVecType ( 0.10128650732346 , 0.10128650732346 ),
//                                               DomVecType ( 0.10128650732346 , 0.79742698535309 ),
//                                               DomVecType ( 0.79742698535309 , 0.10128650732346 )
//                                             };
//     return c[quadpoint];
//   }
// 
//   inline RealType getWeight ( int quadpoint) const {
//     static const RealType w[7] = { 
//                     RealType ( 0.22500000000000 ),
//                     RealType ( 0.13239415278851 ),
//                     RealType ( 0.13239415278851 ),
//                     RealType ( 0.13239415278851 ),
//                     RealType ( 0.12593918054483 ),
//                     RealType ( 0.12593918054483 ),
//                     RealType ( 0.12593918054483 )
//                   };
//     return w[quadpoint];
//   }
// };

  
// DuodecQuadrature, exact for polynomials of degree 6
// template <typename RealType, typename DomVecType>
// class DuodecQuadrature {
// public :
// //   static constexpr unsigned short int _numQuadPoints = 12;
//   static const unsigned short int _numQuadPoints = 12;
// private :
//   
//   static DomVecType _c[_numQuadPoints] = { DomVecType ( 0.24928674517091 , 0.24928674517091 ),
//                                          DomVecType ( 0.24928674517091 , 0.50142650965818 ),
//                                          DomVecType ( 0.50142650965818 , 0.24928674517091 ),
//                                          DomVecType ( 0.06308901449150 , 0.06308901449150 ),
//                                          DomVecType ( 0.06308901449150 , 0.87382197101700 ),
//                                          DomVecType ( 0.87382197101700 , 0.06308901449150 ),
//                                          DomVecType ( 0.31035245103378 , 0.63650249912140 ),
//                                          DomVecType ( 0.63650249912140 , 0.05314504984482 ),
//                                          DomVecType ( 0.05314504984482 , 0.31035245103378 ),
//                                          DomVecType ( 0.63650249912140 , 0.31035245103378 ),
//                                          DomVecType ( 0.31035245103378 , 0.05314504984482 ),
//                                          DomVecType ( 0.05314504984482 , 0.63650249912140 )
//                                        };
//   static RealType _w[_numQuadPoints] = { 
//                      0.11678627572638 ,
//                      0.11678627572638 ,
//                      0.11678627572638 ,
//                      0.05084490637021 ,
//                      0.05084490637021 ,
//                      0.05084490637021 ,
//                      0.08285107561837 ,
//                      0.08285107561837 ,
//                      0.08285107561837 ,
//                      0.08285107561837 ,
//                      0.08285107561837 ,
//                      0.08285107561837 
//                   };
//   
//   
// public:
//   
//   const unsigned short int getNumQuadpoints ( ) const {
//     return _numQuadPoints;
//   }
// 
//   const DomVecType& getRefCoord ( int QuadraturePoint ) const {
//     return _c[QuadraturePoint];
//   }
// 
//   RealType getWeight ( int quadpoint) const {
//     return _w[quadpoint];
//   }
// };


// DuodecQuadrature, exact for polynomials of degree 6
/*template <typename RealType, typename DomVecType>
class DuodecQuadrature {
public :
//   static constexpr unsigned short int _numQuadPoints = 12;
  static const unsigned short int _numQuadPoints = 12;
private :
  
  static const DomVecType _c[_numQuadPoints] = { DomVecType ( 0.24928674517091 , 0.24928674517091 ),
                                         DomVecType ( 0.24928674517091 , 0.50142650965818 ),
                                         DomVecType ( 0.50142650965818 , 0.24928674517091 ),
                                         DomVecType ( 0.06308901449150 , 0.06308901449150 ),
                                         DomVecType ( 0.06308901449150 , 0.87382197101700 ),
                                         DomVecType ( 0.87382197101700 , 0.06308901449150 ),
                                         DomVecType ( 0.31035245103378 , 0.63650249912140 ),
                                         DomVecType ( 0.63650249912140 , 0.05314504984482 ),
                                         DomVecType ( 0.05314504984482 , 0.31035245103378 ),
                                         DomVecType ( 0.63650249912140 , 0.31035245103378 ),
                                         DomVecType ( 0.31035245103378 , 0.05314504984482 ),
                                         DomVecType ( 0.05314504984482 , 0.63650249912140 )
                                       };
  static const RealType _w[_numQuadPoints] = { 
                     0.11678627572638 ,
                     0.11678627572638 ,
                     0.11678627572638 ,
                     0.05084490637021 ,
                     0.05084490637021 ,
                     0.05084490637021 ,
                     0.08285107561837 ,
                     0.08285107561837 ,
                     0.08285107561837 ,
                     0.08285107561837 ,
                     0.08285107561837 ,
                     0.08285107561837 
                  };
  
  
public:
    
  DuodecQuadrature() {}
  
  const int getNumQuadpoints ( ) const {
    return _numQuadPoints;
  }

  const DomVecType& getRefCoord ( int QuadraturePoint ) const {
    return _c[QuadraturePoint];
  }

  RealType getWeight ( int quadpoint) const {
    return _w[quadpoint];
  }
};*/


// SeptiQuadrature, exact for polynomials of degree 5
template <typename RealType, typename DomVecType>
class DuodecQuadrature {
                                   
public:
  enum { numQuadPoints = 12 };

  inline const DomVecType& getRefCoord ( int quadpoint ) const  {
       static const DomVecType c[12] = { DomVecType ( 0.24928674517091 , 0.24928674517091 ),
                                         DomVecType ( 0.24928674517091 , 0.50142650965818 ),
                                         DomVecType ( 0.50142650965818 , 0.24928674517091 ),
                                         DomVecType ( 0.06308901449150 , 0.06308901449150 ),
                                         DomVecType ( 0.06308901449150 , 0.87382197101700 ),
                                         DomVecType ( 0.87382197101700 , 0.06308901449150 ),
                                         DomVecType ( 0.31035245103378 , 0.63650249912140 ),
                                         DomVecType ( 0.63650249912140 , 0.05314504984482 ),
                                         DomVecType ( 0.05314504984482 , 0.31035245103378 ),
                                         DomVecType ( 0.63650249912140 , 0.31035245103378 ),
                                         DomVecType ( 0.31035245103378 , 0.05314504984482 ),
                                         DomVecType ( 0.05314504984482 , 0.63650249912140 )
                                       };
   return c[quadpoint];
  }
  inline RealType getWeight ( int quadpoint) const {
      static const RealType w[12] = { 
                     0.11678627572638 ,
                     0.11678627572638 ,
                     0.11678627572638 ,
                     0.05084490637021 ,
                     0.05084490637021 ,
                     0.05084490637021 ,
                     0.08285107561837 ,
                     0.08285107561837 ,
                     0.08285107561837 ,
                     0.08285107561837 ,
                     0.08285107561837 ,
                     0.08285107561837 
                  };    
    return w[quadpoint];
  }
};




// // TredecQuadrature, exact for polynomials of degree 6
// template <typename RealType>
// class TredecQuadrature {
// public:
//   enum { numQuadPoints = 13 };
// 
//   inline const DomVecType& getRefCoord ( int quadpoint ) const  {
//     static const DomVecType c[13] = {
//                                               DomVecType ( 0.33333333333333 , 0.33333333333333 ),
//                                               DomVecType ( 0.26034596607904 , 0.26034596607904 ),
//                                               DomVecType ( 0.26034596607904 , 0.47930806784192 ),
//                                               DomVecType ( 0.47930806784192 , 0.26034596607904 ),
//                                               DomVecType ( 0.06513010290222 , 0.06513010290222 ),
//                                               DomVecType ( 0.06513010290222 , 0.86973979419557 ),
//                                               DomVecType ( 0.86973979419557 , 0.06513010290222 ),
//                                               DomVecType ( 0.31286549600487 , 0.63844418856981 ),
//                                               DomVecType ( 0.63844418856981 , 0.04869031542532 ),
//                                               DomVecType ( 0.04869031542532 , 0.31286549600487 ),
//                                               DomVecType ( 0.63844418856981 , 0.31286549600487 ),
//                                               DomVecType ( 0.31286549600487 , 0.04869031542532 ),
//                                               DomVecType ( 0.04869031542532 , 0.63844418856981 )
//                                             };
//     return c[quadpoint];
//   }
// 
//   inline RealType getWeight ( int quadpoint) const {
//     static const RealType w[13] = { 
//                     RealType ( -0.14957004446768 ),
//                     RealType ( 0.17561525743321 ),
//                     RealType ( 0.17561525743321 ),
//                     RealType ( 0.17561525743321 ),
//                     RealType ( 0.05334723560884 ),
//                     RealType ( 0.05334723560884 ),
//                     RealType ( 0.05334723560884 ),
//                     RealType ( 0.07711376089026 ),
//                     RealType ( 0.07711376089026 ),
//                     RealType ( 0.07711376089026 ),
//                     RealType ( 0.07711376089026 ),
//                     RealType ( 0.07711376089026 ),
//                     RealType ( 0.07711376089026 )
//                   };
//     return w[quadpoint];
//   }
// };
// 
// /*!

//  */
// // SexdecQuadrature, exact for polynomials of degree 8
// template <typename RealType>
// class SexdecQuadrature {
// public:
//   enum { numQuadPoints = 16 };
// 
//   inline const DomVecType& getRefCoord ( int quadpoint ) const  {
//     static const DomVecType c[16] = {
//     DomVecType ( 0.333333333333333333333333333333333 , 0.333333333333333333333333333333333 ),
//     DomVecType ( 0.459292588292723156028815514494169 , 0.459292588292723156028815514494169 ),
//     DomVecType ( 0.459292588292723156028815514494169 , 0.081414823414553687942368971011661 ),
//     DomVecType ( 0.081414823414553687942368971011661 , 0.459292588292723156028815514494169 ),
//     DomVecType ( 0.170569307751760206622293501491464 , 0.170569307751760206622293501491464 ),
//     DomVecType ( 0.170569307751760206622293501491464 , 0.658861384496479586755412997017071 ),
//     DomVecType ( 0.658861384496479586755412997017071 , 0.170569307751760206622293501491464 ),
//     DomVecType ( 0.050547228317030975458423550596598 , 0.050547228317030975458423550596598 ),
//     DomVecType ( 0.050547228317030975458423550596598 , 0.898905543365938049083152898806802 ),
//     DomVecType ( 0.898905543365938049083152898806802 , 0.050547228317030975458423550596598 ),
//     DomVecType ( 0.263112829634638113421785786284643 , 0.728492392955404281241000379176062 ),
//     DomVecType ( 0.728492392955404281241000379176062 , 0.008394777409957605337213834539294 ),
//     DomVecType ( 0.008394777409957605337213834539294 , 0.263112829634638113421785786284643 ),
//     DomVecType ( 0.728492392955404281241000379176062 , 0.263112829634638113421785786284643 ),
//     DomVecType ( 0.263112829634638113421785786284643 , 0.008394777409957605337213834539294 ),
//     DomVecType ( 0.008394777409957605337213834539294 , 0.728492392955404281241000379176062 )
//     };
//     return c[quadpoint];
//   }
// 
//   inline RealType getWeight ( int quadpoint) const {
//     static const RealType w[16] = { 
//                     RealType ( 0.144315607677787168251091110489064 ),
//                     RealType ( 0.095091634267284624793896104388584 ),
//                     RealType ( 0.095091634267284624793896104388584 ),
//                     RealType ( 0.095091634267284624793896104388584 ),
//                     RealType ( 0.103217370534718250281791550292129 ),
//                     RealType ( 0.103217370534718250281791550292129 ),
//                     RealType ( 0.103217370534718250281791550292129 ),
//                     RealType ( 0.032458497623198080310925928341780 ),
//                     RealType ( 0.032458497623198080310925928341780 ),
//                     RealType ( 0.032458497623198080310925928341780 ),
//                     RealType ( 0.027230314174434994264844690073908 ),
//                     RealType ( 0.027230314174434994264844690073908 ),
//                     RealType ( 0.027230314174434994264844690073908 ),
//                     RealType ( 0.027230314174434994264844690073908 ),
//                     RealType ( 0.027230314174434994264844690073908 ),
//                     RealType ( 0.027230314174434994264844690073908 )
//                   };
//     return w[quadpoint];
//   }
// };
// 
// // GaussDegree8Quadrature, exact for polynomials of degree 10
// template <typename RealType>
// class GaussDegree10Quadrature {
// public:
//   enum { numQuadPoints = 25 };
// 
//   inline const DomVecType& getRefCoord ( int quadpoint ) const  {
//     static const DomVecType c[25] = {
//     DomVecType ( 0.333333333333333 , 0.333333333333333 ),
//     DomVecType ( 0.028844733232685 , 0.485577633383657 ),
//     DomVecType ( 0.485577633383657 , 0.028844733232685 ),
//     DomVecType ( 0.485577633383657 , 0.485577633383657 ),
//     DomVecType ( 0.781036849029926 , 0.109481575485037 ),
//     DomVecType ( 0.109481575485037 , 0.781036849029926 ),
//     DomVecType ( 0.109481575485037 , 0.109481575485037 ),
//     DomVecType ( 0.141707219414880 , 0.307939838764121 ),
//     DomVecType ( 0.141707219414880 , 0.550352941820999 ),
//     DomVecType ( 0.307939838764121 , 0.141707219414880 ),
//     DomVecType ( 0.307939838764121 , 0.550352941820999 ),
//     DomVecType ( 0.550352941820999 , 0.141707219414880 ),
//     DomVecType ( 0.550352941820999 , 0.307939838764121 ),
//     DomVecType ( 0.025003534762686 , 0.246672560639903 ),
//     DomVecType ( 0.025003534762686 , 0.728323904597411 ),
//     DomVecType ( 0.246672560639903 , 0.025003534762686 ),
//     DomVecType ( 0.246672560639903 , 0.728323904597411 ),
//     DomVecType ( 0.728323904597411 , 0.025003534762686 ),
//     DomVecType ( 0.728323904597411 , 0.246672560639903 ),
//     DomVecType ( 0.009540815400299 , 0.066803251012200 ),
//     DomVecType ( 0.009540815400299 , 0.923655933587500 ),
//     DomVecType ( 0.066803251012200 , 0.009540815400299 ),
//     DomVecType ( 0.066803251012200 , 0.923655933587500 ),
//     DomVecType ( 0.923655933587500 , 0.009540815400299 ),
//     DomVecType ( 0.923655933587500 , 0.066803251012200 )
//     };
//     return c[quadpoint];
//   }
// 
//   inline RealType getWeight ( int quadpoint) const {
//     static const RealType w[25] = { 
//                     RealType ( 0.090817990382754 ),
//                     RealType ( 0.036725957756467 ),
//                     RealType ( 0.036725957756467 ),
//                     RealType ( 0.036725957756467 ),
//                     RealType ( 0.045321059435528 ),
//                     RealType ( 0.045321059435528 ),
//                     RealType ( 0.045321059435528 ),
//                     RealType ( 0.072757916845420 ),
//                     RealType ( 0.072757916845420 ),
//                     RealType ( 0.072757916845420 ),
//                     RealType ( 0.072757916845420 ),
//                     RealType ( 0.072757916845420 ),
//                     RealType ( 0.072757916845420 ),
//                     RealType ( 0.028327242531057 ),
//                     RealType ( 0.028327242531057 ),
//                     RealType ( 0.028327242531057 ),
//                     RealType ( 0.028327242531057 ),
//                     RealType ( 0.028327242531057 ),
//                     RealType ( 0.028327242531057 ),
//                     RealType ( 0.009421666963733 ),
//                     RealType ( 0.009421666963733 ),
//                     RealType ( 0.009421666963733 ),
//                     RealType ( 0.009421666963733 ),
//                     RealType ( 0.009421666963733 ),
//                     RealType ( 0.009421666963733 )
//                   };
//     return w[quadpoint];
//   }
// };

}

#endif
