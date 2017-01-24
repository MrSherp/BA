#ifndef FEM_H
#define FEM_H

#include <sys/stat.h>
#include <math.h> 

namespace fe {


/** Class for uniform grids whose extensions do not have to be a power of 2. Moreover the grid can have
 *  different extent in the various coordinate directions
 */

template <qc::Dimension Dim>
class RectangularGrid : public GridStructure {
public:
  typedef RectangularGrid<Dim>          Self;
  typedef Self                          CubicGridType;

protected:
  double _h;
  qc::FastILexMapper< Dim > _indexMapper;

protected:
  void initializeIterators( ) {
    begin_it.getCurrentPosition().set ( 0, 0, 0 );

    if ( Dim == qc::QC_1D ) {
      end_it.getCurrentPosition().set ( 0, 1, 0 );
    } else if ( Dim == qc::QC_2D ) {
      end_it.getCurrentPosition().set ( 0, 0, 1 );
    } else if ( Dim == qc::QC_3D ) {
      end_it.getCurrentPosition().set ( 0, 0, this->getNumZ() - 1 );
    } else {
      throw aol::Exception ( "unsupported dimension", __FILE__, __LINE__ );
    }
  }

public:
  explicit RectangularGrid ( const RectangularGrid &other )
    : GridStructure ( other ),
      _h ( other._h ),
      begin_it ( other.getSize() ), end_it ( other.getSize() ) {
    _indexMapper.resize ( qc::CoordType ( this->getNumX(), this->getNumY(), this->getNumZ() ) );
    initializeIterators( );
  }

  explicit RectangularGrid ( const aol::Vec3<int> &size )
    : GridStructure ( size, Dim ),
      _h ( 1. / aol::Max ( size[0] - 1, size[1] - 1 , size[2] - 1 ) ),
      begin_it ( size ),
      end_it ( size ) {
    _indexMapper.resize ( qc::CoordType ( this->getNumX(), this->getNumY(), this->getNumZ() ) );
    initializeIterators( );
  }

  explicit RectangularGrid ( const qc::GridSize<Dim> &gridSize )
    : GridStructure ( aol::Vec3<int> ( gridSize.getNumX(), ( Dim > QC_1D ? gridSize.getNumY() : 1 ), ( Dim > QC_2D ? gridSize.getNumZ() : 1 ) ), Dim ),
      _h ( 1. / aol::Max ( static_cast<int>( gridSize.getNumX() - 1 ), ( Dim > QC_1D ? gridSize.getNumY() - 1 : 1 ), ( Dim > QC_2D ? gridSize.getNumZ() - 1 : 1 ) ) ),
      begin_it ( aol::Vec3<int> ( gridSize.getNumX(), ( Dim > QC_1D ? gridSize.getNumY() : 1 ), ( Dim > QC_2D ? gridSize.getNumZ() : 1 ) ) ),
      end_it ( aol::Vec3<int> ( gridSize.getNumX(), ( Dim > QC_1D ? gridSize.getNumY() : 1 ), ( Dim > QC_2D ? gridSize.getNumZ() : 1 ) ) ) {
    _indexMapper.resize ( qc::CoordType ( gridSize.getNumX(), ( Dim > QC_1D ? gridSize.getNumY() : 1 ), ( Dim > QC_2D ? gridSize.getNumZ() : 1 ) ) );
    initializeIterators( );
  }

  /**
   * This constructor mimics the behavior of the usual GridDefinition constructor.
   * Used for compatibility reasons.
   *
   * \author Berkels
   */
  explicit RectangularGrid ( const int GridDepth, const Dimension DimOfWorld = qc::QC_2D )
    : GridStructure ( aol::Vec3<int>( (( 1 << GridDepth ) + 1), (( 1 << GridDepth ) + 1), (Dim == QC_2D) ? 1 : (( 1 << GridDepth ) + 1) ), Dim ),
      _h ( 1. / ( aol::Max ( aol::Max ( this->getNumX() - 1, this->getNumY() - 1 ), this->getNumZ() - 1 ) ) ),
      begin_it ( this->getSize() ), // ???????
      end_it ( this->getSize() ) {
    if ( DimOfWorld != Dim )
      throw aol::Exception ( "DimOfWorld != Dim", __FILE__, __LINE__ );
    _indexMapper.resize ( qc::CoordType ( this->getNumX(), this->getNumY(), this->getNumZ() ) );
    initializeIterators( );
  }

  int getNumberOfElements() const {
    if(Dim == QC_3D) {
      return (this->getNumX()-1)*(this->getNumY()-1)*(this->getNumZ()-1);
    } else {
      return (this->getNumX()-1)*(this->getNumY()-1);
    }
  }

  const Self & getCubicGrid() const {
    return *this;
  }

  const Self & getFullGrid() const {
    return *this;
  }

  int getGridDepth() const {
    throw aol::Exception ( "RectangularGrid does not support getGridDepth()!", __FILE__, __LINE__ );
  }

  int getElementIndex ( const Element & /*El*/ ) const {
    throw aol::Exception ( "RectangularGrid does not support getElementIndex()!", __FILE__, __LINE__ );
  }

  // this iterator iterates over the elements of this grid
  struct OldAllElementIterator {
    typedef OldAllElementIterator _Self;
    typedef Element  IteratedType;
    typedef _Self    BeginType;
    typedef _Self    EndType;
    aol::Vec<3, int>  _size;

    OldAllElementIterator() : _size() {}

    explicit OldAllElementIterator ( const aol::Vec<3, int> &size ) : _size ( size )  { }

    qc::Element &getCurrentPosition() {
      return _cur;
    }

    const qc::Element &getCurrentPosition() const {
      return _cur;
    }

    OldAllElementIterator &operator= ( const OldAllElementIterator& Other ) {
      _size = Other._size;
      _cur = Other._cur;
      return *this;
    }

    inline OldAllElementIterator& operator++() {
      _cur.xref() ++;
      if ( _cur.x() == ( _size[0] - 1 ) ) {
        _cur.xref() = 0;
        _cur.yref() ++;
        if ( _cur.y() == ( _size[1] - 1 ) ) {
          _cur.yref() = 0;
          _cur.zref() ++;
        }
      }
      return *this;
    }

    inline OldAllElementIterator operator++ ( int ) {
      OldAllElementIterator copy ( *this );
      ++ ( *this );
      return copy;
    }

    qc::Element& operator*() {
      return _cur;
    }

    const qc::Element& operator*() const {
      return _cur;
    }

    const qc::Element* operator->() const {
      return &_cur;
    }

    inline bool operator!= ( const _Self& Other ) const {
      return ( _cur != Other._cur );
    }

    inline bool operator== ( const _Self& Other ) const {
      return ( _cur == Other._cur );
    }

  protected:
    qc::Element _cur;
    // End of internal class OldAllElementIterator
  };

  OldAllElementIterator begin_it;
  OldAllElementIterator end_it;

  double H() const {
    return _h;
  }

  void setH ( double h ) {
    _h = h;
  }

  //! return width = number of nodes in x direction
  int getWidth() const {
    return ( getNumX() );
  }

  //! return height = number of nodes in y direction
  int getHeight() const {
    return ( getNumY() );
  }

  //! return depth = number of nodes in z direction
  int getDepth() const {
    return ( getNumZ() );
  }

  const qc::FastILexMapper<Dim>& getIndexMapperRef() const {
    return ( _indexMapper );
  }

  //! return grid depth for initializing an Element iterator, returning -1 for rectangularGrids and to be overloaded on derived classes.
  virtual short getElementGridDepth ( ) const {
    return ( -1 );
  }

  //! return the volume fraction of the brick that is represented by this RectangularGrid relative to the unit cube
  double getVolumeFractionOf01Cube ( ) const {
    switch ( Dim ) {
    case QC_2D:
      return ( ( getNumX() - 1 ) * ( getNumY() - 1 ) * aol::Sqr ( _h ) );
    case QC_3D:
      return ( ( getNumX() - 1 ) * ( getNumY() - 1 ) * ( getNumZ() - 1 ) * aol::Cub ( _h ) );
    default:
      throw aol::UnimplementedCodeException ( "qc::RectangularGrid::getVolumeFractionOf01Cube not implemented for dimension != 2, 3", __FILE__, __LINE__ );
    }
  }

  /** Iterator that iterates over all elements in the grid
   *  \author Schwen
   */
  class FullElementIterator : public qc::RectangularIteratorBase< Dim, qc::Element > {
  public:
    explicit FullElementIterator ( const RectangularGrid<Dim> &grid ) : qc::RectangularIteratorBase< Dim, qc::Element > ( qc::Element(), qc::Element() ) {
      this->_lower.set ( 0, 0, 0, grid.getElementGridDepth() );
      this->_current.set ( 0, 0, 0, grid.getElementGridDepth() );
      this->_upper.set ( grid.getSize()[0] - 1, grid.getSize()[1] - 1, grid.getSize()[2] - 1, grid.getElementGridDepth() );
    }

    //! prefix increment operator.
    const qc::Element& operator++ ( ) {
      this->increment();
      return ( this->_current );
    }
  };

  /** Iterator that iterates over all nodes in the grid
   *  \author Schwen
   */
  class FullNodeIterator : public qc::RectangularIterator< Dim > {
  public:
    explicit FullNodeIterator ( const qc::RectangularGrid<Dim> &grid ) : qc::RectangularIterator<Dim> ( qc::CoordType(), qc::CoordType ( grid.getSize() ) ) {
    }
  };

  /** Iterator that iterates over all boundary nodes in the grid
   *  \warning currently uses slow but simple boundary node iterator
   *  \author Schwen
   */
  class FullBoundaryNodeIterator : public qc::RectangularBoundaryIterator< Dim > {
  public:
    explicit FullBoundaryNodeIterator ( const qc::RectangularGrid<Dim> &grid ) : qc::RectangularBoundaryIterator<Dim> ( qc::CoordType(), qc::CoordType( grid.getSize() ) ) {
    }
  };


};


/** Possible replacement for qc::GridDefinition
 *  \todo think about relation to RectangularGrid, move _indexMapper there?
 *  \todo think about placement of code of basis classes (better here ...)
 */
template <qc::Dimension Dim>
class CubicGrid : public RectangularGrid<Dim> {
protected:
  int _depth;

public:
  explicit CubicGrid ( const int Depth ) : RectangularGrid<Dim> ( Depth, Dim ), _depth ( Depth ) {
  }

  //! numX, numY and numZ are the same
  inline int getNumXYZ ( ) const {
    return ( this->getNumX() );
  }

  int getGridDepth() const {
    return ( _depth );
  }

  int getNumberOfBoundaryNodes() const {
    if ( Dim == qc::QC_2D ) {
      return ( 4 * getNumXYZ() - 4 );
    } else if ( Dim == qc::QC_3D ) {
      return ( 6 * ( getNumXYZ() - 1 ) * ( getNumXYZ() - 1 ) + 2 );
    } else {
      throw aol::Exception ( "unsupported dimension", __FILE__, __LINE__ );
      return ( - 1 );
    }
  }

  virtual short getElementGridDepth ( ) const {
    return ( _depth );
  }

};

    
    
    
template <typename RealType, qc::Dimension Dim>
class ConfiguratorBase {
public:

  typedef qc::Element                           ElementType;            //!< use the quoc element here
  typedef qc::GridDefinition::OldAllElementIterator  ElementIteratorType;    //!< use the standard element iterator of qc::GridDefinition
  // typedef UniformGridSparseMatrix<RealType> MatrixType;          //!< for uniform grids, the UniformSparseMatrices are most efficient
  // typedef qc::FastUniformGridMatrix<RealType,Dim>   MatrixType;
  typedef RealType Real;
  typedef qc::GridDefinition                    InitType;               //!< that's the type, that is needed by the constructor of the configurator
  typedef typename qc::BitArray<Dim>            MaskType;

  static const qc::Dimension DimOfWorld = Dim;

  explicit QuocConfiguratorTraitBase ( const InitType &Grid ) :
      _grid ( Grid ) {}

  //! returns the begin iterator of the grid
  const ElementIteratorType &begin( ) const {
    return _grid.begin_it;
  }

  //! returns the end iterator of the grid
  inline const ElementIteratorType &end( ) const {
    return _grid.end_it;
  }

  RealType H ( const qc::Element& ) const {
    return _grid.H();
  }

  const InitType& getInitializer( ) const { return _grid; }

  //! method needed in boundaryIntegration.h
  //! although the name suggests something different this implementation reflects the current behavior
  inline int localOnFaceToLocal ( const qc::Element &/*El*/, const int localIndex ) const {
    return localIndex;
  }

protected:
  const InitType &_grid;   //!< memorize reference to the grid
  //qc::IndexMapper _mapper;         //!< that's the index mapper
};

template < typename RealType,
           qc::Dimension Dim,
           typename QuadType,
           typename _MatrixType = qc::FastUniformGridMatrix< RealType, Dim> >
class QuocConfiguratorTraitMultiLin {};

/**
 * \ingroup FEConfigurator
 */
template <typename _RealType, typename _QuadType, typename _MatrixType >
class QuocConfiguratorTraitMultiLin<_RealType, qc::QC_1D, _QuadType, _MatrixType> : public QuocConfiguratorTraitBase<_RealType, qc::QC_1D> {
  aol::BaseFunctionSetMultiLin<_RealType, qc::QC_1D, _QuadType> _baseFuncSet;
public:
  explicit QuocConfiguratorTraitMultiLin ( const qc::GridDefinition &Grid )
      : QuocConfiguratorTraitBase<_RealType, qc::QC_1D> ( Grid ),  _baseFuncSet ( Grid.H() ), _volEl ( Grid.H() ), _mapper ( Grid ) {}

  typedef aol::Vector<_RealType>       VectorType;
  typedef aol::Vec<1, _RealType>       VecType;
  typedef aol::Vec<1, _RealType>       DomVecType;
  typedef aol::Mat<1, 1, _RealType>    MatType;
  typedef _MatrixType                  MatrixType;
  typedef aol::BaseFunctionSetMultiLin<_RealType, qc::QC_1D, _QuadType> BaseFuncSetType;
  typedef _QuadType                    QuadType;
  typedef _RealType                    RealType;
  typedef qc::ScalarArray<_RealType, qc::QC_1D> ArrayType;
  typedef aol::FullMatrix<RealType>    FullMatrixType;
  typedef qc::QuocConfiguratorTraitMultiLin<_RealType, qc::QC_1D, _QuadType, _MatrixType> FullGridConfiguratorType;

  static const int maxNumLocalDofs = 2;

  static const qc::Dimension Dim = qc::QC_1D;
  static const qc::Dimension DomDim = qc::QC_1D;

  const RealType _volEl;

  inline bool getLocalCoords ( const VecType &Coord, qc::Element &El, DomVecType &LocalCoord ) const {
    return getLocalCoordsRegularRectangularGrid<QuocConfiguratorTraitMultiLin<_RealType, qc::QC_1D, _QuadType, _MatrixType> > ( Coord, this->_grid, El, LocalCoord );
  }

  int maxNumQuadPoints( ) const {
    return _QuadType::numQuadPoints;
  }

  inline void getGlobalCoords ( const qc::Element &El, const DomVecType &LocalCoord, VecType &Coord ) const {
    for ( int c = 0; c < 1; c++ ) {
      Coord[c] = ( static_cast<RealType> ( El[c] ) + LocalCoord[c] ) * this->_grid.H();
    }
  }

  inline int getNumLocalDofs ( const qc::Element & ) const {
    return 2;
  }

  int getNumGlobalDofs( ) const {
    return this->_grid.getNumberOfNodes();
  }

  qc::Element getEmptyElement() const {
    return qc::Element();
  }

  const BaseFuncSetType& getBaseFunctionSet ( const qc::Element& /*El*/ ) const {
    return _baseFuncSet;
  }

  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const qc::Element &El, const int localIndex ) const {
    return _mapper.localToGlobal ( El, localIndex );
  }

  inline void localToGlobal ( const qc::Element &El, const int localIndex0, const int localIndex1, aol::Vec2<int> &glob ) const {
    glob[0] = localToGlobal( El, localIndex0 );
    glob[1] = localToGlobal( El, localIndex1 );
  }

  RealType vol ( const qc::Element & ) const {
    return _volEl;
  }

  RealType H ( const qc::Element & ) const {
    return this->_grid.H();
  }

  //! create a new, clean matrix
  MatrixType* createNewMatrix( ) const {
    MatrixType *mat = new MatrixType ( qc::GridSize<qc::QC_1D>::createFrom ( this->_grid ) );
    // mat->clearRows( );
    return mat;
  }

  const qc::FastILexMapper<Dim>& getIndexMapper() const {
    return _mapper;
  }

  using QuocConfiguratorTraitBase<_RealType, qc::QC_1D>::localOnFaceToLocal;

  template <bool ClipCoord>
  inline bool transformCoord ( const qc::Element &El, const VecType &RefCoord, const VecType &Offset, qc::Element &TransformedEl, VecType &TransformedLocalCoord ) const {
    return qc::transformCoord<QuocConfiguratorTraitMultiLin<_RealType, qc::QC_1D, _QuadType, _MatrixType>, ClipCoord> ( this->_grid, El, RefCoord, Offset, TransformedEl, TransformedLocalCoord );
  }
protected:
  qc::FastILexMapper<Dim> _mapper;         //!< that's the index mapper
};

/**
 * \ingroup FEConfigurator
 */
template <typename _RealType, typename _QuadType, typename _MatrixType >
class QuocConfiguratorTraitMultiLin<_RealType, qc::QC_2D, _QuadType, _MatrixType> : public QuocConfiguratorTraitBase<_RealType, qc::QC_2D> {
  aol::BaseFunctionSetMultiLin<_RealType, qc::QC_2D, _QuadType> _baseFuncSet;
public:
  explicit QuocConfiguratorTraitMultiLin ( const qc::GridDefinition &Grid )
      : QuocConfiguratorTraitBase<_RealType, qc::QC_2D> ( Grid ),  _baseFuncSet ( Grid.H() ), _volEl ( aol::Sqr ( Grid.H() ) ), _mapper ( Grid ) {}

  typedef aol::Vector<_RealType>       VectorType;
  typedef aol::Vec2<_RealType>         VecType;
  typedef aol::Vec2<_RealType>         DomVecType;
  typedef aol::Matrix22<_RealType>     MatType;
  typedef _MatrixType                  MatrixType;
  typedef aol::BaseFunctionSetMultiLin<_RealType, qc::QC_2D, _QuadType> BaseFuncSetType;
  typedef _QuadType                    QuadType;
  typedef _RealType                    RealType;
  typedef qc::ScalarArray<_RealType, qc::QC_2D> ArrayType;
  typedef qc::BitArray<qc::QC_2D>               MaskType;
  typedef aol::FullMatrix<RealType>    FullMatrixType;
  typedef qc::QuocConfiguratorTraitMultiLin<_RealType, qc::QC_2D, _QuadType, _MatrixType> FullGridConfiguratorType;

  static const int maxNumLocalDofs = 4;

  static const qc::Dimension Dim = qc::QC_2D;
  static const qc::Dimension DomDim = qc::QC_2D;
  static const aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE;

  const RealType _volEl;

  inline bool getLocalCoords ( const VecType &Coord, qc::Element &El, DomVecType &LocalCoord ) const {
    return getLocalCoordsRegularRectangularGrid<QuocConfiguratorTraitMultiLin<_RealType, qc::QC_2D, _QuadType, _MatrixType> > ( Coord, this->_grid, El, LocalCoord );
  }


  int maxNumQuadPoints( ) const {
    return _QuadType::numQuadPoints;
  }

  inline void getGlobalCoords ( const qc::Element &El, const DomVecType &LocalCoord, VecType &Coord ) const {
    for ( int c = 0; c < 2; c++ ) {
      Coord[c] = ( static_cast<RealType> ( El[c] ) + LocalCoord[c] ) * this->_grid.H();
    }
  }

  inline int getNumLocalDofs ( const qc::Element & ) const {
    return 4;
  }

  int getNumGlobalDofs( ) const {
    return this->_grid.getNumberOfNodes();
  }

  qc::Element getEmptyElement() const {
    return qc::Element();
  }

  const BaseFuncSetType& getBaseFunctionSet ( const qc::Element& /*El*/ ) const {
    return _baseFuncSet;
  }

  //! get Number of Element ( = Number of ll dof ) (BG)
  int getElementNumber ( const qc::Element & El ) const {
    return _mapper.localToGlobal ( El, 0 );
  }

  //! compute consecutive element numbers (BG)
  int getConsecutiveElementNumber ( const qc::Element & El ) const {
    int elNum = getElementNumber( El );
    return elNum - ( elNum / this->_grid.getNumX() );
  }
  
  //! get element from consecutive number
  void getElementfromConsecutiveNumber ( int number, qc::Element & El ) const {
    number += number / (this->_grid.getNumX()-1);
    _mapper.splitGlobalIndex ( number, El ) ;
    El.setLevel( this->_grid.getGridDepth() );
  }

  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const qc::Element &El, const int localIndex ) const {
    return _mapper.localToGlobal ( El, localIndex );
  }

  //! returns a vec2 of global indices i,j of the dofs with number localIndex_i,j
  inline void localToGlobal ( const qc::Element &El, const int localIndex0, const int localIndex1, aol::Vec2<int> &glob ) const {
    glob[0] = localToGlobal( El, localIndex0 );
    glob[1] = localToGlobal( El, localIndex1 );
  }


  RealType vol ( const qc::Element & ) const {
    return _volEl;
  }

  RealType H ( const qc::Element & ) const {
    return this->_grid.H();
  }

  //! create a new, clean matrix
  MatrixType* createNewMatrix( ) const {
    MatrixType *mat = new MatrixType ( qc::GridSize<qc::QC_2D>::createFrom ( this->_grid ) );
    // mat->clearRows( );
    return mat;
  }

  const qc::FastILexMapper<Dim>& getIndexMapper() const {
    return _mapper;
  }

  using QuocConfiguratorTraitBase<_RealType, qc::QC_2D>::localOnFaceToLocal;

  template <bool ClipCoord>
  inline bool transformCoord ( const qc::Element &El, const VecType &RefCoord, const VecType &Offset, qc::Element &TransformedEl, VecType &TransformedLocalCoord ) const {
    return qc::transformCoord<QuocConfiguratorTraitMultiLin<_RealType, qc::QC_2D, _QuadType, _MatrixType>, ClipCoord> ( this->_grid, El, RefCoord, Offset, TransformedEl, TransformedLocalCoord );
  }
protected:
  qc::FastILexMapper<Dim> _mapper;         //!< that's the index mapper
};


/**
 * \ingroup FEConfigurator
 */
template <typename _RealType, typename _QuadType, typename _MatrixType >
class QuocConfiguratorTraitMultiLin<_RealType, qc::QC_3D, _QuadType, _MatrixType> : public QuocConfiguratorTraitBase<_RealType, qc::QC_3D> {
  aol::BaseFunctionSetMultiLin<_RealType, qc::QC_3D, _QuadType> _baseFuncSet;
public:
  explicit QuocConfiguratorTraitMultiLin ( const qc::GridDefinition &Grid )
    : QuocConfiguratorTraitBase<_RealType, qc::QC_3D> ( Grid ), _baseFuncSet( Grid.H() ), _volEl ( aol::Cub ( Grid.H() ) ), _mapper ( Grid )  {}

  typedef qc::QuocConfiguratorTraitMultiLin<_RealType, qc::QC_3D, _QuadType, _MatrixType> Self;

  typedef aol::Vector<_RealType>       VectorType;
  typedef aol::Vec3<_RealType>         VecType;
  typedef aol::Vec3<_RealType>         DomVecType;
  typedef aol::Matrix33<_RealType>     MatType;
  typedef _MatrixType                  MatrixType;
  typedef aol::BaseFunctionSetMultiLin<_RealType, qc::QC_3D, _QuadType> BaseFuncSetType;
  typedef _QuadType                    QuadType;
  typedef _RealType                    RealType;
  typedef qc::ScalarArray<_RealType, qc::QC_3D> ArrayType;
  typedef qc::BitArray<qc::QC_3D>               MaskType;
  typedef aol::FullMatrix<RealType>    FullMatrixType;
  typedef Self                         FullGridConfiguratorType;

  static const int maxNumLocalDofs = 8;

  static const qc::Dimension Dim = qc::QC_3D;
  static const qc::Dimension DomDim = qc::QC_3D;
  static const aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE;

  const RealType _volEl;

  inline int getNumLocalDofs ( const qc::Element & ) const {
    return 8;
  }

  int getNumGlobalDofs( ) const {
    return this->_grid.getNumberOfNodes();
  }

  int maxNumQuadPoints( ) const {
    return _QuadType::numQuadPoints;
  }

  const BaseFuncSetType& getBaseFunctionSet ( const qc::Element& ) const {
    return _baseFuncSet;
  }

  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const qc::Element &El, const int localIndex ) const {
    return _mapper.localToGlobal ( El, localIndex );
  }

  //! returns a vec2 of global indices i,j of the dofs with number localIndex_i,j
  inline void localToGlobal ( const qc::Element &El, const int localIndex0, const int localIndex1, aol::Vec2<int> &glob ) const {
    glob[0] = localToGlobal( El, localIndex0 );
    glob[1] = localToGlobal( El, localIndex1 );
  }

  inline bool getLocalCoords ( const VecType &Coord, qc::Element &El, DomVecType &LocalCoord ) const {
    return getLocalCoordsRegularRectangularGrid<QuocConfiguratorTraitMultiLin<_RealType, qc::QC_3D, _QuadType, _MatrixType> > ( Coord, this->_grid, El, LocalCoord );
  }

  inline void getGlobalCoords ( const qc::Element &El, const DomVecType &LocalCoord, VecType &Coord ) const {
    for ( int c = 0; c < 3; c++ ) {
      Coord[c] = ( static_cast<RealType> ( El[c] ) + LocalCoord[c] ) * this->_grid.H();
    }
  }

  qc::Element getEmptyElement() const {
    return qc::Element();
  }

  //! create a new, clean matrix
  MatrixType* createNewMatrix( ) const {
    MatrixType *mat = new MatrixType ( qc::GridSize<qc::QC_3D>::createFrom ( this->_grid ) );
    // mat->clearRows( );
    return mat;
  }

  RealType vol ( const qc::Element & ) const {
    return _volEl;
  }

  RealType H ( const qc::Element & ) const {
    return this->_grid.H();
  }

  //! returns consecutive Element Number (ST)
  int getConsecutiveElementNumber (const qc::Element & El) const {
    int ElNum = localToGlobal(El, 0);
    return ElNum - ElNum / this->_grid.getNumX() - ElNum / ( this->_grid.getNumX() * this->_grid.getNumY() ) * ( this->_grid.getNumX() - 1 );
  }

  using QuocConfiguratorTraitBase<_RealType, qc::QC_3D>::localOnFaceToLocal;

  template <bool ClipCoord>
  inline bool transformCoord ( const qc::Element &El, const VecType &RefCoord, const VecType &Offset, qc::Element &TransformedEl, VecType &TransformedLocalCoord ) const {
    return qc::transformCoord<QuocConfiguratorTraitMultiLin<_RealType, qc::QC_3D, _QuadType, _MatrixType>, ClipCoord> ( this->_grid, El, RefCoord, Offset, TransformedEl, TransformedLocalCoord );
  }
protected:
  qc::FastILexMapper<qc::QC_3D> _mapper;         //!< that's the index mapper

};


template <typename RealType, qc::Dimension Dim, typename QuadType, typename _MatrixType = aol::SparseMatrix<RealType> >
class QuocConfiguratorTraitMultiQuad {
};

/**
 * \ingroup FEConfigurator
 */
template <typename _RealType, typename _QuadType, typename _MatrixType >
  class QuocConfiguratorTraitMultiQuad<_RealType, qc::QC_2D, _QuadType, _MatrixType> : public QuocConfiguratorTraitBase<_RealType, qc::QC_2D> {
  aol::BaseFunctionSetMultiQuad<_RealType, qc::QC_2D, _QuadType> _baseFuncSet;
public:
  explicit QuocConfiguratorTraitMultiQuad ( const qc::GridDefinition &Grid )
    : QuocConfiguratorTraitBase<_RealType, qc::QC_2D> ( Grid ), _baseFuncSet( Grid.H() ), _volEl( aol::Sqr( Grid.H() )), _mapper ( Grid ) {}

  typedef qc::QuocConfiguratorTraitMultiQuad<_RealType, qc::QC_2D, _QuadType, _MatrixType> Self;

  typedef aol::Vector<_RealType>       VectorType;
  typedef aol::Vec2<_RealType>         VecType;
  typedef aol::Vec2<_RealType>         DomVecType;
  typedef aol::Matrix22<_RealType>     MatType;
  typedef _MatrixType                  MatrixType;
  typedef aol::BaseFunctionSetMultiQuad<_RealType, qc::QC_2D, _QuadType> BaseFuncSetType;
  typedef _QuadType                    QuadType;
  typedef _RealType                    RealType;
  typedef qc::BitArray<qc::QC_2D>      MaskType;
  typedef ScalarArray<RealType, QC_2D> ArrayType;
  typedef aol::FullMatrix<RealType>    FullMatrixType;

  static const int maxNumLocalDofs = 9;
  static const aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE;

  static const qc::Dimension Dim = qc::QC_2D;
  static const qc::Dimension DomDim = qc::QC_2D;

  const RealType _volEl;

  inline int getNumLocalDofs ( const qc::Element & ) const {
    return 9;
  }

  int getNumGlobalDofs( ) const {
    return aol::Sqr ( ( this->_grid.getWidth() - 1 ) *2 + 1 );
  }

  const BaseFuncSetType& getBaseFunctionSet ( const qc::Element &/*El*/ ) const {
    return _baseFuncSet;
  }

  int maxNumQuadPoints( ) const {
    return _QuadType::numQuadPoints;
  }

  RealType vol ( const qc::Element & ) const {
    return _volEl;
  }

  RealType H ( const qc::Element & ) const {
    return this->_grid.H();
  }

  //! get Number of Element ( = Number of ll dof ) (BG)
  int getElementNumber ( const qc::Element & El ) const {
    return _mapper.localToGlobal ( El, 0 );
  }

  //! compute consecutive element numbers (BG)
  int getConsecutiveElementNumber ( const qc::Element & El ) const {
    int elNum = getElementNumber( El );
    return ( elNum - ( elNum / (2*this->_grid.getWidth()-1) ) )                     // eliminate last node index in each row
           / 2                                                                      // forget about virtual nodes in between
           - ( elNum / (4*this->_grid.getWidth()-2) ) * (this->_grid.getWidth()-1); // subtract intermediate rows of virtual nodes
  }

  //! get element from consecutive number
  void getElementfromConsecutiveNumber ( int number, qc::Element & El ) const {
    // element numbers only refer to real dofs
    number += number / (this->_grid.getNumX()-1);
    // do manual splitting wrt reals dofs
    El.set( number % this->_grid.getNumX(), number / this->_grid.getNumX(), 0, this->_grid.getGridDepth() );
  }

  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const qc::Element &El, const int localIndex ) const {
    return _mapper.localToGlobal ( El, localIndex );
  }

  //! returns a vec2 of global indices i,j of the dofs with number localIndex_i,j
  inline void localToGlobal ( const qc::Element &El, const int localIndex0, const int localIndex1, aol::Vec2<int> &glob ) const {
    glob[0] = localToGlobal( El, localIndex0 );
    glob[1] = localToGlobal( El, localIndex1 );
  }

  inline void getGlobalCoords ( const qc::Element &El, const DomVecType &LocalCoord, VecType &Coord ) const {
    for ( int c = 0; c < 2; c++ ) {
      Coord[c] = ( static_cast<RealType> ( El[c] ) + LocalCoord[c] ) * this->_grid.H();
    }
  }

  inline bool getLocalCoords ( const VecType &Coord, qc::Element &El, DomVecType &LocalCoord ) const {
    return getLocalCoordsRegularRectangularGrid<QuocConfiguratorTraitMultiQuad<_RealType, qc::QC_2D, _QuadType, _MatrixType> > ( Coord, this->_grid, El, LocalCoord );
  }

  //! create a new, clean matrix
  MatrixType* createNewMatrix( ) const {
    MatrixType *mat = new MatrixType ( getNumGlobalDofs( ), getNumGlobalDofs( ) );
    mat->setZero( );
    return mat;
  }

  using QuocConfiguratorTraitBase<_RealType, qc::QC_2D>::localOnFaceToLocal;

protected:
  qc::FastQuadILexMapper<qc::QC_2D> _mapper;         //!< that's the index mapper
};

    
    
    
    
}

#endif