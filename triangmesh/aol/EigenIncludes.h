//TODO in external/eigen
#ifndef __INCLUDEFROMEIGEN_H
#define __INCLUDEFROMEIGEN_H



#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wcast-qual"
#pragma GCC system_header
#endif


#include <Eigen/Core>

#include <Eigen/QR>
#include <Eigen/Eigenvalues>

#include <Eigen/Dense>
#include <Eigen/SparseCore>

// #include <Eigen/SparseQR>
#include <Eigen/LU>
#include <Eigen/UmfPackSupport>
#include <Eigen/CholmodSupport>
#include <Eigen/IterativeLinearSolvers>
// #include <unsupported/Eigen/NonLinearOptimization>

#include <tensor.h>



// #include <iostream>
// #include <Eigen/MPRealSupport>
// #include <Eigen/LU>
// using namespace mpfr;
// using namespace Eigen;
// int main()
// {
//   // set precision to 256 bits (double has only 53 bits)
//   mpreal::set_default_prec(256);
//   // Declare matrix and vector types with multi-precision scalar type
//   typedef Matrix<mpreal,Dynamic,Dynamic>  MatrixXmp;
//   typedef Matrix<mpreal,Dynamic,1>        VectorXmp;
//   MatrixXmp A = MatrixXmp::Random(100,100);
//   VectorXmp b = VectorXmp::Random(100);
//   // Solve Ax=b using LU
//   VectorXmp x = A.lu().solve(b);
//   std::cout << "relative error: " << (A*x - b).norm() / b.norm() << std::endl;
//   return 0;
// }




struct DataTypeContainerDOT {
public :
  typedef double    RealType;
  typedef Eigen::Vector2d GradientType;
  typedef Eigen::VectorXd  VectorType;
  typedef Eigen::VectorXi VectorTypeInt;
  typedef Eigen::MatrixXd FullMatrixType;
  typedef Eigen::SparseMatrix<double, 0, int> SparseMatrixType;
  typedef Eigen::Triplet<RealType> TripletType;
};



struct DataTypeContainerShellFE {
public :
  typedef double    RealType;
  typedef Eigen::Vector2d DomVecType;
  typedef Eigen::Vector3d TangentVecType;
  typedef Eigen::Vector3d Point3DType;
  typedef Eigen::Vector3i Indices3DType;
  
  typedef Eigen::Matrix<double, 2, 2> Matrix22;
  typedef Eigen::Matrix<double, 3, 2> Matrix32;
  typedef Eigen::Matrix<double, 3, 3> Matrix33;
  
  typedef aol::Tensor222<RealType, Matrix22> Tensor222Type;
  typedef aol::Tensor322<RealType, Matrix22, TangentVecType> Tensor322Type;
 
  typedef Eigen::VectorXd  VectorType;
  
  typedef std::vector<bool>  MaskType;
  
  typedef Eigen::MatrixXd FullMatrixType;
//   typedef Eigen::SparseMatrix<double, 0, long int> SparseMatrixType; TODO
  typedef Eigen::SparseMatrix<double, 0, int> SparseMatrixType;
  typedef Eigen::Triplet<RealType> TripletType;
};



namespace shellFE {

enum LinearSolverMethod {
    UmfPackLU = 1,
    CholmodSupernodalLLT = 2,
    ConjugateGradient = 3,
    BiCGSTAB = 4
};

template<typename DataTypeContainer>
void solveLinearSystem ( const typename DataTypeContainer::SparseMatrixType & systemMatrix , 
                         typename DataTypeContainer::VectorType &solution, 
                         const typename DataTypeContainer::VectorType &rhs,
                         const LinearSolverMethod linearSolver = UmfPackLU ) {

    typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;
                             
    switch( linearSolver ){

//       case UmfPackLU:{
//           Eigen::UmfPackLU<SparseMatrixType> solver;
//           solver.compute(systemMatrix);
//           if(solver.info()!=Eigen::Success) {
//               throw std::invalid_argument ( aol::strprintf ( "solver failed in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
//               return;
//           }
//           solution = solver.solve( rhs );
//           break;
//       }
      
      case CholmodSupernodalLLT:{
//           Eigen::CholmodSupernodalLLT<SparseMatrixType> solver;
          Eigen::CholmodDecomposition<SparseMatrixType> solver;
          solver.compute(systemMatrix);
          if(solver.info()!=Eigen::Success) {
              throw std::invalid_argument ( aol::strprintf ( "solver failed in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
              return;
          }
          solution = solver.solve( rhs );         
          break;   
      }
      
      case ConjugateGradient: {
            Eigen::ConjugateGradient<SparseMatrixType, Eigen::Lower|Eigen::Upper> solver;
            solver.compute(systemMatrix);
            solution = solver.solve( rhs );
            //std::cout << "#iterations:     " << cg.iterations() << std::endl;
            //std::cout << "estimated error: " << cg.error()      << std::endl;
            break;
      }
      
      case BiCGSTAB : {
            Eigen::BiCGSTAB<SparseMatrixType> solver;
            solver.compute(systemMatrix);
            solution = solver.solve( rhs );
      }
      
      default:
          throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
          break;
    }
}





template<typename VectorType, typename SparseMatrixType>
void solveLinearSystem ( const SparseMatrixType & systemMatrix , 
                         VectorType &solution, 
                         const VectorType &rhs,
                         const LinearSolverMethod linearSolver = UmfPackLU ) {
                             
    switch( linearSolver ){

//       case UmfPackLU:{
//           Eigen::UmfPackLU<SparseMatrixType> solver;
//           solver.compute(systemMatrix);
//           if(solver.info()!=Eigen::Success) {
//               throw std::invalid_argument ( aol::strprintf ( "solver failed in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
//               return;
//           }
//           solution = solver.solve( rhs );
//           break;
//       }
      
      case CholmodSupernodalLLT:{
//           Eigen::CholmodSupernodalLLT<SparseMatrixType> solver;
          Eigen::CholmodDecomposition<SparseMatrixType> solver;
          solver.compute(systemMatrix);
          if(solver.info()!=Eigen::Success) {
              throw std::invalid_argument ( aol::strprintf ( "solver failed in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
              return;
        }
        solution = solver.solve( rhs );         
        break;   
      }
      
      case ConjugateGradient: {
            Eigen::ConjugateGradient<SparseMatrixType, Eigen::Lower|Eigen::Upper> solver;
            solver.compute(systemMatrix);
            solution = solver.solve( rhs );
            break;
      }
      
      case BiCGSTAB: {
            Eigen::BiCGSTAB<SparseMatrixType> solver;
            solver.compute(systemMatrix);
            solution = solver.solve( rhs );
            break;
      }
      
      default:
          throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
          break;
    }
}



enum NonlinearSolverMethod{ GRADIENTDESCENT, QUASINEWTON, NEWTON }; //TODO not in EigenIncludes

} //end namespace



//Output
template<typename VectorType>
void printVector( const VectorType &vec, const int breakNumber = 5, const std::string name = "", const unsigned precision = 5 ){
    /*Eigen::IOFormat::IOFormat ( 
    int _precision = StreamPrecision, int  _flags = 0, 
    const std::string & _coeffSeparator = " ", 
    const std::string & _rowSeparator = "\n",  const std::string & _rowPrefix = "", const std::string &_rowSuffix = "",
    const std::string &  _matPrefix = "", const std::string & _matSuffix = "" )
    */
    
//     Eigen::IOFormat Fmt( FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]"  );
//     std::cout << m1.format(Fmt) << endl;

  cout << "------------------------------------------------------------------------------------------" << endl;
  cout << "-----------" << name.c_str() << "-----------------------" << endl;

  int iter = 0;
  cout.precision(precision);
  for ( int i = 0; i < vec.size(); ++i ){
      cout << vec[i] << " ; ";
      if( iter == breakNumber - 1 ){
        cout << endl;
        iter = 0;
      }else iter++;
  }
  cout << "------------------------------------------------------------------------------------------" << endl;

}



#endif // __EIGENINCLUDE_H
