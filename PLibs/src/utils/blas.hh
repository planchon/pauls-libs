/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  This program is free software; you can redistribute it and/or modify    |
 |  it under the terms of the GNU General Public License as published by    |
 |  the Free Software Foundation; either version 2, or (at your option)     |
 |  any later version.                                                      |
 |                                                                          |
 |  This program is distributed in the hope that it will be useful,         |
 |  but WITHOUT ANY WARRANTY; without even the implied warranty of          |
 |  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           |
 |  GNU General Public License for more details.                            |
 |                                                                          |
 |  You should have received a copy of the GNU General Public License       |
 |  along with this program; if not, write to the Free Software             |
 |  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               |
 |                                                                          |
 |  Copyright (C) 2008                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                | 
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Meccanica e Strutturale                  |
 |      Universita` degli Studi di Trento                                   |
 |      Via Mesiano 77, I-38050 Trento, Italy                               |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: Blas.hh
///

#ifndef BLAS_HH
#define BLAS_HH

#ifndef F77NAME
  #define F77NAME(A) A##_
#endif

#ifndef ATLASNAME
  #define ATLASNAME(A) ATL_##A
#endif

namespace blas {

  typedef char   character ;
  typedef int    integer ;
  typedef float  single_precision ;
  typedef double double_precision ;

  #ifdef USE_ATLAS

  /*
  //     #    ####### #          #     #####  
  //    # #      #    #         # #   #     # 
  //   #   #     #    #        #   #  #       
  //  #     #    #    #       #     #  #####  
  //  #######    #    #       #######       # 
  //  #     #    #    #       #     # #     # 
  //  #     #    #    ####### #     #  #####
  */

  enum ATLAS_ORDER {AtlasRowMajor=101, AtlasColMajor=102 };
  enum ATLAS_TRANS {AtlasNoTrans=111, AtlasTrans=112,
                    AtlasConjTrans=113, AtlasConj=114};
  enum ATLAS_UPLO  {AtlasUpper=121, AtlasLower=122};
  enum ATLAS_DIAG  {AtlasNonUnit=131, AtlasUnit=132};

  extern "C" {
    void
    ATLASNAME(scopy)( const int N,
                      const float *X, const int incX,
                      float       *Y, const int incY );
    void
    ATLASNAME(dcopy)( const int N,
                      const double *X, const int incX,
                      double       *Y, const int incY );
    void
    ATLASNAME(sswap)( const int N,
                      float *X, const int incX,
                      float *Y, const int incY );
    void
    ATLASNAME(dswap)( const int N,
                      double *X, const int incX,
                      double *Y, const int incY );
    void
    ATLASNAME(sscal)( const int N,
                      const float alpha,
                      float *X, const int incX );
    void
    ATLASNAME(dscal)( const int N,
                      const double alpha,
                      double *X, const int incX );
    void
    ATLASNAME(saxpy)( const int N,
                      const float alpha,
                      const float *X, const int incX,
                      float       *Y, const int incY );
    void
    ATLASNAME(daxpy)( const int N,
                      const double alpha,
                      const double *X, const int incX,
                      double       *Y, const int incY );
    void
    ATLASNAME(sger)( const int M, const int N,
                     const float alpha,
                     const float *X, const int incX,
                     const float *Y, const int incY,
                     float       *A, const int lda );
    void
    ATLASNAME(dger)( const int M, const int N,
                     const double alpha,
                     const double *X, const int incX,
                     const double *Y, const int incY,
                     double       *A, const int lda );
    void
    ATLASNAME(sgemv)( const enum ATLAS_TRANS TransA,
                      const int M, const int N,
                      const float alpha,
                      const float *A, const int lda,
                      const float *X, const int incX,
                      const float beta,
                      float *Y, const int incY );
    void
    ATLASNAME(dgemv)( const enum ATLAS_TRANS TransA,
                      const int M, const int N,
                      const double alpha,
                      const double *A, const int lda,
                      const double *X, const int incX,
                      const double beta,
                      double *Y, const int incY );
    void
    ATLASNAME(sgemm)( const enum ATLAS_TRANS TransA,
                      const enum ATLAS_TRANS TransB,
                      const int M, const int N, const int K,
                      const float alpha,
                      const float *A, const int lda,
                      const float *B, const int ldb,
                      const float beta,
                      float *C, const int ldc );
    void
    ATLASNAME(dgemm)( const enum ATLAS_TRANS TransA,
                      const enum ATLAS_TRANS TransB,
                      const int M, const int N, const int K,
                      const double alpha,
                      const double *A, const int lda,
                      const double *B, const int ldb,
                      const double beta,
                      double *C, const int ldc );
    void
    ATLASNAME(szero)( const int N, float *X, const int incX );
    void
    ATLASNAME(dzero)( const int N, double *X, const int incX );
    float
    ATLASNAME(snrm2)( const int N, const float *X, const int incX );
    double
    ATLASNAME(dnrm2)( const int N, const double *X, const int incX );
    float
    ATLASNAME(sasum)( const int N, const float *X, const int incX );
    double
    ATLASNAME(dasum)( const int N, const double *X, const int incX );
    int
    ATLASNAME(isamax)( const int N, const float *X, const int incX );
    int
    ATLASNAME(idamax)( const int N, const double *X, const int incX );
    float
    ATLASNAME(sdot)( const int N,
                     const float *X, const int incX,
                     const float *Y, const int incY );
    double
    ATLASNAME(ddot)( const int N,
                     const double *X, const int incX,
                     const double *Y, const int incY );
  }

  inline
  ATLAS_TRANS
  trans_atlas(character const TRANS ) {
    switch ( TRANS ) {
    case 'N': return AtlasNoTrans ;
    case 'T': return AtlasTrans ;
    case 'C': return AtlasConjTrans ;
    default:  throw runtime_error("bad traspose mode selection\n") ;
    }
  }

  inline
  ATLAS_UPLO
  uplo_atlas(character const UPLO ) {
    switch ( UPLO ) {
    case 'U': return AtlasUpper ;
    case 'L': return AtlasLower ;
    default:  throw runtime_error("bad UP/LOW mode selection\n") ;
    }
  }

  inline
  ATLAS_DIAG
  diag_atlas(character const DIAG ) {
    switch ( DIAG ) {
    case 'U': return AtlasUnit ;
    case 'N': return AtlasNonUnit ;
    default:  throw runtime_error("bad diagonal mode selection\n") ;
    }
  }

  #else

  /*
  //  #####  #        ##    ####  
  //  #    # #       #  #  #      
  //  #####  #      #    #  ####  
  //  #    # #      ######      # 
  //  #    # #      #    # #    # 
  //  #####  ###### #    #  ####  
  */

  extern "C" {
    void
    F77NAME(scopy)( integer          const * N, 
                    single_precision const   X[], 
                    integer          const * INCX, 
                    single_precision         Y[], 
                    integer          const * INCY ) ;
    void
    F77NAME(dcopy)( integer          const * N, 
                    double_precision const   X[], 
                    integer          const * INCX, 
                    double_precision         Y[], 
                    integer          const * INCY ) ;
    void
    F77NAME(sswap)( integer          const * N, 
                    single_precision         X[], 
                    integer          const * INCX, 
                    single_precision         Y[], 
                    integer          const * INCY ) ;
    void
    F77NAME(dswap)( integer          const * N, 
                    double_precision         X[], 
                    integer          const * INCX, 
                    double_precision         Y[], 
                    integer          const * INCY ) ;
    void
    F77NAME(sscal)( integer          const * N,
                    single_precision const * S,
                    single_precision         X[],
                    integer          const * INCX ) ;
    void
    F77NAME(dscal)( integer          const * N,
                    double_precision const * S,
                    double_precision         X[],
                    integer          const * INCX ) ;
    void
    F77NAME(saxpy)( integer          const * N, 
                    single_precision const * A, 
                    single_precision const   X[], 
                    integer          const * INCX, 
                    single_precision         Y[], 
                    integer          const * INCY ) ;
    void
    F77NAME(daxpy)( integer          const * N, 
                    double_precision const * A, 
                    double_precision const   X[], 
                    integer          const * INCX, 
                    double_precision         Y[], 
                    integer          const * INCY ) ;
    void
    F77NAME(strsv)( character        const   UPLO[], 
                    character        const   TRANS[], 
                    character        const   DIAG[], 
                    integer          const * N, 
                    single_precision const   A[], 
                    integer          const * LDA, 
                    single_precision         X[], 
                    integer          const * INCX ) ;
    void
    F77NAME(dtrsv)( character        const   UPLO[], 
                    character        const   TRANS[], 
                    character        const   DIAG[], 
                    integer          const * N, 
                    double_precision const   A[], 
                    integer          const * LDA, 
                    double_precision         X[], 
                    integer          const * INCX ) ;
    void
    F77NAME(sger)( integer          const * M, 
                   integer          const * N, 
                   single_precision const * ALPHA, 
                   single_precision const   X[], 
                   integer          const * INCX, 
                   single_precision const   Y[], 
                   integer          const * INCY, 
                   single_precision         A[], 
                   integer          const * LDA ) ;
    void
    F77NAME(dger)( integer          const * M, 
                   integer          const * N, 
                   double_precision const * ALPHA, 
                   double_precision const   X[], 
                   integer          const * INCX, 
                   double_precision const   Y[], 
                   integer          const * INCY, 
                   double_precision         A[], 
                   integer          const * LDA ) ;
    void
    F77NAME(sgemv)( character        const   TRANS[], 
                    integer          const * M, 
                    integer          const * N, 
                    single_precision const * ALPHA, 
                    single_precision const   A[], 
                    integer          const * LDA,
                    single_precision const   X[], 
                    integer          const * INCX, 
                    single_precision const * BETA, 
                    single_precision         Y[], 
                    integer          const * INCY ) ;
    void
    F77NAME(dgemv)( character        const   TRANS[], 
                    integer          const * M, 
                    integer          const * N, 
                    double_precision const * ALPHA, 
                    double_precision const   A[], 
                    integer          const * LDA,
                    double_precision const   X[], 
                    integer          const * INCX, 
                    double_precision const * BETA, 
                    double_precision         Y[], 
                    integer          const * INCY ) ;
    void
    F77NAME(sgemm)( character        const   TRANSA[], 
                    character        const   TRANSB[], 
                    integer          const * M, 
                    integer          const * N, 
                    integer          const * K, 
                    single_precision const * ALPHA, 
                    single_precision const   A[], 
                    integer          const * LDA,
                    single_precision const   B[], 
                    integer          const * LDB,
                    single_precision const * BETA, 
                    single_precision         C[], 
                    integer          const * LDC ) ;
    void
    F77NAME(dgemm)( character        const   TRANSA[], 
                    character        const   TRANSB[], 
                    integer          const * M, 
                    integer          const * N, 
                    integer          const * K, 
                    double_precision const * ALPHA, 
                    double_precision const   A[], 
                    integer          const * LDA,
                    double_precision const   B[], 
                    integer          const * LDB,
                    double_precision const * BETA, 
                    double_precision         C[], 
                    integer          const * LDC ) ;
    single_precision
    F77NAME(snrm2)( integer          const * N,
                    single_precision const   X[],
                    integer          const * INCX ) ;

    double_precision
    F77NAME(dnrm2)( integer          const * N,
                    double_precision const   X[],
                    integer          const * INCX ) ;
    single_precision
    F77NAME(sasum)( integer          const * N,
                    single_precision const   X[],
                    integer          const * INCX ) ;
    double_precision
    F77NAME(dasum)( integer          const * N,
                    double_precision const   X[],
                    integer          const * INCX ) ;
    integer
    F77NAME(isamax)( integer          const * N,
                     single_precision const   X[],
                     integer          const * INCX ) ;
    integer
    F77NAME(idamax)( integer          const * N,
                     double_precision const   X[],
                     integer          const * INCX ) ;
    single_precision
    F77NAME(sdot)( integer          const * N, 
                   single_precision const   SX[], 
                   integer          const * INCX, 
                   single_precision const   SY[], 
                   integer          const * INCY ) ;

    double_precision
    F77NAME(ddot)( integer          const * N, 
                   double_precision const   SX[], 
                   integer          const * INCX, 
                   double_precision const   SY[], 
                   integer          const * INCY ) ;
  }

  #endif

  /* ____ ____ ___  _   _ 
  // |    |  | |__]  \_/  
  // |___ |__| |      |   
  */
  inline 
  void
  copy( integer          const N, 
        single_precision const X[], 
        integer          const INCX, 
        single_precision       Y[], 
        integer          const INCY )
  #ifdef USE_ATLAS
  { ATLASNAME(scopy)( N, X, INCX, Y, INCY ) ; }
  #else
  { F77NAME(scopy)( &N, X, &INCX, Y, &INCY ) ; }
  #endif

  inline
  void
  copy( integer          const N, 
        double_precision const X[], 
        integer          const INCX, 
        double_precision       Y[], 
        integer          const INCY )
  #ifdef USE_ATLAS
  { ATLASNAME(dcopy)( N, X, INCX, Y, INCY ) ; }
  #else
  { F77NAME(dcopy)( &N, X, &INCX, Y, &INCY ) ; }
  #endif

  /*   __            _  
  //  (_ \    / /\  |_) 
  //  __) \/\/ /--\ |   
  */
  inline
  void
  swap( integer          const N, 
        single_precision       X[], 
        integer          const INCX, 
        single_precision       Y[], 
        integer          const INCY )
  #ifdef USE_ATLAS
  { ATLASNAME(sswap)( N, X, INCX, Y, INCY ) ; }
  #else
  { F77NAME(sswap)( &N, X, &INCX, Y, &INCY ) ; }
  #endif

  inline
  void
  swap( integer          const N, 
        double_precision       X[], 
        integer          const INCX, 
        double_precision       Y[], 
        integer          const INCY )
  #ifdef USE_ATLAS
  { ATLASNAME(dswap)( N, X, INCX, Y, INCY ) ; }
  #else
  { F77NAME(dswap)( &N, X, &INCX, Y, &INCY ) ; }
  #endif

  /*  ____ ____ ____ _    
  //  [__  |    |__| |    
  //  ___] |___ |  | |___ 
  */
  inline
  void
  scal( integer          const N,
        single_precision const S,
        single_precision       X[],
        integer          const INCX )
  #ifdef USE_ATLAS
  { ATLASNAME(sscal)(N, S, X, INCX) ; }
  #else
  { F77NAME(sscal)(&N, &S, X, &INCX) ; }
  #endif

  inline
  void
  scal( integer          const N,
        double_precision const S,
        double_precision       X[],
        integer          const INCX )
  #ifdef USE_ATLAS
  { ATLASNAME(dscal)(N, S, X, INCX) ; }
  #else
  { F77NAME(dscal)(&N, &S, X, &INCX) ; }
  #endif

  /*  ____ _  _ ___  _   _ 
  //  |__|  \/  |__]  \_/  
  //  |  | _/\_ |      |   
  */
  inline
  void
  axpy( integer          const N, 
        single_precision const A, 
        single_precision const X[], 
        integer          const INCX, 
        single_precision       Y[], 
        integer          const INCY )
  #ifdef USE_ATLAS
  { ATLASNAME(saxpy)( N, A, X, INCX, Y, INCY ) ; }
  #else
  { F77NAME(saxpy)( &N, &A, X, &INCX, Y, &INCY ) ; }
  #endif

  inline
  void
  axpy( integer          const N, 
        double_precision const A, 
        double_precision const X[], 
        integer          const INCX, 
        double_precision       Y[], 
        integer          const INCY )
  #ifdef USE_ATLAS
  { ATLASNAME(daxpy)( N, A, X, INCX, Y, INCY ) ; }
  #else
  { F77NAME(daxpy)( &N, &A, X, &INCX, Y, &INCY ) ; }
  #endif

  /*  __  _  _   _  
  //   / |_ |_) / \ 
  //  /_ |_ | \ \_/
  */
  inline
  void
  zero( integer          const N,
        single_precision       X[],
        integer          const INCX )
  #ifdef USE_ATLAS
  { ATLASNAME(szero)( N, X, INCX ) ; }
  #else
  { single_precision z  = 0 ;
    integer          iz = 0 ;
    F77NAME(scopy)( &N, &z, &iz, X, &INCX ) ; }
  #endif

  inline
  void
  zero( integer          const N,
        double_precision       X[],
        integer          const INCX )
  #ifdef USE_ATLAS
  { ATLASNAME(dzero)( N, X, INCX ) ; }
  #else
  { double_precision z  = 0 ;
    integer          iz = 0 ;
    F77NAME(dcopy)( &N, &z, &iz, X, &INCX ) ; }
  #endif

  /*   __  _  _  
  //  /__ |_ |_) 
  //  \_| |_ | \ 
  */
  inline
  void
  ger( integer          const M, 
       integer          const N, 
       single_precision const ALPHA, 
       single_precision const X[], 
       integer          const INCX, 
       single_precision const Y[], 
       integer          const INCY, 
       single_precision       A[], 
       integer          const LDA )
  #ifdef USE_ATLAS
  { ATLASNAME(sger)( M, N, ALPHA, X, INCX, Y, INCY, A, LDA ) ; }
  #else
  { F77NAME(sger)( &M, &N, &ALPHA, X, &INCX, Y, &INCY, A, &LDA ) ; }
  #endif

  inline
  void
  ger( integer          const M, 
       integer          const N, 
       double_precision const ALPHA, 
       double_precision const X[], 
       integer          const INCX, 
       double_precision const Y[], 
       integer          const INCY, 
       double_precision       A[], 
       integer          const LDA )
  #ifdef USE_ATLAS
  { ATLASNAME(dger)( M, N, ALPHA, X, INCX, Y, INCY, A, LDA ) ; }
  #else
  { F77NAME(dger)( &M, &N, &ALPHA, X, &INCX, Y, &INCY, A, &LDA ) ; }
  #endif

  /*   __  _           
  //  /__ |_ |\/| \  / 
  //  \_| |_ |  |  \/  
  */
  inline
  void
  gemv( character        const TRANS[], 
        integer          const M, 
        integer          const N, 
        single_precision const ALPHA, 
        single_precision const A[], 
        integer          const LDA,
        single_precision const X[], 
        integer          const INCX, 
        single_precision const BETA, 
        single_precision       Y[], 
        integer          const INCY )
  #ifdef USE_ATLAS
  { ATLASNAME(sgemv)( trans_atlas(TRANS[0]), M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY ) ; }
  #else
  { F77NAME(sgemv)( TRANS, &M, &N, &ALPHA, A, &LDA, X, &INCX, &BETA, Y, &INCY ) ; }
  #endif

  inline
  void
  gemv( character        const TRANS[], 
        integer          const M, 
        integer          const N, 
        double_precision const ALPHA, 
        double_precision const A[], 
        integer          const LDA,
        double_precision const X[], 
        integer          const INCX, 
        double_precision const BETA, 
        double_precision       Y[], 
        integer          const INCY )
  #ifdef USE_ATLAS
  { ATLASNAME(dgemv)( trans_atlas(TRANS[0]), M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY ) ; }
  #else
  { F77NAME(dgemv)( TRANS, &M, &N, &ALPHA, A, &LDA, X, &INCX, &BETA, Y, &INCY ) ; }
  #endif

  /*   ___ _   __     
  //    | |_) (_ \  / 
  //    | | \ __) \/  
  */
  inline
  void
  trsv( character        const UPLO[], 
        character        const TRANS[], 
        character        const DIAG[],
        integer          const N, 
        single_precision const A[], 
        integer          const LDA, 
        single_precision       X[], 
        integer          const INCX )
  #ifdef USE_ATLAS
  { ATLASNAME(strsv)( uplo_atlas(UPLO[0]),
                      trans_atlas(TRANS[0]),
                      diag_atlas(DIAG[0]),
                      N, A, LDA, X, INCX ) ; }
  #else
  { F77NAME(strsv)( UPLO, TRANS, DIAG, &N, A, &LDA, X, &INCX ) ; }
  #endif

  inline
  void
  trsv( character        const UPLO[],
        character        const TRANS[],
        character        const DIAG[],
        integer          const N, 
        double_precision const A[], 
        integer          const LDA, 
        double_precision       X[], 
        integer          const INCX )
  #ifdef USE_ATLAS
  { ATLASNAME(dtrsv)( uplo_atlas(UPLO[0]),
                      trans_atlas(TRANS[0]),
                      diag_atlas(DIAG[0]),
                      N, A, LDA, X, INCX ) ; }
  #else
  { F77NAME(dtrsv)( UPLO, TRANS, DIAG, &N, A, &LDA, X, &INCX ) ; }
  #endif

  /*   __  _           
  //  /__ |_ |\/| |\/| 
  //  \_| |_ |  | |  | 
  */
  inline
  void
  gemm( character        const TRANSA[], 
        character        const TRANSB[], 
        integer          const M, 
        integer          const N, 
        integer          const K, 
        single_precision const ALPHA, 
        single_precision const A[], 
        integer          const LDA,
        single_precision const B[], 
        integer          const LDB,
        single_precision const BETA, 
        single_precision       C[], 
        integer          const LDC )
  #ifdef USE_ATLAS
  { ATLASNAME(sgemm)( trans_atlas(TRANSA[0]), trans_atlas(TRANSB[0]),
                      M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC ) ; }
  #else
  { F77NAME(sgemm)( TRANSA, TRANSB,
                    &M, &N, &K,
                    &ALPHA, A, &LDA,
                    B, &LDB,
                    &BETA, C, &LDC ) ; }
  #endif

  inline
  void
  gemm( character        const TRANSA[],
        character        const TRANSB[], 
        integer          const M, 
        integer          const N, 
        integer          const K, 
        double_precision const ALPHA,
        double_precision const A[], 
        integer          const LDA,
        double_precision const B[], 
        integer          const LDB,
        double_precision const BETA, 
        double_precision       C[], 
        integer          const LDC )
  #ifdef USE_ATLAS
  { ATLASNAME(dgemm)( trans_atlas(TRANSA[0]), trans_atlas(TRANSB[0]),
                      M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC ) ; }
  #else
  { F77NAME(dgemm)( TRANSA, TRANSB,
                    &M, &N, &K,
                    &ALPHA, A, &LDA,
                    B, &LDB,
                    &BETA, C, &LDC ) ; }
  #endif

  /*        _       _  
  //  |\ | |_) |\/|  ) 
  //  | \| | \ |  | /_ 
  */
  inline
  single_precision
  nrm2( integer          const N,
        single_precision const X[],
        integer          const INCX )
  #ifdef USE_ATLAS
  { return ATLASNAME(snrm2)( N, X, INCX ) ; }
  #else
  { return F77NAME(snrm2)( &N, X, &INCX ) ; }
  #endif

  inline
  double_precision
  nrm2( integer          const N,
        double_precision const X[],
        integer          const INCX )
  #ifdef USE_ATLAS
  { return ATLASNAME(dnrm2)( N, X, INCX ) ; }
  #else
  { return F77NAME(dnrm2)( &N, X, &INCX ) ; }
  #endif

  /*       __        
  //   /\ (_ | ||\/| 
  //  /--\__)|_||  | 
  */
  inline
  single_precision
  asum(integer          const N,
       single_precision const X[],
       integer          const INCX)
  #ifdef USE_ATLAS
  { return ATLASNAME(sasum)( N, X, INCX ) ; }
  #else
  { return F77NAME(sasum)( &N, X, &INCX ) ; }
  #endif
  
  inline
  double_precision
  asum(integer          const N,
       double_precision const X[],
       integer          const INCX)
  #ifdef USE_ATLAS
  { return ATLASNAME(dasum)( N, X, INCX ) ; }
  #else
  { return F77NAME(dasum)( &N, X, &INCX ) ; }
  #endif

  /*        
  //    /\  |\/|  /\  \/ 
  //   /--\ |  | /--\ /\ 
  */
  inline
  single_precision
  amax(integer          const N,
       single_precision const X[],
       integer          const INCX)
  #ifdef USE_ATLAS
  { return X[ATLASNAME(isamax)( N, X, INCX )-1] ; }
  #else
  { return X[F77NAME(isamax)( &N, X, &INCX )-1] ; }
  #endif
  
  inline
  double_precision
  amax(integer          const N,
       double_precision const X[],
       integer          const INCX)
  #ifdef USE_ATLAS
  { return X[ATLASNAME(idamax)( N, X, INCX )-1] ; }
  #else
  { return X[F77NAME(idamax)( &N, X, &INCX )-1] ; }
  #endif

  /*   _   _ ___ 
  //  | \ / \ |  
  //  |_/ \_/ |  
  */
  inline
  single_precision
  dot(integer          const N, 
      single_precision const SX[], 
      integer          const INCX, 
      single_precision const SY[], 
      integer          const INCY)
  #ifdef USE_ATLAS
  { return ATLASNAME(sdot)( N, SX, INCX, SY, INCY ) ; }
  #else
  { return F77NAME(sdot)( &N, SX, &INCX, SY, &INCY ) ; }
  #endif
  
  inline
  double_precision
  dot(integer          const N, 
      double_precision const SX[], 
      integer          const INCX, 
      double_precision const SY[], 
      integer          const INCY) 
  #ifdef USE_ATLAS
  { return ATLASNAME(ddot)( N, SX, INCX, SY, INCY ) ; }
  #else
  { return F77NAME(ddot)( &N, SX, &INCX, SY, &INCY ) ; }
  #endif

} // end namespace blas

#endif

///
/// eof: Blas.hh
///

