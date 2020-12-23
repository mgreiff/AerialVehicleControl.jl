/*******************************************************************************
* Copyright (C) Marcus Greiff 2020
*
* @file cont_matrix_math.h
* @author Marcus Greiff
* @date June 2020
*******************************************************************************/
#include "cont_matrix_math.h"

#if defined(LINK_LAPACK)
/* load dgemm as an external function from LAPACK */
extern void dgemm_(
  char * transa,
  char * transb,
  int * m,
  int * n,
  int * k,
  double * alpha,
  double * A,
  int * lda,
  double * B,
  int * ldb,
  double * beta,
  double * C,
  int * ldc
);

/* load dpotrf as an external function from LAPACK */
extern void dpotrf_(
  char * uplo,
  int * n,
  double * A,
  int * lda,
  int * info
);

/* load dposv as an external function from LAPACK */
extern void dposv_(
  char* uplo,
  int* n,
  int* nrhs,
  double* a,
  int* lda,
  double* b,
  int* ldb,
  int* info
);
#endif

int matrix_double_addition(
  matrix_double_t *Amat,
  matrix_double_t *Bmat,
  matrix_double_t *Cmat
  ) {
  int i;
  /* Input data check */
  if ((Amat->numRows != Bmat->numRows) ||
      (Bmat->numRows != Cmat->numRows) ||
      (Amat->numCols != Bmat->numCols) ||
      (Bmat->numCols != Cmat->numCols)){
    return 0;
  }
  /* Perform matrix addition */
  for (i = 0; i < Amat->numRows*Amat->numCols; i++) Cmat->pData[i] = Amat->pData[i] + Bmat->pData[i];
  return 1;
}

int matrix_double_addition_inplace(
  matrix_double_t *Amat,
  matrix_double_t *Bmat
  ) {
  int i;
  /* Input data check */
  if ((Amat->numRows != Bmat->numRows) ||
      (Amat->numCols != Bmat->numCols)){
    return 0;
  }
  /* Perform matrix addition */
  for (i = 0; i < Amat->numRows*Amat->numCols; i++) Amat->pData[i] += Bmat->pData[i];
  return 1;
}

int matrix_double_subtraction(
  matrix_double_t *Amat,
  matrix_double_t *Bmat,
  matrix_double_t *Cmat
  ) {
  int i;
  /* Input data check */
  if ((Amat->numRows != Bmat->numRows) ||
      (Bmat->numRows != Cmat->numRows) ||
      (Amat->numCols != Bmat->numCols) ||
      (Bmat->numCols != Cmat->numCols)){
    return 0;
  }
  /* Perform matrix addition */
  for (i = 0; i < Amat->numRows*Amat->numCols; i++) Cmat->pData[i] = Amat->pData[i] - Bmat->pData[i];
  return 1;
}

int matrix_double_subtraction_inplace(
  matrix_double_t *Amat,
  matrix_double_t *Bmat
  ) {
  int i;
  /* Input data check */
  if ((Amat->numRows != Bmat->numRows) ||
      (Amat->numCols != Bmat->numCols)){
    return 0;
  }
  /* Perform matrix addition */
  for (i = 0; i < Amat->numRows*Amat->numCols; i++) Amat->pData[i] -= Bmat->pData[i];
  return 1;
}

int matrix_double_multiplication(
  matrix_double_t *Amat,
  matrix_double_t *Bmat,
  matrix_double_t *Cmat
){
#if defined(LINK_LAPACK)
  double alpha = 1.0, beta = 0.0;
#else
  int i, j, k;
#endif

  if ((Amat->numRows != Cmat->numRows) ||
      (Bmat->numCols != Cmat->numCols) ||
      (Amat->numCols != Bmat->numRows)){
    return 0;
  }
#if defined(LINK_LAPACK)
  dgemm_(
    "N",
    "N",
    &Amat->numRows,
    &Bmat->numCols,
    &Amat->numCols,
    &alpha,
    Amat->pData,
    &Amat->numRows,
    Bmat->pData,
    &Bmat->numRows,
    &beta,
    Cmat->pData,
    &Cmat->numRows
  );
#else
  matrix_zero(Cmat);
  for (i = 0; i < Amat->numRows; i++ ){
    for (k = 0; k < Amat->numCols; k++ ){
      for (j = 0; j < Bmat->numCols; j++ ){
        matrix_set(Cmat, i, j, matrix_get(Amat, i, k) * matrix_get(Bmat, k, j) + matrix_get(Cmat, i, j));
      }
    }
  }
#endif
  return 1;
}

int matrix_double_transposition(
  matrix_double_t *Amat,
  matrix_double_t *ATmat
) {
  int i, j;
  if ((Amat->numRows != ATmat->numCols) ||
      (Amat->numCols != ATmat->numRows)){
    return 0;
  }
  for (i = 0; i < Amat->numRows; i++ ){
    for (j = 0; j < Amat->numCols; j++ ){
      ATmat->pData[i*Amat->numCols + j] = Amat->pData[j*Amat->numRows + i];
    }
  }
  return 1;
}

int matrix_double_cholesky_inplace(
  matrix_double_t *Amat
) {
#if defined(LINK_LAPACK)
  int i, j;
#endif
  int info = 0;
  if ( Amat->numRows != Amat->numCols ) return 0;

#if defined(LINK_LAPACK)
  dpotrf_(
    "L",
    &Amat->numRows,
    Amat->pData,
    &Amat->numRows,
    &info
  );
  for (i = 0; i < Amat->numRows; i++ ){
    for (j = i + 1; j < Amat->numRows; j++ ) Amat->pData[i + j * Amat->numRows] = 0.0;
  }
#else
  TRACE(5,("The in-place cholesky decomposition has not yet been implemented dependency free", info));
#endif
  if (0 == info){
    return 1;
  } else {
    if (0 < info) TRACE(5,("The %i-th argument in dpotrf had an illegal value\n", -info));
    if (0 > info) TRACE(5,("The leading minor of order %i is not positive definite\n", info));
    return 0;
  }
}

int matrix_double_cholesky(
  matrix_double_t *Amat,
  matrix_double_t *Lmat
) {
  int i, status;
  if ((Amat->numRows != Lmat->numRows) ||
      (Amat->numCols != Lmat->numCols)){
    return 0;
  }
  for (i = 0; i < Amat->numRows * Amat->numCols; i++) Lmat->pData[i] = Amat->pData[i];
  status = matrix_double_cholesky_inplace(Lmat);
  return status;
}

int matrix_double_solve_posdef(
  matrix_double_t *Amat,
  matrix_double_t *Bmat
) {
  int info = 0;
  if ((Amat->numRows != Amat->numCols) ||
      (Amat->numRows != Bmat->numRows) ||
      (Bmat->numCols < 1)){
    return 0;
  }

#if defined(LINK_LAPACK)
  dposv_(
    "L",
    &Amat->numCols,
    &Bmat->numCols,
    Amat->pData,
    &Amat->numRows,
    Bmat->pData,
    &Bmat->numRows,
    &info
  );
#else
  TRACE(5,("The PSD linear system solver has not het been implemented dependency free", info));
#endif

  if (0 == info){
    return 1;
  } else {
    if (0 < info) TRACE(5,("The %i-th argument in dposv had an illegal value\n", -info));
    if (0 > info) TRACE(5,("The leading minor of order %i is not positive definite\n", info));
    return 0;
  }
}

void mat_add(matrix_double_t *Amat, matrix_double_t *Bmat, matrix_double_t *Cmat)
{ assert(1 == matrix_double_addition(Amat, Bmat, Cmat)); }

void mat_add_inplace(matrix_double_t *Amat, matrix_double_t *Bmat)
{ assert(1 == matrix_double_addition_inplace(Amat, Bmat)); }

void mat_sub(matrix_double_t *Amat, matrix_double_t *Bmat, matrix_double_t *Cmat)
{ assert(1 == matrix_double_subtraction(Amat, Bmat, Cmat)); }

void mat_sub_inplace(matrix_double_t *Amat, matrix_double_t *Bmat)
{ assert(1 == matrix_double_subtraction_inplace(Amat, Bmat)); }

void mat_mul(matrix_double_t *Amat, matrix_double_t *Bmat, matrix_double_t *Cmat)
{ assert(1 == matrix_double_multiplication(Amat, Bmat, Cmat)); }

void mat_trans(matrix_double_t *Amat, matrix_double_t *ATmat)
{ assert(1 == matrix_double_transposition(Amat, ATmat)); }

void mat_chol(matrix_double_t *Amat, matrix_double_t *Lmat)
{ assert(1 == matrix_double_cholesky(Amat, Lmat)); }

void mat_chol_inplace(matrix_double_t *Amat)
{ assert(1 == matrix_double_cholesky_inplace(Amat)); }

void mat_sol(matrix_double_t *Amat, matrix_double_t *Bmat)
{ assert(1 == matrix_double_solve_posdef(Amat, Bmat)); }
