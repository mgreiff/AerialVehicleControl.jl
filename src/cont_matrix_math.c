/*******************************************************************************
* Copyright (C) Marcus Greiff 2020
*
* @file cont_matrix_math.h
* @author Marcus Greiff
* @date June 2020
*******************************************************************************/
#include "cont_matrix_math.h"

#if defined(LINK_LAPACK)

/* Load the necessary functions from LAPACK */

extern void dgemm_(char * transa, char * transb, int * m, int * n, int * k,
                   double * alpha, double * A, int * lda, double * B, int * ldb,
                   double * beta, double * C, int * ldc);

extern void dpotrf_(char * uplo, int * n, double * A, int * lda, int * info);

extern void dposv_(char* uplo, int* n, int* nrhs, double* a, int* lda,
                   double* b, int* ldb, int* info);

extern void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda,
                   double* w, double* work, int* lwork, int* info);
#endif

int matrix_double_addition(
  matrix_double_t *Amat,
  matrix_double_t *Bmat,
  matrix_double_t *Cmat
){
  int i;
  /* Input data check */
  if ((Amat->numRows != Bmat->numRows) ||
      (Bmat->numRows != Cmat->numRows) ||
      (Amat->numCols != Bmat->numCols) ||
      (Bmat->numCols != Cmat->numCols)) return 0;
  /* Perform matrix addition */
  for (i = 0; i < Amat->numRows*Amat->numCols; i++) Cmat->pData[i] = Amat->pData[i] + Bmat->pData[i];
  return 1;
}

int matrix_double_addition_inplace(
  matrix_double_t *Amat,
  matrix_double_t *Bmat
){
  int i;
  /* Input data check */
  if ((Amat->numRows != Bmat->numRows) ||
      (Amat->numCols != Bmat->numCols)) return 0;
  /* Perform matrix addition */
  for (i = 0; i < Amat->numRows*Amat->numCols; i++) Amat->pData[i] += Bmat->pData[i];
  return 1;
}

int matrix_double_subtraction(
  matrix_double_t *Amat,
  matrix_double_t *Bmat,
  matrix_double_t *Cmat
){
  int i;
  /* Input data check */
  if ((Amat->numRows != Bmat->numRows) ||
      (Bmat->numRows != Cmat->numRows) ||
      (Amat->numCols != Bmat->numCols) ||
      (Bmat->numCols != Cmat->numCols)) return 0;
  /* Perform matrix addition */
  for (i = 0; i < Amat->numRows*Amat->numCols; i++) Cmat->pData[i] = Amat->pData[i] - Bmat->pData[i];
  return 1;
}

int matrix_double_subtraction_inplace(
  matrix_double_t *Amat,
  matrix_double_t *Bmat
){
  int i;
  /* Input data check */
  if ((Amat->numRows != Bmat->numRows) ||
      (Amat->numCols != Bmat->numCols)) return 0;
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
      (Amat->numCols != Bmat->numRows)) return 0;
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
  /* TODO: Make this more efficient when not using LAPACK*/
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
){
  int i, j;
  if ((Amat->numRows != ATmat->numCols) ||
      (Amat->numCols != ATmat->numRows)) return 0;
  for (i = 0; i < Amat->numRows; i++ ){
    for (j = 0; j < Amat->numCols; j++ ){
      ATmat->pData[i*Amat->numCols + j] = Amat->pData[j*Amat->numRows + i];
    }
  }
  return 1;
}

int matrix_double_solve_posdef(
  matrix_double_t *Amat,
  matrix_double_t *Bmat
){
  int info = 0;
  if ((Amat->numRows != Amat->numCols) ||
      (Amat->numRows != Bmat->numRows) ||
      (Bmat->numCols < 1)) return 0;

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
  TRACE(5,("The PSD linear system solver has not yet been implemented dependency free", info));
#endif

  if (0 == info){
    return 1;
  } else {
    if (0 < info) TRACE(5,("The %i-th argument in dposv had an illegal value\n", -info));
    if (0 > info) TRACE(5,("The leading minor of order %i is not positive definite\n", info));
    return 0;
  }
}

int matrix_double_symmetric_real_eigenvalues(
  matrix_double_t *Amat,
  matrix_double_t *eigVals
){
  int info, lwork;
  double* work;
  double wkopt;
  matrix_double_t tmpm;

  /* Check that the input dimensionlity is correct */
  if ((Amat->numRows != Amat->numCols) ||
      (eigVals->numRows != Amat->numCols) ||
      (eigVals->numCols != 1)) {
    TRACE(5,("Input dimensionality error\n"));
    return 0;
  }
  /* Check that the matrix is symmetric
  if (0 == isSymmetric(Amat)) return 0;*/

  matrix_allocate(&tmpm, Amat->numRows, Amat->numCols);
  matrix_copy(Amat, &tmpm);

  /* Query and allocate the optimal workspace */
  lwork = -1;
  dsyev_(
      "N",
      "L",
      &tmpm.numCols,
      tmpm.pData,
      &tmpm.numCols,
      eigVals->pData,
      &wkopt,
      &lwork,
      &info
    );

  lwork = (int)wkopt;
  if (0 >= lwork) {
    TRACE(5,("Could not compute the optimal work size\n"));
    return 0;
  }
  assert(lwork > 0);

  work  = (double*)malloc( lwork*sizeof(double) );

  /* Solve eigenproblem */
  dsyev_(
    "N",
    "L",
    &tmpm.numCols,
    tmpm.pData,
    &tmpm.numCols,
    eigVals->pData,
    work,
    &lwork,
    &info
  );

  free(work);
  free(tmpm.pData);

  if (0 == info){
    return 1;
  } else {
    if (0 < info) TRACE(5,("The %i-th argument in dposv had an illegal value\n", -info));
    if (0 > info) TRACE(5,("The the algorithm failed to converge; %i off-diagonal elements of an intermediate tridiagonal form did not converge to zero\n", info));
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

void mat_sol(matrix_double_t *Amat, matrix_double_t *Bmat)
{ assert(1 == matrix_double_solve_posdef(Amat, Bmat)); }

void mat_eigvals(matrix_double_t *Amat, matrix_double_t *Bmat)
{ assert(1 == matrix_double_symmetric_real_eigenvalues(Amat, Bmat)); }
