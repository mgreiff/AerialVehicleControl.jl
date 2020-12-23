/*******************************************************************************
* Copyright (C) Marcus Greiff 2020
*
* @file cont_matrix_math.h
* @author Marcus Greiff
* @date June 2020
*******************************************************************************/
#ifndef __CONT_MATRIX_MATH_H__
#define __CONT_MATRIX_MATH_H__

#include "cont_types.h"

/***************************************************************************//**
* @brief Matrix addition, A + B = C
*******************************************************************************/
int matrix_double_addition(
  matrix_double_t *Amat,
  matrix_double_t *Bmat,
  matrix_double_t *Cmat
);

/***************************************************************************//**
* @brief Matrix subtraction, A - B = C
*******************************************************************************/
int matrix_double_subtraction(
  matrix_double_t *Amat,
  matrix_double_t *Bmat,
  matrix_double_t *Cmat
);

/***************************************************************************//**
* @brief Matrix multiplication, A * B = C
*******************************************************************************/
int matrix_double_multiplication(
  matrix_double_t *Amat,
  matrix_double_t *Bmat,
  matrix_double_t *Cmat
);

/***************************************************************************//**
* @brief Matrix transposition, A^T = A
*******************************************************************************/
int matrix_double_transposition(
  matrix_double_t *Amat,
  matrix_double_t *ATmat
);

/***************************************************************************//**
* @brief Inplace cholesky decomposition A = L * L^T, writes over A with L
*******************************************************************************/
int matrix_double_cholesky_inplace(
  matrix_double_t *Amat
);

/***************************************************************************//**
* @brief Inplace cholesky decomposition A = L * L^T, writes L to Lmat
*******************************************************************************/
int matrix_double_cholesky(
  matrix_double_t *Amat,
  matrix_double_t *Lmat
);

/***************************************************************************//**
* @brief Solve the linear system AX=B, writing X into B, for PSD real A
*******************************************************************************/
int matrix_double_solve_posdef(
  matrix_double_t *Amat,
  matrix_double_t *Bmat
);

/***************************************************************************//**
* @brief Compute, without changing A, the eigenvalues of a real symmetric A to B
*******************************************************************************/
int matrix_double_symmetric_real_eigenvalues(
  matrix_double_t *Amat,
  matrix_double_t *Bmat
);

/*******************************************************************************
* External handles for the various operations, raising assertions and calling
* whatever library that is linked in - if any
*******************************************************************************/
void mat_add(matrix_double_t *Amat, matrix_double_t *Bmat, matrix_double_t *Cmat);
void mat_add_inplace(matrix_double_t *Amat, matrix_double_t *Bmat);
void mat_sub(matrix_double_t *Amat, matrix_double_t *Bmat, matrix_double_t *Cmat);
void mat_sub_inplace(matrix_double_t *Amat, matrix_double_t *Bmat);
void mat_mul(matrix_double_t *Amat, matrix_double_t *Bmat, matrix_double_t *Cmat);
void mat_trans(matrix_double_t *Amat, matrix_double_t *ATmat);
void mat_chol(matrix_double_t *Amat, matrix_double_t *Lmat);
void mat_chol_inplace(matrix_double_t *Amat);
void mat_sol(matrix_double_t *Amat, matrix_double_t *Bmat);
void mat_eigvals(matrix_double_t *Amat, matrix_double_t *Bmat);

#endif /* __CONT_MATRIX_MATH_H__ */
