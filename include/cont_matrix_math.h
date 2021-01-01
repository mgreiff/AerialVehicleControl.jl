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
* @brief Matrix addition, C <-- A + B
*******************************************************************************/
int matrix_double_addition(
  matrix_double_t *Amat,
  matrix_double_t *Bmat,
  matrix_double_t *Cmat
);

/***************************************************************************//**
* @brief In-place matrix addition, A <-- A + B
*******************************************************************************/
int matrix_double_addition_inplace(
  matrix_double_t *Amat,
  matrix_double_t *Bmat
);

/***************************************************************************//**
* @brief Matrix subtraction, C <-- A - B
*******************************************************************************/
int matrix_double_subtraction(
  matrix_double_t *Amat,
  matrix_double_t *Bmat,
  matrix_double_t *Cmat
);

/***************************************************************************//**
* @brief In-place matrix subtraction, A <-- A - B
*******************************************************************************/
int matrix_double_subtraction_inplace(
  matrix_double_t *Amat,
  matrix_double_t *Bmat
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
void mat_sol(matrix_double_t *Amat, matrix_double_t *Bmat);
void mat_eigvals(matrix_double_t *Amat, matrix_double_t *Bmat);

#endif /* __CONT_MATRIX_MATH_H__ */
