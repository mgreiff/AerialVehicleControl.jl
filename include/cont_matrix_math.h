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
* @param[in] Amat - Input matrix
* @param[in] Bmat - Input matrix
* @param[out] Cmat - Output matrix sum
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int matrix_double_addition(
  matrix_double_t *Amat,
  matrix_double_t *Bmat,
  matrix_double_t *Cmat
);

/***************************************************************************//**
* @brief In-place matrix addition, A <-- A + B
* @param[out] Amat - Overwrites the data in Amat with the sum
* @param[in] Bmat - Input matrix
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int matrix_double_addition_inplace(
  matrix_double_t *Amat,
  matrix_double_t *Bmat
);

/***************************************************************************//**
* @brief Matrix subtraction, C <-- A - B
* @param[in] Amat - Input matrix
* @param[in] Bmat - Input matrix
* @param[out] Cmat - Output matrix difference
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int matrix_double_subtraction(
  matrix_double_t *Amat,
  matrix_double_t *Bmat,
  matrix_double_t *Cmat
);

/***************************************************************************//**
* @brief In-place subtraction addition, A <-- A - B
* @param[out] Amat - Overwrites the data in Amat with the difference
* @param[in] Bmat - Input matrix
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int matrix_double_subtraction_inplace(
  matrix_double_t *Amat,
  matrix_double_t *Bmat
);

/***************************************************************************//**
* @brief Matrix multiplication, C <-- A * B
* @param[in] Amat - Input matrix
* @param[in] Bmat - Input matrix
* @param[out] Cmat - Output product matrix
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int matrix_double_multiplication(
  matrix_double_t *Amat,
  matrix_double_t *Bmat,
  matrix_double_t *Cmat
);

/***************************************************************************//**
* @brief Matrix transposition, A^T <-- transpose(A)
* @param[in] Amat - Input matrix
* @param[out] ATmat - Output matrix transpose
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int matrix_double_transposition(
  matrix_double_t *Amat,
  matrix_double_t *ATmat
);

/***************************************************************************//**
* @brief Solve the linear system AX=B, writing X into B, for PSD real A
*
* Note: This destroys and overwrites whatever is in A.
*
* @param[in] Amat - Input matrix
* @param[out] Bmat - Output solution to the linear system
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int matrix_double_solve_posdef(
  matrix_double_t *Amat,
  matrix_double_t *Bmat
);

/***************************************************************************//**
* @brief Compute, without changing A, the eigenvalues of a real symmetric A to B
*
* Note: This does not overwrite whatever is in A
*
* @param[in] Amat - Input matrix
* @param[out] Bmat - Output real-valued eigenvalues of A
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int matrix_double_symmetric_real_eigenvalues(
  matrix_double_t *Amat,
  matrix_double_t *Bmat
);

/*******************************************************************************
* External handles for the various operations, raising assertions and calling
* whatever library that is linked in. This is effectively the interface that is
* used in the controllers. Raises assertions if use dincorrectly.
*******************************************************************************/
/***************************************************************************//**
* @brief Matrix addition, C <-- A + B, asserts that this is done correctly
*******************************************************************************/
void mat_add(matrix_double_t *Amat, matrix_double_t *Bmat, matrix_double_t *Cmat);
/***************************************************************************//**
* @brief In-place matrix addition, A <-- A + B, asserts that this is done correctly
*******************************************************************************/
void mat_add_inplace(matrix_double_t *Amat, matrix_double_t *Bmat);
/***************************************************************************//**
* @brief Matrix subtraction, C <-- A - B, and asserts that this is done correctly
*******************************************************************************/
void mat_sub(matrix_double_t *Amat, matrix_double_t *Bmat, matrix_double_t *Cmat);
/***************************************************************************//**
* @brief In-place subtraction addition, A <-- A - B, asserts that this is done correctly
*******************************************************************************/
void mat_sub_inplace(matrix_double_t *Amat, matrix_double_t *Bmat);
/***************************************************************************//**
* @brief and asserts that this is done correctly
*******************************************************************************/
void mat_mul(matrix_double_t *Amat, matrix_double_t *Bmat, matrix_double_t *Cmat);
/***************************************************************************//**
* @brief Matrix transposition, A^T <-- transpose(A), asserts that this is done correctly
*******************************************************************************/
void mat_trans(matrix_double_t *Amat, matrix_double_t *ATmat);
/***************************************************************************//**
* @brief Compute, without changing A, the eigenvalues of a real symmetric A to B, asserts that this is done correctly
*******************************************************************************/
void mat_sol(matrix_double_t *Amat, matrix_double_t *Bmat);
/***************************************************************************//**
* @brief Computethe eigenvalues of a real symmetric A to B, asserts that this is done correctly
*******************************************************************************/
void mat_eigvals(matrix_double_t *Amat, matrix_double_t *Bmat);

#endif /* __CONT_MATRIX_MATH_H__ */
