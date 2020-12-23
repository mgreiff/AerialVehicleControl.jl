/*******************************************************************************
* Copyright (C) Marcus Greiff 2020
*
* @file cont_math.h
* @author Marcus Greiff
* @date June 2020
*******************************************************************************/
#ifndef __CONT_MATH_H__
#define __CONT_MATH_H__

#include "cont_types.h"
#include "cont_matrix_math.h"
#include "math.h"

/***************************************************************************//**
* @brief Conjugate an element of SU(2)
* @param[in] qin - Pointer to a four-dimensional input quaternion
* @param[out] qout - Pointer to a four-dimensional output quaternion
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int cont_SU2_conjugate(matrix_double_t * qin, matrix_double_t * qout);

/***************************************************************************//**
* @brief Product of two elements of SU(2)
* @param[in] q - Pointer to an input element of SU(2) (four-dimensional vector)
* @param[in] p - Pointer to an input element of SU(2) (four-dimensional vector)
* @param[out] out - Pointer to an output quaternion product (four-dimensional vector)
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int cont_SU2_product( matrix_double_t * q, matrix_double_t * p, matrix_double_t * out);

/***************************************************************************//**
* @brief The vee map on SU(2)
* @param[in] in - Pointer to an input element of the lie-agbera of SU(2)
* @param[out] out - Pointer to a three-dimensional output (one parameter group)
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int cont_SU2_vee(matrix_double_t * in, matrix_double_t * out);

/***************************************************************************//**
* @brief The hat map on SU(2)
* @param[in] in - Pointer to a three-dimensional output (one parameter group)
* @param[out] out - Pointer to an input element of the lie-agbera of SU(2)
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int cont_SU2_hat(matrix_double_t * in, matrix_double_t * out);

/***************************************************************************//**
* @brief The Exp map on SU(2)
* @param[in] in - Pointer to a three-dimensional vector (one parameter group)
* @param[out] out - Pointer to an element of SU(2) (here a 4-dimenisonal vector)
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int cont_SU2_Exp(matrix_double_t * in, matrix_double_t * out);

/***************************************************************************//**
* @brief Computes a triple product on SU(2)
*
* The inputs a, b, and c all represent angles in radians, and the triple product
* is computed by taking the exponential map with respect to these angles as
*
*     qa = Exp([a, 0, 0])
*     qb = Exp([0, b, 0])
*     qc = Exp([0, 0, c])
*
* and then computing their consequtive product
*
*     out = qa * qb * qc
*
* @param[in] in - Pointer to a three-dimensional vector (one parameter group)
* @param[out] out - Pointer to an element of SU(2) (here a 4-dimenisonal vector)
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int cont_SU2_triple( double a, double b, double c, matrix_double_t * out);

/***************************************************************************//**
* @brief Evaluate distance Gamma(q1, q2) in [0,2] defined on SU2
*
* This distance is defined as Gamma(X1, X2) = (1/2)Tr(I - X1'*X2) in [0,2]
*
* @param[in] q1 - Pointer to a four-dimensional input quaternion
* @param[in] q2 - Pointer to a four-dimensional output quaternion
* @return distance - The evaluated distance on SU2
*******************************************************************************/
double cont_SU2_distance(matrix_double_t * q1, matrix_double_t * q2);

/***************************************************************************//**
* @brief Embedding relating SU(2) to SO(3)
* @param[in] q - Pointer to a four-dimensional input quaternion
* @param[out] R - Pointer to an output 3x3 rotation matrix
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int cont_quat_2_SO3(matrix_double_t * q, matrix_double_t * R);

/***************************************************************************//**
* @brief The hat operator on SO3
* @param[in] in - Pointer to a three-dimensional vector (one parameter group)
* @param[out] out - Pointer to an element of the 3x3 lie algebra of SO(3)
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int cont_SO3_hat(matrix_double_t * in, matrix_double_t * out);

/***************************************************************************//**
* @brief The hat operator on SO3
* @param[in] in - Pointer to an element of the 3x3 lie algebra of SO(3)
* @param[out] out - Pointer to a three-dimensional vector (one parameter group)
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int cont_SO3_vee(matrix_double_t * in, matrix_double_t * out);

/***************************************************************************//**
* @brief The Exp map on SO3
* @param[in] u - Pointer to a three-dimensional vector (one parameter group)
* @param[out] R - Pointer to an output 3x3 rotation matrix
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int cont_SO3_Exp(matrix_double_t * u, matrix_double_t * R);

/***************************************************************************//**
* @brief The Log map on SO3
* @param[in] R - Pointer to an output 3x3 rotation matrix
* @param[out] u - Pointer to a three-dimensional vector (one parameter group)
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int cont_SO3_Log(matrix_double_t * R, matrix_double_t * u );

/***************************************************************************//**
* @brief The distance Phi(R1, R2) defined on SO3
*
* This distance is defined as Phi(R1, R2) = (1/2)Tr(I - R1'*R2) in [0,2]
*
* @param[in] R1 - Pointer to an output 3x3 rotation matrix
* @param[in] R2 - Pointer to an output 3x3 rotation matrix
* @return distance - The evaluated distance on SO(3)
*******************************************************************************/
double cont_SO3_distance(matrix_double_t * R1, matrix_double_t * R2);

/***************************************************************************//**
* @brief Evaluates the sinc function as f(x) = sin(x)/x
* @param[in] x - the value at which the function is to be evaluated
* @return output - the output of the function f(x)
*******************************************************************************/
double cont_sinc(double x);

/***************************************************************************//**
* @brief Evaluates the dot product of two vectors
 * @param[in] vecA - Input vector
 * @param[in] vecB - Input vector
 * @return output - dot product
*******************************************************************************/
double cont_dot_product(matrix_double_t * vecA, matrix_double_t * vecB);

/***************************************************************************//**
* @brief  Cross product operation, To be migrated to the cont_math stack
*
* @param[in] inAm  - Input vector A in R^3
* @param[in] inBm  - Input vector B in R^3
* @param[out] outm - Output as the cross-product of A and B, out=AxB, in R^3
* @return void
*******************************************************************************/
void cont_cross_product( matrix_double_t * inAm, matrix_double_t * inBm, matrix_double_t * outm);

/***************************************************************************//**
* @brief Takes the sign of the input, returning +1 if the input is zero.
* @param[in] in - Input argument
* @return output - Output sign in {-1.0, 1.0}
*******************************************************************************/
double cont_sign_func(double in);

/***************************************************************************//**
* @brief Normalizes the vector
* @param[in/out] vec - Input vector with either 1 row, or 1 column
* @return output - 1 if successful, 0 otherwise
*******************************************************************************/
int cont_normalize(matrix_double_t * vec);

#endif /* __CONT_MATH.H__ */
