/*******************************************************************************
* Copyright (C) Marcus Greiff 2020
*
* @file cont_types.h
* @author Marcus Greiff
* @date June 2020
*******************************************************************************/
#ifndef __CONT_TYPES_H__
#define __CONT_TYPES_H__

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#define CONT_MASS 0.1
#define CONT_GRAVACC 9.81

/** \cond INTERNAL */
/* Trace level and print options - nominally set to zero to disable outputs   */
#ifndef TRACE_LEVEL
#define TRACE_LEVEL 5
#endif
#ifndef TRACE_LINES
#define TRACE_LINES 1
#endif

/***************************************************************************//**
* @brief Custom debug print function
*******************************************************************************/
void dbg_printf(const char *fmt, ...);
/** \endcond */

/***************************************************************************//**
* @brief Trace print macro - can output both to console and file
* @param[in] level - The trace level
* @param[in] x - The arguments which are to be printed
*******************************************************************************/
#define TRACE(level, x) do {\
if (level <= TRACE_LEVEL) {\
if (1 == TRACE_LINES) {\
fprintf(stderr, "\033[0;34m--- %s - %d ---\033[0m\n", __FILE__, __LINE__);\
} dbg_printf x;}} while (0)


/***************************************************************************//**
* @brief Matrix object used for all matrix manipulation
*******************************************************************************/
typedef struct matrix_double_s {
  int numRows;          /**< Number of columns                                */
  int numCols;          /**< Number of columns                                */
  double *pData;        /**< Pointer to continuous chuck of memory            */
} matrix_double_t;

/***************************************************************************//**
* @brief Dynamical states assumed known in the attitude FSF on S(3) or SU(2)
*******************************************************************************/
typedef struct dyn_state_qw_s {
  double quaternion[4]; /**< Quaternion attitude                              */
  double omega[3];      /**< Attitude rates                                   */
} dyn_state_qw_t;

/***************************************************************************//**
* @brief Complete state of the attitude FSF on SO(3) or SU(2)
*******************************************************************************/
typedef struct con_state_qw_fsf_s {
  int status;           /**< Controller status - 0 for idle 1 for active      */
  double thrust;        /**< Controller thrust  (for actuation)               */
  double torque[3];     /**< Controller torques (for actuation)               */
  double inertia[9];    /**< Controller inertia                               */
  double gain_kR;       /**< Gain relating to the attitude error              */
  double gain_kc;       /**< Gain relating to the cross-terms                 */
  double gain_kw;       /**< Gain relating to the attitude rate error         */
  double param_a;       /**< Constant related to the PWM-thrust map           */
  double param_b;       /**< Constant related to the PWM-thrust map           */
  double param_c;       /**< Constant related to rotor motor torque           */
  double param_d;       /**< Distance of the rotors to the center of mass     */
  double dist_Psi;      /**< Attitude error on SO(3) (for analysis)           */
  double dist_Gamma;    /**< Attitude error on SU(2) (for analysis)           */
  double dist_lyapunov; /**< Lyapunov function (for analysis)                 */
} con_state_qw_fsf_t;

/***************************************************************************//**
* @brief Complete state of the attitude FOF on SO(3) or SU(2)
*******************************************************************************/
typedef struct con_state_qw_fof_s {
  int status;                  /**< Controller status, 0 idle, 1 for active   */
  double torque[3];            /**< Controller torques (for actuation)        */
  double inertia[9];           /**< Controller inertia                        */
  double invinertia[9];        /**< Controller inertia inverse                */
  double gain_Kw[9];           /**< Gain relating to the attitude error       */
  double gain_ki[3];           /**< Gain relating to the attitude rate error  */
  double gain_Cw[9];           /**< Gain relating to the cross-terms          */
  double gain_cw;              /**< Gain relating to the attitude rate error  */
  double gain_cR;              /**< Gain relating to the attitude rate error  */
  double param_a;              /**< Constant related to the PWM-thrust map    */
  double param_b;              /**< Constant related to the PWM-thrust map    */
  double param_c;              /**< Constant related to rotor motor torque    */
  double param_d;              /**< Distance of the rotors to the c.o.m.      */
  double measuredGyrorates[3]; /**< Meas. 0   (gyroscopic rates) y0           */
  double measuredDirections[9];/**< Meas. 1-3 (dirs., acc./mag.) [y1, y2, y3] */
  double globalDirections[9];  /**< Global directions            [v1, v2, v3] */
} con_state_qw_fof_t;

/***************************************************************************//**
* @brief Reference signal structure for the attitude FSF on S(3) or SU(2)
*******************************************************************************/
typedef struct ref_state_qw_s {
  double quaternion[4]; /**< Reference quaternion attitude                    */
  double omega[3];      /**< Reference attitude rates (body frame)            */
  double alpha[3];      /**< Reference attitude acceleratoins (body frame)    */
  double thrust;        /**< Smoothened commanded lateral thrust              */
} ref_state_qw_t;

/* Matrix utility functions */
void matrix_print( matrix_double_t *mat, char * variableName);
double matrix_get( matrix_double_t *mat, int row, int column );
void matrix_set( matrix_double_t *mat, int row, int column, double value );
void matrix_zero( matrix_double_t *matrix );
void matrix_identity( matrix_double_t *matrix );
void matrix_allocate( matrix_double_t *matrix, int numRows, int numCols );
void matrix_define( matrix_double_t *matrix, int numRows, int numCols, double *data );
void matrix_copy(matrix_double_t *Amat, matrix_double_t *Bmat);

#endif /* __CONT_TYPES_H__ */
