/*******************************************************************************
* Copyright (C) Marcus Greiff 2020
*
* @file cont_attitude_utils.h
* @author Marcus Greiff
* @date June 2020
*******************************************************************************/
#ifndef __CONT_ATTITUDE_UTILS_H__
#define __CONT_ATTITUDE_UTILS_H__

#include "cont_types.h"
#include "cont_math.h"
#include "cont_matrix_math.h"

/***************************************************************************//**
* @brief Assert that the SO(3) controller is feasible
*
* This function should be run any time that the controller parameters are
* changed, but may also be called on each time-step as a precauiton.
*
* @param[in] controller - State of the controller
* @return status        - 1 if feasible, 0 otherwise.
*******************************************************************************/
int assert_attitude_FSF_SO3(
  con_state_qw_fsf_t * controller
);

/***************************************************************************//**
* @brief Assert that the SU(2) controller is feasible
*
* This function should be run any time that the controller parameters are
* changed, but may also be called on each time-step as a precauiton.
*
* @param[in] controller - State of the controller
* @return status        - 1 if feasible, 0 otherwise.
*******************************************************************************/
int assert_attitude_FSF_SU2(
  con_state_qw_fsf_t * controller
);

/***************************************************************************//**
* @brief Assert that the SO(3) controller is feasible
*
* Here, the maximum and minimum eigenvalues of the inertia matrx can potentially
* be computed on tsartup if the inertia matrix is constant in time
*
* @param[in] controller - State of the controller
* @param[out] *minJ     - Output minimum eigenvalue of the inertia matrix
* @param[out] *maxJ     - Output maximum eigenvalue of the inertia matrix
* @return status        - 1 if feasible, 0 otherwise.
*******************************************************************************/
int assert_attitude_FSF_parameters(
  con_state_qw_fsf_t * controller,
  double *minJ,
  double *maxJ
);

#endif /* __CONT_ATTITUDE_UTILS_H__ */
