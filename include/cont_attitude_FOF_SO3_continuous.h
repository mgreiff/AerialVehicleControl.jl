/*******************************************************************************
* Copyright (C) Marcus Greiff 2020
*
* @file cont_attitude_FOF_SO3_continuous.h
* @author Marcus Greiff
* @date June 2020
*******************************************************************************/
#ifndef __CONT_ATTITUDE_FOF_SO3_CONTINUOUS_H__
#define __CONT_ATTITUDE_FOF_SO3_CONTINUOUS_H__

#include "cont_math.h"

/***************************************************************************//**
* @brief  Update using the continuous filtered outpot feedback (FOF) on SO(3)
*
* This function takes the states of the reference dynamics in "reference"
* structure, as well as the current states of the system in the "state"
* structure, and the controller settings/memory in the "controller" structure,
* and updates the torque fields in the controller for actuation, as well as the
* observer states and relevant distances based on this information.
*
* @param[in] reference       - State of the reference dynamics
* @param[in] state           - State of the system dynamics
* @param[in, out] controller - State of the controller, updated by the call
* @param[in] dt              - Time since the previous controller update [s]
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int update_attitude_FOF_SO3_continuous(
  ref_state_qw_t     * reference,
  dyn_state_qw_t     * state,
  con_state_qw_fof_t * controller,
  double               dt
);

/***************************************************************************//**
* @brief  Help function used to update the measurements in the controller struct
*
* This is only required when calling the function from Julia, where the structs
* are immutable, but for safety reasons, it is best to set the measurements
* using this function instead of direct manipulation of the controller struct
* memory.
*
* @param[in] controller - The Controlle robject
* @param[in] y0m        - Gyroscopic measurement
* @param[in] y1m        - Directional measurement, corresponding to v1
* @param[in] y2m        - Directional measurement, corresponding to v2
* @param[in] y3m        - Directional measurement, corresponding to v3
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int update_attitude_FOF_SO3_continuous_measurements(
  con_state_qw_fof_t * controller,
  matrix_double_t * y0m,
  matrix_double_t * y1m,
  matrix_double_t * y2m,
  matrix_double_t * y3m
);

/***************************************************************************//**
* @brief  Help function used to compute the innovation terms in the controller
*
* The gains are stored in an input array on the form gain = [k1,...,kN], and the
* directions are similarly stored as YA = [yA1,...,YAN], and YA = [yB1,...,YBN],
* and the function evaluate sthe weighted cross-product sum on the form.
*
*     out = sum_{i = 1}^N k_i * S(yAi) * yBi
*
* where S(*) denotes the hat map on SO(3), while checking for dimensionality
* errors and verifying that the gains are postrictly positive.
*
* @param[in] YAm   - Input directions (3xN) dimensional
* @param[in] YBm   - Input directions (3xN) dimensional
* @param[in] gainm - Controller gains (1xN) dimensional
* @param[out] outm - Function output,
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int attitude_FOF_SO3_continuous_cross_terms(
  matrix_double_t * YAm,
  matrix_double_t * YBm,
  matrix_double_t * gainm,
  matrix_double_t * outm
);

/***************************************************************************//**
* @brief  Help function to simulate the attitude observer dynamics
*
*   d(Q)/dt = Q * [ deltaQ / 2 ]_{SU(2)}^{\^}
*   d(W)/dt = deltaW
*
* For simplicity, take a first order euler method and let
*
*   Q(t + dt) = Q(t) * expm(dt * [ deltaQ(t) / 2 ]_{SU(2)}^{\^})
*   W(t + dt) = W(t) + dt * deltaW(t)
*
* @param[in] Qm      - Attitude configured on SU(2)
* @param[in] Wm      - Attitude rates configured on R^3
* @param[in] deltaQm - Time derivative of the signal Qm (one parameter group, R^3)
* @param[in] deltaWm - Time derivative of the signal Wm (on R^3)
* @param[in] dt      - Propagates the continuous time dynamics dt [s] ahead
*                      a first order euler method. This should eventually be
*                      replaced by a Crouch-Grossman method.
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int attitude_FOF_SO3_continuous_simulate_observer(
  matrix_double_t * Qm,
  matrix_double_t * Wm,
  matrix_double_t * deltaQm,
  matrix_double_t * deltaWm,
  double dt
);

#endif /* __CONT_ATTITUDE_FOF_SO3_CONTINUOUS_H__ */
