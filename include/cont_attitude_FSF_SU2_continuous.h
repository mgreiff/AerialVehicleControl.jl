/*******************************************************************************
* Copyright (C) Marcus Greiff 2020
*
* @file cont_attitude_FSF_SU2_continuous.h
* @author Marcus Greiff
* @date June 2020
*******************************************************************************/
#ifndef __CONT_ATTITUDE_FSF_SU2_CONTINUOUS_H__
#define __CONT_ATTITUDE_FSF_SU2_CONTINUOUS_H__

#include "cont_types.h"
#include "cont_math.h"
#include "cont_matrix_math.h"
#include "cont_attitude_utils.h"

/***************************************************************************//**
* @brief  Update using the continuous full state feedback (FSF) on SU(2)
*
* This function takes the states of the reference dynamics in ''reference''
* structure, as well sa the current states of the system in the "state"
* structure, and the controller settings/memory in the "controller" structure,
* and updates the torque fields in the controller, as well as the various
* distances based on this information. The updateis done using the continuous
* feedback law n SU(2).
*
* @param[in] reference - State of the reference dynamics
* @param[in] state - State of the system dynamics
* @param[in, out] controller - State of the controller, updated by the call
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int update_attitude_FSF_SU2_continuous(
  ref_state_qw_t     * reference,
  dyn_state_qw_t     * state,
  con_state_qw_fsf_t * controller
);

/***************************************************************************//**
* @brief  Assert that the controller tuning is feasible in the FSF controller
*
* This function takes the controller object and outputs an integer 1 if the
* controller parameters are feasible, and a 0 if not.
*
* @param[in] controller - State of the controller
* @param[in] phi - Upper bound on the initial attitude error phi \in (0, 2)
* @return status - 1 if feasible, 0 otherwise
*******************************************************************************/
int assert_attitude_FSF_SU2_continuous(
  con_state_qw_fsf_t * controller
);

#endif /* __CONT_ATTITUDE_FSF_SU2_CONTINUOUS_H__ */
