/*******************************************************************************
* Copyright (C) Marcus Greiff 2020
*
* @file cont_attitude_FSF_SO3_continuous.h
* @author Marcus Greiff
* @date June 2020
*******************************************************************************/
#ifndef __CONT_ATTITUDE_FSF_SO3_ROBUST_H__
#define __CONT_ATTITUDE_FSF_SO3_ROBUST_H__

#include "cont_attitude_utils.h"

/***************************************************************************//**
* @brief  Update using the robust full state feedback (FSF) on SO(3)
*
* This function takes the states of the reference dynamics in "reference"
* structure, as well as the current states of the system in the "state"
* structure, and the controller settings/memory in the "controller" structure,
* and updates the torque fields in the controller, as well as the various
* distances based on this information. The update is done using the robust
* feedback law defined with distances on SO(3).
*
* @param[in] reference - State of the reference dynamics
* @param[in] state - State of the system dynamics
* @param[in, out] controller - State of the controller, updated by the call
* @return status - 1 if successful, 0 otherwise
*******************************************************************************/
int update_attitude_FSF_SO3_robust(
  ref_state_qw_t     * reference,
  dyn_state_qw_t     * state,
  con_state_qw_fsf_t * controller
);

#endif /* __CONT_ATTITUDE_FSF_SO3_ROBUST_H__ */
