/*******************************************************************************
* Copyright (C) Marcus Greiff 2020
*
* @file cont_power_distribution.h
* @author Marcus Greiff
* @date June 2020
*******************************************************************************/
#ifndef __CONT_POWER_DISTRIBUTION_H__
#define __CONT_POWER_DISTRIBUTION_H__

#include "cont_math.h"
#include "cont_matrix_math.h"

/***************************************************************************//**
* @brief  Compute the desired rotor thrusts.

* This function assumes a linear relationship between the rotor thrusts and the
* rigid-body forces and torques as
*
*     |f |   | +1, +1, +1, +1| |f1|
*     |t1|   |  0, -d,  0, +d| |f2|
*     |t2| = | +d,  0, -d,  0| |f3|
*     |t3|   | -c, +c, -c, +c| |f4|
*
* This map is subsequently inverted, as
*
*     |f1|         | 1,        0,  2/d, -1/c| |f |
*     |f2|         | 1, -2/d,        0, +1/c| |t1|
*     |f3| = (1/4) | 1,        0, -2/d, -1/c| |t2|
*     |f4|         | 1, +2/d,        0, +1/c| |t3|
*
* These forces should be positive at all times, and are therefore
* converted to a duty-cycle for each of the rotors using a polynomial
* approximation of the duty cycle to thrust map, identified offline. In general,
* this map can be well approximated with second order polynomials which are zero
* at the origin. Consequently, we assume that
*
*     fi = h(di) = a * di + b * di^2
*
* and the inverse map is then taken as
*
*     di = g(fi) =  -(a/(2*b)) + sqrt((a/(2*b))^2 + fi/b) > 0 for all fi > 0
*
* @param[in]  controller - The current controller state
* @param[out] dutyCycles - The corresponding cycles thfor the rotor control
* @return status - 1 if successful, 0 otherwise.
*******************************************************************************/
int compute_power_distribution(
  con_state_qw_fsf_t * controller,
  matrix_double_t * dutyCycles
);

#endif /* __CONT_POWER_DISTRIBUTION_H__ */
