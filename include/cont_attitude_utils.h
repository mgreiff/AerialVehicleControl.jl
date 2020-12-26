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
* @brief TODO: Write docstring
*
*
*
*
*******************************************************************************/
int assert_attitude_FSF_SU2(
  con_state_qw_fsf_t * controller
);

int assert_attitude_FSF_parameters(
  con_state_qw_fsf_t * controller,
  double *minJ,
  double *maxJ
);

#endif /* __CONT_ATTITUDE_UTILS_H__ */
