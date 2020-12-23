/*******************************************************************************
* Copyright (C) Marcus Greiff 2020
*
* @file cont_main.h
* @author Marcus Greiff
* @date June 2020
*******************************************************************************/
#ifndef __CONT_MAIN_H__
#define __CONT_MAIN_H__

#include "cont_types.h"
#include "cont_attitude_reference_generator.h"
#include "cont_power_distribution.h"
#include "cont_attitude_FSF_SO3_continuous.h"
#include "cont_attitude_FSF_SU2_continuous.h"
#include "cont_attitude_FSF_SU2_discontinuous.h"
#include "cont_attitude_FOF_SO3_continuous.h"

int main();

void example_attitude_FSF(void);

void example_attitude_FOF(void);

#endif /* __CONT_MAIN_H__ */
