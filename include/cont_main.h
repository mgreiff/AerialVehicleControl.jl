/*******************************************************************************
* Copyright (C) Marcus Greiff 2020
*
* @file cont_main.h
* @author Marcus Greiff
* @date June 2020
*******************************************************************************/
#ifndef __CONT_MAIN_H__
#define __CONT_MAIN_H__

#include "cont_attitude_reference_generator.h"
#include "cont_power_distribution.h"
#include "cont_attitude_FSF_SO3_continuous.h"
#include "cont_attitude_FSF_SO3_robust.h"
#include "cont_attitude_FSF_SU2_continuous.h"
#include "cont_attitude_FSF_SU2_discontinuous.h"
#include "cont_attitude_FSF_SU2_robust.h"
#include "cont_attitude_FOF_SO3_continuous.h"

/***************************************************************************//**
* @brief Main program, facilitating memory checks and analysis of the code
*
* This function calls the two attitude examples, in order to check that all
* of the controllers run and to test their performance and memory usage.
*******************************************************************************/
int main(int argc, char *argv[]);

/***************************************************************************//**
* @brief Example of how the attitude FSF controllers can be implemented
*******************************************************************************/
void example_attitude_FSF(void);

/***************************************************************************//**
* @brief Example of how the attitude FOF controllers can be implemented
*******************************************************************************/
void example_attitude_FOF(void);

#endif /* __CONT_MAIN_H__ */
