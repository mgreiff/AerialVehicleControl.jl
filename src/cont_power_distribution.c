/*******************************************************************************
* Copyright (C) Marcus Greiff 2020
*
* @file cont_power_distribution.c
* @author Marcus Greiff
* @date June 2020
*******************************************************************************/
#include "cont_power_distribution.h"

int compute_power_distribution(
  con_state_qw_fsf_t * controller,
  matrix_double_t * dutyCycles
){
  int i;
  double a = controller->param_a;
  double b = controller->param_b;
  double c = controller->param_c;
  double dx = controller->param_d;
  double dy = controller->param_d;
  double adtb;
  double f[4];

  /* Input error check*/
  if (controller->thrust < 0.0) {
    TRACE(5, ("The target thrust cannot be negative in the power distribution\n"));
    return 0;
  }
  if (dutyCycles->numCols != 1 || dutyCycles->numRows != 4) {
    TRACE(5, ("The dutyCycle vector is not of the correct dimensions\n"));
    return 0;
  }
  if ((a <= 0.0) || (b <= 0.0) || (c <= 0.0)) {
    TRACE(5, ("The input parameters must be positive\n"));
    return 0;
  }

  /* Compute the individual rotor forces */
  /*
  f[0] = controller->thrust +tid*controller->torque[1] -controller->torque[2]/c;
  f[1] = controller->thrust -tid*controller->torque[0] +controller->torque[2]/c;
  f[2] = controller->thrust -tid*controller->torque[1] -controller->torque[2]/c;
  f[3] = controller->thrust +tid*controller->torque[0] +controller->torque[2]/c;
  */

  /*
  for the razor drone, with

  A = [   1,   1,  1,   1]
      [ -dy, -dy, dy,  dy]
      [ -dx,  dx, dx, -dx]
      [   c,  -c,  c,  -c]

  we instead get

      inv(A) = [ 1/4, -1/(4*dy), -1/(4*dx),  1/(4*c)]
               [ 1/4, -1/(4*dy),  1/(4*dx), -1/(4*c)]
               [ 1/4,  1/(4*dy),  1/(4*dx),  1/(4*c)]
               [ 1/4,  1/(4*dy), -1/(4*dx), -1/(4*c)]
  */
  /* Todo remove hard coded coefficients
  dx =  0.09/2.0;
  dy = 0.152/2.0;*/

  f[0] = controller->thrust -(1/dy)*controller->torque[0] -(1/dx)*controller->torque[1] +(1/c)*controller->torque[2];
  f[1] = controller->thrust -(1/dy)*controller->torque[0] +(1/dx)*controller->torque[1] -(1/c)*controller->torque[2];
  f[2] = controller->thrust +(1/dy)*controller->torque[0] +(1/dx)*controller->torque[1] +(1/c)*controller->torque[2];
  f[3] = controller->thrust +(1/dy)*controller->torque[0] -(1/dx)*controller->torque[1] -(1/c)*controller->torque[2];
  for (i = 0; i < 4; i++) f[i] /= 4.0;

  /* TODO make this saturation continuous using a hyperbolic tan function*/
  for (i = 0; i < 4; i++){
    if(f[i] < 0.0){
      TRACE(5, ("Got a negative force, f[%i]=%f, rounding to zero.\nThis should not generally happen - your controller might be too agressively tuned.\n", i, f[i]));
      f[i] = 0.0;
    }
  }

  /* Precomputation of parameters*/
  adtb = a/(2.0*b);

  /* Compute the duty cycles */
  for (i = 0; i < 4; i++ ) dutyCycles->pData[i] = - adtb + sqrt(adtb*adtb + f[i]/b);

  /* Assert that the duty cycles are positive*/
  for (i = 0; i < 4; i++ ){
    if (dutyCycles->pData[i]  < 0.0){
      TRACE(5, ("The duty cycles must be positive\n"));
      return 0;
    }
    if (dutyCycles->pData[i] > 1.0){
      dutyCycles->pData[i] = 1.0;
      TRACE(5, ("Warning, attemting to set a duty cycle > 1 - should never happen\n"));
    }
  }
  return 1;
}
