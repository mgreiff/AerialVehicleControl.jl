/*******************************************************************************
* Copyright (C) Marcus Greiff 2020
*
* @file cont_attitude_utils.c
* @author Marcus Greiff
* @date June 2020
*******************************************************************************/
#include "cont_attitude_utils.h"

int assert_attitude_FSF_SO3(
  con_state_qw_fsf_t * controller
){
  double maxJ, minJ, minVal, temp;
  double kR = controller->gain_kR;
  double kc = controller->gain_kc;
  double kw = controller->gain_kw;

  /* Check that the parameters have the correct sign, and evaluate eig(J) */
  if (0 == assert_attitude_FSF_parameters(controller, &minJ, &maxJ)) return 0;

  /* Check that kc is set sufficiently small w.r.t. kX and kw */
  minVal = kw;
  temp   = sqrt(kR * minJ);
  if (temp < minVal) minVal = temp;
  temp   = 4.0 * kw * kR * minJ * minJ;
  temp  /= (maxJ * kw * kw + 4.0 * minJ * minJ * kR);
  if (temp < minVal) minVal = temp;
  if (kc  > minVal){
    TRACE(5, ("The controller parameters do not result in a feasible controller\n"));
    return 0;
  }
  return 1;
}

int assert_attitude_FSF_SU2(
  con_state_qw_fsf_t * controller
){
  double maxJ, minJ, minVal, temp;
  double kR = controller->gain_kR;
  double kc = controller->gain_kc;
  double kw = controller->gain_kw;

  /* Check that the parameters have the correct sign, and evaluate eig(J) */
  if (0 == assert_attitude_FSF_parameters(controller, &minJ, &maxJ)) return 0;

  /* Check that kc is set sufficiently small w.r.t. kX and kw */
  minVal = 4.0 * kw;
  temp   = 2.0 * sqrt(kw * minJ);
  if (temp < minVal) minVal = temp;
  temp   = 4.0 * kw * kR * minJ * minJ;
  temp  /= (maxJ * kw * kw + minJ * minJ * kR);
  if (temp < minVal) minVal = temp;
  if (kc  > minVal){
    TRACE(5, ("The controller parameters do not result in a feasible controller\n"));
    return 0;
  }
  return 1;
}

int assert_attitude_FSF_parameters(
  con_state_qw_fsf_t * controller,
  double *minJ,
  double *maxJ
){
  double tol = 0.000001;
  double kR = controller->gain_kR;
  double kc = controller->gain_kc;
  double kw = controller->gain_kw;
  matrix_double_t Jm, Em;

  matrix_define(&Jm,  3, 3, controller->inertia);
  matrix_allocate(&Em,  3, 1);

  /* Check for symmetry */
  if ((tol < abs(matrix_get(&Jm, 0, 1) - matrix_get(&Jm, 1, 0))) ||
      (tol < abs(matrix_get(&Jm, 0, 2) - matrix_get(&Jm, 2, 0))) ||
      (tol < abs(matrix_get(&Jm, 1, 2) - matrix_get(&Jm, 2, 1)))) {
    TRACE(5, ("The inertia matrix must be symmetric\n"));
    free(Em.pData);
    return 0;
  }
  /* Compute the maximum and minimum eigenvalues of the inertia matrix */
  mat_eigvals(&Jm, &Em);
  *minJ = Em.pData[0];
  *maxJ = Em.pData[2];

  if (*minJ < -tol){
    TRACE(5, ("The inertia matrix must be positive definite\n"));
    free(Em.pData);
    return 0;
  }
  if ((kR <= 0) ||
      (kc <= 0) ||
      (kw <= 0) ||
      (controller->param_a <= 0) ||
      (controller->param_b <= 0) ||
      (controller->param_c <= 0) ||
      (controller->param_d <= 0)){
    TRACE(5, ("The controller parameters must be positive\n"));
    free(Em.pData);
    return 0;
  }
  free(Em.pData);
  return 1;
}
