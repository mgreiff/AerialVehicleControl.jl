/*******************************************************************************
* Copyright (C) Marcus Greiff 2020
*
* @file cont_attitude_FSF_SU2_discontinuous.c
* @author Marcus Greiff
* @date June 2020
*******************************************************************************/
#include "cont_attitude_FSF_SU2_discontinuous.h"

int update_attitude_FSF_SU2_discontinuous(
  ref_state_qw_t * reference,
  dyn_state_qw_t * state,
  con_state_qw_fsf_t * controller
){
  int i;
  matrix_double_t Xm, Xrm, Xem, XeCm, eXm, Wm, Wrm, ewm, Arm;
  matrix_double_t tmp41Am, tmp41Bm, tmp41Cm, XeCmhatWrmXem;
  matrix_double_t tmp31Am, tmp31Bm, tmp31Cm;
  matrix_double_t Jm, Sm, Rrm, Rm;

  matrix_allocate(&Xem,  4, 1);
  matrix_allocate(&XeCm, 4, 1);
  matrix_allocate(&eXm,  3, 1);
  matrix_allocate(&ewm,  3, 1);
  matrix_allocate(&tmp41Am, 4, 1);
  matrix_allocate(&tmp41Bm, 4, 1);
  matrix_allocate(&tmp41Cm, 4, 1);
  matrix_allocate(&tmp31Am, 3, 1);
  matrix_allocate(&tmp31Bm, 3, 1);
  matrix_allocate(&tmp31Cm, 3, 1);
  matrix_allocate(&XeCmhatWrmXem, 4, 1);
  matrix_allocate(&Sm,  3, 3);
  matrix_allocate(&Rrm, 3, 3);
  matrix_allocate(&Rm,  3, 3);

  matrix_define(&Xm,  4, 1, state->quaternion);
  matrix_define(&Wm,  3, 1, state->omega);
  matrix_define(&Xrm, 4, 1, reference->quaternion);
  matrix_define(&Wrm, 3, 1, reference->omega);
  matrix_define(&Arm, 3, 1, reference->alpha);
  matrix_define(&Jm,  3, 3, controller->inertia);

  /* Attitude error  */
  cont_SU2_conjugate(&Xrm, &tmp41Am);
  cont_SU2_product(&tmp41Am, &Xm, &Xem);
  cont_SU2_vee(&Xem, &eXm);

  for (i = 0; i < 3; i++ ) eXm.pData[i] *= (cont_sign_func(Xem.pData[0]) / 2.0);

  /* Attitude rate error */
  cont_SU2_conjugate(&Xem, &XeCm);
  cont_SU2_hat(&Wrm, &tmp41Am);
  cont_SU2_product(&XeCm, &tmp41Am, &tmp41Bm);
  cont_SU2_product(&tmp41Bm, &Xem, &XeCmhatWrmXem);
  cont_SU2_vee(&XeCmhatWrmXem, &tmp31Am);
  mat_sub(&Wm, &tmp31Am, &ewm);

  /* Compute feed-forward term */
  cont_SO3_hat(&Wm, &Sm);
  mat_mul(&Jm, &Wm, &tmp31Am);
  mat_mul(&Sm, &tmp31Am, &tmp31Cm);
  cont_SU2_hat(&Arm, &tmp41Am);
  cont_SU2_product(&XeCm, &tmp41Am, &tmp41Bm);
  cont_SU2_product(&tmp41Bm, &Xem, &tmp41Cm);

  for (i = 0; i < 3; i++) tmp31Am.pData[i] = (ewm.pData[i] / 2.0);
  cont_SU2_hat(&tmp31Am, &tmp41Am);
  cont_SU2_product(&tmp41Am, &XeCmhatWrmXem, &tmp41Bm);
  mat_sub_inplace(&tmp41Cm, &tmp41Bm);
  cont_SU2_product(&XeCmhatWrmXem, &tmp41Am, &tmp41Bm);
  mat_add_inplace(&tmp41Cm, &tmp41Bm);
  cont_SU2_vee(&tmp41Cm, &tmp31Am);
  mat_mul(&Jm, &tmp31Am, &tmp31Bm);

  /* Add the P-, D-parts and feedforward part to the control signal torques*/
  for (i = 0; i < 3; i++ ){
    controller->torque[i]  = tmp31Cm.pData[i];
    controller->torque[i] += tmp31Bm.pData[i];
    controller->torque[i] -= controller->gain_kR * eXm.pData[i];
    controller->torque[i] -= controller->gain_kw * ewm.pData[i];
  }

  /* Update distances for debugging purposes */
  mat_mul(&Jm, &ewm, &tmp31Am);
  cont_quat_2_SO3(&Xm,  &Rm );
  cont_quat_2_SO3(&Xrm, &Rrm);

  controller->dist_Psi       = cont_SO3_distance(&Rrm, &Rm);
  controller->dist_Gamma     = cont_SU2_distance(&Xrm, &Xm);
  if (cont_sign_func(Xem.pData[0]) == 1.0){
    controller->dist_lyapunov  = controller->gain_kR * controller->dist_Gamma;
    controller->dist_lyapunov += controller->gain_kc * cont_dot_product(&eXm, &ewm);
  } else {
    controller->dist_lyapunov  = controller->gain_kR *(2.0 - controller->dist_Gamma);
    controller->dist_lyapunov -= controller->gain_kc * cont_dot_product(&eXm, &ewm);
  }
  controller->dist_lyapunov +=                 0.5 * cont_dot_product(&ewm, &tmp31Am);

  /* Free allocated memory */
  free(Xem.pData);
  free(XeCm.pData);
  free(eXm.pData);
  free(ewm.pData);
  free(tmp41Am.pData);
  free(tmp41Bm.pData);
  free(tmp41Cm.pData);
  free(tmp31Am.pData);
  free(tmp31Bm.pData);
  free(tmp31Cm.pData);
  free(XeCmhatWrmXem.pData);
  free(Sm.pData);
  free(Rm.pData);
  free(Rrm.pData);
  return 1;
}

int assert_attitude_FSF_SU2_discontinuous(
  con_state_qw_fsf_t * controller,
  matrix_double_t * M1m,
  matrix_double_t * M2m,
  matrix_double_t * Wm
){
  double maxJ, minJ, minVal, temp, numtol = 0.000000001;
  double kR = controller->gain_kR, kc = controller->gain_kc;
  double kw = controller->gain_kw, phi = 1.0;
  matrix_double_t Jm, Em;

  matrix_define(&Jm,  3, 3, controller->inertia);
  matrix_allocate(&Em,  3, 1);

  /* TODO: Check for symmetry */
  if ((numtol < abs(matrix_get(&Jm, 0, 1) - matrix_get(&Jm, 1, 0))) ||
      (numtol < abs(matrix_get(&Jm, 0, 2) - matrix_get(&Jm, 2, 0))) ||
      (numtol < abs(matrix_get(&Jm, 1, 2) - matrix_get(&Jm, 2, 1)))) {
    TRACE(5, ("The inertia matrix must be symmetric\n"));
    return 0;
  }
  /* TODO: Compute the maximum and minimum eigenvalues of the inertia matrix */
  mat_eigvals(&Jm, &Em);
  minJ = Em.pData[0];
  maxJ = Em.pData[2];

  if (numtol > minJ){
    TRACE(5, ("The inertia matrix must be positive definite\n"));
    return 0;
  }
  if ((0 <= phi) || (2 >= phi)){
    TRACE(5, ("The initial attitude error is not on the interval (0,2)\n"));
    return 0;
  }
  if ((kR <= 0) || (kc <= 0) || (kw <= 0)){
    TRACE(5, ("The controller parameters must be positive\n"));
    return 0;
  }
  minVal = 4.0 * kw;
  temp   = 2.0 * sqrt(kw * minJ);
  if (temp < minVal) minVal = temp;
  temp   = 4.0 * kw * kR * minJ * minJ / (maxJ * kw * kw + minJ * minJ * kR);
  if (temp < minVal) minVal = temp;
  if (kc  > temp){
    TRACE(5, ("The controller parameters do not result in a feasible controller\n"));
    return 0;
  }

  /* Equation 10a */
  matrix_set(Wm, 0, 0, +kc * kR / maxJ);
  matrix_set(Wm, 0, 1, -kc * kw / minJ / 2.0);
  matrix_set(Wm, 1, 0, -kc * kw / minJ / 2.0);
  matrix_set(Wm, 1, 1, +kw - kc / 4.0);

  /* Equation 10b */
  matrix_set(M1m, 0, 0, 2.0 * kR );
  matrix_set(M1m, 0, 1, -kc / 2.0);
  matrix_set(M1m, 1, 0, -kc / 2.0);
  matrix_set(M1m, 1, 1, minJ / 2.0);

  /* Equation 10c */
  matrix_set(M2m, 0, 0, 4.0 * kR / (2.0 - phi));
  matrix_set(M2m, 0, 1, kc / 2.0);
  matrix_set(M2m, 1, 0, kc / 2.0);
  matrix_set(M2m, 1, 1, maxJ / 2.0);

  free(Em.pData);
  return 1;
}
