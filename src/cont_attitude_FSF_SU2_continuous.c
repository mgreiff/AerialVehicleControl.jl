/*******************************************************************************
* Copyright (C) Marcus Greiff 2020
*
* @file cont_attitude_FSF_SU2_continuous.c
* @author Marcus Greiff
* @date June 2020
*******************************************************************************/
#include "cont_attitude_FSF_SU2_continuous.h"

int update_attitude_FSF_SU2_continuous(
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
  for (i = 0; i < 3; i++ ) eXm.pData[i] /= 2.0;

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
  controller->dist_lyapunov  = controller->gain_kR * controller->dist_Gamma;
  controller->dist_lyapunov += controller->gain_kc * cont_dot_product(&eXm, &ewm);
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
