/*******************************************************************************
* Copyright (C) Marcus Greiff 2020
*
* @file cont_attitude_FSF_SO3_continuous.c
* @author Marcus Greiff
* @date June 2020
*******************************************************************************/
#include "cont_attitude_FSF_SO3_continuous.h"

int update_attitude_FSF_SO3_continuous(
  ref_state_qw_t * reference,
  dyn_state_qw_t * state,
  con_state_qw_fsf_t * controller
){
  int i;
  matrix_double_t Rm, Rrm, RrTm, Rem, ReTm;
  matrix_double_t Sm;
  matrix_double_t Ewm, ERm, tmp31Am, tmp31Bm, tmp31Cm;
  matrix_double_t Qm, Qrm, Wm, Wrm, Arm, Tm, Jm;

  /* Assert controller tuning feasibility */
  if (0 == assert_attitude_FSF_SO3(controller)) return 0;

  matrix_allocate(&Rm,      3, 3);
  matrix_allocate(&Rrm,     3, 3);
  matrix_allocate(&RrTm,    3, 3);
  matrix_allocate(&Rem,     3, 3);
  matrix_allocate(&ReTm,    3, 3);
  matrix_allocate(&Sm,      3, 3);
  matrix_allocate(&Ewm,     3, 1);
  matrix_allocate(&ERm,     3, 1);
  matrix_allocate(&tmp31Am, 3, 1);
  matrix_allocate(&tmp31Bm, 3, 1);
  matrix_allocate(&tmp31Cm, 3, 1);

  matrix_define(&Qm,  4, 1, state->quaternion);
  matrix_define(&Wm,  3, 1, state->omega);
  matrix_define(&Qrm, 4, 1, reference->quaternion);
  matrix_define(&Wrm, 3, 1, reference->omega);
  matrix_define(&Arm, 3, 1, reference->alpha);
  matrix_define(&Tm,  3, 1, controller->torque);
  matrix_define(&Jm,  3, 3, controller->inertia);

  /* Compute the attitude error  */
  cont_quat_2_SO3(&Qm,  &Rm );
  cont_quat_2_SO3(&Qrm, &Rrm);
  mat_trans(&Rrm, &RrTm);
  mat_mul(&RrTm, &Rm, &Rem);
  mat_trans(&Rem, &ReTm);
  mat_sub(&Rem, &ReTm, &Sm);
  cont_SO3_vee(&Sm, &ERm);
  for (i = 0; i < 3; i++) ERm.pData[i] /= 2.0;

  /* Compute the attitude rate error*/
  mat_mul(&ReTm, &Wrm, &tmp31Am);
  mat_sub(&Wm, &tmp31Am, &Ewm);

  /* Compute the feed-forward term  */
  cont_SO3_hat(&Wm, &Sm);
  mat_mul(&Jm, &Wm, &tmp31Am);
  mat_mul(&Sm, &tmp31Am, &tmp31Cm);
  mat_mul(&ReTm, &Wrm, &tmp31Am);
  mat_mul(&Sm, &tmp31Am, &tmp31Bm);
  mat_mul(&ReTm, &Arm, &tmp31Am);
  mat_sub_inplace(&tmp31Am, &tmp31Bm);
  mat_mul(&Jm, &tmp31Am, &tmp31Bm);

    /* Add the P-, D-parts and feedforward part to the control signal torques*/
  for (i = 0; i < 3; i++ ){
    controller->torque[i]  = tmp31Cm.pData[i];
    controller->torque[i] += tmp31Bm.pData[i];
    controller->torque[i] -= controller->gain_kR * ERm.pData[i];
    controller->torque[i] -= controller->gain_kw * Ewm.pData[i];
  }

  /* Update distances for debugging purposes */
  mat_mul(&Jm, &Ewm, &tmp31Am);
  controller->dist_Psi       = cont_SO3_distance(&Rrm, &Rm);
  controller->dist_Gamma     = cont_SU2_distance(&Qrm, &Qm);
  controller->dist_lyapunov  = controller->gain_kR * controller->dist_Psi;
  controller->dist_lyapunov += controller->gain_kc * cont_dot_product(&ERm, &Ewm);
  controller->dist_lyapunov +=                 0.5 * cont_dot_product(&Ewm, &tmp31Am);

  /* Free all allocated memory */
  free(Rm.pData);
  free(Rrm.pData);
  free(RrTm.pData);
  free(Rem.pData);
  free(ReTm.pData);
  free(Sm.pData);
  free(Ewm.pData);
  free(ERm.pData);
  free(tmp31Am.pData);
  free(tmp31Bm.pData);
  free(tmp31Cm.pData);

  return 1;
}
