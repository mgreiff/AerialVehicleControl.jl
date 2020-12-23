/*******************************************************************************
* Copyright (C) Marcus Greiff 2020
*
* @file cont_attitude_FOF_SO3_continuous.c
* @author Marcus Greiff
* @date June 2020
*******************************************************************************/
#include "cont_attitude_FOF_SO3_continuous.h"

int update_attitude_FOF_SO3_continuous_measurements(
  con_state_qw_fof_t * controller,
  matrix_double_t * y0m,
  matrix_double_t * y1m,
  matrix_double_t * y2m,
  matrix_double_t * y3m
){
  int i;
  matrix_double_t Wm, Ym;

  if ((y0m->numCols != 1) ||
      (y1m->numCols != 1) ||
      (y1m->numCols != 1) ||
      (y1m->numCols != 1) ||
      (y1m->numRows != 3) ||
      (y1m->numRows != 3) ||
      (y1m->numRows != 3) ||
      (y1m->numRows != 3)) {
    TRACE(5, ("Wrong input dimensions\n"));
    return 0;
  }

  matrix_define(&Wm, 3, 1, controller->measuredGyrorates);
  matrix_define(&Ym, 3, 3, controller->measuredDirections);

  for (i = 0; i < 3; i++){
    matrix_set(&Wm, i, 0, y0m->pData[i]);
    matrix_set(&Ym, i, 0, y1m->pData[i]);
    matrix_set(&Ym, i, 1, y2m->pData[i]);
    matrix_set(&Ym, i, 2, y3m->pData[i]);
  }
  return 1;
}

int update_attitude_FOF_SO3_continuous(
  ref_state_qw_t     * reference,
  dyn_state_qw_t     * state,
  con_state_qw_fof_t * controller,
  double               dt
){

  int i;

  matrix_double_t Qrm, Wrm, Arm;
  matrix_double_t Qhatm, Whatm;
  matrix_double_t Tm, Jm, invJm, gain_ki_m, gain_Kwm, gain_Cwm;

  matrix_double_t Rrm, Rhatm;
  matrix_double_t Vm, Ym, Wm, Yrm, Yhatm;

  matrix_double_t Whatem, Wem, Wtildem;

  matrix_double_t DRm, DWm;

  matrix_double_t vecAm, vecBm, vecCm, vecDm;

  matrix_double_t tmp33m, tmp31Am, tmp31Bm;

  /* Define the reference trajectory terms */
  matrix_define(&Qrm, 4, 1, reference->quaternion);
  matrix_define(&Wrm, 3, 1, reference->omega);
  matrix_define(&Arm, 3, 1, reference->alpha);

  /* Compute the estimator state (here stored in the state object) */
  matrix_define(&Qhatm,  4, 1, state->quaternion);
  matrix_define(&Whatm,  3, 1, state->omega);

  /* Extract the controller parameters */
  matrix_define(&Tm,  3, 1, controller->torque);
  matrix_define(&Jm,  3, 3, controller->inertia);
  matrix_define(&invJm,  3, 3, controller->invinertia);
  matrix_define(&gain_ki_m, 1, 3, controller->gain_ki);
  matrix_define(&gain_Kwm,  3, 3, controller->gain_Kw);
  matrix_define(&gain_Cwm,  3, 3, controller->gain_Cw);

  /* Define the directions and the measured directions */
  matrix_define(&Vm, 3, 3, controller->globalDirections);
  matrix_define(&Ym, 3, 3, controller->measuredDirections);
  matrix_define(&Wm, 3, 1, controller->measuredGyrorates);

  /* Allocate the necessary memory */
  matrix_allocate(&Rrm,   3, 3);
  matrix_allocate(&Rhatm, 3, 3);
  matrix_allocate(&Yrm,   3, 3);
  matrix_allocate(&Yhatm, 3, 3);

  matrix_allocate(&DRm, 3, 1);
  matrix_allocate(&DWm, 3, 1);

  matrix_allocate(&vecAm, 3, 1);
  matrix_allocate(&vecBm, 3, 1);
  matrix_allocate(&vecCm, 3, 1);
  matrix_allocate(&vecDm, 3, 1);

  matrix_allocate(&tmp33m,  3, 3);
  matrix_allocate(&tmp31Am, 3, 1);
  matrix_allocate(&tmp31Bm, 3, 1);

  /* Convert quaternions to elements on SO(3) */
  cont_quat_2_SO3(&Qrm,   &Rrm);
  cont_quat_2_SO3(&Qhatm, &Rhatm );

  /* Compute the directions rotated by the refernce and estimate rotations */
  mat_trans(&Rrm, &tmp33m);
  mat_mul(&tmp33m, &Vm, &Yrm);
  mat_trans(&Rhatm, &tmp33m);
  mat_mul(&tmp33m, &Vm, &Yhatm);

  /* Compute the error terms {whate = wr - what; we = wr - w; wtilde = what - w; } */
  matrix_allocate(&Whatem,  3, 1);
  matrix_allocate(&Wem,     3, 1);
  matrix_allocate(&Wtildem, 3, 1);
  mat_sub(&Wrm,   &Whatm, &Whatem );
  mat_sub(&Wrm,   &Wm   , &Wem    );
  mat_sub(&Whatm, &Wm   , &Wtildem);

  /*****************************************************************************
   * Compute the control signal torques in four distict components, as
   *
   *   tau = +taur ...                      =A
   *         +S(J*whate)*wr...              =B
   *         +Kw*whate...                   =C
   *         +k1*S(Rr'*v1)*(Rhat'*v1)...
   *         +k2*S(Rr'*v2)*(Rhat'*v2)...    =D
   *         +k3*S(Rr'*v3)*(Rhat'*v3)
   ****************************************************************************/

  /****** Computation of A ******/
  /* vecA = J*dwr               */
  mat_mul(&Jm, &Arm, &vecAm);
  /* tmpA = J*wr                */
  mat_mul(&Jm, &Wrm, &tmp31Am);
  /* tmpB = S(J*wr)*wr          */
  cont_cross_product(&tmp31Am, &Wrm, &tmp31Bm);
  /* vecA = vecA - S(J*wr)*wr   */
  mat_sub_inplace(&vecAm, &tmp31Bm);

  /****** Computation of B ******/
  /* tmpA = J*whate             */
  mat_mul(&Jm, &Whatem, &tmp31Am);
  /* vecB = S(J*whate)*wr       */
  cont_cross_product(&tmp31Am, &Wrm, &vecBm);

  /****** Computation of C ******/
  /* vecC = Kw*whate)           */
  mat_mul(&gain_Kwm, &Whatem, &vecCm);

  /****** Computation of D ******/
  assert(1 == attitude_FOF_SO3_continuous_cross_terms(&Yrm, &Yhatm, &gain_ki_m, &vecDm));

  /* Form the control signal torque as tau = vecA + vecB + vecC + vecD        */
  for ( i = 0; i < 3; i++ ) Tm.pData[i] = vecAm.pData[i] + vecBm.pData[i] + vecCm.pData[i] + vecDm.pData[i];

  /*****************************************************************************
  * Compute the estimator rotation innovation terms, as
  *
  * DR = -cR * (+k1 * S(Rhat'*v1)*(Rr'*v1 + R'*v1)...
  *             +k2 * S(Rhat'*v2)*(Rr'*v2 + R'*v2)...
  *             +k3 * S(Rhat'*v3)*(Rr'*v3 + R'*v3))
  *****************************************************************************/
  mat_add(&Yrm, &Ym, &tmp33m);
  assert(1 == attitude_FOF_SO3_continuous_cross_terms(&Yhatm, &tmp33m, &gain_ki_m, &DRm));
  for (i = 0; i < 3; i++ ) DRm.pData[i] *= -controller->gain_cR;

  /*****************************************************************************
  * Compute the estimator rotation innovation terms, as
  *
  * DW = -cw * J * S(wr)*we...
  *      -cw * Kw * we...
  *      -Cw * wtilde;
  *****************************************************************************/
  matrix_zero(&DWm);
  /* tmpA = S(wr)*we */
  cont_cross_product(&Wrm, &Wem, &tmp31Am);
  /* DW   = J*S(wr)*we */
  mat_mul(&Jm, &tmp31Am, &DWm);
  /* tmpA = Kw * we  */
  mat_mul(&gain_Kwm, &Wem, &tmp31Am);
  /* DW   = J*S(wr)*we + Kw * we */
  mat_add_inplace(&DWm, &tmp31Am);
  /* DW   = -cw * (J*S(wr)*we + Kw * we) */
  for (i = 0; i < 3; i++ ) DWm.pData[i] *= -controller->gain_cw;
  /* tmpA = Cw * wtilde */
  mat_mul(&gain_Cwm, &Wtildem, &tmp31Am);
  /* DW   = -cw * (J*S(wr)*we + Kw * we) - Cw * wtilde */
  mat_sub_inplace(&DWm, &tmp31Am);

  /*****************************************************************************
  * Simulate the observer one time-step using a first order euler method, with
  *
  * Rhatdot = Rhat*S(w + DR);
  * whatdot = (J)\(S(J*w)*w + tau + DW);
  *
  * Although, the simulation of the attitude will be done on SU(2).
  *****************************************************************************/
  /* tmpB = J*w */
  mat_mul(&Jm, &Wm, &tmp31Bm);
  /* tmpA = S(J*w)*w */
  cont_cross_product(&tmp31Bm, &Wm, &tmp31Am);
  /* tmpA = S(J*w)*w + tau */
  mat_add_inplace(&tmp31Am, &Tm);
  /* tmpA = S(J*w)*w + tau + DW */
  mat_add_inplace(&tmp31Am, &DWm);
  /* tmpB = J\(S(J*w)*w + tau + DW) */
  mat_mul(&invJm, &tmp31Am, &tmp31Bm);
  /* tmpA = w + DR */
  mat_add(&Wm, &DRm, &tmp31Am);
  /* Simulate the dynamics one time-step and project the attitude onto SU(2)  */
  assert(1 == attitude_FOF_SO3_continuous_simulate_observer(&Qhatm, &Whatm, &tmp31Am, &tmp31Bm, dt));

  /* Free the allocated memory */
  free(Rrm.pData);
  free(Rhatm.pData);
  free(Yrm.pData);
  free(Yhatm.pData);

  free(DRm.pData);
  free(DWm.pData);

  free(Whatem.pData);
  free(Wem.pData);
  free(Wtildem.pData);

  free(tmp33m.pData);
  free(tmp31Am.pData);
  free(tmp31Bm.pData);

  free(vecAm.pData);
  free(vecBm.pData);
  free(vecCm.pData);
  free(vecDm.pData);

  return 1;
}

int attitude_FOF_SO3_continuous_simulate_observer(
  matrix_double_t * Qm,
  matrix_double_t * Wm,
  matrix_double_t * deltaQm,
  matrix_double_t * deltaWm,
  double dt
){

  int i;
  matrix_double_t tmp31m, tmp41Am, tmp41Bm;

  if ((Qm->numRows != 4) ||
      (Qm->numCols != 1) ||
      (Wm->numRows != 3) ||
      (Wm->numCols != 1) ||
      (deltaQm->numRows != 3) ||
      (deltaQm->numCols != 1) ||
      (deltaWm->numRows != 3) ||
      (deltaWm->numCols != 1) ||
      (dt < 0)) {
    TRACE(5, ("Bad inputs to the FOF attitude observer simulator\n"));
    return 0;
  }

  matrix_allocate(&tmp31m,  3, 1);
  matrix_allocate(&tmp41Am, 4, 1);
  matrix_allocate(&tmp41Bm, 4, 1);

  /* Simulate the attitude observer */
  for (i = 0; i < 3; i++ ) tmp31m.pData[i] = dt * deltaQm->pData[i];
  if(0 == cont_SU2_Exp(&tmp31m, &tmp41Am)) return 0;
  if(0 == cont_SU2_product(Qm, &tmp41Am, &tmp41Bm)) return 0;
  for (i = 0; i < 4; i++ ) Qm->pData[i] = tmp41Bm.pData[i];

  /* Normalize the resulting quaternion (should not be needed) given the
     implemeted numerical integration scheme, but is done anyway */
  cont_normalize(Qm);

  /* Simulate the rate observer */
  for (i = 0; i < 3; i++ ) Wm->pData[i] += dt * deltaWm->pData[i];

  /* free allocated memory */
  free(tmp31m.pData);
  free(tmp41Am.pData);
  free(tmp41Bm.pData);
  return 1;
}

int attitude_FOF_SO3_continuous_cross_terms(
  matrix_double_t * YAm,
  matrix_double_t * YBm,
  matrix_double_t * gainm,
  matrix_double_t * outm
){
  matrix_double_t tmp31Am, tmp31Bm, tmp31Cm;
  int i, j;

  if ((YAm->numRows != 3) || (YBm->numRows != 3) || (gainm->numRows != 1)) {
    TRACE(5, ("The directional matrices must have three rows, and the gains must have a single row\n"));
    return 0;
  }
  if ((YAm->numCols < 3) || (YBm->numCols < 3)) {
    TRACE(5, ("The directional matrices must have at least three columns\n"));
    return 0;
  }
  if ((YAm->numCols != YBm->numCols) || (YAm->numCols != gainm->numCols)) {
    TRACE(5, ("The directional matrices must have the same number of columns\n"));
    return 0;
  }

  /* Input error check */
  matrix_allocate(&tmp31Am, 3, 1);
  matrix_allocate(&tmp31Bm, 3, 1);
  matrix_allocate(&tmp31Cm, 3, 1);

  matrix_zero(outm);
  for ( i = 0; i < YAm->numCols; i++ ){
    /* Extract vectors and muliply the second vector with the associated gain */
    for ( j = 0; j < 3; j++ ){
      matrix_set(&tmp31Am, j, 0, matrix_get(YAm, j, i));
      matrix_set(&tmp31Bm, j, 0, gainm->pData[i]*matrix_get(YBm, j, i));
    }
    /* Take the cross product and add to the output */
    cont_cross_product(&tmp31Am, &tmp31Bm, &tmp31Cm);
    mat_add_inplace(outm, &tmp31Cm);
  }

  free(tmp31Am.pData);
  free(tmp31Bm.pData);
  free(tmp31Cm.pData);

  return 1;
}
