/*******************************************************************************
* Copyright (C) Marcus Greiff 2020
*
* @file cont_attitude_reference_generator.c
* @author Marcus Greiff
* @date June 2020
*******************************************************************************/
#include "cont_attitude_reference_generator.h"

int update_attitude_references(
  matrix_double_t * commands,
  matrix_double_t * mem,
  ref_state_qw_t * reference,
  double h
){
  int i, j;
  double tmp[3], Bm[3];
  matrix_double_t tmp31m, tmp41Am, tmp41Bm, tmp41Cm, tmp91m, Xm, Am, Qrm, Wrm, Arm, dQrm, ddQrm;

  double p   = 4.0;
  double psq = p*p;
  double pcb = psq*p;
  double hp  = h*p;
  double hsq = h*h;
  double a = exp(-h * p);

  matrix_allocate(&Am,      3, 3);
  matrix_allocate(&tmp41Am, 4, 1);
  matrix_allocate(&tmp41Bm, 4, 1);
  matrix_allocate(&tmp41Cm, 4, 1);
  matrix_allocate(&dQrm,    4, 1);
  matrix_allocate(&ddQrm,   4, 1);

  matrix_set(&Am, 0, 0,  (a*(hsq*psq + 2.0*hp + 2.0))/2.0);
  matrix_set(&Am, 0, 1,                    h*a*(hp + 1.0));
  matrix_set(&Am, 0, 2,                       (hsq*a)/2.0);
  matrix_set(&Am, 1, 0,                  -(hp*hp*p*a)/2.0);
  matrix_set(&Am, 1, 1,          a*(- hsq*psq + hp + 1.0));
  matrix_set(&Am, 1, 2,             -(h*a*(hp - 2.0))/2.0);
  matrix_set(&Am, 2, 0,         (hp*psq*a*(hp - 2.0))/2.0);
  matrix_set(&Am, 2, 1,                  h*psq*a*(hp - 3));
  matrix_set(&Am, 2, 2, (a*(hsq*psq - 4.0*h*p + 2.0))/2.0);
  Bm[0] = 1.0 - (hsq*psq*a)/2.0 - hp*a - a;
  Bm[1] =                  (hsq*pcb*a)/2.0;
  Bm[2] =       -(hp*psq*a*(hp - 2.0))/2.0;

  for (i = 0; i < 4; i++){
    matrix_define(&Xm,     3, 1, &mem->pData[i*3]);
    matrix_define(&tmp31m, 3, 1, tmp);
    mat_mul(&Am, &Xm, &tmp31m);
    for (j = 0; j < 3; j++ ) mem->pData[i*3 + j] = tmp[j] + Bm[j] * commands->pData[i];
  }

  matrix_define(&Qrm, 4, 1, reference->quaternion);
  matrix_define(&Wrm, 3, 1, reference->omega);
  matrix_define(&Arm, 3, 1, reference->alpha);

  matrix_define(&tmp91m, 9, 1, &mem->pData[3]);

  /* Compute the reference quaternion */
  cont_SU2_triple_derivatives(&tmp91m, 0, 0, 0, &Qrm);

  /* Compute the reference quaternion time terivative */
  cont_SU2_triple_derivatives(&tmp91m, 1, 0, 0, &dQrm);
  cont_SU2_triple_derivatives(&tmp91m, 0, 1, 0, &tmp41Am);
  mat_add_inplace(&dQrm, &tmp41Am);
  cont_SU2_triple_derivatives(&tmp91m, 0, 0, 1, &tmp41Am);
  mat_add_inplace(&dQrm, &tmp41Am);

  /* Compute the reference quaternion second time derivative */
  cont_SU2_triple_derivatives(&tmp91m, 2, 0, 0, &ddQrm);
  cont_SU2_triple_derivatives(&tmp91m, 0, 2, 0, &tmp41Am);
  mat_add_inplace(&ddQrm, &tmp41Am);
  cont_SU2_triple_derivatives(&tmp91m, 0, 0, 2, &tmp41Am);
  mat_add_inplace(&ddQrm, &tmp41Am);
  cont_SU2_triple_derivatives(&tmp91m, 1, 1, 0, &tmp41Am);
  mat_add_inplace(&ddQrm, &tmp41Am);
  mat_add_inplace(&ddQrm, &tmp41Am);
  cont_SU2_triple_derivatives(&tmp91m, 1, 0, 1, &tmp41Am);
  mat_add_inplace(&ddQrm, &tmp41Am);
  mat_add_inplace(&ddQrm, &tmp41Am);
  cont_SU2_triple_derivatives(&tmp91m, 0, 1, 1, &tmp41Am);
  mat_add_inplace(&ddQrm, &tmp41Am);
  mat_add_inplace(&ddQrm, &tmp41Am);

  /* Ang. rate ref. w = 2 * [conj(q) * dq]^{vee}  */
  cont_SU2_conjugate(&Qrm, &tmp41Am);
  cont_SU2_product(&tmp41Am, &dQrm, &tmp41Bm);
  for (i = 0; i < 3; i++ ) Wrm.pData[i] = 2.0*tmp41Bm.pData[i+1];

  /* Ang. acc. ref. a = [conj(q) * (2*ddq - dq * [w]^{hat})]^{vee} */
  cont_SU2_hat(&Wrm, &tmp41Am);
  cont_SU2_product(&dQrm, &tmp41Am, &tmp41Bm);
  for (i = 0; i < 4; i++ ) tmp41Am.pData[i] = 2.0*ddQrm.pData[i] - tmp41Bm.pData[i];
  cont_SU2_conjugate(&Qrm, &tmp41Bm);
  cont_SU2_product(&tmp41Bm, &tmp41Am, &tmp41Cm);
  cont_SU2_vee(&tmp41Cm, &Arm);

  /* Set the reference thrust to the filtered thrust */
  reference->thrust = CONT_MASS * CONT_GRAVACC * (mem->pData[0] / 2.0 + 0.5);

  /* Free malloced memory */
  free(Am.pData);
  free(tmp41Am.pData);
  free(tmp41Bm.pData);
  free(tmp41Cm.pData);
  free(dQrm.pData);
  free(ddQrm.pData);
  return 1;
}

/* TODO refactor this */
int cont_SU2_triple_derivatives(
  matrix_double_t * in,
  int derivA,
  int derivB,
  int derivC,
  matrix_double_t * out
){
  double ca, cb, cc, sa, sb, sc;
  matrix_double_t tmp31m;

  if((in->numCols != 1) ||
     (in->numRows != 9) ||
     (out->numCols != 1) ||
     (out->numRows != 4)){
    TRACE(5, ("Bad input dimensions\n"));
    return 0;
  }

  matrix_define(&tmp31m, 3, 1, &in->pData[0]);
  cont_get_cossin(&tmp31m, derivA, &ca, &sa);
  matrix_define(&tmp31m, 3, 1, &in->pData[3]);
  cont_get_cossin(&tmp31m, derivB, &cb, &sb);
  matrix_define(&tmp31m, 3, 1, &in->pData[6]);
  cont_get_cossin(&tmp31m, derivC, &cc, &sc);

  out->pData[0] = ca*cb*cc - sa*sb*sc;
  out->pData[1] = cb*cc*sa + ca*sb*sc;
  out->pData[2] = ca*cc*sb - cb*sa*sc;
  out->pData[3] = ca*cb*sc + cc*sa*sb;

  return 1;
}

/* TODO refactor this */
int cont_get_cossin(matrix_double_t * input, int deriv, double * ca, double * sa) {
  double a   = input->pData[0] / 2.0;
  double da  = input->pData[1] / 2.0;
  double dda = input->pData[2] / 2.0;

  *ca = 0.0;
  *sa = 0.0;
  if (deriv == 0) {
    *ca = cos(a);
    *sa = sin(a);
  } else if (deriv == 1){
    *ca = -sin(a)*da;
    *sa = +cos(a)*da;
  } else if (deriv == 2){
    *ca = -cos(a)*da*da - sin(a)*dda;
    *sa = -sin(a)*da*da + cos(a)*dda ;
  } else {
    TRACE(5, ("Bad input derivative order - must be {0,1,2}\n"));
    return 0;
  }
  return 1;
}
