/*******************************************************************************
* Copyright (C) Marcus Greiff 2020
*
* @file cont_math.c
* @author Marcus Greiff
* @date June 2020
*******************************************************************************/
#include "cont_math.h"
/* Maps pertaining to the SU2 manifold */
int cont_SU2_conjugate(
  matrix_double_t * qin,
  matrix_double_t * qout
){
  int i;
  /* Input error check */
  if((qin->numCols  != 1) ||
     (qin->numRows  != 4) ||
     (qout->numCols != 1) ||
     (qout->numRows != 4)){
    TRACE(5, ("Bad input dimensions\n"));
    return 0;
  }
  qout->pData[0] = qin->pData[0];
  for (i = 1; i < 4; i++) qout->pData[i] = -qin->pData[i];
  return 1;
}

int cont_SU2_product(
  matrix_double_t * q,
  matrix_double_t * p,
  matrix_double_t * out
){
  if((q->numCols != 1) ||
     (q->numRows != 4) ||
     (p->numCols != 1) ||
     (p->numRows != 4) ||
     (out->numCols != 1) ||
     (out->numRows != 4)){
    TRACE(5, ("Bad input dimensions\n"));
    return 0;
  }
  out->pData[0]  =   p->pData[0]*q->pData[0] - p->pData[1]*q->pData[1];
  out->pData[0] += - p->pData[2]*q->pData[2] - p->pData[3]*q->pData[3];
  out->pData[1]  =   p->pData[0]*q->pData[1] + p->pData[1]*q->pData[0];
  out->pData[1] += - p->pData[2]*q->pData[3] + p->pData[3]*q->pData[2];
  out->pData[2]  =   p->pData[0]*q->pData[2] + p->pData[2]*q->pData[0];
  out->pData[2] +=   p->pData[1]*q->pData[3] - p->pData[3]*q->pData[1];
  out->pData[3]  =   p->pData[0]*q->pData[3] - p->pData[1]*q->pData[2];
  out->pData[3] +=   p->pData[2]*q->pData[1] + p->pData[3]*q->pData[0];
  return 1;
}

int cont_SU2_vee(
  matrix_double_t * in,
  matrix_double_t * out
){
  int i;
  if((in->numCols != 1) ||
     (in->numRows != 4) ||
     (out->numCols != 1) ||
     (out->numRows != 3)){
    TRACE(5, ("Bad input dimensions\n"));
    return 0;
  }
  for (i = 0; i < 3; i++) out->pData[i] = in->pData[i+1];
  return 1;
}

int cont_SU2_hat(
  matrix_double_t * in,
  matrix_double_t * out
){
  int i;
  if((in->numCols != 1) ||
     (in->numRows != 3) ||
     (out->numCols != 1) ||
     (out->numRows != 4)){
    TRACE(5, ("Bad input dimensions\n"));
    return 0;
  }
  out->pData[0] = 0.0;
  for (i = 0; i < 3; i++) out->pData[i+1] = in->pData[i];
  return 1;
}

int cont_SU2_Exp(
  matrix_double_t * in,
  matrix_double_t * out
){
  int i;
  double theta = 0.0;
  double coeff_im;

  if((in->numCols != 1) ||
     (in->numRows != 3) ||
     (out->numCols != 1) ||
     (out->numRows != 4)){
    TRACE(5, ("Bad input dimensions\n"));
    return 0;
  }

  for (i = 0; i < 3; i++) theta += in->pData[i]*in->pData[i];
  theta = sqrt(theta);
  coeff_im = cont_sinc(theta / 2.0);
  out->pData[0] = cos(theta / 2.0);
  for (i = 0; i < 3; i++) out->pData[i+1] = coeff_im * in->pData[i] / 2.0;
  return 1;
}

int cont_SU2_triple(
  double a,
  double b,
  double c,
  matrix_double_t * out
){
  double cosad2, cosbd2, coscd2, sinad2, sinbd2, sincd2;
  if((out->numCols != 1) ||
     (out->numRows != 4)){
    TRACE(5, ("Bad input dimensions\n"));
    return 0;
  }

  cosad2 = cos(a/2.0);
  cosbd2 = cos(b/2.0);
  coscd2 = cos(c/2.0);
  sinad2 = sin(a/2.0);
  sinbd2 = sin(b/2.0);
  sincd2 = sin(c/2.0);

  out->pData[0] = cosad2*cosbd2*coscd2 - sinad2*sinbd2*sincd2;
  out->pData[1] = cosbd2*coscd2*sinad2 + cosad2*sinbd2*sincd2;
  out->pData[2] = cosad2*coscd2*sinbd2 - cosbd2*sinad2*sincd2;
  out->pData[3] = cosad2*cosbd2*sincd2 + coscd2*sinad2*sinbd2;

  return 1;
}

double cont_SU2_distance(
  matrix_double_t * q1,
  matrix_double_t * q2
){
  int i;
  double out = 1.0;
  /* Input error check */
  if((q1->numCols  != 1) ||
     (q1->numRows  != 4) ||
     (q2->numCols != 1) ||
     (q2->numRows != 4)){
    TRACE(5, ("Bad input dimensions\n"));
    return 0;
  }
  for (i = 0; i < 4; i++ ) out -= (q1->pData[i] * q2->pData[i]);
  return out;
}

/* Embedding relating SU(2) to SO(3) */
int cont_quat_2_SO3(
  matrix_double_t * q,
  matrix_double_t * R
){
  double a, b, c, d, asq, bsq, csq, dsq;

  /* Input error check */
  if((q->numCols != 1) ||
     (q->numRows != 4) ||
     (R->numCols != 3) ||
     (R->numRows != 3)){
    TRACE(5, ("Bad input dimensions\n"));
    return 0;
  }
  a = q->pData[0]; b = q->pData[1]; c = q->pData[2]; d = q->pData[3];
  asq = a * a; bsq = b * b; csq = c * c; dsq = d * d;
  matrix_set(R, 0, 0, asq+bsq-csq-dsq);
  matrix_set(R, 1, 0, 2.0*(b*c+a*d));
  matrix_set(R, 2, 0, 2.0*(b*d-a*c));
  matrix_set(R, 0, 1, 2.0*(b*c-a*d));
  matrix_set(R, 1, 1, asq-bsq+csq-dsq);
  matrix_set(R, 2, 1, 2.0*(c*d+a*b));
  matrix_set(R, 0, 2, 2.0*(b*d+a*c));
  matrix_set(R, 1, 2, 2.0*(c*d-a*b));
  matrix_set(R, 2, 2, asq-bsq-csq+dsq);
  return 1;
}

/* Maps pertaining to the SO3 manifold */
int cont_SO3_hat(
  matrix_double_t * in,
  matrix_double_t * out
){
  if((in->numCols != 1) ||
     (in->numRows != 3) ||
     (out->numCols != 3) ||
     (out->numRows != 3)){
    TRACE(5, ("Bad input dimensions\n"));
    return 0;
  }
  matrix_zero(out);
  matrix_set(out, 0, 1, -in->pData[2]);
  matrix_set(out, 0, 2, +in->pData[1]);
  matrix_set(out, 1, 0, +in->pData[2]);
  matrix_set(out, 1, 2, -in->pData[0]);
  matrix_set(out, 2, 0, -in->pData[1]);
  matrix_set(out, 2, 1, +in->pData[0]);
  return 1;
}

int cont_SO3_vee(
  matrix_double_t * in,
  matrix_double_t * out
){
  if((in->numCols != 3) ||
     (in->numRows != 3) ||
     (out->numCols != 1) ||
     (out->numRows != 3)){
    TRACE(5, ("Bad input dimensions\n"));
    return 0;
  }
  out->pData[0] = matrix_get(in, 2, 1);
  out->pData[1] = matrix_get(in, 0, 2);
  out->pData[2] = matrix_get(in, 1, 0);
  return 1;
}

int cont_SO3_Exp(
  matrix_double_t * in,
  matrix_double_t * out
){
  int i;
  matrix_double_t Am, Bm;
  double coeff_a, coeff_b, theta = 0.0;

  if((in->numCols  != 1) ||
     (in->numRows  != 3) ||
     (out->numCols != 3) ||
     (out->numRows != 3)){
    TRACE(5, ("Bad input dimensions\n"));
    return 0;
  }

  matrix_allocate(&Am, 3, 3);
  matrix_allocate(&Bm, 3, 3);

  for (i = 0; i < 3; i++ ) theta += pow(in->pData[i], 2.0);
  theta = sqrt(theta);

  /* Compute coefficients using the sinc function */
  coeff_a = cont_sinc(theta);
  coeff_b = cont_sinc(theta / 2.0);
  coeff_b = coeff_b*coeff_b / 2.0;

  /* Compute the skew symetric matrix and its power */
  cont_SO3_hat(in, &Am);
  mat_mul(&Am, &Am, &Bm);

  /* out = I + sinc(theta) * S(in) + sinc(theta/2)^2/2 * S(in) * S(in) */
  matrix_identity(out);
  for (i = 0; i < 9;  i++ ) out->pData[i] += (coeff_a * Am.pData[i] + coeff_b * Bm.pData[i]);

  /* Free malloced memory */
  free(Am.pData);
  free(Bm.pData);
  return 1;
}

int cont_SO3_Log(
  matrix_double_t * in,
  matrix_double_t * out
){
  int i;
  matrix_double_t RTm;
  double coeff, theta = 0;

  if((in->numCols  != 3) ||
     (in->numRows  != 3) ||
     (out->numCols != 1) ||
     (out->numRows != 3)){
    TRACE(5, ("Bad input dimensions\n"));
    return 0;
  }
  /* TODO this is currently not safe, as the sinc funciton may dip down to 0
     yielding a division by zero. This  can be fixed at a later time, currently
     not done here in order to keep everything clean in the initial implementation*/
  matrix_allocate(&RTm, 3, 3);
  mat_trans(in, &RTm);
  mat_sub_inplace(&RTm, in);  /* RT = -(R - R'); */

  for (i = 0; i < 3; i++ ) theta += matrix_get(in, i, i);
  theta /=2;
  theta -=0.5;
  theta = acos(theta);
  coeff = 0.5 / cont_sinc(theta);

  cont_SO3_vee(&RTm, out);

  /* out = [R - R']^v / (2*sinc(theta)), theta = acos(trace(R) -  1) / 2 */
  for (i = 0; i < 3; i++) out->pData[i] *= -coeff;
  free(RTm.pData);

  return 1;
}

double cont_SO3_distance(
  matrix_double_t * R1,
  matrix_double_t * R2
){

  matrix_double_t R1T, R1TR2;
  double psi;

  /* Input error check */
  if((R1->numCols != 3) ||
     (R1->numRows != 3) ||
     (R2->numCols != 3) ||
     (R2->numRows != 3)){
    TRACE(5, ("Bad input dimensions\n"));
    assert(1 == 0);
  }

  matrix_allocate(   &R1T, 3, 3);
  matrix_allocate( &R1TR2, 3, 3);
  mat_trans(R1, &R1T);
  mat_mul(&R1T, R2,  &R1TR2);

  /* Compute Psi = (1/2)*trace(I - R1'*R2) */
  psi = (3.0 - R1TR2.pData[0] - R1TR2.pData[4] -  R1TR2.pData[8]) / 2.0;

  free(R1T.pData);
  free(R1TR2.pData);

  return psi;
}

double cont_sinc(
  double x
){
  if (x == 0.0) {
    return 1.0;
  } else {
    return sin(x) / x;
  }
}

double cont_dot_product(
  matrix_double_t * vecA,
  matrix_double_t * vecB
){
  int i;
  double output = 0.0;
  if((vecA->numCols != 1) ||
     (vecA->numRows != vecB->numRows) ||
     (vecB->numCols != 1)){
    TRACE(5, ("Bad input dimensions\n"));
    assert(1 == 0);
  }
  for (i = 0; i < vecA->numRows; i++) output += (vecA->pData[i]*vecB->pData[i]);
  return output;
}

void cont_cross_product(
  matrix_double_t * inAm,
  matrix_double_t * inBm,
  matrix_double_t * outm
){
  outm->pData[0] = inAm->pData[1]*inBm->pData[2] - inAm->pData[2]*inBm->pData[1];
  outm->pData[1] = inAm->pData[2]*inBm->pData[0] - inAm->pData[0]*inBm->pData[2];
  outm->pData[2] = inAm->pData[0]*inBm->pData[1] - inAm->pData[1]*inBm->pData[0];
  return;
}

double cont_sign_func(
  double in
){
  if (in >= 0.0) return 1.0;
  return -1.0;
}


int cont_normalize(matrix_double_t * vec) {
  int N, i;
  double val = 0.0;

  if ((vec->numRows != 1) &&  (vec->numCols != 1)) return 0;
  if (vec->numRows == 1){
    if (vec->numCols < 1) return 0;
    N = vec->numCols;
  }
  if (vec->numCols == 1){
    if (vec->numRows < 1) return 0;
    N = vec->numRows;
  }

  for ( i = 0; i < N; i++) val += (vec->pData[i]*vec->pData[i]);
  if (val <= 0.0) return 0;
  val = sqrt(val);
  for ( i = 0; i < N; i++) vec->pData[i] /= val;
  return 1;
}
