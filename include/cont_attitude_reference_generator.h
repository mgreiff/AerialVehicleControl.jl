/*******************************************************************************
* Copyright (C) Marcus Greiff 2020
*
* @file cont_attitude_reference_generator.h
* @author Marcus Greiff
* @date June 2020
*******************************************************************************/
#ifndef __CONT_ATTITUDE_REFERENCE_GENERATOR_H__
#define __CONT_ATTITUDE_REFERENCE_GENERATOR_H__

#include "cont_math.h"

/***************************************************************************//**
* @brief  Update the reference trajectory based on input commands
*
* This function computes a set of references of the attitude control from a set
* of commands in anlges r, t and s. These angles are filtered through a
* critically damped second order system to yield smooth reference signals. Let
*
*     q_{p} = exp([p, 0, 0]^hat/2) = [ cos(p/2), sin(p/2), 0, 0]'
*     q_{t} = exp([0, t, 0]^hat/2) = [ cos(t/2), 0, sin(t/2), 0]'            (1)
*     q_{s} = exp([0, 0, s]^hat/2) = [ cos(s/2), 0, 0, sin(t/2)]'
*
* and take
*
*     q   = q_{p} * q_{t} * q_{s}                                            (2)
*
* where then
*
*     dq  = + dq_{p} *  q_{t}  *  q_{s}
*           +  q_{p} * dq_{t}  *  q_{s}                                      (3)
*           +  q_{p} *  q_{t}  * dq_{s}
*
* and
*
*     ddq = +     ddq_{p} *   q_{t}  *   q_{s}
*           +       q_{p} * ddq_{t}  *   q_{s}
*           +       q_{p} *   q_{t}  * ddq_{s}                               (4)
*           + 2 * (dq_{p} *  dq_{t}  *   q_{s})
*           + 2 * (dq_{p} *   q_{t}  *  dq_{s})
*           + 2 * ( q_{p} *  dq_{t}  *  dq_{s})
*
* From these objects, the attitude rate can be computed through the kinematics
*
*     dq = q * [w]^{hat} / 2                                                 (5)
*
* The attitude rates can be computed as
*
*     w = 2 * [conj(q) * dq]^{vee}                                           (6)
*
* And since
*
*     ddq = (dq * [w]^{hat} + q * [a]^{hat}) / 2,                            (7)
*
* the atittude accelerations can be computed through
*
*     a = [conj(q) * (2*ddq - dq * [w]^{hat})]^{vee},                        (8)
*
* For the filtering of the inputs, we consider a second order critically damped
* system with a double pole in -p, characterized by the LTI system
*
*     dx(t) = A*x(t) + B*u(t)                                                (9)
*
* where
*
*     A = [    0,     1,   0]    B = [   0]
*         [    0,     0,   1]        [   0]
*         [ -p^3, -3p^2, -3p]        [ p^3]
*
* when discretized usign zero order hold at a time-step of h, we obtain
*
*     x(k+h) = A(h;p)*x(k) + B(h;p)*u(k),                                   (10)
*
* were
*
*     A(h;p) = [ (exp(-h*p)*(h^2*p^2 + 2*h*p + 2))/2,           h*exp(-h*p)*(h*p + 1),                   (h^2*exp(-h*p))/2]
*              [              -(h^2*p^3*exp(-h*p))/2, exp(-h*p)*(- h^2*p^2 + h*p + 1),          -(h*exp(-h*p)*(h*p - 2))/2]
*              [       (h*p^3*exp(-h*p)*(h*p - 2))/2,       h*p^2*exp(-h*p)*(h*p - 3), (exp(-h*p)*(h^2*p^2 - 4*h*p + 2))/2]
*
*     B(h;p) = [ -(exp(-h*p)*(2*h*p - 2*exp(h*p) + h^2*p^2 + 2))/2]
*              [                             (h^2*p^3*exp(-h*p))/2]
*              [                    -(h*p^3*exp(-h*p)*(h*p - 2))/2]
*
* For this system, the positional states y = [1, 0, 0]x are bound to the domain
* of the input u, and the magnitude of the velocities and accelerations are
* given by the pole location p, being greater for faster dynamics (larger p).
*
* In this implementation, we take the signals r s t and f to be configured on
* rectangular intervals r in [-r0, r0], use the second order dynamics in (9) to
* generate smooth reference trajectotries in the filter memory "mem", as
*
*     mem = (r, dr, ddr, s, ds, dds, t, dt, ddt, d, df, ddf)
*
* which are subsequently converted to a refeence quaternion and attitude rates
* through equations (1)-(8). The commended force is given in the global
* z-direction such that the default (r, s, t, f) = (0, 0, 0, 0) results in a
* stable hovering state, and that an excitation of the angles at a force command
* of f = 0 will keep the system at roughly the same height.
*
* @param[in] reference - The desired normalized thrust on an interval [-1,1]
* @param[in] commands  - The desired commands in normalized angles and thrust
* @param[in] dt        - The time-step between updates in the reference gen.
* @return status       - 1 if successful, 0 otherwise.
*******************************************************************************/
int update_attitude_references(
  matrix_double_t * commands,
  matrix_double_t * filterMemory,
  ref_state_qw_t * reference,
  double h
);

/* TODO refactor and test this */
int cont_SU2_triple_derivatives(
  matrix_double_t * in,
  int derivA,
  int derivB,
  int derivC,
  matrix_double_t * out
);

/* TODO refactor and test this */
int cont_get_cossin(matrix_double_t * input, int deriv, double * ca, double * sa);
#endif /* __CONT_ATTITUDE_REFERENCE_GENERATOR_H__ */
