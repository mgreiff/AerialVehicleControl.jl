#include "cont_main.h"

/*******************************************************************************
 * This is an example of how the attitude controllers could be implemented
 * used mainly for checking of memory leaks using valgrind.
 ******************************************************************************/
int main(int argc, char *argv[]) {
  /* Run the FSF attitude cotroller in closed loop for 100 time-steps */
  example_attitude_FSF();
  /* Run the FOF attitude cotroller in closed loop for 100 time-steps */
  example_attitude_FOF();
  return 0;
}

void example_attitude_FSF(void){

  int i;
  double dt;
  matrix_double_t Jm, commandm, filtermemm, dutycyclem;
  dyn_state_qw_t * state;
  con_state_qw_fsf_t * controller;
  ref_state_qw_t * reference;

  printf("Allocate memory for the reference, state and controller structs...\n");
  state      = (dyn_state_qw_t*)     calloc(1, sizeof(dyn_state_qw_t));
  controller = (con_state_qw_fsf_t*) calloc(1, sizeof(con_state_qw_fsf_t));
  reference  = (ref_state_qw_t*)     calloc(1, sizeof(ref_state_qw_t));

  matrix_allocate(&commandm,    4, 1);
  matrix_allocate(&filtermemm, 12, 1);
  matrix_allocate(&dutycyclem,  4, 1);

  printf("Set some controller parameters, define dynamics and reference\n");
  /* Initialize filter reference smoother state, commands and set a timestep  */
  for (i = 0; i < 4  ; i++) commandm.pData[i]   = 0.0;
  for (i = 0; i < 12 ; i++) filtermemm.pData[i] = 0.0;
  dt = 0.1;

  controller->status   = 1;    /* Controller status - 0 for idle 1 for active */
  controller->gain_kR  = 4.0;  /* Gain relating to the attitude error         */
  controller->gain_kc  = 0.1;  /* Gain relating to the cross-terms            */
  controller->gain_kw  = 2.0;  /* Gain relating to the attitude rate error    */
  controller->gain_eps = 0.1;  /* Gain relating to the cross-terms            */
  controller->gain_L   = 2.0;  /* Gain relating to the attitude rate error    */
  controller->param_a  = 0.35; /* Constant related to the PWM-thrust map      */
  controller->param_b  = 0.26; /* Constant related to the PWM-thrust map      */
  controller->param_c  = 0.1;  /* Constant related to rotor motor torque      */
  controller->param_d  = 0.2;  /* Distance from rotors to the center of mass  */

  /* Controller inertia */
  matrix_define(&Jm, 3, 3, controller->inertia);
  matrix_set(&Jm, 0, 0, 1.0);
  matrix_set(&Jm, 0, 2, 0.2);
  matrix_set(&Jm, 0, 1, 0.1);
  matrix_set(&Jm, 1, 0, 0.1);
  matrix_set(&Jm, 1, 1, 1.0);
  matrix_set(&Jm, 1, 2, 0.3);
  matrix_set(&Jm, 2, 0, 0.2);
  matrix_set(&Jm, 2, 1, 0.3);
  matrix_set(&Jm, 2, 2, 1.0);
  /*matrix_print(&Jm, "Defined inertia");*/

  for (i = 0; i < 100; i++){
    /* 1. Update the system state in "state" and update user commands in "comamndm"*/
    /* state->quaternion  = ... */
    /* state->omega       = ... */
    /* commandm->pData[0] = ... */ /* lateral thrust force reference */
    /* commandm->pData[1] = ... */ /* pitch ref (rotation about the e1 axis)*/
    /* commandm->pData[2] = ... */ /* roll ref (rotation about the e2 axis)*/
    /* commandm->pData[3] = ... */ /* yaw ref (rotation about the e3 axis)*/
    /* dt                 = ... */ /* time-step since last update */

    /* 2. Update the reference trajectory */
    assert(1==update_attitude_references(&commandm, &filtermemm, reference, dt));

    /* 3. Compute the rigid-body forces and torques in the controller, we can
    use any of the controllers below, but all are called for debugging purposes
    to check that there are no memory leaks                                   */
    assert(1==update_attitude_FSF_SO3_continuous(reference, state, controller));
    assert(1==update_attitude_FSF_SU2_continuous(reference, state, controller));
    assert(1==update_attitude_FSF_SU2_discontinuous(reference, state, controller));

    /* 4. Compute the desired thrust - now just setting the reference directly
    here the force could be projected so as to be constant in the z-direction if
    unchanged by the uer. */
    controller->thrust = reference->thrust;

    /* 5. Externalize the moments and forces in the controller to duty cycles */
    assert(1 == compute_power_distribution(controller, &dutycyclem));
  }

  printf("Free allocated memory...\n");
  free(state);
  free(controller);
  free(reference);

  free(commandm.pData);
  free(filtermemm.pData);
  free(dutycyclem.pData);

}

void example_attitude_FOF(void){

  int i;
  double dt;
  matrix_double_t Jm, iJm, Kwm, Cwm, Vim, y0m, y1m, y2m, y3m, commandm, filtermemm, dutycyclem;
  dyn_state_qw_t * state;
  con_state_qw_fof_t * controller;
  ref_state_qw_t * reference;

  printf("Allocate memory for the reference, state and controller structs...\n");
  state      = (dyn_state_qw_t*)     calloc(1, sizeof(dyn_state_qw_t));
  controller = (con_state_qw_fof_t*) calloc(1, sizeof(con_state_qw_fof_t));
  reference  = (ref_state_qw_t*)     calloc(1, sizeof(ref_state_qw_t));

  matrix_allocate(&commandm,    4, 1);
  matrix_allocate(&filtermemm, 12, 1);
  matrix_allocate(&dutycyclem,  4, 1);

  printf("Set some controller parameters, define dynamics and reference\n");
  /* Initialize filter reference smoother state, commands and set a timestep  */
  for (i = 0; i < 4  ; i++) commandm.pData[i]   = 0.0;
  for (i = 0; i < 12 ; i++) filtermemm.pData[i] = 0.0;
  dt = 0.1;

  /* Scalar controller parameters */
  controller->status = 1;      /* Controller status - 0 for idle 1 for active */
  controller->gain_cw = 1.0;   /* Gain relating to the attitude rate error    */
  controller->gain_cR = 1.0;   /* Gain relating to the attitude rate error    */
  controller->param_a = 0.35;  /* Constant related to the PWM-thrust map      */
  controller->param_b = 0.26;  /* Constant related to the PWM-thrust map      */
  controller->param_c = 0.1;   /* Constant related to rotor motor torque      */
  controller->param_d = 0.2;   /* Distance from rotors to the center of mass  */

  /* Controller inertia */
  matrix_define(&Jm, 3, 3, controller->inertia);
  matrix_zero(&Jm);
  matrix_set(&Jm, 0, 0, 0.1);
  matrix_set(&Jm, 1, 1, 0.1);
  matrix_set(&Jm, 2, 2, 0.1);

  /* Controller inertia inverse */
  matrix_define(&iJm, 3, 3, controller->invinertia);
  matrix_zero(&iJm);
  matrix_set(&iJm, 0, 0, 1.0/0.1);
  matrix_set(&iJm, 1, 1, 1.0/0.1);
  matrix_set(&iJm, 2, 2, 1.0/0.1);

  /* Controller gain Kw */
  matrix_define(&Kwm, 3, 3, controller->gain_Kw);
  matrix_identity(&Kwm);

  /* Controller gain ki */
  for ( i = 0; i < 3; i++ ) controller->gain_ki[i] = 1.0;

  /* Controller gain Cw */
  matrix_define(&Cwm, 3, 3, controller->gain_Cw);
  matrix_identity(&Cwm);

  /* Controller global directions Cw */
  matrix_define(&Vim, 3, 3, controller->globalDirections);
  matrix_identity(&Vim);

  /* Matrices pointing to the measurements */
  matrix_allocate(&y0m, 3, 1);
  matrix_allocate(&y1m, 3, 1);
  matrix_allocate(&y2m, 3, 1);
  matrix_allocate(&y3m, 3, 1);

  for (i = 0; i < 100; i++){
    /* 1. Update user commands in "comamndm"*/
    /* commandm->pData[0] = ... */ /* lateral thrust force reference */
    /* commandm->pData[1] = ... */ /* pitch ref (rotation about the e1 axis)*/
    /* commandm->pData[2] = ... */ /* roll ref (rotation about the e2 axis)*/
    /* commandm->pData[3] = ... */ /* yaw ref (rotation about the e3 axis)*/
    /* dt                 = ... */ /* time-step since last update */

    /* 2. Update the reference trajectory*/
    assert(1==update_attitude_references(&commandm, &filtermemm, reference, dt));

    /* 3. Update the measurements in the controller struct - gyro rates in y0,
          and the accelerometer directions in y1, with mg- directions in y2
          and a virtual measurements taken as the cross-product of the two    */
    matrix_set(&y0m, 0, 0, +0.1);
    matrix_set(&y0m, 1, 0, -0.1);
    matrix_set(&y0m, 2, 0, +0.2);
    matrix_set(&y1m, 0, 0, +0.3);
    matrix_set(&y1m, 1, 0, -0.1);
    matrix_set(&y1m, 2, 0, -0.1);
    matrix_set(&y2m, 0, 0, -0.1);
    matrix_set(&y2m, 1, 0, -0.1);
    matrix_set(&y2m, 2, 0, +0.3);
    cont_normalize(&y1m);                 /* Normalize the dir. meas.         */
    cont_normalize(&y2m);                 /* Normalize the dir. meas.         */
    cont_cross_product(&y1m, &y2m, &y3m); /* Construct virtual measurement    */
    cont_normalize(&y3m);                 /* Normalize virtual measurement    */
    assert(1==update_attitude_FOF_SO3_continuous_measurements(controller, &y0m, &y1m, &y2m, &y3m));

    /* 4. Compute the rigid-body forces and torques in the controller, we can
    use any of the controllers below, but all are called for debugging purposes
    to check that there are no memory leaks                                   */
    assert(1==update_attitude_FOF_SO3_continuous(reference, state, controller, dt));

    /* 5. Compute the desired thrust - now just setting the reference directly
    here the force could be projected so as to be constant in the z-direction if
    unchanged by the uer. */

    /* TODO fix the rotation of the thrust
    controller->thrust = reference->thrust;
    */

    /* 6. Externalize the moments and forces in the controller to duty cycles
    assert(1 == compute_power_distribution(controller, &dutycyclem));
    TODO, write the power distribution controller independent
    */
  }

  printf("Free allocated memory...\n");
  free(state);
  free(controller);
  free(reference);

  free(commandm.pData);
  free(filtermemm.pData);
  free(dutycyclem.pData);
  free(y0m.pData);
  free(y1m.pData);
  free(y2m.pData);
  free(y3m.pData);

}
