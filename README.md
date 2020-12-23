@mainpage
This software stack serves as
a framework for the analysis and generation of controllers pertaining to
quad-rotor dynamics. The controllers are implemented in C89 and can wrap
[LAPACK](http://www.netlib.org/lapack/) and
[BLAS](http://www.netlib.org/blas/), or the
[CMSIS DSP ARM math library](http://www.keil.com/pack/doc/CMSIS/DSP/html/index.html),
should these be available - if not, the
code runs without any dependencies. The general idea is to call these
controllers through Julia, to (i) test all of the components of the controllers
using unit testing, (ii) quantitatively study their closed loop
behavior in simulation, and (iii) tune them and offline in a data-driven manner.

Importantly, once a controller has been chosen for a given problem, the exact
code used in the simulations can be run directly on the application, being
platform independent C-code. Thus the project does not aim to solve a single
control problem, but rather serves as an analysis tool and implementation aid
for any given UAV control problem.


Install Julia: export PATH="$PATH:/path/to/<Julia directory>/bin"
Install lapack and BLAS: sudo apt-get install libblas-dev liblapack-dev
Install Doxygen:

### Table of contents
* @ref Testing
* @ref Analysis
* @ref Tuning

### Implemented controllers
The library consists of several controllers, distinguished by begin either
full-state-feedback controllers (FSF) or filtered output feedback controllers
(FOF). Furthermore, the controllers are either continuous or discontinuous. A
complete summary of the controllers, with their current implementation status
is given by the table below.

| Controller          | Configuration manifold | Implemented in C | Tested in Julia | Example in Julia |
|---------------------|------------------------|------------------|-----------------|------------------|
| FSF (continuous)    | SO(3)                  | \ref cont_attitude_FSF_SO3_continuous.c "Yes" | \ref test_FSF_continuous_SO3.jl "Yes" | \ref example_cont_attitude_FSF_SO3_continuous "Yes" |
| FSF (continuous)    | SU(2)                  | \ref cont_attitude_FSF_SU2_continuous.c "Yes" | \ref test_FSF_continuous_SU2.jl "Yes" | \ref example_cont_attitude_FSF_SU2_continuous "Yes" |
| FSF (discontinuous) | SU(2)                  | \ref cont_attitude_FSF_SU2_discontinuous.c "Yes" | \ref test_FSF_discontinuous_SU2.jl "Yes" | \ref example_cont_attitude_FSF_SU2_discontinuous "Yes" |
| FSF (continuous)    | SO(3) x R^3            | No               | No              | No               |
| FSF (continuous)    | SU(2) x R^3            | No               | No              | No               |
| FSF (discontinuous) | SU(2) x R^3            | No               | No              | No               |
| FOF (continuous)    | SO(3) x SO(3)          | \ref cont_attitude_FOF_SO3_continuous.c "Yes" | \ref test_FOF_continuous_SO3.jl "Yes" | \ref example_cont_attitude_FOF_SO3_continuous "Yes" |
| FOF (continuous)    | SU(2) x SU(2)          | No               | No              | No               |
| FOF (discontinuous) | SU(2) x SU(2)          | No               | No              | No               |
| FOF (continuous)    | SO(3) x SO(3) x R^6    | No               | No              | No               |
| FOF (continuous)    | SU(2) x SU(2) x R^6    | No               | No              | No               |
| FOF (discontinuous) | SU(2) x SU(2) x R^6    | No               | No              | No               |

### Implemented utilities
In addition to the controllers, basic functionality such as the computation of
duty cycles from rigid-body forces and torques, as well as reference generators
utilizing differential flatness are all provided and also tested in Julia.

| Utility             | Configuration manifold | Implemented in C | Tested in Julia | Example in Julia |
|---------------------|------------------------|------------------|-----------------|------------------|
| Power Distribution  | -                      | \ref cont_power_distribution.c "Yes" | Yes | No |
| Flatness eq.        | SU(2) x R^3            | No               | No              | No               |
| Reference generator | SU(2)                  | \ref cont_attitude_reference_generator.c "Yes"               | No              | \ref example_cont_attitude_reference_generator "Yes" |
| Reference generator | SU(2) x R^3            | No               | No              | No               |

@page Testing
### Testing
The C code is tested using Julia 1.4.2 (although any version >1.0 should work),
and all tests can be run by executing ``include("runtests.jl")`` in the Julia
REPL. This compiles the cont-gen stack and calls it through Libdl using the help
functions in ``/util``, and subsequently verifies that each function in the
stack works as intended. From the base directory, running

```
cd tests && julia -e 'include("runtests.jl")'
```

launches Julia and produces an output similar to the figure below

```
Test Summary:                                         | Pass  Fail  Total
All tests                                             | 3896     4   3900
  Math library tests                                  |  328          328
  Attitude controllers                                |  250          250
  Utility functions                                   |    8     4     12
    Test the power distribution                       |    8     4     12
      Good inputs                                     |    1     4      5
      Bad inputs (negative thrust)                    |    1            1
      Bad inputs (negative parameter value)           |    4            4
      Bad inputs (wrong duty cycle row dimensions)    |    1            1
      Bad inputs (wrong duty cycle column dimensions) |    1            1
  Matrix math functions                               | 3310         3310
ERROR: LoadError: Some tests did not pass: 3896 passed, 4 failed, 0 errored, 0 broken.
in expression starting at /home/mgreiff/Desktop/cont-gen/tests/runtests.jl:15
```

In addition, examples of C-implementations with complete control loops are given
in the cont_main stack, and by first compiling the stack and then running the
debug command from the bash shell, Valgrind will be run on the compiled stack to
check for memory leaks. From the base directory, running

```
cd tests && make clean && make && make debug
```

produces an output similar to that below

```
Allocate memory for the reference, state and controller structs...
Set some controller parameters, define dynamics and reference
Free allocated memory...
Allocate memory for the reference, state and controller structs...
Set some controller parameters, define dynamics and reference
Free allocated memory...
==13480==
==13480== HEAP SUMMARY:
==13480==     in use at exit: 0 bytes in 0 blocks
==13480==   total heap usage: 8,238 allocs, 8,238 frees, 351,976 bytes allocated
==13480==
==13480== All heap blocks were freed -- no leaks are possible
==13480==
==13480== For counts of detected and suppressed errors, rerun with: -v
==13480== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)
```

@page "Analysis"
### Analysis
The C-code can be used in the loop by calls to the relevant controller functions
using
[``ccall()``](https://docs.julialang.org/en/v1/manual/calling-c-and-fortran-code/)
and the differential equation solvers in
[``DifferentialEquations.jl``](https://docs.sciml.ai/stable/). The former
permits the wrapping of individual
controllers, such as the call to the FSF attitude feedback on SO(3) below
```
function update_control( R::ref_state_qw_t, S::dyn_state_qw_t, C::con_state_qw_fsf_t)
    status = ccall((:update_control, CONT_LIB_NAME),
        Cint,
        (Ref{ref_state_qw_t}, Ref{dyn_state_qw_t}, Ref{con_state_qw_fsf_t},),
        R,
        S,
        C)
    return status;
end
```
This function can subsequently be called in the loop when simulating the NLTV
systems using [``DifferentialEquations.jl``](https://docs.sciml.ai/stable/), as
```
x0, tspan, C = initialize_example()
prob         = ODEProblem(odefun!, x0, tspan, C);
sol          = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
```
A set of such examples implemented in this way can be found in ``/examples``,
and can be run out of the box. This permits a study of the controller
performance, and also facilitates the plotting of signals computed internally in
the C-code, such as the Lyapunov function associated with each controller. Some
simulation examples of the code is given below as follows

* @ref example_cont_attitude_FSF_SO3_continuous
* @ref example_cont_attitude_FSF_SU2_continuous
* @ref example_cont_attitude_FSF_SU2_discontinuous
* @ref example_cont_attitude_FOF_SO3_continuous
* @ref example_cont_attitude_reference_generator

@section example_cont_attitude_FSF_SO3_continuous Continuous attitude FSF on SO(3)

Here we run the full state feedback (FSF) continuous SO(3) controller defined in
cont_attitude_FSF_SO3_continuous.c, steering the SU(2)-configured attitude
dynamics. This example can be reproduced by running
``include("example_FSF_attitude.jl")``, makes use of the
infrastructure provided by the [``DifferentialEquations.jl``](https://docs.sciml.ai/stable/) stack, and permits
the generation of plots such as the one below. Note that we may converge to an
error on SU(2) which is  \f$\pm I\f$, both representing the same
element on SO(3). Hence, \f$\Gamma(X_r, X)\rightarrow \{0\lor 2\}\f$,
\f$\bar{\Gamma}(X_r, X)\rightarrow \{0\lor 2\}\f$, and
\f$\Psi(R_r, R)\rightarrow 0\f$. Furthermore, with an appropriate tuning, the
Lyapunov function will be monotonically decreasing.

\image html attitude_dynamics_cont_SO3_states.png "States and controls when calling the continuous FSF attitude controller on SO(3) in Julia" width=500px
\image html attitude_dynamics_cont_SO3_errors.png "Distances and errors when calling the continuous FSF attitude controller on SO(3) in Julia" width=500px

@section example_cont_attitude_FSF_SU2_continuous Continuous attitude FSF on SU(2)

Here we run the full state feedback (FSF) continuous SU(2) controller defined in
cont_attitude_FSF_SU2_continuous.c, steering the SU(2)-configured attitude
dynamics. This example can be reproduced by running
``include("example_FSF_attitude.jl")``, makes use of the
infrastructure provided by the [``DifferentialEquations.jl``](https://docs.sciml.ai/stable/) stack, and permits
the generation of plots such as the one below. Here we may converge to an
error on SU(2) \f$+I\f$. Hence, \f$\Gamma(X_r, X)\rightarrow 0\f$,
\f$\bar{\Gamma}(X_r, X)\rightarrow 2 \f$, and
\f$\Psi(R_r, R)\rightarrow 0\f$. Furthermore, with an appropriate tuning, the
Lyapunov function will be monotonically decreasing.

\image html attitude_dynamics_cont_SU2_states.png "States and controls when calling the continuous FSF attitude controller on SO(3) in Julia" width=500px
\image html attitude_dynamics_cont_SU2_errors.png "Distances and errors when calling the continuous FSF attitude controller on SO(3) in Julia" width=500px

@section example_cont_attitude_FSF_SU2_discontinuous Discontinuous attitude FSF on SU(2)

Here we run the full state feedback (FSF) discontinuous SU(2) controller defined in
cont_attitude_FSF_SU2_discontinuous.c, steering the SU(2)-configured attitude
dynamics. This example can be reproduced by running
``include("example_FSF_attitude.jl")``, makes use of the
infrastructure provided by the [``DifferentialEquations.jl``](https://docs.sciml.ai/stable/) stack, and permits
the generation of plots such as the one below. Here we may converge to an
error on SU(2) \f$\pm I\f$, boh representinhg the same attitude on SO(3). Hence,
\f$\Gamma(X_r, X)\rightarrow \{0\lor 2\}\f$,
\f$\bar{\Gamma}(X_r, X)\rightarrow \{0\lor 2\} \f$, and
\f$\Psi(R_r, R)\rightarrow 0\f$. Furthermore, with an appropriate tuning, the
Lyapunov function will be monotonically decreasing.

\image html attitude_dynamics_disc_SU2_states.png "States and controls when calling the continuous FSF attitude controller on SO(3) in Julia" width=500px
\image html attitude_dynamics_disc_SU2_errors.png "Distances and errors when calling the continuous FSF attitude controller on SO(3) in Julia" width=500px

@section example_cont_attitude_FOF_SO3_continuous Continuous attitude FOF on SO(3)

Here we run the filtered output feedback (FOF) continuous SO(3) controller
defined in cont_attitude_FOF_SO3_continuous.c, steering the SU(2)-configured
attitude dynamics \f$\{X, \omega\}\f$, to a reference \f$\{X_r, \omega_r\}\f$
by estimating \f$\{\hat{X}, \hat{\omega}\}\f$ without direct knowledge of the
system states.

This is a radically different controller to the previous FSF controller, as we
now provide directional and gyroscopic measurements, which are used to form the
innovation terms of an estimator. The estimator state is integrated on each
time-step in the C-implementation using a CG-method, and as such, we need to use
single-step integration methods in Julia to simulate it.

 This example can be reproduced by running
``include("example_FOF_attitude.jl")``, makes use of the
infrastructure provided by the [``DifferentialEquations.jl``](https://docs.sciml.ai/stable/) stack, and permits
the generation of plots such as the one below. Here we may converge to an
error on SU(2) of \f$\pm I\f$, both in the observer and controller error. Hence,
in some simulations, the quaternion trajectories will mirror each-other as time
goes to infinity. This behaviour is to be expected - as the FOF is designed on
SO(3), the attitude errors \f$\pm I\in SU(2)\f$ corresponds to the element
\f$I\in SO(3)\f$ on SO(3). Furthermore, with any tuning, the Lyapunov function
will be strictly decreasing.

\image html attitude_dynamics_cont_FOF_SO3_states.png "States and controls when calling the continuous FOF attitude controller on SO(3) in Julia" width=500px
\image html attitude_dynamics_cont_FOF_SO3_errors.png "Distances and errors when calling the continuous FOF attitude controller on SO(3) in Julia" width=500px

@section example_cont_attitude_reference_generator The reference generation on SU(2)
The attitude reference generator is a utility that takes noisy and discontinuous
user commands in terms of normalized angular commands and a normalized thrust,
all defined on on the interval [-1,1], and filters these commands to generate
continuous first and second derivatives of the signal through a critically damped
third-order system with a pole in -p on the form

\f[
  G(s) = \frac{p^3}{s^3 + 3ps^2 + 3p^3s + p^3}.
\f]

These signals are subsequently multiplied by a positive constant, defining the
size of the interval on which the signals are defined. For the angles, it would
be reasonable to bound the angles to something like
\f$\gamma_\phi = \gamma_\theta = 0.4\f$ [rad], and the yaw can be set much
larger, at something like \f$\gamma_\psi =2\pi\f$ [rad], and the force should be
centered around the stable hovering force \f$\gamma_{f,1} = mg\f$ with the interval
confined to something like \f$\gamma_{f,2} = mg/2\f$.

\f{eqnarray*}{
     y_{\phi}(s) &=& G(s)\phi_c(s) \qquad \normalfont{(normalized\;filtered\;pitch\;command)},\\
     y_{\phi}(s) &=& G(s)\phi_c(s) \qquad (normalized\;filtered\;roll\;command),\\
     y_{\theta}(s) &=& G(s)\theta_c(s) \qquad (normalized\;filtered\;yaw\;command),\\
     y_{f}(s) &=& G(s)f_c(s) \qquad (normalized\;filtered\;thust\;command).\\
\f}

For the angles, it would be reasonable to bound the angles to something like
\f$\gamma_\phi = \gamma_\theta = 0.4\f$ [rad], and the yaw can be set much
larger, at something like \f$\gamma_\psi =2\pi\f$ [rad], and the force should be
centered around the stable hovering force \f$\gamma_{f,1} = mg\f$ with the interval
confined to something like \f$\gamma_{f,2} = mg/2\f$.

\f{eqnarray*}{
     {\phi}(t) &=& \gamma_{\phi}y_{\phi}(t)\in[-\gamma_{\phi}, \gamma_{\phi}] \quad\qquad\qquad\qquad\qquad\qquad (filtered\;pitch\;command),\\
     {\phi}(t) &=& \gamma_{\theta}y_{\theta}(t)\in[-\gamma_{\theta}, \gamma_{\theta}] \quad\;\;\qquad\qquad\qquad\qquad\qquad (filtered\;roll\;command),\\
     {\theta}(t)) &=& \gamma_{\psi}y_{\psi}(t)\in[-\gamma_{\psi}, \gamma_{\psi}] \;\;\;\qquad\qquad\qquad\qquad\qquad (filtered\;yaw\;command),\\
     {f}(t) &=& \gamma_{f,1} + \gamma_{f,2}y_f(t)\in[\gamma_{f,1}-\gamma_{f,2}  \gamma_{f,1}+\gamma_{f,2}] \qquad (filtered\;thust\;command).\\
\f}

The filtered commands are subsequently expanded into a feasible reference
trajectory in a quaternion \f$q(t)\f$, and attitude rate in the body frame
\f$\omega(t)\f$ and the attitude accelerations \f$\alpha(t)\f$, such that

\f{eqnarray*}{
  q_r(t) &=& Exp([(\phi, 0, 0)/2]_{SU(2)}^\land)\odot Exp([(0, \theta, 0)/2]_{SU(2)}^\land)\odot Exp([(0, 0, \psi)/2]_{SU(2)}^\land)\\
  \dot{q}_r(t) &=& \frac{1}{2}{q}_r(t) \odot [\omega]^\land_{SU(2)},\\
  \dot{\omega}_r(t) &=& \alpha_r(t)
\f}

Below is a simulation example which can be reproduced by running
``include("example_referece_attitude.jl")``, which enters some discontinuous
time-varying commands in \f$(\phi_c(t), \theta_c(t), \psi_c(t), f_c(t))\f$, and
subsequently visualizes the filter memory in time, along with the generated
continuous reference trajectory \f$(q_r(t), \omega_r(t), \alpha_r(t))\f$. To,
verify that the reference trajectory indeed satisfies the attitude dynamics,
the signals \f$\omega_r(t)\f$ (generated in C) are plotted against the signal
\f$2[q_r^*(t)\odot \dot{q}_r(t)]_{SU(2)}^{\lor}\f$ where \f$\dot{q}_r(t)\f$ is
evaluated numerically from \f$q_r(t)\f$ (evaluated in Julia). Similarly, the
signals \f$\alpha_r(t)\f$ (generated in C) is plotted against the numerically
differentiated signals \f$\omega_r(t)\f$ (evaluated in Julia).

This shows that the implementation has been done correctly, and the reference
generator can be used together with any of the attitude controller and the power
distribution to close the loop from commanded inputs to generated PWM duty
cycles.

\image html example_reference_attitude_filtered.png "Normalized filtered commanded inputs" width=500px
\image html example_reference_attitude_refs.png "Generated references based on the commanded inputs with the gamma-constants defining the signal intervals all set to 1" width=500px

@page Tuning
The tuning algorithms are currently implemented in Matlab, but will be migrated
to the cont-gen stack at a later time. This includes robust worst-case tuning by
BMIs, and data-driven tuning algorithms using the NLopts stack.
