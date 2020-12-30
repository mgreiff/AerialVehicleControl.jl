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
for any given aerial vehicle control problem.

### Table of contents
* @ref Installation
* @ref Testing
* @ref Examples
* @ref Tuning

### Implemented controllers
The library consists of several controllers, distinguished by begin either
full-state-feedback controllers (FSF) or filtered output feedback controllers
(FOF). Furthermore, the controllers are either continuous, discontinuous, or
robust. A complete summary of the controllers, with their current implementation
status is given by the table below.

| Controller          | Configuration manifold | Implemented in C | Tested in Julia | Example in Julia |
|---------------------|------------------------|------------------|-----------------|------------------|
| FSF (continuous)    | SO(3)                  | \ref cont_attitude_FSF_SO3_continuous.c "Yes" | \ref Testing "Yes" | \ref example_cont_attitude_FSF_SO3_continuous "Yes" |
| FSF (robust)        | SO(3)                  | \ref cont_attitude_FSF_SO3_robust.c "Yes" | \ref Testing "Yes" | \ref example_cont_attitude_FSF_SO3_robust "Yes" |
| FSF (continuous)    | SU(2)                  | \ref cont_attitude_FSF_SU2_continuous.c "Yes" | \ref Testing "Yes" | \ref example_cont_attitude_FSF_SU2_continuous "Yes" |
| FSF (discontinuous) | SU(2)                  | \ref cont_attitude_FSF_SU2_discontinuous.c "Yes" | \ref Testing "Yes" | \ref example_cont_attitude_FSF_SU2_discontinuous "Yes" |
| FSF (robust)        | SU(2)                  | \ref cont_attitude_FSF_SU2_robust.c "Yes" | \ref Testing "Yes" | \ref example_cont_attitude_FSF_SU2_robust "Yes" |
| FSF (continuous)    | SO(3) x R^3            | No               | No              | No               |
| FSF (continuous)    | SU(2) x R^3            | No               | No              | No               |
| FSF (discontinuous) | SU(2) x R^3            | No               | No              | No               |
| FOF (continuous)    | SO(3) x SO(3)          | \ref cont_attitude_FOF_SO3_continuous.c "Yes" | \ref Testing "Yes" | \ref example_cont_attitude_FOF_SO3_continuous "Yes" |
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

@page Installation
To run the code, you will need to install
* [Julia](https://julialang.org/downloads/platform/)
* [GCC](https://gcc.gnu.org/)
* [LAPACK](http://www.netlib.org/lapack/) and [BLAS](http://www.netlib.org/blas/)
* [Doxygen](https://www.doxygen.nl/download.html)

The last point is only required if you wish to regenerate the docs (see the
above links for platform specific installation instructions). On an Ubuntu
20.04 Distribution, and in the order that they were mentioned, this can be done
by executing the following
```
sudo apt install julia
sudo apt install build-essential
sudo apt-get install libblas-dev liblapack-dev
sudo apt-get install doxygen
sudo apt-get install graphviz
```
Next, clone the repository and enter the Julia REPL by
```
git clone https://github.com/mgreiff/AerialVehicleControl.jl
cd AerialVehicleControl.jl
julia
```
and activate the module it by running the following command in the Julia REPL
```
]activate .
```
You can now run the tests and examples in \ref Testing and \ref Examples respectively.

@page Testing
The C code is tested using Julia 1.5.1 (although any version >1.0 should work),
and all tests can be run by executing ``include("runtests.jl")`` in the Julia
REPL. This compiles the C-code and calls it through Libdl using ``ccall``,
and subsequently verifies that each function in the C-stack works as intended.

From the base directory, all tests can be run after activating the package, by

```
]activate .
include("test/runtests.jl")
```

This compiles the C-code as a shared library and proceeds to call individual functions, testing them against the same functions implemented in Julia. You can also run the tests individually to see exactly how specific functions have been implemented and their intended usage. Running all the tests produces the output below.

```
Testing math functions...
Testing the SO(3) maps...
Testing the SU(2) maps...
Testing FSF utilities...
Testing FSF continuous on SO(3)...
Testing FSF robust on SO(3)...
Testing FSF continuous on SU(2)...
Testing FSF discontinuous on SU(2)...
Testing FSF robust on SU(2)...
Testing FOF continuous on SO(3)...
Testing the power distribution...
Matrix addition (inplace)...
Matrix addition...
Matrix cholesky decpomosition (inplace)...
Matrix cholesky decpomosition...
Matrix PSD solver...
Matrix multiplication...
Matrix PSD solver...
Matrix subtraction (inplace)...
Matrix subtraction...
Matrix transposition...
Test Summary: | Pass  Total
All tests     | 4123   4123
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

@page "Examples"
The C-code can be used in the loop by calls to the relevant controller functions
using
[``ccall()``](https://docs.julialang.org/en/v1/manual/calling-c-and-fortran-code/)
and the differential equation solvers in
[``DifferentialEquations.jl``](https://docs.sciml.ai/stable/). The former
permits the wrapping of individual
controllers, and a set of examples implemented in this way can be found in
``/examples``, and can be run out of the box. This permits a study of the
controller performance, and also facilitates the plotting of signals computed
internally in the C-code, such as the Lyapunov function associated with each
controller and its theoretical upper and lower bounds.

@section example_FSF_attitude FSF Attitude control
All of the attitude controllers are implemented so as to take a reference
(\ref ref_state_qw_t R) which is updated by the user and includes the reference
system state: a reference attitude (as a quaternion), the reference attitude
rates, and the reference attitude accelerations. This structure can be filled
from user commands by means of filtering using the
\ref update_attitude_references() utility function. The controllers also take
take a state object (\ref dyn_state_qw_t S) containing the measured or estimated
system state: an attitude (as a quaternion) and attitude rates. The FSF attitude
controllers are all defined by a set of parameters and controller memory
(\ref con_state_qw_fsf_t C), and in the attitude control example, found in
``examples/example_FSF_attitude.jl``, these structures are (i) first initialized
(with a nominal tuning in C and randomized initial conditions in R and S), then
(ii) used to define a function handle for numerical integration, (iii) solved
using [``DifferentialEquations.jl``](https://docs.sciml.ai/stable/) and finally
(iv) the solution is visualized using the ``solution_analysis()`` function.
In this example, you can define the controller that is to be used by setting
the controllerType variable

* controllerType=1 \f$\to\f$ @ref example_cont_attitude_FSF_SO3_continuous
* controllerType=2 \f$\to\f$ @ref example_cont_attitude_FSF_SO3_robust
* controllerType=3 \f$\to\f$ @ref example_cont_attitude_FSF_SU2_continuous
* controllerType=4 \f$\to\f$ @ref example_cont_attitude_FSF_SU2_discontinuous
* controllerType=5 \f$\to\f$ @ref example_cont_attitude_FSF_SU2_robust

You can also choose if you wish to include a disturbance or not. The controllers
will be tuned to guarantee feasibility, and the example warns if the
initialization is such that the error dynamics are initialized outside of the
domain of exponential attraction. This is done using the ``get_info()`` function
implemented in the ``solution_analysis()`` function, whose output in terms of
analysis and plotting is shown in the following examples, with notation
following \cite mgreiff2020robust.

@subsection example_cont_attitude_FSF_SO3_continuous Continuous control on SO(3)

Here we run the continuous FSF on SO(3) defined in
cont_attitude_FSF_SO3_continuous.c (see \cite lee2010geometric for the original
reference). The controller is run without any disturbance, and both the
system state and the reference are internally represented by quaternions. As
such, we may converge to an error on SU(2) which is \f$\pm I\f$, both
representing the same element on SO(3). Hence, \f$\Gamma(X_r, X)\rightarrow \{0\lor 2\}\f$,
\f$\bar{\Gamma}(X_r, X)\rightarrow \{0\lor 2\}\f$, and
\f$\Psi(R_r, R)\rightarrow 0\f$ in both of these cases. Furthermore, with an
appropriate tuning, the Lyapunov function will be monotonically decreasing. In
this particular example, the ``get_info()`` function returns the following
```
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~ Simulation done with the continuous FSF controller on SO(3) (called from C) ~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* The initial attitude error is: Psi(Rr(t0), R(t0)) = 1.37
* Any V(t0)/kR = 1.40 < phi < 2 can be used in the stability proof. We let phi = 1.40.
* Worst case decay rate of the Lyapunov function : 0.04120564135100125
```
and the control system states and errors are visualized below.

\image html continuous_SO3_FSF_states.png "States and controls when calling the continuous FSF attitude controller on SO(3)" width=500px
\image html continuous_SO3_FSF_errors.png "Errors and attitude distances when calling the continuous FSF attitude controller on SO(3)" width=500px
\image html continuous_SO3_FSF_analysis.png "Lyapunov function with theoretical bounds for continuous FSF attitude controller on SO(3)" width=500px

Note that \f$\mathcal{V}(t)\f$ decays down to the precision set in the ODE
solver, here an absolute and relative tolerance of \f$10^{-8}\f$, so the
violation of the theoretical bounds is expected as \f$\mathcal{V}(t)\f$ becomes
small, and arises due to numerical errors in the ODE solver.


@subsection example_cont_attitude_FSF_SO3_robust Robust control on SO(3)

Here we run the robust FSF on SO(3) defined in cont_attitude_FSF_SO3_robust.c
(see \cite lee2013nonlinear for the original reference). The controller is run
with a sinusoidal disturbance, which is upper bound in the \f$l_2\f$-norm by
\f$\sup_{t \geq 0}\|d(t)\|_2<L=1.0\f$, and uses the same seed as the previous
example for the generation of the initial conditions.

The robust controller should behave much better in the face of the applied
load disturbances than the continuous controller on SO(3), and for this
particular example, the Lyapunov function settles around
\f$\mathcal{V}(t)\approx 10^{-3}\f$ for the robust controller, whereas
\f$\mathcal{V}(t)\approx 0.5\f$ if running the continuous controller in
\ref example_cont_attitude_FSF_SO3_continuous with the same disturbance.

Again, both the system state and the reference are internally represented by
quaternions. As such, we may converge to an error on SU(2) which is \f$\pm I\f$, both
representing the same element on SO(3). Hence, \f$\Gamma(X_r, X)\rightarrow \{0\lor 2\}\f$,
\f$\bar{\Gamma}(X_r, X)\rightarrow \{0\lor 2\}\f$, and
\f$\Psi(R_r, R)\rightarrow 0\f$ in both of these cases. Furthermore, with an
appropriate tuning, the Lyapunov function will be monotonically decreasing. In
this particular example, the ``get_info()`` function returns the following
```
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~ Simulation done with the robust FSF controller on SO(3) (called from C) ~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* The initial attitude error is: Psi(Rr(t0), R(t0)) = 1.37
* Any V(t0)/kR = 1.40 < phi < 2 can be used in the stability proof. We let phi = 1.40.
* Worst case decay rate of the Lyapunov function : 0.0412
* Upper bound on the allowed epsilon given phi   : 0.00595
* Uniform ultimate bound K*eps, where K is       : 141.06
```
and the control system states and errors are visualized below

\image html robust_SO3_FSF_states.png "States, controls and disturbance when calling the robust FSF attitude controller on SO(3)" width=500px
\image html robust_SO3_FSF_errors.png "Errors and attitude distances when calling the robust FSF attitude controller on SO(3)" width=500px
\image html robust_SO3_FSF_analysis.png "Lyapunov function and normed errors with theoretical bounds for robust FSF attitude controller on SO(3)" width=500px


@subsection example_cont_attitude_FSF_SU2_continuous Continuous control on SU(2)

Here we run the continuous FSF on SU(2) defined in
cont_attitude_FSF_SU2_continuous.c (see \cite mgreiff2020robust). The controller
is run without any disturbance, and both the
system state and the reference are internally represented by quaternions. Due to the
1:1 correspondence between quaternions and elements of SU(2), we always converge
to an error element on SU(2) which is \f$X_e\to I\f$. Hence \f$\Gamma(X_r, X)\rightarrow 0\f$,
\f$\bar{\Gamma}(X_r, X)\rightarrow 2\f$, and \f$\Psi(R_r, R)\rightarrow 0\f$. Furthermore, with an
appropriate tuning, the Lyapunov function will be monotonically decreasing. In
this particular example, the ``get_info()`` function returns the following
```
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~ Simulation done with the continuous FSF controller on SU(2) (called from C) ~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* The initial attitude error is: Gamma(Xr(t0), X(t0)) = 1.56
* Warning: This is relatively large, requiring very small initial attitude rate errors
* Any V(t0)/kR = 1.57 < phi < 2 can be used in the stability proof. We let phi = 1.57
* Worst case decay rate of the Lyapunov function : 0.014
```
\image html continuous_SU2_FSF_states.png "States, controls and disturbance when calling the continuous FSF attitude controller on SU(2)" width=500px
\image html continuous_SU2_FSF_errors.png "Errors and attitude distances when calling the continuous FSF attitude controller on SU(2)" width=500px
\image html continuous_SU2_FSF_analysis.png "Lyapunov function and normed errors with theoretical bounds for continuous FSF attitude controller on SU(2)" width=500px

Note that \f$\mathcal{V}(t)\f$ decays down to the precision set in the ODE
solver, here an absolute and relative tolerance of \f$10^{-8}\f$, so the
violation of the theoretical bounds is expected as \f$\mathcal{V}(t)\f$ becomes
small, and arises due to numerical errors in the ODE solver.

@subsection example_cont_attitude_FSF_SU2_discontinuous Discontinuous control on SU(2)

Here we run the discontinuous FSF on SU(2) defined in
cont_attitude_FSF_SU2_discontinuous.c (see \cite mgreiff2020robust). The controller
is run without any disturbance, and both the
system state and the reference are internally represented by quaternions. This
example is run with the same initialization as the previous with the continuous
control on SU(2) to showcase their difference. Here
the SU(2) manifold is partitioned into two regions, on which the attitude error
will converge to \f$X_e\to I\f$ or \f$X_e\to I\f$ respectively. Hence
\f$\Gamma(X_r, X)\rightarrow \{0\lor 2\}\f$,
\f$\bar{\Gamma}(X_r, X)\rightarrow \{2\lor 0\}\f$, and
\f$\Psi(R_r, R)\rightarrow 0\f$.
Furthermore, with an appropriate tuning, the Lyapunov function will be monotonically decreasing. In
this particular example, the ``get_info()`` function returns the following
```
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~ Simulation done with the discontinuous FSF controller on SU(2) (called from C) ~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* The initial attitude error is: Gamma(Xr(t0), X(t0)) = 1.56
* Warning: This is relatively large, requiring very small initial attitude rate errors
* Any V(t0)/kR = 1.0 < phi < 2 can be used in the stability proof. We let phi = 1.0.
* Worst case decay rate of the Lyapunov function : 0.033
```
In general, if \f$\Gamma(X_r, X)>1\f$ the decay
of the error will be more rapid in the discontinuous controller. This is seen in
the worst case decay rates, being 0.032 with the discontinuous controller and
0.014 with the previous example, and this difference can also be seen in the
plots below.

\image html discontinuous_SU2_FSF_states.png "States, controls and disturbance when calling the discontinuous FSF attitude controller on SU(2)" width=500px
\image html discontinuous_SU2_FSF_errors.png "Errors and attitude distances when calling the discontinuous FSF attitude controller on SU(2)" width=500px
\image html discontinuous_SU2_FSF_analysis.png "Lyapunov function and normed errors with theoretical bounds for discontinuous FSF attitude controller on SU(2)" width=500px

@subsection example_cont_attitude_FSF_SU2_robust Robust control on SU(2)
Here we run the robust FSF on SU(2) defined in
cont_attitude_FSF_SU2_robust.c (see \cite mgreiff2020robust). The controller
is run with a sinusoidal disturbance, which is upper bound in the \f$l_2\f$-norm by
\f$\sup_{t \geq 0}\|d(t)\|_2<L=1.0\f$, just as with the simulation of the robust
controller on SO(3). Due to the
1:1 correspondence between quaternions and elements of SU(2), we always converge
to an error element on SU(2) which is \f$X_e\to I\f$. Hence \f$\Gamma(X_r, X)\rightarrow 0\f$,
\f$\bar{\Gamma}(X_r, X)\rightarrow 2\f$, and \f$\Psi(R_r, R)\rightarrow 0\f$. Furthermore, with an
appropriate tuning, the Lyapunov function will be monotonically decreasing. In
this particular example, the ``get_info()`` function returns the following
```
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~ Simulation done with the robust FSF controller on SU(2) (called from C) ~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* The initial attitude error is: Gamma(Xr(t0), X(t0)) = 1.56
* Warning: This is relatively large, requiring very small initial attitude rate errors
* Any V(t0)/kR = 1.5747146320898264 < phi < 2 can be used in the stability proof. We let phi = 1.57.
* Worst case decay rate of the Lyapunov function : 0.014
* Upper bound on the allowed epsilon given phi   : 0.0016
* Uniform ultimate bound K*eps, where K is       : 417.78
```

\image html robust_SU2_FSF_states.png "States, controls and disturbance when calling the robust FSF attitude controller on SU(2)" width=500px
\image html robust_SU2_FSF_errors.png "Errors and attitude distances when calling the robust FSF attitude controller on SU(2)" width=500px
\image html robust_SU2_FSF_analysis.png "Lyapunov function and normed errors with theoretical bounds for robust FSF attitude controller on SU(2)" width=500px

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
