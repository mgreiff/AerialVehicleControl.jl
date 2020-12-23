################################################################################
# Copyright (C) Marcus Greiff 2020
#
# @file example_FOF_observer_sim.jl
# @author Marcus Greiff
# @date July 2020
#
# This example simply demonstrates the integration of the observer dynamics in
# the FOF attitude feedback. This is a bit special, as it does not adhere to
# the regular attitude kinematics, but involves a set of innovation terms. In
# addition, the numerical integration needs to be done in such a way that the
# attitude remains on SU(2), which is the internal representation of the filter
# memory in the numerical integration. As such, the function which is
# exemplified here considers the simulation of the dynamics
#
#     d(Q(t))/dt = Q(t) * [ deltaQ(t) / 2 ]_{SU(2)}^{\land}
#     d(W(t))/dt = deltaW(t)
#
# one step ahead with a time-step dt [s]. This is done by evaluating
#
#     Q(t + dt) = Q(t) * exp([ deltaQ(t) / 2 ]_{SU(2)}^{\land})
#     W(t + dt) = W(t) + dt * deltaW(t)
#
# whih is nontrivial to implement in C, hence verified extrnally by comparison
# to a simulation of the same dynamics using the Tsit5 method in the
# DifferentialEquations.jl stack.
#
# This example computes these to separate solutions and subsequently plots them
# on-top of eachother, showing that we compute the exact same solution in C as
# we do in Julia, and verifying that there are no discernable bugs in the
#
#     attitude_FOF_SO3_continuous_simulate_observer()
#
# function located in cont_attitude_FOF_SO3_continuous.c.
#
################################################################################
using DifferentialEquations
using LinearAlgebra
using LaTeXStrings
using AerialVehicleControl

showPlot = 1;     # Set to 0 for no plot, or any number in (1,2)
savePlot = false; # Set to true if the plot is to be saved

# Recompile the C-code
recompile()

# Wrapper for the calling of the control update law from C
function update_estimator(
    Qm::matrix_double_t,
    Wm::matrix_double_t,
    deltaQm::matrix_double_t,
    deltaWm::matrix_double_t,
    dt::Cdouble)
    status = ccall((:attitude_FOF_SO3_continuous_simulate_observer, AerialVehicleControl.CONT_LIB_PATH),
        Cint,
        (Ref{matrix_double_t},
         Ref{matrix_double_t},
         Ref{matrix_double_t},
         Ref{matrix_double_t},
         Cdouble,),
        Qm,
        Wm,
        deltaQm,
        deltaWm,
        dt)
    return status;
end

function odefun!(dx,x,param,t)
    # Extract the states
    q  = x[1:4];
    w  = x[5:7];

    we      = zeros(4,1);
    dw      = zeros(3,1);
    we[2:4] = w;
    dq      = quaternion_left_product(q) * we / 2.0;
    dw[1]   = cos(t)+1;
    dw[2]   = sin(2*t);
    dw[3]   = sin(t)/2;


    # Externalize state derivatives
    dx[1:4]   = dq;
    dx[5:7]   = dw;
end

function initialize_example()

    # Initial states of the physical system
    q0 = rand(4,1).-0.5;
    q0 = q0./norm(q0);
    w0 = 1*randn(3,1);

    x0    = vcat(q0, w0)
    tspan = (0.0,20.0)
    return x0, tspan
end

# Compute solution using DifferentialEquations.jl
x0, tspan    = initialize_example()
prob         = ODEProblem(odefun!, x0, tspan, []);
sol          = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)

# Compute solution using the implementation in C by the first order Euler method
# implemented in the C-file
time  = collect(0:0.01:20)
Qtraj = zeros(4, length(time));
Wtraj = zeros(3, length(time));
Qtraj[:,1] = x0[1:4,1];
Wtraj[:,1] = x0[5:7,1];

# Run a constant time-step simulation using the update_estimator funciton
# implemented in the FOF attitude controller
for i = 2:length(time);

    dw      = zeros(3,1);
    dw[1]   = cos(time[i-1])+1;
    dw[2]   = sin(2*time[i-1]);
    dw[3]   = sin(time[i-1])/2;

    Q,  _,  Qm = get_matrix_and_copy(Qtraj[:,i-1]);
    W,  _,  Wm = get_matrix_and_copy(Wtraj[:,i-1]);
    dQ, _, dQm = get_matrix_and_copy(Wtraj[:,i-1]);
    dW, _, dWm = get_matrix_and_copy(dw);
    dt = Cdouble(time[i] - time[i - 1]);

    update_estimator(Qm, Wm, dQm, dWm, dt)

    Qtraj[:,i] = Q;
    Wtraj[:,i] = W;

end

# Using the Tsit5() method in Julia as a benchmark,  the solutions in computed
# in DifferentialEquations.jl can be compared to the solution computed using
# the numerical CG-intergration cheme implemented in the FOF controller.
using Plots
pq = plot(sol, vars = (0,1:4), color=:black, linewidth=2, xaxis="Time (t)", yaxis="Attitude [.]", label=[L"q(t)\;\;(computed\;in\;Julia)" nothing nothing nothing])
plot!(time, Qtraj',  color=:red, style=:dash,  linewidth=2, xaxis="Time (t)", yaxis="Attitude [.]", label=[L"q(t)\;\;(computed\;in\;C)" nothing nothing nothing])
pw = plot(sol, vars = (0,5:7), color=:black, linewidth=2, xaxis="Time (t)",yaxis="Attitude rate [rad/s]", label=[L"\omega(t)\;\;(computed\;in\;Julia)" nothing nothing])
plot!(time, Wtraj', color=:red, style=:dash, linewidth=2, xaxis="Time (t)",yaxis="Attitude rate [rad/s]", label=[L"\omega(t)\;\;(computed\;in\;C)" nothing nothing])
plot(pq, pw, layout=(2,1), size=(1000,750))
