################################################################################
# Copyright (C) Marcus Greiff 2020
#
# @file example_FOF_attitude.jl
# @author Marcus Greiff
# @date July 2020
#
# This example demonstrates the filtered output feedback for the attitude
# dynamics with errors related to distances on SO(3), and a closed loop
# simulation, visualizing the state trajectories along with the observer states
# and the state trajectories. SImilar to the FSF attitude examples, it also
# plots the associated Lyapunov function and distances as computed in the
# controller. Contrary to the FSF controllers, this FOF controller has memory,
# by virtue of the integration of the observer dynamcis - consequently, we can
# only call the function once pwe time-step, which means the ay multistep-method
# has to be avoided in the simulation of the closed loop system in Julia.
# Consequently, the system is simulated with a constant step Euler method
# instead of Tsit5().
#
# Furthermore, the method does not have access to the true system states, so
# The lyapunov function and errors have to be evaluated externally in Julia,
# instead of in C. So contrary to the FSF feedback controllers, the estimated
# state (updated by the FOF controller) is stored in the State structure,
# instead of the actual system state.
#
################################################################################
using DifferentialEquations
using LinearAlgebra
using LaTeXStrings
using AerialVehicleControl

showPlot = [1,2];     # Set to 0 for no plot, or any number in (1,2)
savePlot = false; # Set to true if the plot is to be saved
namePlot = "continuous_SO3"

# Recompile the C-code
recompile()

# Wrapper for the calling of the FOF controller from C
function update_control(R::ref_state_qw_t, S::dyn_state_qw_t, C::con_state_qw_fof_t, dt::Cdouble)
    status = ccall((:update_attitude_FOF_SO3_continuous, AerialVehicleControl.CONT_LIB_PATH),
        Cint,
        (Ref{ref_state_qw_t}, Ref{dyn_state_qw_t}, Ref{con_state_qw_fof_t}, Cdouble,),
        R,
        S,
        C,
        dt)
    return status;
end

# Wrapper for the updating of the measurements in the controller structure
function update_measurements(
    C::con_state_qw_fof_t,
    y0m::matrix_double_t,
    y1m::matrix_double_t,
    y2m::matrix_double_t,
    y3m::matrix_double_t
)
    status = ccall((:update_attitude_FOF_SO3_continuous_measurements, AerialVehicleControl.CONT_LIB_PATH),
        Cint,
        (Ref{con_state_qw_fof_t}, Ref{matrix_double_t}, Ref{matrix_double_t},Ref{matrix_double_t}, Ref{matrix_double_t},),
        C,
        y0m,
        y1m,
        y2m,
        y3m)
    return status;
end

# Right-hand side in the ODE used for numerical integration. Here, we have
#
# x[1:4]   = q     (quaternion, integrated in Julia)
# x[5:7]   = w     (att. rates, integrated in Julia)
# x[8:11]  = qr    (quaternion reference, integrated in Julia)
# x[12:14] = wr    (att. rates reference, integrated in Julia)
# x[15:18] = qhat  (quaternion estimate, integrated in C using CG)
# x[19:21] = what  (att. rates estimate, integrated in C using CG)
# x[22:24] = taur  (torques reference, computed in Julia)
# x[25:27] = tau   (torques actuation, computed in C)
# x[28]    = prevt (time at the previous call of the OF-controller)
# x[29]    = lyap  (Lyapunov function, computed in Julia)
# x[30]    = Psi1  (Estimate attitude error, computed in Julia)
# x[31]    = Psi2  (Tracking attitude error, computed in Julia)
#
function odefun!(dx,x,params,t)
    # Extract the states
    q     = x[1:4];
    w     = x[5:7];
    qr    = x[8:11];
    wr    = x[12:14];
    tprev = x[28];

    S = params[1];
    C = params[2];

    # Extract the parameters
    J = ntuple_2_mat(C.inertia, 3, 3)

    # Reference torques driving the reference dynamics
    taur = [sin(t), sin(2*t), sin(3*t)];

    # Reference dynamics
    dqr, dwr = attitude_dynamics(qr, wr, J, taur)

    # Create state structs for ccall
    R = ref_state_qw_t(hcat(qr), hcat(wr), hcat(dwr), 0.1)

    # Compute control signal through a ccall to the FSF configured on SO3
    if tprev != t
        timestep = Cdouble(t - tprev);
        #println("Computed new control action at a time t=$t, dt = $timestep")

        # Update the measurements in the controller object
        v1 = [C.globalDirections[1][i] for i = 1:3];
        v2 = [C.globalDirections[2][i] for i = 1:3];
        v3 = [C.globalDirections[3][i] for i = 1:3];
        y0 = hcat(w);
        y1 = hcat(quat_2_SO3(hcat(q))'*v1);
        y2 = hcat(quat_2_SO3(hcat(q))'*v2);
        y3 = hcat(quat_2_SO3(hcat(q))'*v3);
        y0m = matrix_double_t(y0);
        y1m = matrix_double_t(y1);
        y2m = matrix_double_t(y2);
        y3m = matrix_double_t(y3);
        @assert(1 == update_measurements(C, y0m, y1m, y2m, y3m));

        # Update the controller
        @assert(1 == update_control(R, S, C, timestep));
    end


    # Actuate the dynamics with the computed control signal
    tau = [C.torques[i] for i in 1:3]
    # Physical dynamics
    dq , dw  = attitude_dynamics( q,  w, J, tau )

    # Externalize state derivatives (the observer dynamics are simulated in the FOF controller)
    dx[1:4]   = dq;
    dx[5:7]   = dw;
    dx[8:11]  = dqr;
    dx[12:14] = dwr;

    # Normalize the quaternions which are simulated by DifferentialEquations.jl
    x[1:4]    = q  / norm(q );
    x[8:11]   = qr / norm(qr);

    # Output the estimates
    qhat = [S.quaternion[i] for i = 1:4];
    what = [S.attrates[i] for i = 1:3];
    x[15:18] = qhat;
    x[19:21] = what;

    # Externalize other signals for visualization
    x[22:24]  = taur;
    x[25:27]  = tau;
    x[28]     = t;   # Set the current time so that we can compute the dt

    # Evaluate distances and the Lyapunov fucntion
    Rmat   = quat_2_SO3(hcat(q))
    Rhat   = quat_2_SO3(hcat(qhat))
    Rr     = quat_2_SO3(hcat(qr))
    Re     = Rr*Rmat';
    Rtilde = Rhat*Rmat';
    we     = wr   - w;
    wtilde = what - w;

    v1 = [C.globalDirections[1][i] for i = 1:3];
    v2 = [C.globalDirections[2][i] for i = 1:3];
    v3 = [C.globalDirections[3][i] for i = 1:3];

    erre1 = Re * Rtilde' * v1 - v1;
    erre2 = Re * Rtilde' * v2 - v2;
    erre3 = Re * Rtilde' * v3 - v3;

    errt1 = Rtilde' * v1 - v1;
    errt2 = Rtilde' * v2 - v2;
    errt3 = Rtilde' * v3 - v3;

    x[29] = (+(
               +erre1'*erre1*C.gain_ki[1]
               +erre2'*erre2*C.gain_ki[2]
               +erre3'*erre3*C.gain_ki[3]
             ) / 2.0
             +we'*J*we
             +(
               +errt1'*errt1*C.gain_ki[1]
               +errt2'*errt2*C.gain_ki[2]
               +errt3'*errt3*C.gain_ki[3]
             ) / 2.0
             +(1.0/(2.0 * C.gain_cw)) * wtilde'*J*wtilde)

    x[30] = tr(Matrix{Float64}(I, 3, 3) - Rhat'*Rmat) / 2.0
    x[31] = tr(Matrix{Float64}(I, 3, 3) - Rr'*Rmat) / 2.0
end

function initialize_example()
    # Randomize the initial system states
    q0 = randn(4,1); q0 = q0./norm(q0); w0 = 2*randn(3,1);

    # Randomize an initial reference of states
    qr0 = randn(4,1); qr0 = qr0./norm(qr0); wr0 = 0*rand(3,1); ar0 = 0.0*rand(3,1); fr0 = 0.5;
    R = ref_state_qw_t(hcat(qr0), hcat(wr0), hcat(ar0), fr0)

    # Randomize a set of observer states
    qhat0 = randn(4,1); qhat0 = qhat0./norm(qhat0); what0 = 2*randn(3,1);
    S = dyn_state_qw_t(hcat(qhat0), hcat(what0))

    # Randomize the controller struct
    J  = rand_PSD(3); J  = J ./norm(J );
    Kw = rand_PSD(3); Kw = Kw./norm(Kw); Kw = 2*Kw;
    ki = (rand(1,3) .+ 0.5); ki = ki./norm(ki); ki = 2*ki;
    Cw = rand_PSD(3); Cw = Cw./norm(Cw); Cw = 2*Cw;
    cw, cR = 2.0, 2.0;
    a, b, c, d = 1.0, 2.0, 3.0, 4.0;
    Y0 = rand(3,1) .- 0.5;
    Yi = rand(3,3) .- 0.5;
    Vi = rand(3,3) .- 0.5;
    Vi = Matrix{Float64}(I, 3, 3);
    # Normalize the directions in Vi
    for i = 1:3
        Vi[:,i] = Vi[:,i] / norm(Vi[:,i])
    end

    # Form controller object
    C  = con_state_qw_fof_t(J, Kw, ki, Cw, cw, cR, a, b, c, d, Y0, Yi, Vi);

    # Initial states for the debugging and logging of signals
    output = zeros(10,1);

    x0    = vcat(q0, w0, qr0, wr0, qhat0, what0, output)
    tspan = (0.0,20.0)
    return x0, tspan, R, S, C
end

# Compute the solution
x0, tspan, _, S, C = initialize_example()
prob               = ODEProblem(odefun!, x0, tspan, [S, C]);
sol                = solve(prob, Euler(), adaptive=false, dt = 1e-4)

# Visualization of the results
if 1 in showPlot
    using Plots
    pqqr = plot(sol, vars = (0,1:4),   color=:black, style=:solid, linewidth=2, xaxis="Time (t)",yaxis="Attitude [.]", label=[L"q(t)" nothing nothing nothing])
    plot!(      sol, vars = (0,8:11),  color=:blue,  style=:dash,  linewidth=2, xaxis="Time (t)",yaxis="Attitude [.]", label=[L"q_r(t)" nothing nothing nothing])
    plot!(      sol, vars = (0,15:18), color=:red,   style=:dot,   linewidth=2, xaxis="Time (t)",yaxis="Attitude [.]", label=[L"\hat{q}(t)" nothing nothing nothing])
    pwwr = plot(sol, vars = (0,5:7),   color=:black, style=:solid, linewidth=2, xaxis="Time (t)",yaxis="Attitude rate [rad/s]", label=[L"\omega(t)" nothing nothing])
    plot!(      sol, vars = (0,12:14), color=:blue,  style=:dash,  linewidth=2, xaxis="Time (t)",yaxis="Attitude rate [rad/s]", label=[L"\omega_{r}(t)" nothing nothing])
    plot!(      sol, vars = (0,19:21), color=:red,   style=:dot,   linewidth=2, xaxis="Time (t)",yaxis="Attitude rate [rad/s]", label=[L"\hat{\omega}(t)" nothing nothing])
    pttr = plot(sol, vars = (0,22:24), color=:black, style=:solid, linewidth=2, xaxis="Time (t)",yaxis="Torque [Nm]",           label=[L"\tau_r(t)" nothing nothing])
    plot!(      sol, vars = (0,25:27), color=:blue,  style=:solid, linewidth=2, xaxis="Time (t)",yaxis="Torque [Nm]",           label=[L"\tau(t)" nothing nothing])
    plot(pqqr, pwwr, pttr, layout=(3,1), size=(1000,750))
    gui()
    if savePlot
        savefig(joinpath(AerialVehicleControl.CONT_LIB_DOCS, "$(namePlot)_FOF_states.png"))
    end
end

if 2 in showPlot
    using Plots
    plotattd = plot(sol, vars = (0,30), color=:black, linewidth=2, xaxis="Time (t)",yaxis="Attitude error",    label=L"\Psi(\hat{R}(t), R(t))")
    plot!(          sol, vars = (0,31), color=:blue,  linewidth=2, xaxis="Time (t)",yaxis="Attitude error",    label=L"\Gamma(R_r(t), R(t))")
    plotlyap    = plot(sol, vars = (0,29), color=:black, linewidth=2, xaxis="Time (t)",yaxis="Lyapunov function", label=L"\mathcal{V}(t)")
    plotlyaplog = plot(sol.t, log10.(abs.(sol[29,:])), color=:black, linewidth=2, xaxis="Time (t)",yaxis="Lyapunov function", label=L"log_{10}(\mathcal{V}(t))")
    plot(plotattd, plotlyap, plotlyaplog, layout=(3,1), size=(1000,500))
    gui()
    if savePlot
        savefig(joinpath(AerialVehicleControl.CONT_LIB_DOCS, "$(namePlot)_FOF_errors.png"))
    end
end
