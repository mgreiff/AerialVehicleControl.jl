################################################################################
# Copyright (C) Marcus Greiff 2020
#
# @file example_FSF_attitude.jl
# @author Marcus Greiff
# @date July 2020
#
# This example contains the code necessary to call the FSF attitude controllers
# in the loop for simulation and verification by the DifferentialEquations.jl
# stack. By wrapping the functions update_attitude_FSF_*_* un the update_control
# fucntion, this is subsequently called in the odefun!() for simulation Tsit5().
#
# Currently, any of the functions below can be used
#
#     update_attitude_FSF_SO3_continuous()
#     update_attitude_FSF_SU2_continuous()
#     update_attitude_FSF_SU2_discontinuous()
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
function update_control( R::ref_state_qw_t, S::dyn_state_qw_t, C::con_state_qw_fsf_t)
    # Here we can change the called function to any of
    # 1. update_attitude_FSF_SO3_continuous
    # 2. update_attitude_FSF_SU2_continuous
    # 3. update_attitude_FSF_SU2_discontinuous
    status = ccall((:update_attitude_FSF_SO3_continuous, AerialVehicleControl.CONT_LIB_PATH),
        Cint,
        (Ref{ref_state_qw_t}, Ref{dyn_state_qw_t}, Ref{con_state_qw_fsf_t},),
        R,
        S,
        C)
    return status;
end

function odefun!(dx,x,C,t)
    # Extract the states
    q  = x[1:4];
    w  = x[5:7];
    qr = x[8:11];
    wr = x[12:14];

    # Extract the parameters
    J = ntuple_2_mat(C.inertia, 3, 3)

    # Reference torques driving the reference dynamics
    taur = [sin(t), sin(2*t), sin(3*t)];

    # Reference dynamics
    dqr, dwr = attitude_dynamics(qr, wr, J, taur)

    # Create state structs for ccall
    R = ref_state_qw_t(hcat(qr), hcat(wr), hcat(dwr), 0.5)
    S = dyn_state_qw_t(hcat(q ), hcat(w ))

    # Compute control signal through a ccall to the FSF configured on SO3
    status = update_control(R, S, C)

    # Actuate the dynamics with the computed control signal
    tau = [C.torques[i] for i in 1:3]
    # Physical dynamics
    dq , dw  = attitude_dynamics( q,  w, J, tau )

    # Externalize state derivatives
    dx[1:4]   = dq;
    dx[5:7]   = dw;
    dx[8:11]  = dqr;
    dx[12:14] = dwr;

    # Normalize the quaternions
    x[1:4]    = q  / norm(q );
    x[8:11]   = qr / norm(qr);

    # Externalize other signals for visualization
    x[15:17]  = tau;
    x[18:20]  = taur;
    x[21]     = C.dist_Psi;
    x[22]     = C.dist_Gamma;
    x[23]     = C.dist_lyapunov;
    x[24]     = norm(q);
    x[25]     = norm(qr);
end

function initialize_example()
    # Define the controller parameters
    J  = rand_PSD(3);
    J  = J./opnorm(J);
    kR, kc, kw = 5.0, 1.0, 2.0;
    a, b, c, d = 1.0, 1.0, 1.0, 1.0;
    C  = con_state_qw_fsf_t(J, kR, kc, kw, a, b, c, d)

    # Initial states of the physical system
    q0 = rand(4,1).-0.5;
    q0 = q0./norm(q0);
    w0 = 4*rand(3,1);

    # Initial states of the reference system
    qr0 = rand(4,1).-0.5;
    qr0 = qr0./norm(qr0);
    wr0 = zeros(3,1);

    # Initial states for the debugging signals
    output = zeros(11,1);
    output[10,1] = 1.0;
    output[11,1] = 1.0;

    x0    = vcat(q0, w0, qr0, wr0, output)
    tspan = (0.0,20.0)
    return x0, tspan, C
end

# Compute solution
x0, tspan, C = initialize_example()
prob         = ODEProblem(odefun!, x0, tspan, C);
sol          = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)

using Plots
if showPlot == 1
    pqqr = plot(sol, vars = (0,1:4),   color=:black, linewidth=2, xaxis="Time (t)",yaxis="Attitude [.]",          label=[L"q(t)" nothing nothing nothing])
    plot!(      sol, vars = (0,8:11),  color=:red,   linewidth=2, xaxis="Time (t)",yaxis="Attitude [.]",          label=[L"q_r(t)" nothing nothing nothing])
    pwwr = plot(sol, vars = (0,5:7),   color=:black, linewidth=2, xaxis="Time (t)",yaxis="Attitude rate [rad/s]", label=[L"\omega(t)" nothing nothing])
    plot!(      sol, vars = (0,12:14), color=:red,   linewidth=2, xaxis="Time (t)",yaxis="Attitude rate [rad/s]", label=[L"\omega_{r}(t)" nothing nothing])
    pttr = plot(sol, vars = (0,15:17), color=:black, linewidth=2, xaxis="Time (t)",yaxis="Torque [Nm]",           label=[L"\tau(t)" nothing nothing])
    plot!(      sol, vars = (0,18:20), color=:red,   linewidth=2, xaxis="Time (t)",yaxis="Torque [Nm]",           label=[L"\tau_{r}(t)" nothing nothing])
    plot(pqqr, pwwr, pttr, layout=(3,1), size=(1000,750))
    if savePlot
        savefig("../docs/images/attitude_dynamics_disc_SU2_states.png")
    end
end

if showPlot == 2
    plotattd = plot(sol, vars = (0,21), color=:black, linewidth=2, xaxis="Time (t)",yaxis="Attitude error",    label=L"\Psi(R_r(t), R(t))")
    plot!(          sol, vars = (0,22), color=:red, linewidth=2, xaxis="Time (t)",yaxis="Attitude error",    label=L"\Gamma(X_r(t), X(t))")
    plotlyap    = plot(sol, vars = (0,23), color=:black, linewidth=2, xaxis="Time (t)",yaxis="Lyapunov function", label=L"\mathcal{V}(t)")
    plotlyaplog = plot(sol.t, log10.(abs.(sol[23,:])), color=:black, linewidth=2, xaxis="Time (t)",yaxis="Lyapunov function", label=L"log_{10}(\mathcal{V}(t))")
    plot(plotattd, plotlyap, plotlyaplog, layout=(3,1), size=(1000,500))
    if savePlot
        savefig("../docs/images/attitude_dynamics_disc_SU2_errors.png")
    end
end
