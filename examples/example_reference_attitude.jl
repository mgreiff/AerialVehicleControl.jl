################################################################################
# Copyright (C) Marcus Greiff 2020
#
# @file example_reference_attitude.jl
# @author Marcus Greiff
# @date July 2020
#
# Example of the reference attitude generation functionality, from commanded
# inputs (hete steps) to feasible quaternion, attitude rate, and attitude
# acceleration trajectories, which correspond to a physical attitude converging
# with a set speed to the commanded inputs.
#
################################################################################

using LinearAlgebra
using LaTeXStrings
using AerialVehicleControl

showPlot = 1;     # Set to 0 for no plot, or any number in (1,2)
savePlot = false; # Set to true if the plot is to be saved

# Recompile the C-code
recompile()

T = 1000;
time  = zeros(1,T);
traj  = zeros(12,T);
comm  = zeros(4,T);
qtraj = zeros(4,T);
wtraj = zeros(3,T);
atraj = zeros(3,T);
ftraj = zeros(1,T);

c, cc, cm = get_matrix_and_copy(rand(4, 1));
t, tc, tm = get_matrix_and_copy(zeros(12,1));

qr0 = rand(4,1);
qr0 = qr0 ./ norm(qr0);
wr0 = rand(3,1).-0.5;
ar0 = rand(3,1).-0.5;
t0  = 0.0;
R = ref_state_qw_t(qr0, wr0, ar0, t0);


# Run a a simulation at a fixed time-step of 0.01, continuously updating the
# references by the update_attitude_references() function given discrete
# commanded step inputs
dt = 0.01;
for i = 1:T;
    time[1,i] = dt*i;
    c[1,1] = (i > 500) ? +0.3 : -0.5
    c[2,1] = (i > 300) ? -0.3 : -0.5
    c[3,1] = (i > 600) ? +0.7 : -0.1
    c[4,1] = (i > 800) ? +0.9 : +0.5
    comm[:,i] = c[:,1];

    status = ccall((:update_attitude_references, AerialVehicleControl.CONT_LIB_PATH),
        Cint,
        (Ref{matrix_double_t}, Ref{matrix_double_t}, Ref{ref_state_qw_t}, Cdouble,),
        cm,
        tm,
        R,
        dt)

    traj[:,i] = t[:,1];
    for j = 1:4
        qtraj[j,i] = R.quaternion[j];
    end
    for j = 1:3
        wtraj[j,i] = R.attrates[j];
        atraj[j,i] = R.attaccel[j];
    end
    ftraj[1,i] = R.thrust;
end

# Visualize the filtered command trajectories

if showPlot == 1
    using Plots
    plot1 = plot(time[1,:], comm[1,:], linestyle = :dash, color=:black, linewidth=4, xaxis="Time (t)", yaxis="Commanded normalized angle", label=L"f_c(t)")
    plot!(time[1,:], traj[1,:], color=:red, linewidth=2, label=L"y_f(t)")
    plot!(time[1,:], traj[2,:], linestyle = :dash, color=:red, linewidth=2, label=L"(d/dt){y}_f(t)")
    plot!(time[1,:], traj[3,:], linestyle = :dot,  color=:red, linewidth=2, label=L"(d^2/dt^2){y}_f(t)")
    plot2 = plot(time[1,:], comm[2,:], linestyle = :dash, color=:black, linewidth=4, xaxis="Time (t)", yaxis="Commanded normalized angle", label=L"\phi_c(t)")
    plot!(time[1,:], traj[4,:], color=:red, linewidth=2, label=L"y_\phi(t)")
    plot!(time[1,:], traj[5,:], linestyle = :dash, color=:red, linewidth=2, label=L"(d/dt){y}_\phi(t)")
    plot!(time[1,:], traj[6,:], linestyle = :dot,  color=:red, linewidth=2, label=L"(d^2/dt^2){y}_\phi(t)")
    plot3 = plot(time[1,:], comm[3,:], linestyle = :dash, color=:black, linewidth=4, xaxis="Time (t)", yaxis="Commanded normalized angle", label=L"\phi_c(t)")
    plot!(time[1,:], traj[7,:], color=:red, linewidth=2, label=L"y_\psi(t)")
    plot!(time[1,:], traj[8,:], linestyle = :dash, color=:red, linewidth=2, label=L"(d/dt){y}_\theta(t)")
    plot!(time[1,:], traj[9,:], linestyle = :dot,  color=:red, linewidth=2, label=L"(d^2/dt^2){y}_\theta(t)")
    plot4 = plot(time[1,:], comm[4,:], linestyle = :dash, color=:black, linewidth=4, xaxis="Time (t)",yaxis="Commanded normalized force", label=L"\theta_c(t)")
    plot!(time[1,:], traj[10,:], color=:red, linewidth=2, label=L"y_\psi(t)")
    plot!(time[1,:], traj[11,:], linestyle = :dash, color=:red, linewidth=2, label=L"(d/dt){y}_\psi(t)")
    plot!(time[1,:], traj[12,:], linestyle = :dot,  color=:red, linewidth=2, label=L"(d^2/dt^2){y}_\psi(t)")
    plot(plot1, plot2, plot3, plot4, layout=(2,2), size=(1000,750));
    #savefig("../docs/images/example_reference_attitude_filtered.png")
end

# Visualize the resulting attitude reference trajectories
if showPlot == 2
    using Plots
    dq       = diff(qtraj, dims=2)./dt;
    w_approx = zeros(3,T - 1);
    for i = 1:T-1
        qcon = qtraj[:,i];
        qcon[2:4,1] = -qcon[2:4,1];
        prod = quaternion_left_product(qcon) * dq[:,i];
        w_approx[:,i] = 2*prod[2:4,1];
    end

    plotq = plot(time', qtraj', linewidth=2, linestyle = [:solid :dash :dashdot :dot], color=[:black :red :red :red], xaxis="Time (t)",yaxis="Reference attitude", label=[L"q_1(t)" L"q_2(t)" L"q_3(t)" L"q_4(t)"])
    plotw = plot(time', wtraj', linewidth=2, linestyle = [:dash :dashdot :dot], color=:red, xaxis="Time (t)",yaxis="Reference att. rates", label=[L"\omega_1(t)" L"\omega_2(t)" L"\omega_3(t)"])
    plot!(time[1:1:end-1], w_approx', color=:black, label=[L"2 [q^*(t)\otimes ({q}(t) - {q}(t-h))/h]_{SU(2)}^{hat}" nothing nothing])
    plota = plot(time', atraj', linewidth=2, linestyle = [:dash :dashdot :dot], color=:red, xaxis="Time (t)",yaxis="Reference att. accels", label=[L"\alpha_1(t)" L"\alpha_2(t)" L"\alpha_3(t)"])
    plot!(time[1:1:end-1], diff(wtraj, dims=2)'./dt, color=:black, label=[L"(\omega(t) - \omega(t-h))/h" nothing nothing])
    plotf = plot(time', ftraj', linewidth=2, linestyle = [:solid], color=:black, xaxis="Time (t)",yaxis="Reference force", label=L"f(t)")
    plot(plotq, plotw, plota, plotf, layout=(4,1), size=(1000,750));
    #savefig("../docs/images/example_reference_attitude_refs.png")
end
