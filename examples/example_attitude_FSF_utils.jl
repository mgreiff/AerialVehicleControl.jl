################################################################################
# Copyright (C) Marcus Greiff 2020
#
# @file example_attitude_FSF_utils.jl
# @author Marcus Greiff
# @date July 2020
#
# Utility functions for initialization and plotting in the FSF example.
#
################################################################################
function initialize_FSF_attitude_example(controllerType::Integer)

    # Define the controller parameters
    J  = rand_PSD(3);
    J  = J./opnorm(J);
    kR = 5.0;
    kw = 2.0;

    # Compute a feasible cross gain kc
    if controllerType in [1,2]
        # In case of an SO(3) configured controller (continuous or robust)
        maxkc = maximal_feasible_kc_SO3(kR, kw, J);
    else
        error("Not yet implemented")
    end
    kc = 0.1 * maxkc;

    # Define the controller state object
    C  = con_state_qw_fsf_t(J, kR, kc, kw)

    # Initial states of the reference system
    qr0 = rand(4,1).-0.5;
    qr0 = qr0./norm(qr0);
    wr0 = zeros(3,1);

    # Initial states of the physical system
    q0 = rand(4,1).-0.5;
    q0 = q0./norm(q0);
    w0 = 0.4*rand(3,1);

    # Select a feasible epsilon gain, with the uptimate bound characterize dby alpha
    if controllerType in [1,2]
        # In case of an SO(3) configured controller (continuous or robust)
        maxeps, _, _, _, _, _ = maximal_feasible_eps_SO3(C, quat_2_SO3(q0), quat_2_SO3(qr0), w0, wr0);
    else
        # In case of an SU(2) configured controller (continuous, discontinuous, or robust)
        error("Not yet implemented")
    end
    alpha = 1e-2;
    C.gain_eps = alpha*maxeps;

    # Initial states for the debugging signals
    output = zeros(14,1);
    output[13,1] = 1.0;
    output[14,1] = 1.0;

    x0    = vcat(q0, w0, qr0, wr0, output)
    tspan = (0.0,20.0)
    return x0, tspan, C
end

function solution_analysis(
    sol::ODESolution,
    C::con_state_qw_fsf_t,
    controllerType::Integer;
    showPlot=[1,2,3],
    savePlot=false,
    namePlot="unnamed"
)

    # Posterior solution analysis
    t     = sol.t;
    data  = zeros(length(t),length(sol.u[1]));
    for i = 1:length(t); data[i,:] = sol.u[i]; end

    q     = data[:,1:4]';
    w     = data[:,5:7]';
    qr    = data[:,8:11]';
    wr    = data[:,12:14]';
    tau   = data[:,15:17]';
    taur  = data[:,18:20]';
    psi   = data[:,21]';
    gamma = data[:,22]';
    dist  = data[:,24:26]';
    lyap  = data[:,23]';
    lyapd = numerical_differentiation(lyap, t);

    # In case we are dealine with an SO(3) controller
    if controllerType in [1,2]
        R0  = quat_2_SO3(reshape(q[:,1],  (4,1)));
        Rr0 = quat_2_SO3(reshape(qr[:,1], (4,1)));
        w0  = reshape(w[:,1],  (3,1));
        wr0 = reshape(wr[:,1], (3,1));

        M1, M2, W, K     = get_info_SO3(C, R0, Rr0, w0, wr0);
        eR, ew, z, znorm = get_errors_SO3(q, qr, w, wr);
        Vmin, Vmax, dV   = get_bounds_SO3(M1, M2, W, z);
        zbound           = ones(length(t),1).*K.*C.gain_eps;
    end

    # In case we are dealine with an SU(2) controller
    if controllerType in [3,4,5]
        error("Not yet implemented")
    end

    if 1 in showPlot
        pqqr = plot( t, q',    color=:black, linewidth=2, xaxis="Time (t)",yaxis="Attitude [.]",          label=[L"q(t)" nothing nothing nothing])
               plot!(t, qr',   color=:red,   linewidth=2, xaxis="Time (t)",yaxis="Attitude [.]",          label=[L"q_r(t)" nothing nothing nothing])
        pwwr = plot( t, w',    color=:black, linewidth=2, xaxis="Time (t)",yaxis="Attitude rate [rad/s]", label=[L"\omega(t)" nothing nothing])
               plot!(t, wr',   color=:red,   linewidth=2, xaxis="Time (t)",yaxis="Attitude rate [rad/s]", label=[L"\omega_{r}(t)" nothing nothing])
        pttr = plot( t, tau',  color=:black, linewidth=2, xaxis="Time (t)",yaxis="Torque [Nm]",           label=[L"\tau(t)" nothing nothing])
               plot!(t, taur', color=:red,   linewidth=2, xaxis="Time (t)",yaxis="Torque [Nm]",           label=[L"\tau_{r}(t)" nothing nothing])
        pdis = plot( t, dist', color=:black, linewidth=2, xaxis="Time (t)",yaxis="Load disturbance [Nm]", label=[L"d(t)" nothing nothing])
        plot(pqqr, pwwr, pttr, pdis, layout=(4,1), size=(1000,1000))
        gui()
        if savePlot
            savefig(joinpath(AerialVehicleControl.CONT_LIB_DOCS, "$(namePlot)_states.png"))
        end
    end

    if 2 in showPlot
        plotattd    = plot( t, psi',   color=:black, linewidth=2, xaxis="Time (t)",yaxis="Attitude error",    label=L"\Psi(R_r(t), R(t))")
                      plot!(t, gamma', color=:red,   linewidth=2, xaxis="Time (t)",yaxis="Attitude error",    label=L"\Gamma(X_r(t), X(t))")
        plotlyap    = plot( t, lyap',               color=:black, linewidth=2, xaxis="Time (t)",yaxis="Lyapunov function", label=L"\mathcal{V}(t)")
        plotlyaplog = plot( t, log10.(abs.(lyap')), color=:black, linewidth=2, xaxis="Time (t)",yaxis="Lyapunov function", label=L"log_{10}(\mathcal{V}(t))")
        plot(plotattd, plotlyap, plotlyaplog, layout=(3,1), size=(1000,750))
        gui()
        if savePlot
            savefig(joinpath(AerialVehicleControl.CONT_LIB_DOCS, "$(namePlot)_errors.png"))
        end
    end

    if 3 in showPlot
        ubstyle    = :dash;
        lbstyle    = :dot;
        sigcolor   = :black;
        boundcolor = :red;
        boundalpha = 0.1;
        fillcolor  = :black;

        bounds     = (Vmin[:].*0, Vmax[:]-Vmin[:]);
        plotlyap   = plot(t,  Vmin[:], ribbon=bounds, color=boundcolor, linewidth=2, style=lbstyle, fillcolor=fillcolor, fillalpha=boundalpha, label=L"z^{\top}M_1z")
                     plot!(t, Vmax[:],                color=boundcolor, linewidth=2, style=ubstyle,                                            label=L"z^{\top}M_2z")
                     plot!(t, lyap[:],                color=sigcolor,   linewidth=2, xaxis="Time (t)",yaxis=L"\mathcal{V}(t)",                 label=L"$\mathcal{V}(t)$")

        logVmin = log10.(abs.(Vmin[:]))
        logVmax = log10.(abs.(Vmax[:]))
        logLyap = log10.(abs.(lyap[:]))
        bounds  = (logVmin.*0, logVmax - logVmin)
        plotlyaplog = plot(t,  logVmin, ribbon=bounds, color=boundcolor, linewidth=2, style=lbstyle, fillcolor=fillcolor, fillalpha=boundalpha, label=L"$log_{10}(z^{\top}M_1z)$")
                      plot!(t, logVmax,                color=boundcolor, linewidth=2, style=ubstyle,                                            label=L"$log_{10}(z^{\top}M_2z)$")
                      plot!(t, logLyap,                color=sigcolor,   linewidth=2, xaxis="Time (t)",yaxis=L"log_{10}(\mathcal{V}(t))",       label=L"$log_{10}(\mathcal{V}(t))$")

        ind           = 3;
        tshort        = t[ind:end];
        dVshort       = dV[1,ind:end];
        lyapshort     = lyapd[1,ind:end];
        lbval         = 1.1*minimum([dVshort; lyapshort])
        bounds        = (dVshort .- lbval, 0 .* dVshort)
        plotlyapderiv = plot(tshort,    dVshort, ribbon=bounds, color=boundcolor, linewidth=2, style=ubstyle, fillcolor=fillcolor, fillalpha=boundalpha,   label=L"-z^{\top}Wz")
                        plot!(tshort, lyapshort,                color=sigcolor,   linewidth=2, xaxis="Time (t)",yaxis=L"\frac{d}{dt}\mathcal{V}(t)",       label=L"$\dot{\mathcal{V}}(t)$")

        log10zbound   = log10.(abs.(zbound[:]))
        log10znorm    = log10.(abs.(znorm[:]))
        lbval         = minimum([log10zbound; log10znorm]) - 1.0;
        bounds        = (log10zbound .- lbval, 0 .*log10zbound);
        plotznorm     = plot(t,  log10zbound, ribbon=bounds, color=boundcolor, linewidth=2, style=ubstyle, fillcolor=fillcolor, fillalpha=boundalpha, label=L"$log_{10}(\lambda_M(M_2)(\lambda_m(M_1)\lambda_m(W))^{-1}\epsilon)$")
                        plot!(t, log10znorm,                 color=sigcolor,   linewidth=2, xaxis="Time (t)",yaxis=L"$log_{10}(\|\|z(t)\|\|_2^2)$", label=L"$log_{10}(\|\|z(t)\|\|_2^2)$")

        plot(plotlyap, plotlyaplog, plotlyapderiv, plotznorm, layout=(4,1), size=(1000,750))
        gui()
        if savePlot
            savefig(joinpath(AerialVehicleControl.CONT_LIB_DOCS, "$(namePlot)_analysis.png"))
        end
    end
end
