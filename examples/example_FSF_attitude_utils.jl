################################################################################
# Copyright (C) Marcus Greiff 2020
#
# @file example_FSF_attitude_utils.jl
# @author Marcus Greiff
# @date July 2020
#
# Utility functions for initialization and plotting in the FSF example.
#
################################################################################

function get_errors_SO3(qvec, qrvec, wvec, wrvec)
    NT = size(qvec, 2);
    eRvec = zeros(3, NT);
    ewvec = zeros(3, NT);
    zvec  = zeros(2, NT);
    znorm = zeros(1, NT);
    for t = 1:NT
        R  = quat_2_SO3(reshape(qvec[:,t],  (4,1)))
        Rr = quat_2_SO3(reshape(qrvec[:,t], (4,1)))
        w  = wvec[:,t];
        wr = wrvec[:,t];

        eRvec[:,t]  = SO3_vee(Rr'*R - R'*Rr) / 2.0;
        ewvec[:,t]  = w -  R'*Rr * wr;
        zvec[1,t]   = norm(eRvec[:,t])
        zvec[2,t]   = norm(ewvec[:,t])
        znorm[1,t]  = zvec[:,t]'*zvec[:,t]

    end
    return eRvec, ewvec, zvec, znorm
end

function get_errors_SU2(qvec, qrvec, wvec, wrvec)
    NT = size(qvec, 2);
    eXvec = zeros(3, NT);
    ewvec = zeros(3, NT);
    zvec  = zeros(2, NT);
    znorm = zeros(1, NT);
    for t = 1:NT
        X  = quat_2_SU2(reshape(qvec[:,t],  (4,1)))
        Xr = quat_2_SU2(reshape(qrvec[:,t], (4,1)))
        w  = wvec[:,t];
        wr = wrvec[:,t];

        Xe = Xr'*X;

        eXvec[:,t]  = SU2_vee(Xe - tr(Xe)*Matrix{Float64}(I, 2, 2)) / 2.0;
        ewvec[:,t]  = w -  SU2_vee(Xe'*SU2_hat(wr)*Xe);
        zvec[1,t]   = norm(eXvec[:,t])
        zvec[2,t]   = norm(ewvec[:,t])
        znorm[1,t]  = zvec[:,t]'*zvec[:,t]

    end
    return eXvec, ewvec, zvec, znorm
end

function get_bounds(M1, M2, W, zvec)
    NT = size(zvec, 2);
    Vmin = zeros(1, NT);
    Vmax = zeros(1, NT);
    dV   = zeros(1, NT);
    for t = 1:NT
        z = zvec[:,t];
        Vmin[1,t] = +z'*M1*z;
        Vmax[1,t] = +z'*M2*z;
        dV[1,t]   = -z'*W*z;
    end
    return Vmin, Vmax, dV
end

function maximal_feasible_eps_SU2(
    C::con_state_qw_fsf_t,
    X0::Array{Complex{T},2},
    Xr0::Array{Complex{T},2},
    w0::Array{T,2},
    wr0::Array{T,2};
    isContinuous=true
) where {T<:Number}

    # Extract parameters
    kX, kc, kw  = C.gain_kR, C.gain_kc, C.gain_kw;
    J           = ntuple_2_mat(C.inertia, 3, 3);

    if isContinuous
        # Compute maximum permissible phi
        Gamma     = abs(0.5 * tr(Matrix{Float64}(I, 2, 2) - Xr0'*X0));
        Xe = Xr0'*X0;
        eX = SU2_vee(Xe - tr(Xe)*Matrix{Float64}(I, 2, 2)) / 2.0
        ew = w0 -  SU2_vee(Xe'*SU2_hat(wr0)*Xe)
        attpart   = kX * Gamma;
        crosspart = kc * ew' * eX;
        ratepart  = ew'*J*ew / 2.0;
        Lyap0     = attpart + crosspart[1,1] + ratepart[1,1];
        phi       = Lyap0/kX;
    else
        # phi is fixed at 1.0 in the discontinuous case
        phi = 1.0
    end

    M1, M2, W   = get_M1M2W_SU2(C, phi)

    M1max, M1min = maximum(eigvals(M1)), minimum(eigvals(M1));
    M2max, M2min = maximum(eigvals(M2)), minimum(eigvals(M2));
    Wmax,  Wmin  = maximum(eigvals(W)),  minimum(eigvals(W));

    mindecayrate = Wmin/M2max;
    maxeps       = ((M1min*Wmin) / M2max)*phi*(2.0 - phi);
    uniboundgain = (M2max / (M1min*Wmin));

    return maxeps, mindecayrate, uniboundgain, phi, M1, M2, W
end

function get_info(
    C::con_state_qw_fsf_t,
    controllerType::Integer,
    q0::Array{T,2},
    qr0::Array{T,2},
    w0::Array{T,2},
    wr0::Array{T,2}
) where {T<:Number}

    controllerStrings = [
        "continuous FSF controller on SO(3) (called from C)",
        "robust FSF controller on SO(3) (called from C)",
        "continuous FSF controller on SU(2) (called from C)",
        "discontinuous FSF controller on SU(2) (called from C)",
        "robust FSF controller on SU(2) (called from C)"
    ]

    nstr = length(controllerStrings[controllerType]);
    println("\n$('~'^(34+nstr))")
    println("~~~ Simulation done with the $(controllerStrings[controllerType]) ~~~")
    println("$('~'^(34+nstr))")

    if controllerType in [1,2]
        R0  = quat_2_SO3(q0);
        Rr0 = quat_2_SO3(qr0);

        phi0 = tr(I-Rr0'*R0)/2.0;
        println("* The initial attitude error is: Psi(Rr(t0), R(t0)) = $(phi0)")
        if (phi0 > 1.5)
            println("* Warning: This is relatively large, requiring very small initial attitude rate errors")
        end
        maxeps, mindecayrate, uniboundgain, phi, M1, M2, W = maximal_feasible_eps_SO3(C, R0, Rr0, w0, wr0)
    end
    if controllerType in [3,4,5]
        X0  = quat_2_SU2(q0);
        Xr0 = quat_2_SU2(qr0);

        phi0 = abs(0.5 * tr(Matrix{Float64}(I, 2, 2) - Xr0'*X0));
        println("* The initial attitude error is: Gamma(Xr(t0), X(t0)) = $(phi0)")
        if (phi0 > 1.5)
            println("* Warning: This is relatively large, requiring very small initial attitude rate errors")
        end
        if controllerType in [3,5]
            maxeps, mindecayrate, uniboundgain, phi, M1, M2, W = maximal_feasible_eps_SU2(C, X0, Xr0, w0, wr0; isContinuous=true)
        end
        if controllerType in [4]
            maxeps, mindecayrate, uniboundgain, phi, M1, M2, W = maximal_feasible_eps_SU2(C, X0, Xr0, w0, wr0; isContinuous=false)
        end
    end

    println("* Any V(t0)/kR = $(phi) < phi < 2 can be used in the stability proof. We let phi = $(phi).")
    if (phi > 2);
        println("* Waring. The system may move outside of the domain of exponential attraction given phi")
        println("* Lyap0/kR = $(phi), please retune the controller, or lower the initial rate error")
    end

    println("* Worst case decay rate of the Lyapunov function : $(mindecayrate)")
    if controllerType in [2,5]
        println("* Upper bound on the allowed epsilon given phi   : $(maxeps)")
        println("* Uniform ultimate bound K*eps, where K is       : $(uniboundgain)")
    end
    return M1, M2, W, uniboundgain
end

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
        maxkc = maximal_feasible_kc_SU2(kR, kw, J); # here we let kR=kX for simplicity
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
    w0 = rand(3,1).-0.5;

    # Select a feasible epsilon gain, with the uptimate bound characterize by alpha
    alpha = 1e-2;
    if controllerType == 2
        # In case of an SO(3) configured controller (only applicable in the robust case)
        maxeps, _, _, _, _, _ = maximal_feasible_eps_SO3(C, quat_2_SO3(q0), quat_2_SO3(qr0), w0, wr0);
        C.gain_eps = alpha*maxeps;
    end
    if controllerType == 5
        # In case of an SU(2) configured controller (only applicable in the robust case)
        maxeps, _, _, _, _, _ = maximal_feasible_eps_SU2(C, quat_2_SU2(q0), quat_2_SU2(qr0), w0, wr0);
        C.gain_eps = alpha*maxeps;
    end

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
    namePlot="unnamed",
    showDisturbance=false
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

    # The first element of certain logged variables (computed in the C-code)
    # are not set at t0, hence we let tis value be the one at the next time step
    # to tidy up and simplify the plotting.
    psi[:,1], gamma[:,1] = psi[:,2], gamma[:,2];
    tau[:,1] = tau[:,2];
    taur[:,1] = taur[:,2];
    dist[:,1] = dist[:,2];
    lyap[:,1] = lyap[:,2];


    q0  = reshape(q[:,1],  (4,1));
    qr0 = reshape(qr[:,1], (4,1));
    w0  = reshape(w[:,1],  (3,1));
    wr0 = reshape(wr[:,1], (3,1));

    M1, M2, W, K = get_info(C, controllerType, q0, qr0, w0, wr0);

    # In case we are dealine with an SO(3) controller (continuous, or robust)
    if controllerType in [1,2]
        eAtt, eW, z, znorm = get_errors_SO3(q, qr, w, wr);
    end
    # In case we are dealine with an SU(2) controller
    if controllerType in [3,4,5]
        eAtt, eW, z, znorm = get_errors_SU2(q, qr, w, wr)
    end

    Vmin, Vmax, dV   = get_bounds(M1, M2, W, z);
    zbound           = ones(length(t),1).*K.*C.gain_eps;

    if 1 in showPlot
        pqqr = plot( t, q',    color=:black, linewidth=2, xaxis="Time (t)",yaxis="Attitude [.]",          label=[L"q(t)" nothing nothing nothing])
               plot!(t, qr',   color=:red,   linewidth=2, xaxis="Time (t)",yaxis="Attitude [.]",          label=[L"q_r(t)" nothing nothing nothing])
        pwwr = plot( t, w',    color=:black, linewidth=2, xaxis="Time (t)",yaxis="Attitude rate [rad/s]", label=[L"\omega(t)" nothing nothing])
               plot!(t, wr',   color=:red,   linewidth=2, xaxis="Time (t)",yaxis="Attitude rate [rad/s]", label=[L"\omega_{r}(t)" nothing nothing])
        pttr = plot( t, tau',  color=:black, linewidth=2, xaxis="Time (t)",yaxis="Torque [Nm]",           label=[L"\tau(t)" nothing nothing])
               plot!(t, taur', color=:red,   linewidth=2, xaxis="Time (t)",yaxis="Torque [Nm]",           label=[L"\tau_{r}(t)" nothing nothing])
        pdis = plot( t, dist', color=:black, linewidth=2, xaxis="Time (t)",yaxis="Load disturbance [Nm]", label=[L"d(t)" nothing nothing])

        if showDisturbance
            plot(pqqr, pwwr, pttr, pdis, layout=(4,1), size=(1000,1000))
        else
            plot(pqqr, pwwr, pttr, layout=(3,1), size=(1000,750))
        end

        gui()
        if savePlot
            savefig(joinpath(AerialVehicleControl.CONT_LIB_DOCS, "$(namePlot)_states.png"))
        end
    end

    if 2 in showPlot

        attitudeErrors = [L"e_{R}(t)", L"e_{R}(t)", L"e_{X}(t)", L"e_{X}^{\pm}(t)", L"e_{X}(t)"]
        plotattd    = plot( t, psi',        color=:black,             linewidth=2,                            label=L"\Psi(R_r(t), R(t))")
                      plot!(t, gamma',      color=:red,               linewidth=2,                            label=L"\Gamma(X_r(t), X(t))")
                      plot!(t, 2 .- gamma', color=:red,  style=:dash, linewidth=2, yaxis="Attitude distance", label=L"\bar{\Gamma}(X_r(t), X(t))")
        plotatterr  = plot( t, eAtt', color=:black, linewidth=2,                   yaxis="Attitude error",    label=[attitudeErrors[controllerType] nothing nothing])
        plotrateerr = plot( t, eW',   color=:black, linewidth=2, xaxis="Time (t)", yaxis="Rate error",        label=[L"e_{\omega}(t)" nothing nothing])
        plot(plotattd, plotatterr,  plotrateerr, layout=(3,1), size=(1000,750))
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
                        plot!(tshort, lyapshort,                color=sigcolor,   linewidth=2, xaxis="Time (t)",yaxis=L"\frac{d}{dt}\mathcal{V}(t)",       label=L"$({d}/{dt}){\mathcal{V}}(t)$")

        log10zbound   = zbound[:]; #log10.(abs.(zbound[:]))
        log10znorm    = znorm[:];  #log10.(abs.(znorm[:]))
        lbval         = minimum([log10zbound; log10znorm]) - 1.0;
        bounds        = (log10zbound .- lbval, 0 .*log10zbound);
        plotznorm     = plot(t,  log10zbound, ribbon=bounds, color=boundcolor, linewidth=2, style=ubstyle, fillcolor=fillcolor, fillalpha=boundalpha, label=L"$\lambda_M(M_2)(\lambda_m(M_1)\lambda_m(W))^{-1}\epsilon$")
                        plot!(t, log10znorm,                 color=sigcolor,   linewidth=2, xaxis="Time (t)",yaxis=L"$\|\|z(t)\|\|_2^2$", label=L"$\|\|z(t)\|\|_2^2$")

        if showDisturbance
            plot(plotlyap, plotlyaplog, plotlyapderiv, plotznorm, layout=(4,1), size=(1000,1000))
        else
            plot(plotlyap, plotlyaplog, plotlyapderiv, layout=(3,1), size=(1000,750))
        end
        gui()
        if savePlot
            savefig(joinpath(AerialVehicleControl.CONT_LIB_DOCS, "$(namePlot)_analysis.png"))
        end
    end
end
