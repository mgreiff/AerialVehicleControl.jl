################################################################################
# Copyright (C) Marcus Greiff 2020
#
# @file example_attitude_FSF.jl
# @author Marcus Greiff
# @date July 2020
#
# This example contains the code necessary to call the FSF attitude controllers
# in the loop for simulation and verification by the DifferentialEquations.jl
# stack. By wrapping the functions update_attitude_FSF_*_* un the update_control
# function, this is subsequently called in the odefun!() for simulation Tsit5().
#
# Currently, any of the functions below can be used in this example
#
#     1. update_attitude_FSF_SO3_continuous()
#     2. update_attitude_FSF_SO3_robust()
#     3. update_attitude_FSF_SU2_continuous()
#     4. update_attitude_FSF_SU2_discontinuous()
#     5. update_attitude_FSF_SU2_robust()
#
# and you can choose which controller to use by setting the controllerType flag
# to the number before each of the functions. For instance, if you would like to
# simulate and analyze the continuous FSF on SU(2), let controllerType=3.
#
# Note thet you can also define disturbances to test the controller. These can
# be written to the function attitude_disturbance(), which takes the maximum
# allowed gain in the l2-norm that you have permitted in the controller tuning
# as well as time as arguments. However, you could also include feedback in the
# states here (provided that you smoothly saturate these). If no disturbance is
# to be applied, the function attitude_disturbance() should simply return an
# array of zeros.
#
# Finally, you can use the solution_analysis() function to analyze the computed
# solution. This function takes the ODS solution and generates any combination
# of three plots, with
#
#     1. State trajectories
#     2. Error trajectories
#     3. Analysis plots with theoretical bounds
#
# The kwarg "showPlot" is used to selevtively generate plots, the "savePlot"
# flag indicates if the plots should be save to file, and if so, they are saved
# to /docs/images/ under the filename specified by the kwarg "namePlot". For
# instance, if you wish to save the state and analysis plots, you would write
#
#     solution_analysis(...; showPlot=[1,3], savePlot=true, namePlot="myplot")
#
################################################################################
using DifferentialEquations
using LinearAlgebra
using LaTeXStrings
using AerialVehicleControl
using Plots
import Random

# Set seed, include util functions for example, and recompile the C-code
include("example_attitude_FSF_utils.jl")
Random.seed!(0)
recompile()

# Set the controller type
controllerType = 1;

# Disturbance (must be smaller than L in the 2-norm) used in odefun!()
function attitude_disturbance(L::Float64, t::Float64)
    useDisturbance = false
    if useDisturbance
        d = [sin(t), cos(t), sin(3*t)];
        d = d./norm(d);
        d = 0.99*d * L;
    else
        d = zeros(3,);
    end
    return d
end

# Wrapper for the calling of the hosen update law from the C-code
function update_FSF_attitude_control(
    R::ref_state_qw_t,
    S::dyn_state_qw_t,
    C::con_state_qw_fsf_t,
    controllerType::Integer
)
    @assert(controllerType in 1:5, "Underfined controllerType=$(controllerType)")
    if controllerType == 1
        return ccall((:update_attitude_FSF_SO3_continuous, AerialVehicleControl.CONT_LIB_PATH),
            Cint, (Ref{ref_state_qw_t}, Ref{dyn_state_qw_t}, Ref{con_state_qw_fsf_t},), R, S, C)
    end
    if controllerType == 2
        return ccall((:update_attitude_FSF_SO3_robust, AerialVehicleControl.CONT_LIB_PATH),
            Cint, (Ref{ref_state_qw_t}, Ref{dyn_state_qw_t}, Ref{con_state_qw_fsf_t},), R, S, C)
    end
    if controllerType == 3
        return ccall((:update_attitude_FSF_SU2_continuous, AerialVehicleControl.CONT_LIB_PATH),
            Cint, (Ref{ref_state_qw_t}, Ref{dyn_state_qw_t}, Ref{con_state_qw_fsf_t},), R, S, C)
    end
    if controllerType == 4
        return ccall((:update_attitude_FSF_SU2_discontinuous, AerialVehicleControl.CONT_LIB_PATH),
            Cint, (Ref{ref_state_qw_t}, Ref{dyn_state_qw_t}, Ref{con_state_qw_fsf_t},), R, S, C)
    end
    if controllerType == 5
        return ccall((:update_attitude_FSF_SU2_robust, AerialVehicleControl.CONT_LIB_PATH),
            Cint, (Ref{ref_state_qw_t}, Ref{dyn_state_qw_t}, Ref{con_state_qw_fsf_t},), R, S, C)
    end
end

# Function handle for the numerical integration
function odefun!(dx,x,params,t)

    # Get the controller object and the chosen controller type
    C, controllerType = params;

    # Extract the states
    q  = x[1:4];
    w  = x[5:7];
    qr = x[8:11];
    wr = x[12:14];

    # Extract the parameters
    J = ntuple_2_mat(C.inertia, 3, 3)
    L = C.gain_L;

    # Reference torques driving the reference dynamics
    taur = [sin(t), sin(2*t), sin(3*t)];

    # Reference dynamics
    dqr, dwr = attitude_dynamics(qr, wr, J, taur)

    # Create state structs for ccall
    R = ref_state_qw_t(hcat(qr), hcat(wr), hcat(dwr), 0.5)
    S = dyn_state_qw_t(hcat(q ), hcat(w ))

    # Compute control signal through a ccall to the FSF configured on SO3
    status = update_FSF_attitude_control(R, S, C, controllerType)

    # Actuate the dynamics with the computed control signal
    tau  = [C.torques[i] for i in 1:3]

    # Realization of the additive load disturbance (if applicable)
    disturbance = attitude_disturbance(L, t);
    @assert(norm(disturbance) <= L, "Disturbance exceeds L in l2-norm");

    # Physical dynamics
    dq , dw  = attitude_dynamics( q,  w, J, tau; delta=disturbance )

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
    x[24:26]  = disturbance;
    x[27]     = norm(q);
    x[28]     = norm(qr);
end

# Compute solution
x0, tspan, C = initialize_FSF_attitude_example(controllerType)
parameters   = (C, controllerType)
prob         = ODEProblem(odefun!, x0, tspan, parameters);
sol          = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-21)

# Visualize solutioncontrollerType
solution_analysis(sol, C, controllerType; showPlot=3)
