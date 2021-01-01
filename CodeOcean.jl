################################################################################
# Copyright (C) Marcus Greiff 2020
#
# @file CodeOcean.jl
# @author Marcus Greiff
# @date January 2021
#
# This file simply one of the relevant Julia examples for the LCSS-publication,
# it runs all of the implemented examples pertaining to the FSF attitude
# control, and saves the plots to the ../results directory.
################################################################################

# This is needed to make GR work in CodeOcean
ENV["GKSwstype"] = "100";
# Load and activate the AerialVehicleControl.jl package
using Pkg; Pkg.activate(".");

################################################################################
############################# Run the tests ####################################
################################################################################
println("Running tests...")
include("test/runtests.jl")

################################################################################
################### Run the FSF controller Feedback examples ###################
################################################################################
fileNames = ["continuous_SO3", "robust_SO3", "continuous_SU2", "discontinuous_SU2", "robust_SU2"]
applyDist = [false, true, false, false, true];

for i = 1:5
    # Simulation with a sought controller type
    global controllerType  = i;

    # Plotting settings
    global showDisturbance = applyDist[i];
    global showPlot        = [1,2,3];
    global savePlot        = true;
    global namePlot        = fileNames[i];

    # Simulation setting
    global useDisturbance = showDisturbance;

    # Run example
    include("examples/example_FSF_attitude.jl")
end
