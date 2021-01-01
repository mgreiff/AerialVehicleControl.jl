"""
Copyright Marcus Greiff 2020

The Aerial Vehicle Control package
"""
module AerialVehicleControl

export @verbose,
    recompile,
    get_matrix_and_copy,
    skew,
    SO3_hat,
    SO3_vee,
    SU2_hat,
    SU2_vee,
    quat_2_SU2,
    SU2_2_quat,
    quat_2_SO3,
    quaternion_left_product,
    attitude_dynamics,
    ntuple_2_mat,
    rand_PSD,
    allFieldsEqual,
    matrix_double_t,
    ref_state_qw_t,
    dyn_state_qw_t,
    con_state_qw_fsf_t,
    con_state_qw_fof_t,
    maximal_feasible_kc_SO3,
    maximal_feasible_kc_SU2,
    maximal_feasible_eps_SO3,
    get_M1M2W_SO3,
    get_M1M2W_SU2,
    numerical_differentiation

#using Plots
using LaTeXStrings
using LinearAlgebra
using Libdl
using Test

const CONT_LIB_BASE = dirname(@__DIR__)
const CONT_LIB_SRC  = joinpath(CONT_LIB_BASE, "src")
const CONT_LIB_UTIL = joinpath(CONT_LIB_BASE, "utils")
const CONT_LIB_PATH = joinpath(CONT_LIB_SRC,  "tests")
const CONT_LIB_DOCS = joinpath(CONT_LIB_BASE, joinpath("docs", "images"))

# Compilation and recompilation of the C-code
include(joinpath(CONT_LIB_UTIL, "cont_compile.jl"))
# The unmutable Julia types mirroring this in C
include(joinpath(CONT_LIB_UTIL, "cont_types.jl"))
# The help functions used in the tests and examples
include(joinpath(CONT_LIB_UTIL, "cont_help_functions.jl"))

end
