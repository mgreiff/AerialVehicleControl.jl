"""
Copyright Marcus Greiff 2020

Full test suite for the aerial vehicle controller generation
"""

using LinearAlgebra, Test, AerialVehicleControl

recompile()        # Recompile the C-code
runAllTests = true # If defined, will not recompile the C-code

@testset "All tests" begin
    @testset "Math library tests" begin
        include("controller_math/test_math_functions.jl")
        include("controller_math/test_SO3_maps.jl")
        include("controller_math/test_SU2_maps.jl")
    end
    @testset "Attitude controllers" begin
        include("controller_implementations/test_FSF_continuous_SO3.jl")
        include("controller_implementations/test_FSF_continuous_SU2.jl")
        include("controller_implementations/test_FSF_discontinuous_SU2.jl")
        include("controller_implementations/test_FOF_continuous_SO3.jl")
    end
    @testset "Utility functions" begin
        include("controller_utils/test_power_distribution.jl")
    end
    @testset "Matrix math functions" begin
        include("matrix_math/test_matrix_addition_inplace.jl")
        include("matrix_math/test_matrix_addition.jl")
        include("matrix_math/test_matrix_cholesky_inplace.jl")
        include("matrix_math/test_matrix_cholesky.jl")
        include("matrix_math/test_matrix_eigvals.jl")
        include("matrix_math/test_matrix_multiplication.jl")
        include("matrix_math/test_matrix_solve_posdef.jl")
        include("matrix_math/test_matrix_subtraction_inplace.jl")
        include("matrix_math/test_matrix_subtraction.jl")
        include("matrix_math/test_matrix_transposition.jl")
    end
end
