"""
Copyright Marcus Greiff 2020

Tests for the matrix_double_cholesky_inplace() function
"""

# Include modules and recompile if running test standalone
if !isdefined(Main, :runAllTests)
    using Test, AerialVehicleControl, LinearAlgebra
    recompile();
end

# Testing fixture
function fixture(m)
    A = rand(m,m)
    A = A + A'
    A = A + (abs(minimum(eigvals(A))) + 0.01) * Matrix{Float64}(I, m, m)
    Ac = deepcopy(A)
    Am = matrix_double_t(A)
    return A, Ac, Am
end

# Compact call for the matrix_double_addition function
function test_call_matrix_double_cholesky_inplace(Am)
    status = ccall((:matrix_double_cholesky_inplace, AerialVehicleControl.CONT_LIB_PATH),
        Cint,
        (Ref{matrix_double_t},),
        Am)
    return status
end

@verbose(1, "Matrix cholesky decpomosition (inplace)...")
@testset "Matrix cholesky inplace" begin
    @testset "Good inputs" begin
        for m = 1:10
            # Fixture
            A, Ac, Am = fixture(m)
            # Compute cholesky decomposition in Julia
            julia_L = cholesky(A).L + zeros(size(A))
            # Call
            status = test_call_matrix_double_cholesky_inplace(Am)
            # Test
            @test isapprox(status, 1)
            @test !isapprox(A, Ac)
            @test isapprox(julia_L*julia_L',Ac)
            @test isapprox(A*A',Ac)
        end
    end
end