"""
Copyright Marcus Greiff 2020

Tests for the matrix_double_symmetric_real_eigenvalues() function
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
    Ac = deepcopy(A)
    Am = matrix_double_t(A)
    B  = zeros(m,1);
    Bc = deepcopy(B)
    Bm = matrix_double_t(B)
    return A, Ac, Am, B, Bc, Bm
end

# Compact call for the matrix_double_symmetric_real_eigenvalues function
function test_call_matrix_double_symmetric_real_eigenvalues(Am, Bm)
    status = ccall((:matrix_double_symmetric_real_eigenvalues, AerialVehicleControl.CONT_LIB_PATH),
        Cint,
        (Ref{matrix_double_t},Ref{matrix_double_t}),
        Am,
        Bm)
    return status
end

@verbose(1, "Matrix PSD solver...")
@testset "Matrix solve PD linear systems" begin
    @testset "Good inputs" begin
        for m = 1:10
            # Fixture
            A, Ac, Am, B, Bc, Bm = fixture(m)
            # Solve system in Julia
            X = eigvals(A);
            # Call
            status = test_call_matrix_double_symmetric_real_eigenvalues(Am, Bm)
            # Test
            @test isapprox(status, 1)
            @test isapprox(A, Ac)
            @test !isapprox(B, Bc)
            @test isapprox(X, B)
        end
    end
end
