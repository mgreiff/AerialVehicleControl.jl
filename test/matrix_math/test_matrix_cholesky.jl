"""
Copyright Marcus Greiff 2020

Tests for the matrix_double_cholesky() function
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
    L = rand(m,m)
    Lc = deepcopy(L)
    Lm = matrix_double_t(L)
    return A, Ac, Am, L, Lc, Lm
end

# Compact call for the matrix_double_addition function
function test_call_matrix_double_cholesky(Am, Lm)
    status = ccall((:matrix_double_cholesky, AerialVehicleControl.CONT_LIB_PATH),
        Cint,
        (Ref{matrix_double_t},Ref{matrix_double_t}),
        Am,
        Lm)
    return status
end

@verbose(1, "Testing matrix cholesky decpomosition...")
@testset "Matrix cholesky inplace" begin
    @testset "Good inputs" begin
        for m = 1:10
            # Fixture
            A, Ac, Am, L, Lc, Lm = fixture(m)
            # Compute cholesky decomposition in Julia
            julia_L = cholesky(A).L + zeros(size(A))
            # Call
            status = test_call_matrix_double_cholesky(Am, Lm)
            # Test
            @test isapprox(status, 1)
            @test isapprox(A, Ac)
            @test !isapprox(L, Lc)
            @test isapprox(julia_L*julia_L',Ac)
            @test isapprox(L*L',Ac)
        end
    end
end
