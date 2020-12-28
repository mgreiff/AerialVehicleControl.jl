"""
Copyright Marcus Greiff 2020

Tests for the matrix_solve_posdef() function
"""

# Include modules and recompile if running test standalone
if !isdefined(Main, :runAllTests)
    using Test, AerialVehicleControl, LinearAlgebra
    recompile();
end

# Testing fixture
function fixture(m,n)
    A = rand(m,m)
    A = A + A'
    A = A + (abs(minimum(eigvals(A))) + 0.01) * Matrix{Float64}(I, m, m)
    Ac = deepcopy(A)
    Am = matrix_double_t(A)
    B = rand(m,n)
    Bc = deepcopy(B)
    Bm = matrix_double_t(B)
    return A, Ac, Am, B, Bc, Bm
end

# Compact call for the matrix_double_addition function
function test_call_matrix_double_solve_posdef(Am, Bm)
    status = ccall((:matrix_double_solve_posdef, AerialVehicleControl.CONT_LIB_PATH),
        Cint,
        (Ref{matrix_double_t},Ref{matrix_double_t}),
        Am,
        Bm)
    return status
end

@verbose(1, "Testing matrix PSD solver...")
@testset "Matrix solve PD linear systems" begin
    @testset "Good inputs" begin
        for m = 1:10
            for n = 1:10
                # Fixture
                A, Ac, Am, B, Bc, Bm = fixture(m, n)
                # Solve system in Julia
                X = A\B;
                # Call
                status = test_call_matrix_double_solve_posdef(Am, Bm)
                # Test
                @test isapprox(status, 1)
                @test !isapprox(A, Ac)
                @test !isapprox(B, Bc)
                @test isapprox(X,B)
                @test isapprox(Ac*B,Bc)
            end
        end
    end
end
