"""
Copyright Marcus Greiff 2020

Tests for the matrix_transposition() function
"""

# Include modules and recompile if running test standalone
if !isdefined(Main, :runAllTests)
    using Test, AerialVehicleControl, LinearAlgebra
    recompile();
end

# Testing fixture
function fixture(mA,nA,mB,nB)
    A  = rand(Cdouble,mA,nA)
    Ac = deepcopy(A)
    Am = matrix_double_t(A)
    B  = rand(Cdouble,mB,nB)
    Bc = deepcopy(B)
    Bm = matrix_double_t(B)
    return A, Ac, Am, B, Bc, Bm
end

# Compact call for the matrix_double_addition function
function test_call_matrix_double_transposition(Am, Bm)
    status = ccall((:matrix_double_transposition, AerialVehicleControl.CONT_LIB_PATH),
        Cint,
        (Ref{matrix_double_t}, Ref{matrix_double_t}),
        Am,
        Bm)
    return status
end

@verbose(1, "Matrix transposition...")
@testset "Matrix transposition" begin
    @testset "Good inputs" begin
        for m = 1:10
            for n = 1:10
                # Fixture
                mA, nA, mB, nB = m, n, n, m
                A, Ac, Am, B, Bc, Bm = fixture(mA,nA,mB,nB)
                # Call
                status = test_call_matrix_double_transposition(Am, Bm)
                # Test
                @test isapprox(status, 1)
                @test isapprox(A, Ac)
                @test !isapprox(B, Bc)
                @test isapprox(B, A')
            end
        end
    end
end
