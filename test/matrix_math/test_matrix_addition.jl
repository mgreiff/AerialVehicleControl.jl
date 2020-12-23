"""
Copyright Marcus Greiff 2020

Tests for the matrix_double_addition() function
"""

# Include modules and recompile if running test standalone
if !isdefined(Main, :runAllTests)
    using Test, AerialVehicleControl, LinearAlgebra
    recompile();
end

# Testing fixture
function fixture(m,n)
    A  = rand(Cdouble,m,n)
    Ac = deepcopy(A)
    Am = matrix_double_t(A)
    B  = rand(Cdouble,m,n)
    Bc = deepcopy(B)
    Bm = matrix_double_t(B)
    C  = rand(Cdouble,m,n)
    Cc = deepcopy(C)
    Cm = matrix_double_t(C)
    return A, Ac, Am, B, Bc, Bm, C, Cc, Cm
end

# Compact call for the matrix_double_addition function
function test_call_matrix_double_addition(Am, Bm, Cm)
    status = ccall((:matrix_double_addition, AerialVehicleControl.CONT_LIB_PATH),
        Cint,
        (Ref{matrix_double_t}, Ref{matrix_double_t}, Ref{matrix_double_t}),
        Am,
        Bm,
        Cm)
    return status
end

@verbose(1, "Matrix addition...")
@testset "Matrix addition" begin
    @testset "Good inputs" begin
        for m = 1:10
            for n = 1:10
                # Fixture
                A, Ac, Am, B, Bc, Bm, C, Cc, Cm = fixture(m,n)
                # Call
                status = test_call_matrix_double_addition(Am, Bm, Cm)
                # Test
                @test isapprox(status, 1)
                @test isapprox(A, Ac)
                @test isapprox(B, Bc)
                @test !isapprox(C, Cc)
                @test isapprox(C, A + B)
            end
        end
    end
    @testset "Bad inputs" begin
        # Fixture
        A, Ac, Am, B, Bc, Bm, C, Cc, Cm = fixture(5,5)
        Am.numRows = 0
        # Call
        status = test_call_matrix_double_addition(Am, Bm, Cm)
        # Test
        @test isapprox(status, 0)

        # Fixture
        A, Ac, Am, B, Bc, Bm, C, Cc, Cm = fixture(5,5)
        Am.numCols = 0
        # Call
        status = test_call_matrix_double_addition(Am, Bm, Cm)
        # Test
        @test isapprox(status, 0)

        # Fixture
        A, Ac, Am, B, Bc, Bm, C, Cc, Cm = fixture(5,5)
        Bm.numRows = 0
        # Call
        status = test_call_matrix_double_addition(Am, Bm, Cm)
        # Test
        @test isapprox(status, 0)

        # Fixture
        A, Ac, Am, B, Bc, Bm, C, Cc, Cm = fixture(5,5)
        Bm.numCols = 0
        # Call
        status = test_call_matrix_double_addition(Am, Bm, Cm)
        # Test
        @test isapprox(status, 0)

        # Fixture
        A, Ac, Am, B, Bc, Bm, C, Cc, Cm = fixture(5,5)
        Cm.numRows = 0
        # Call
        status = test_call_matrix_double_addition(Am, Bm, Cm)
        # Test
        @test isapprox(status, 0)

        # Fixture
        A, Ac, Am, B, Bc, Bm, C, Cc, Cm = fixture(5,5)
        Cm.numCols = 0
        # Call
        status = test_call_matrix_double_addition(Am, Bm, Cm)
        # Test
        @test isapprox(status, 0)
    end
end
