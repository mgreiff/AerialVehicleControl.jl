"""
Copyright Marcus Greiff 2020

Tests for the matrix_double_multiplication() function
"""

# Include modules and recompile if running test standalone
if !isdefined(Main, :runAllTests)
    using Test, AerialVehicleControl, LinearAlgebra
    recompile();
end

# Testing fixture
function fixture(mA,nA,mB,nB,mC,nC)
    A  = rand(Cdouble,mA,nA)
    Ac = deepcopy(A)
    Am = matrix_double_t(A)
    B  = rand(Cdouble,mB,nB)
    Bc = deepcopy(B)
    Bm = matrix_double_t(B)
    C  = rand(Cdouble,mC,nC)
    Cc = deepcopy(C)
    Cm = matrix_double_t(C)
    return A, Ac, Am, B, Bc, Bm, C, Cc, Cm
end

# Compact call for the matrix_double_addition function with NN
#function test_call_matrix_double_multiplication(Am, Bm, Cm)
#    status = ccall((:matrix_double_multiplication, AerialVehicleControl.CONT_LIB_PATH),
#        Cint,
#        (Ref{matrix_double_t}, Ref{matrix_double_t}, Ref{matrix_double_t}),
#        Am,
#        Bm,
#        Cm)
#    return status
#end

# Compact call for the matrix_double_addition function with NT
function test_call_matrix_double_multiplication(Am, Bm, Cm)
    status = ccall((:matrix_double_multiplication, AerialVehicleControl.CONT_LIB_PATH),
        Cint,
        (Ref{matrix_double_t}, Ref{matrix_double_t}, Ref{matrix_double_t}),
        Am,
        Bm,
        Cm)
    return status
end

@verbose(1, "Matrix multiplication...")
@testset "Matrix multiplication" begin
    @testset "Good inputs NN" begin
        for m = 1:10
            for n = 1:10
                # Fixture
                A, Ac, Am, B, Bc, Bm, C, Cc, Cm = fixture(m,n,n,m,m,m)
                # Call
                status = test_call_matrix_double_multiplication(Am, Bm, Cm)
                # Test
                @test isapprox(status, 1)
                @test isapprox(A, Ac)
                @test isapprox(B, Bc)
                @test !isapprox(C, Cc)
                @test isapprox(C, A * B)
            end
        end
    end
end
