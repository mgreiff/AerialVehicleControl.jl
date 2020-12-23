"""
Copyright Marcus Greiff 2020

Tests of the attitude conversions and distances pertaining to SO(3)
"""

# Include modules and recompile if running test standalone
if !isdefined(Main, :runAllTests)
    using Test, AerialVehicleControl, LinearAlgebra
    recompile();
end

@verbose(1, "Testing the SO(3) maps...")
@testset "Attitude maps" begin
    @testset "Test the SO(3) hat map" begin
        @testset "Good inputs" begin
            # Fixture
            A, Ac, Am = get_matrix_and_copy(rand(3,3))
            u, uc, um = get_matrix_and_copy(rand(3,1))
            # Test
            status = ccall((:cont_SO3_hat, AerialVehicleControl.CONT_LIB_PATH),
                Cint,
                (Ref{matrix_double_t}, Ref{matrix_double_t}),
                um,
                Am)
            # Assert
            @test isapprox(status,1);
            @test isapprox(u, uc);
            @test !isapprox(A, Ac);
            @test isapprox(A, skew(u));
        end
        @testset "Bad inputs" begin
            # Fixture
            A, Ac, Am = get_matrix_and_copy(rand(3,4))
            u, uc, um = get_matrix_and_copy(rand(3,1))
            # Test
            status = ccall((:cont_SO3_hat, AerialVehicleControl.CONT_LIB_PATH),
                Cint,
                (Ref{matrix_double_t}, Ref{matrix_double_t}),
                um,
                Am)
            # Assert
            @test isapprox(status, 0);
        end
    end
    @testset "Test the SO(3) vee map" begin
        @testset "Good inputs" begin
            # Fixture
            u               = rand(3,1);
            S, Sc, Sm       = get_matrix_and_copy(skew(u))
            out, outc, outm = get_matrix_and_copy(rand(3,1))
            # Test
            status = ccall((:cont_SO3_vee, AerialVehicleControl.CONT_LIB_PATH),
                Cint,
                (Ref{matrix_double_t}, Ref{matrix_double_t}),
                Sm,
                outm)
            # Assert
            @test isapprox(status,1);
            @test isapprox(S, Sc);
            @test !isapprox(out, outc);
            @test isapprox(u, out);
        end
        @testset "Bad inputs" begin
            # Fixture
            u               = rand(3,1);
            S, Sc, Sm       = get_matrix_and_copy(skew(u))
            out, outc, outm = get_matrix_and_copy(rand(3,2))
            # Test
            status = ccall((:cont_SO3_vee, AerialVehicleControl.CONT_LIB_PATH),
                Cint,
                (Ref{matrix_double_t}, Ref{matrix_double_t}),
                Sm,
                outm)
            # Assert
            @test isapprox(status, 0);
        end
    end
    @testset "Test the SO(3) Psi distance" begin
        @testset "Good inputs" begin
            for i = 1:10
                # Fixture
                R1, R1c, R1m  = get_matrix_and_copy(exp(skew(rand(3,1))))
                R2, R2c, R2m  = get_matrix_and_copy(exp(skew(rand(3,1))))
                # Test
                value = ccall((:cont_SO3_distance, AerialVehicleControl.CONT_LIB_PATH),
                    Cdouble,
                    (Ref{matrix_double_t}, Ref{matrix_double_t}),
                    R1m,
                    R2m)
                # Assert
                @test isapprox(value, tr(I-R1'R2)/2.0);
                @test isapprox(R1, R1c)
                @test isapprox(R2, R2c)
            end
        end
        # TODO write bad input tests
    end
    @testset "Test the SO(3) exponential map" begin
        @testset "Good inputs" begin
            for i = 1:5
                # Fixture
                u, uc, um = get_matrix_and_copy(rand(3,1));
                R, Rc, Rm = get_matrix_and_copy(rand(3,3));
                # Test
                status = ccall((:cont_SO3_Exp, AerialVehicleControl.CONT_LIB_PATH),
                    Cint,
                    (Ref{matrix_double_t}, Ref{matrix_double_t}),
                    um,
                    Rm)
                # Assert
                @test isapprox(status,1);
                @test isapprox(u, uc);
                @test !isapprox(R, Rc);
                @test isapprox(R, exp(skew(uc)));
            end
        end
        # TODO write bad input tests
    end
    @testset "Test the SO(3) logarithmic map" begin
        @testset "Good inputs" begin
            for i = 1:5
                # Fixture
                u, uc, um = get_matrix_and_copy(rand(3,1));
                R, Rc, Rm = get_matrix_and_copy(exp(skew(u)));
                out, outc, outm = get_matrix_and_copy(rand(3,1));
                # Test
                status = ccall((:cont_SO3_Log, AerialVehicleControl.CONT_LIB_PATH),
                    Cint,
                    (Ref{matrix_double_t}, Ref{matrix_double_t}),
                    Rm,
                    outm)

                theta = acos((tr(R) - 1)/2)
                # Assert
                @test isapprox(status,1);
                # We should have changed out, but not R or u
                @test isapprox(R, Rc);
                @test isapprox(u, uc);
                @test !isapprox(out, outc);
                # We chould be able to recover u from R
                @test isapprox(u, out);
            end
        end
        # TODO write bad input tests
    end
end
