"""
Copyright Marcus Greiff 2020

Tests of the attitude conversions and distances pertaining to SU(2)
"""

# Include modules and recompile if running test standalone
if !isdefined(Main, :runAllTests)
    using Test, AerialVehicleControl, LinearAlgebra
    recompile();
end

@verbose(1, "Testing the SU(2) maps...")
@testset "Attitude maps on SU(2)" begin
    @testset "Test the conjugate operation" begin
        @testset "Good inputs" begin
            # Fixture
            q, qc, qm = get_matrix_and_copy(rand(4,1))
            p, pc, pm = get_matrix_and_copy(rand(4,1))
            # Test
            status = ccall((:cont_SU2_conjugate, AerialVehicleControl.CONT_LIB_PATH),
                Cint,
                (Ref{matrix_double_t}, Ref{matrix_double_t}),
                qm,
                pm)
            # Assert
            @test isapprox(status,1);
            @test isapprox(q, qc);
            @test !isapprox(p, pc);
            @test isapprox(q[1], p[1]);
            @test isapprox(q[2:4], -p[2:4]);
        end
        @testset "Bad inputs" begin
            # Fixture
            q, qc, qm = get_matrix_and_copy(rand(3,1))
            p, pc, pm = get_matrix_and_copy(rand(4,1))
            # Test
            status = ccall((:cont_SU2_conjugate, AerialVehicleControl.CONT_LIB_PATH),
                Cint,
                (Ref{matrix_double_t}, Ref{matrix_double_t}),
                qm,
                pm)
            # Assert
            @test isapprox(status, 0);
        end
    end
    @testset "Test the product operation" begin
        @testset "Good inputs" begin
            # Fixture
            q, qc, qm = get_matrix_and_copy(rand(4,1))
            p, pc, pm = get_matrix_and_copy(rand(4,1))
            o, oc, om = get_matrix_and_copy(rand(4,1))
            qpprod    = quaternion_left_product(q)*p;

            # Test
            status = ccall((:cont_SU2_product, AerialVehicleControl.CONT_LIB_PATH),
                Cint,
                (Ref{matrix_double_t}, Ref{matrix_double_t}, Ref{matrix_double_t}),
                qm,
                pm,
                om)

            # Assert
            @test  isapprox(status,1);
            @test  isapprox(q, qc);
            @test  isapprox(p, pc);
            @test !isapprox(o, oc);
            @test  isapprox(o, qpprod);
        end
        @testset "Bad inputs" begin
            # Fixture
            q, qc, qm = get_matrix_and_copy(rand(4,1))
            p, pc, pm = get_matrix_and_copy(rand(4,1))
            o, oc, om = get_matrix_and_copy(rand(3,1))
            qpprod    = quaternion_left_product(q)*p;

            # Test
            status = ccall((:cont_SU2_product, AerialVehicleControl.CONT_LIB_PATH),
                Cint,
                (Ref{matrix_double_t}, Ref{matrix_double_t}, Ref{matrix_double_t}),
                qm,
                pm,
                om)
            # Assert
            @test isapprox(status, 0);
        end
    end
    @testset "Test the vee operation" begin
        @testset "Good inputs" begin
            # Fixture
            q, qc, qm = get_matrix_and_copy(rand(4,1))
            o, oc, om = get_matrix_and_copy(rand(3,1))
            X         = quat_2_SU2(q);
            output    = SU2_vee(X);

            # Test
            status = ccall((:cont_SU2_vee, AerialVehicleControl.CONT_LIB_PATH),
                Cint,
                (Ref{matrix_double_t}, Ref{matrix_double_t}),
                qm,
                om)

            # Assert
            @test  isapprox(status,1);
            @test  isapprox(q, qc);
            @test !isapprox(o, oc);
            @test  isapprox(o, output);
        end
        # Todo write tests for bad inputs
    end
    @testset "Test the hat operation" begin
        @testset "Good inputs" begin
            # Fixture
            u, uc, um = get_matrix_and_copy(rand(3,1))
            o, oc, om = get_matrix_and_copy(rand(4,1))
            output    = SU2_hat(u);

            # Test
            status = ccall((:cont_SU2_hat, AerialVehicleControl.CONT_LIB_PATH),
                Cint,
                (Ref{matrix_double_t}, Ref{matrix_double_t}),
                um,
                om)

            # Assert
            @test  isapprox(status,1);
            @test  isapprox(u, uc);
            @test !isapprox(o, oc);
            @test  isapprox(quat_2_SU2(o), output);
        end
        # Todo write tests for bad inputs
    end
    @testset "Test the Exponential map on SU(2)" begin
        @testset "Good inputs" begin
            # Fixture
            u, uc, um = get_matrix_and_copy(rand(3,1))
            q, qc, qm = get_matrix_and_copy(rand(4,1))

            R_J = exp(SO3_hat(u))

            # Test
            status = ccall((:cont_SU2_Exp, AerialVehicleControl.CONT_LIB_PATH),
                Cint,
                (Ref{matrix_double_t}, Ref{matrix_double_t}),
                um,
                qm)

            R_C = quat_2_SO3(q);

            # Assert
            @test  isapprox(status, 1);
            @test  isapprox(u, uc);
            @test !isapprox(q, qc);
            @test  isapprox(R_J, R_C);
        end
        # Todo write tests for bad inputs
    end
    @testset "Test the triple product function on SU(2)" begin
        @testset "Good inputs" begin
            # Fixture
            a = 2*pi*rand()-pi;
            b = 2*pi*rand()-pi;
            c = 2*pi*rand()-pi;

            q, qc, qm = get_matrix_and_copy(rand(4,1))

            R_A = exp(SO3_hat([a;0;0]))
            R_B = exp(SO3_hat([0;b;0]))
            R_C = exp(SO3_hat([0;0;c]))
            R_J = R_A * R_B * R_C;

            # Test
            status = ccall((:cont_SU2_triple, AerialVehicleControl.CONT_LIB_PATH),
                Cint,
                (Cdouble, Cdouble, Cdouble, Ref{matrix_double_t}),
                a,
                b,
                c,
                qm)

            R_C = quat_2_SO3(q);

            # Assert
            @test  isapprox(status, 1);
            @test !isapprox(q, qc);
            @test  isapprox(R_J, R_C);
        end
        # Todo write tests for bad inputs
    end
    @testset "Test the distance on SU(2)" begin
        @testset "Good inputs" begin
            # Fixture
            q = rand(4,1); q = q/norm(q);
            p = rand(4,1); p = p/norm(p);
            Xq = quat_2_SU2(q);
            Xp = quat_2_SU2(p);
            Gamma = real(tr(I - Xq'*Xp) / 2.0)

            q, qc, qm = get_matrix_and_copy(q)
            p, pc, pm = get_matrix_and_copy(p)

            # Test
            distance = ccall((:cont_SU2_distance, AerialVehicleControl.CONT_LIB_PATH),
                Cdouble,
                (Ref{matrix_double_t}, Ref{matrix_double_t}),
                qm,
                pm)

            # Assert
            @test isapprox(q, qc);
            @test isapprox(p, pc);
            @test isapprox(distance, Gamma);
        end
        # Todo write tests for bad inputs
    end
end
