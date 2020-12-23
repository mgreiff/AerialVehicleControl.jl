"""
Copyright Marcus Greiff 2020

Tests of the attitude conversions general math functions
"""

# Include modules and recompile if running test standalone
if !isdefined(Main, :runAllTests)
    using Test, AerialVehicleControl, LinearAlgebra
    recompile();
end


# Testing fixture
function fixture()
    # Generate a random normalized Quaternion
    q = rand(4,1); q = q./norm(q);
    # Generate a random rotation matrix
    R = exp(skew(rand(3,1)));
    return q, R
end

@verbose(1, "Testing math functions...")
@testset "General math functions" begin
    @testset "Test the dot product function" begin
        @testset "Good inputs" begin
            for i = 2:5
                # Fixture
                v1, v1c, v1m = get_matrix_and_copy(rand(i,1));
                v2, v2c, v2m = get_matrix_and_copy(rand(i,1));
                prod = v1'*v2;

                # Test
                output = ccall((:cont_dot_product, AerialVehicleControl.CONT_LIB_PATH),
                    Cdouble,
                    (Ref{matrix_double_t}, Ref{matrix_double_t}),
                    v1m,
                    v2m)

                # Assert
                @test isapprox(output, prod[1,1]);
                @test isapprox(v1, v1c);
                @test isapprox(v2, v2c);
            end
        end
        # TODO write bad input tests
    end
    @testset "Test the cross-product function" begin
        for i = 1:5
            # Fixture
            yA, yAc, yAm = get_matrix_and_copy(rand(3,1));
            yB, yBc, yBm = get_matrix_and_copy(rand(3,1));
            o,  oc,  om  = get_matrix_and_copy(rand(3,1));

            out_J      = zeros(3,1);
            out_J[:,1] = cross(yA[:,1], yB[:,1]);

            # Test
            ccall((:cont_cross_product, AerialVehicleControl.CONT_LIB_PATH),
                Cvoid,
                (Ref{matrix_double_t}, Ref{matrix_double_t}, Ref{matrix_double_t},),
                yAm,
                yBm,
                om)

            # Test that only the output has changed
            @test isapprox(yA, yAc);
            @test isapprox(yB, yBc);
            @test !isapprox(o,  oc);

            # Test that the correct values are computed
            @test isapprox(o,  out_J);

        end
    end
    @testset "Test the sinc function" begin
        @testset "Good inputs" begin
            for i = 1:10
                # Fixture
                input = 10.0* rand() - 5.0;

                # Test
                output = ccall((:cont_sinc, AerialVehicleControl.CONT_LIB_PATH),
                    Cdouble,
                    (Ref{Cdouble},),
                    input)

                # Assert
                @test isapprox(sin(input)/input, output);
            end
        end
    end
    @testset "Test the sign function" begin
        @testset "Good inputs" begin
            for i = 1:10
                # Fixture
                input = 10.0* rand() - 5.0;
                if i == 1;
                    input = 0.0;
                end
                sign_in = (input<0.0) ? -1.0 : 1.0;

                # Test
                output = ccall((:cont_sign_func, AerialVehicleControl.CONT_LIB_PATH),
                    Cdouble,
                    (Ref{Cdouble},),
                    input)

                # Assert
                @test isapprox(sign_in, output);
            end
        end
    end
    @testset "Test the conversion or quaternion to R in SO(3)" begin
        @testset "Check that R is in SO(3)" begin
            for i = 1:10
                # Fixture
                q, _ = fixture();
                q, qc, qm = get_matrix_and_copy(q);
                R, Rc, Rm = get_matrix_and_copy(rand(3,3));

                # Test
                status = ccall((:cont_quat_2_SO3, AerialVehicleControl.CONT_LIB_PATH),
                    Cint,
                    (Ref{matrix_double_t}, Ref{matrix_double_t}),
                    qm,
                    Rm)

                # Assert
                @test  isapprox(status, 1);
                @test  isapprox(q, qc);
                @test !isapprox(R, Rc);
                @test  isapprox(R'*R, Matrix{Float64}(I, 3, 3));
                @test  isapprox(det(R), 1.0);
            end
        end
        @testset "Check that the rotation is corresponds to the quaternion definition" begin
            for i = 1:10
                # Fixture
                theta = 2*pi*rand();
                u      = rand(3,1);
                u      = u ./ norm(u);
                q      = rand(4,1);
                q[1]   = cos(theta/2.0);
                q[2:4] = sin(theta/2.0)*u;

                Rtrue = exp(skew(u*theta));

                q, qc, qm = get_matrix_and_copy(q);
                R, Rc, Rm = get_matrix_and_copy(rand(3,3));

                # Test
                status = ccall((:cont_quat_2_SO3, AerialVehicleControl.CONT_LIB_PATH),
                    Cint,
                    (Ref{matrix_double_t}, Ref{matrix_double_t}),
                    qm,
                    Rm)

                # Assert
                @test  isapprox(status, 1);
                @test  isapprox(q, qc);
                @test !isapprox(R, Rc);
                @test  isapprox(R, Rtrue);
            end
        end
    end
    @testset "Test the vector normalization function" begin
        @testset "Column vector" begin
            for N = 1:5
                v, vc, vm = get_matrix_and_copy(rand(N,1));

                v_normalized = v ./ norm(v);

                status = ccall((:cont_normalize, AerialVehicleControl.CONT_LIB_PATH),
                    Cint,
                    (Ref{matrix_double_t},),
                    vm)

                @test isapprox(status, 1)
                for i = 1:N
                    @test  isapprox(v_normalized[i,1], v[i,1])
                    @test !isapprox(vc[i,1],           v[i,1])
                end
            end
        end
        @testset "Row vector" begin
            for N = 1:5
                v, vc, vm = get_matrix_and_copy(rand(1,N));

                v_normalized = v ./ norm(v);

                status = ccall((:cont_normalize, AerialVehicleControl.CONT_LIB_PATH),
                    Cint,
                    (Ref{matrix_double_t},),
                    vm)

                @test isapprox(status, 1)
                for i = 1:N
                    @test  isapprox(v_normalized[1,i], v[1,i])
                    @test !isapprox(vc[1,i],           v[1,i])
                end
            end
        end
        @testset "Bad inputs" begin
            v, vc, vm = get_matrix_and_copy(rand(2,2));

            v_normalized = v ./ norm(v);

            status = ccall((:cont_normalize, AerialVehicleControl.CONT_LIB_PATH),
                Cint,
                (Ref{matrix_double_t},),
                vm)

            @test isapprox(status, 0)
        end
    end
end
