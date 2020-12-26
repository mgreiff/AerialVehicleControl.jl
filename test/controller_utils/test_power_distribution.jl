"""
Copyright Marcus Greiff 2020

Tests for the power distribution functionality
"""

# Include modules and recompile if running test standalone
if !isdefined(Main, :runAllTests)
    using Test, AerialVehicleControl, LinearAlgebra
    recompile();
end

# Testing fixture
function fixture()
    J  = rand(3,3); J  = J + J'; J  = J + I*maximum(abs.(eigvals(J)));
    kR, kc, kw = 1.0, 1.0, 1.0;
    a, b, c, d = 0.35, 0.26, 1.0+rand(), 1.0+rand();
    C  = con_state_qw_fsf_t(J, kR, kc, kw; a=a, b=b, c=c, d=d);
    C.thrust = 0.5;
    t = rand(3,1);
    C.torques = ntuple(i->Cdouble(t[i,1]), 3)
    Cc = deepcopy(C);
    d, dc, dm = get_matrix_and_copy(rand(4,1))
    return C, Cc, t, d, dc, dm;
end

@verbose(1, "Testing the power distribution...")
@testset "Test the power distribution" begin
    @testset "Good inputs" begin
        # Fixture
        C, Cc, t, d, dc, dm = fixture();

        apar = C.param_a;
        bpar = C.param_b;
        cpar = C.param_c;
        dpar = C.param_d;
        T = C.thrust;

        A = [+1    +1    +1    +1   ;
             +0    -dpar +0    +dpar;
             +dpar +0    -dpar +0   ;
             -cpar +cpar -cpar +cpar]

        f_julia = A \ ([T; t]);

        f_julia[f_julia .<0.0] .= 0.0

        d_julia = -(apar/(2.0*bpar)) .+ sqrt.((apar/(2*bpar)).^2 .+ f_julia./bpar);

        # Test
        status = ccall((:compute_power_distribution, AerialVehicleControl.CONT_LIB_PATH),
            Cint,
            (Ref{con_state_qw_fsf_t}, Ref{matrix_double_t}),
            C,
            dm)

        # Assert
        @test isapprox(status, 1);
        for i = 1:4
            @test isapprox(d_julia[i,1], d[i,1]);
        end
    end
    @testset "Bad inputs (negative thrust)" begin
        # Fixture
        C, Cc, t, d, dc, dm = fixture();
        C.thrust = -C.thrust;
        Cc = deepcopy(C);

        # Test
        status = ccall((:compute_power_distribution, AerialVehicleControl.CONT_LIB_PATH),
            Cint,
            (Ref{con_state_qw_fsf_t}, Ref{matrix_double_t}),
            C,
            dm)

        # Assert
        @test isapprox(status, 0);
    end
    @testset "Bad inputs (negative parameter value)" begin
        for i = 1:4
            # Fixture
            C, Cc, t, d, dc, dm = fixture();
            C.param_a = -C.param_a;
            Cc = deepcopy(C);

            # Test
            status = ccall((:compute_power_distribution, AerialVehicleControl.CONT_LIB_PATH),
                Cint,
                (Ref{con_state_qw_fsf_t}, Ref{matrix_double_t}),
                C,
                dm)

            # Assert
            @test isapprox(status, 0);
        end
    end
    @testset "Bad inputs (wrong duty cycle row dimensions)" begin
        # Fixture
        C, Cc, t, _, _, _ = fixture();
        d, dc, dm = get_matrix_and_copy(rand(3,1));

        # Test
        status = ccall((:compute_power_distribution, AerialVehicleControl.CONT_LIB_PATH),
            Cint,
            (Ref{con_state_qw_fsf_t}, Ref{matrix_double_t}),
            C,
            dm)
        # Assert
        @test isapprox(status, 0);
    end
    @testset "Bad inputs (wrong duty cycle column dimensions)" begin
        # Fixture
        C, Cc, t, _, _, _ = fixture();
        d, dc, dm = get_matrix_and_copy(rand(4,2));

        # Test
        status = ccall((:compute_power_distribution, AerialVehicleControl.CONT_LIB_PATH),
            Cint,
            (Ref{con_state_qw_fsf_t}, Ref{matrix_double_t}),
            C,
            dm)
        # Assert
        @test isapprox(status, 0);
    end
end
