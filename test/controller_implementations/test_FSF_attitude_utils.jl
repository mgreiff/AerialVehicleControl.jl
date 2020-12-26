"""
Copyright Marcus Greiff 2020

Tests of the discontinuous attitude FSF controller on SU(2)
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
    a, b, c, d = 1.0, 2.0, 3.0, 4.0;
    return J, kR, kc, kw, a, b, c, d
end

@verbose(1, "Testing FSF utilities on SU(2)...")
@testset "Test the continuous FSF controller on SU(2)" begin
    @testset "Feasibility check - Non-symmetric inertia" begin
        # Fixture
        J, kX, kc, kw, a, b, c, d = fixture()
        J[1,2], J[2,1] = 2.0, 1.0;
        C      = con_state_qw_fsf_t(J, kX, kc, kw, a, b, c, d);
        C_copy = deepcopy(C)

        # Test
        status = ccall((:assert_attitude_FSF_SU2, AerialVehicleControl.CONT_LIB_PATH),
            Cint,
            (Ref{con_state_qw_fsf_t},),
            C)

        @test isapprox(status, 0);
        allFieldsEqual(C, C_copy, true);
    end

    @testset "Feasibility check - Non-PSD inertia" begin
        # Fixture
        J, kX, kc, kw, a, b, c, d = fixture()
        J  = J - I*1.01*minimum(abs.(eigvals(J)));
        C      = con_state_qw_fsf_t(J, kX, kc, kw, a, b, c, d);
        C_copy = deepcopy(C)

        # Test
        status = ccall((:assert_attitude_FSF_SU2, AerialVehicleControl.CONT_LIB_PATH),
            Cint,
            (Ref{con_state_qw_fsf_t},),
            C)

        @test isapprox(status, 0);
        allFieldsEqual(C, C_copy, true);
    end
    @testset "Feasibility check - Controller parameters - kc large" begin
        # Fixture
        J, kX, _, kw, a, b, c, d = fixture()
        kc     = 1.01 * maximal_feasible_kc_SU2(kX, kw, J);
        C      = con_state_qw_fsf_t(J, kX, kc, kw, a, b, c, d);
        C_copy = deepcopy(C)

        # Test
        status = ccall((:assert_attitude_FSF_SU2, AerialVehicleControl.CONT_LIB_PATH),
            Cint,
            (Ref{con_state_qw_fsf_t},),
            C)

        # Assert
        @test isapprox(status, 0);
        allFieldsEqual(C, C_copy, true);
    end
    @testset "Feasibility check - Controller parameters - kc small" begin
        # Fixture
        J, kX, _, kw, a, b, c, d = fixture()
        kc     = 0.99 * maximal_feasible_kc_SU2(kX, kw, J);
        C      = con_state_qw_fsf_t(J, kX, kc, kw, a, b, c, d);
        C_copy = deepcopy(C)

        # Test
        status = ccall((:assert_attitude_FSF_SU2, AerialVehicleControl.CONT_LIB_PATH),
            Cint,
            (Ref{con_state_qw_fsf_t},),
            C)

        # Assert
        @test isapprox(status, 1);
        allFieldsEqual(C, C_copy, true);
    end
end
