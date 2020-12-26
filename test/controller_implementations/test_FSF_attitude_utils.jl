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
    return J, kR, kc, kw
end

@verbose(1, "Testing FSF utilities...")
@testset "Test FSF utilities" begin
    @testset "SU2 - Feasibility check - Non-symmetric inertia" begin
        # Fixture
        J, kX, kc, kw = fixture()
        J[1,2], J[2,1] = 2.0, 1.0;
        C      = con_state_qw_fsf_t(J, kX, kc, kw);
        C_copy = deepcopy(C)

        # Test
        status = ccall((:assert_attitude_FSF_SU2, AerialVehicleControl.CONT_LIB_PATH),
            Cint,
            (Ref{con_state_qw_fsf_t},),
            C)

        @test isapprox(status, 0);
        allFieldsEqual(C, C_copy, true);
    end

    @testset "SU2 - Feasibility check - Non-PSD inertia" begin
        # Fixture
        J, kX, kc, kw = fixture()
        J  = J - I*1.01*minimum(abs.(eigvals(J)));
        C      = con_state_qw_fsf_t(J, kX, kc, kw);
        C_copy = deepcopy(C)

        # Test
        status = ccall((:assert_attitude_FSF_SU2, AerialVehicleControl.CONT_LIB_PATH),
            Cint,
            (Ref{con_state_qw_fsf_t},),
            C)

        @test isapprox(status, 0);
        allFieldsEqual(C, C_copy, true);
    end
    @testset "SU2 - Feasibility check - Controller parameters - kc large" begin
        # Fixture
        J, kX, _, kw = fixture()
        kc     = 1.01 * maximal_feasible_kc_SU2(kX, kw, J);
        C      = con_state_qw_fsf_t(J, kX, kc, kw);
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
    @testset "SU2 - Feasibility check - Controller parameters - kc small" begin
        # Fixture
        J, kX, _, kw = fixture()
        kc     = 0.99 * maximal_feasible_kc_SU2(kX, kw, J);
        C      = con_state_qw_fsf_t(J, kX, kc, kw);
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
    @testset "SO3 - Feasibility check - Non-symmetric inertia" begin
        # Fixture
        J, kR, kc, kw = fixture()
        J[1,2], J[2,1] = 2.0, 1.0;
        C      = con_state_qw_fsf_t(J, kR, kc, kw);
        C_copy = deepcopy(C)

        # Test
        status = ccall((:assert_attitude_FSF_SO3, AerialVehicleControl.CONT_LIB_PATH),
            Cint,
            (Ref{con_state_qw_fsf_t},),
            C)

        @test isapprox(status, 0);
        allFieldsEqual(C, C_copy, true);
    end

    @testset "SO3 - Feasibility check - Non-PSD inertia" begin
        # Fixture
        J, kR, kc, kw = fixture()
        J  = J - I*1.01*minimum(abs.(eigvals(J)));
        C      = con_state_qw_fsf_t(J, kR, kc, kw);
        C_copy = deepcopy(C)

        # Test
        status = ccall((:assert_attitude_FSF_SO3, AerialVehicleControl.CONT_LIB_PATH),
            Cint,
            (Ref{con_state_qw_fsf_t},),
            C)

        @test isapprox(status, 0);
        allFieldsEqual(C, C_copy, true);
    end
    @testset "SO3 - Feasibility check - Controller parameters - kc large" begin
        # Fixture
        J, kR, _, kw = fixture()
        kc     = 1.01 * maximal_feasible_kc_SO3(kR, kw, J);
        C      = con_state_qw_fsf_t(J, kR, kc, kw);
        C_copy = deepcopy(C)

        # Test
        status = ccall((:assert_attitude_FSF_SO3, AerialVehicleControl.CONT_LIB_PATH),
            Cint,
            (Ref{con_state_qw_fsf_t},),
            C)

        # Assert
        @test isapprox(status, 0);
        allFieldsEqual(C, C_copy, true);
    end
    @testset "SO3 - Feasibility check - Controller parameters - kc small" begin
        # Fixture
        J, kR, _, kw = fixture()
        kc     = 0.99 * maximal_feasible_kc_SO3(kR, kw, J);
        C      = con_state_qw_fsf_t(J, kR, kc, kw);
        C_copy = deepcopy(C)

        # Test
        status = ccall((:assert_attitude_FSF_SO3, AerialVehicleControl.CONT_LIB_PATH),
            Cint,
            (Ref{con_state_qw_fsf_t},),
            C)

        # Assert
        @test isapprox(status, 1);
        allFieldsEqual(C, C_copy, true);
    end
end
