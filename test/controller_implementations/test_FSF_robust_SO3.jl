"""
Copyright Marcus Greiff 2020

Tests of the robust attitude FSF controller on SO(3)
"""

# Include modules and recompile if running test standalone
if !isdefined(Main, :runAllTests)
    using Test, AerialVehicleControl, LinearAlgebra
    recompile();
end

# Testing fixture
function fixture()

    qref = rand(4,1); qref = qref./norm(qref); wref = rand(3,1); aref = rand(3,1); fref = 0.5;
    R = ref_state_qw_t(hcat(qref), hcat(wref), hcat(aref), fref)

    q = rand(4,1); q = q./norm(q); w = rand(3,1);
    S = dyn_state_qw_t(hcat(q), hcat(w))

    J   = rand(3,3); J  = J + J'; J  = J + I*maximum(abs.(eigvals(J)));
    kR  = 1.0;
    kw  = 1.0;
    kc  = 0.5 * maximal_feasible_kc_SO3(kR, kw, J);
    eps = 0.001;
    L   = 1.0;
    C   = con_state_qw_fsf_t(J, kR, kc, kw; eps=eps, L=L);

    return R, S, C, qref, wref, aref, q, w, kR, kc, kw, eps, L, J
end

@verbose(1, "Testing robust FSF on SO(3)...")
@testset "Test the robust FSF controller on SO(3)" begin
    @testset "Check that R is in SO(3)" begin
        # Fixture
        Rstruct, Sstruct, Cstruct, qref, wref, aref, q, w, kR, kc, kw, eps, L, J = fixture()
        Rstruct_copy = deepcopy(Rstruct)
        Sstruct_copy = deepcopy(Sstruct)
        Cstruct_copy = deepcopy(Cstruct)

        # Compute what the controls/distances should be based on the current memory
        R   = quat_2_SO3(q );
        Rr  = quat_2_SO3(qref);
        eR  = SO3_vee(Rr'*R - R'*Rr) / 2.0;
        ew  = w -  R'*Rr * wref;
        eA  = ew + kc * (J \ eR);
        mu  = -(L^2 / (L * norm(eA) + eps)) * eA;
        ff  = SO3_hat(w) * (J * w) - J*(SO3_hat(w)*R'*Rr*wref -  R'*Rr*aref);

        tau = -kR * eR - kw * ew + ff + mu;
        Psi   = 0.5 * tr(Matrix{Float64}(I, 3, 3) - Rr'*R);
        Gamma = 0.5 * tr(Matrix{Float64}(I, 2, 2) - quat_2_SU2(qref)'*quat_2_SU2(q));
        attpart  = kR * Psi;
        crosspart= kc * ew'*eR;
        ratepart = ew'*J*ew/2.0;
        Lyap     = attpart + crosspart[1,1] + ratepart[1,1];

        # Test
        status = ccall((:update_attitude_FSF_SO3_robust, AerialVehicleControl.CONT_LIB_PATH),
            Cint,
            (Ref{ref_state_qw_t}, Ref{dyn_state_qw_t}, Ref{con_state_qw_fsf_t},),
            Rstruct,
            Sstruct,
            Cstruct)

        # Test that the computations passed
        @test isapprox(status, 1);

        # Check that all of the values in the controller memory that should remain constant are unchanged
        allFieldsEqual(Sstruct, Sstruct_copy);
        allFieldsEqual(Rstruct, Rstruct_copy);
        allFieldsEqual(Cstruct, Cstruct_copy, false);

        # Check that the distances have changed
        @test !isapprox(Cstruct.dist_Psi,       Cstruct_copy.dist_Psi );
        @test !isapprox(Cstruct.dist_Gamma,     Cstruct_copy.dist_Gamma );
        @test !isapprox(Cstruct.dist_lyapunov,  Cstruct_copy.dist_lyapunov );

        # Check that the values of the controller structure that should be
        # updated is updated, and that these values corresponds to the above
        # computed torques and distances
        for i = 1:3
            @test isapprox(Cstruct.torques[i], tau[i,1]);
        end
        @test isapprox(Cstruct.dist_Psi,      Psi );
        @test isapprox(Cstruct.dist_Gamma,    Gamma );
        @test isapprox(Cstruct.dist_lyapunov, Lyap  );

    end
end
