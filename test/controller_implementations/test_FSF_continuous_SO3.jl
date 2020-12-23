"""
Copyright Marcus Greiff 2020

Tests of the continuous attitude FSF controller on SO(3)
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

    J  = rand(3,3); J  = J + J'; J  = J + I*maximum(abs.(eigvals(J)));
    kR, kc, kw = 1.0, 1.0, 1.0;
    a, b, c, d = 1.0, 2.0, 3.0, 4.0;
    C  = con_state_qw_fsf_t(J, kR, kc, kw, a, b, c, d);

    return R, S, C, qref, wref, aref, q, w, kR, kc, kw, J
end

@verbose(1, "Testing FSF continuous on SO(3)...")
@testset "Test the continuous FSF controller on SO(3)" begin
    @testset "Check that R is in SO(3)" begin
        # Fixture
        Rstruct, Sstruct, Cstruct, qref, wref, aref, q, w, kR, kc, kw, J = fixture()
        Rstruct_copy = deepcopy(Rstruct)
        Sstruct_copy = deepcopy(Sstruct)
        Cstruct_copy = deepcopy(Cstruct)

        # Compute what the controls/distances should be based on the current memory
        R   = quat_2_SO3(q );
        Rr  = quat_2_SO3(qref);
        eR  = SO3_vee(Rr'*R - R'*Rr) / 2.0;
        ew  = w -  R'*Rr * wref;
        ff  = SO3_hat(w) * (J * w) - J*(SO3_hat(w)*R'*Rr*wref -  R'*Rr*aref);

        tau = -kR * eR - kw * ew + ff;
        Psi   = 0.5 * tr(Matrix{Float64}(I, 3, 3) - Rr'*R);
        Gamma = 0.5 * tr(Matrix{Float64}(I, 2, 2) - quat_2_SU2(qref)'*quat_2_SU2(q));
        attpart  = kR * Psi;
        crosspart= kc * ew'*eR;
        ratepart = ew'*J*ew/2.0;
        Lyap     = attpart + crosspart[1,1] + ratepart[1,1];

        # Test
        status = ccall((:update_attitude_FSF_SO3_continuous, AerialVehicleControl.CONT_LIB_PATH),
            Cint,
            (Ref{ref_state_qw_t}, Ref{dyn_state_qw_t}, Ref{con_state_qw_fsf_t},),
            Rstruct,
            Sstruct,
            Cstruct)

        @test  isapprox(status, 1);

        # Check that the values of the structs R and S that should remain constant are unchanged
        for i = 1:4
            @test  isapprox(Rstruct.quaternion[i], Rstruct_copy.quaternion[i]);
            @test  isapprox(Sstruct.quaternion[i], Sstruct_copy.quaternion[i]);
        end
        for i = 1:3
            @test  isapprox(Rstruct.attrates[i], Rstruct_copy.attrates[i]);
            @test  isapprox(Rstruct.attaccel[i], Rstruct_copy.attaccel[i]);
            @test  isapprox(Sstruct.attrates[i], Sstruct_copy.attrates[i]);
        end

        # Check that all of the values in the controller memory that should remain constant are unchanged
        @test isapprox(Cstruct.gain_kR, Cstruct_copy.gain_kR);
        @test isapprox(Cstruct.gain_kc, Cstruct_copy.gain_kc);
        @test isapprox(Cstruct.gain_kw, Cstruct_copy.gain_kw);
        @test isapprox(Cstruct.param_a, Cstruct_copy.param_a);
        @test isapprox(Cstruct.param_b, Cstruct_copy.param_b);
        @test isapprox(Cstruct.param_c, Cstruct_copy.param_c);
        @test isapprox(Cstruct.param_d, Cstruct_copy.param_d);
        @test isapprox(Cstruct.thrust, Cstruct_copy.thrust);
        for i = 1:3
            for j = 1:3
                @test isapprox(Cstruct.inertia[i][j], Cstruct_copy.inertia[i][j]);
            end
        end

        # Check that the control signal torques have changed
        for i = 1:3
            @test !isapprox(Cstruct.torques[i],  Cstruct_copy.torques[i] );
        end

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
