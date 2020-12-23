"""
Copyright Marcus Greiff 2020

Tests of the continuous attitude FOF controller on SO(3)
"""

# Include modules and recompile if running test standalone
if !isdefined(Main, :runAllTests)
    using Test, AerialVehicleControl, LinearAlgebra
    recompile();
end

# Testing fixture
function fixture()
    # Randomize a set of reference states
    qref = rand(4,1); qref = qref./norm(qref); wref = rand(3,1); aref = rand(3,1); fref = 0.5;
    R = ref_state_qw_t(hcat(qref), hcat(wref), hcat(aref), fref)
    # Randomize a set of states
    q = rand(4,1); q = q./norm(q); w = rand(3,1);
    S = dyn_state_qw_t(hcat(q), hcat(w))
    # Randomize the controller struct
    J  = rand_PSD(3);
    Kw = rand_PSD(3);
    ki = rand(1,3);
    Cw = rand_PSD(3);
    cw, cR = 1.0, 1.0;
    a, b, c, d = 1.0, 2.0, 3.0, 4.0;
    Y0 = rand(3,1);
    Yi = rand(3,3);
    Vi = rand(3,3);
    # Form controller object
    C  = con_state_qw_fof_t(J, Kw, ki, Cw, cw, cR, a, b, c, d, Y0, Yi, Vi);

    return R, S, C, qref, wref, aref, q, w, J, Kw, ki, Cw, cw, cR, Y0, Yi, Vi;
end


@verbose(1, "Testing FOF continuous on SO(3)...")
@testset "Test the continuous FOF controller on SO(3)" begin
    @testset "Test the crosproduct sum function" begin
        # Fixture
        N = 5;
        YA, YAc, YAm = get_matrix_and_copy(randn(3,N));
        YB, YBc, YBm = get_matrix_and_copy(randn(3,N));
        g,  gc,  gm  = get_matrix_and_copy(rand(1,N).+0.1);
        o,  oc,  om  = get_matrix_and_copy(rand(3,1));


        out_J = zeros(3,1);
        for i = 1:N
            out_J[:,1] = out_J[:,1] + g[1,i] * cross(YA[:,i], YB[:,i]);
        end

        # Test
        status = ccall((:attitude_FOF_SO3_continuous_cross_terms, AerialVehicleControl.CONT_LIB_PATH),
            Cint,
            (Ref{matrix_double_t}, Ref{matrix_double_t}, Ref{matrix_double_t}, Ref{matrix_double_t},),
            YAm,
            YBm,
            gm,
            om)

        @test isapprox(status, 1);

        # Assert that none of the input arguments have changed
        @test isapprox(YA, YAc);
        @test isapprox(YB, YBc);
        @test isapprox(g,  gc );
        # Assert that the output has changed
        @test !isapprox(o,  oc );
        # Assert that the output is computed correctly
        @test isapprox(o,  out_J);
    end
    @testset "Test computation of the control signals" begin
        # Fixture
        R, S, C, qr, wr, ar, q, w, J, Kw, ki, Cw, cw, cR, Y0, Yi, Vi = fixture();
        Rc = deepcopy(R);
        Sc = deepcopy(S);
        Cc = deepcopy(C);
        dt = Cdouble(0.1);

        # Compute the controls in Julia based on the current memory
        qhat  = zeros(4,1);
        qhat[1,1]  = S.quaternion[1];
        qhat[2,1]  = S.quaternion[2];
        qhat[3,1]  = S.quaternion[3];
        qhat[4,1]  = S.quaternion[4];
        what  = zeros(3,1);
        what[1,1]  = S.attrates[1];
        what[2,1]  = S.attrates[2];
        what[3,1]  = S.attrates[3];
        k1, k2, k3 = ki[1], ki[2], ki[3];

        taur       = J*ar - skew(J*wr)*wr;
        Rr         = quat_2_SO3(qr);
        Rhat       = quat_2_SO3(qhat);
        Rmat       = quat_2_SO3(q);
        v1, v2, v3 = Vi[:,1], Vi[:,2], Vi[:,3];

        whate  = wr - what;
        we     = wr - Y0;
        wtilde = what - Y0;

        vecA   = taur;
        vecB   = skew(J*whate)*wr;
        vecC   = Kw*whate;
        vecD   = (+k1*skew(Rr'*v1)*(Rhat'*v1)
                  +k2*skew(Rr'*v2)*(Rhat'*v2)
                  +k3*skew(Rr'*v3)*(Rhat'*v3));
        tau_J  = vecA + vecB + vecC + vecD;


        deltaR = -cR * (+k1 * skew(Rhat'*v1)*(Rr'*v1 + Yi[:,1])
                        +k2 * skew(Rhat'*v2)*(Rr'*v2 + Yi[:,2])
                        +k3 * skew(Rhat'*v3)*(Rr'*v3 + Yi[:,3]));

        deltaW = (-cw * J * skew(wr)*we
                  -cw * Kw * we
                  -Cw * wtilde);

        Rhatp1_J = Rmat * exp(skew( dt * (Y0 + deltaR)));
        whatp1_J = w + dt * inv(J) * (skew(J*Y0)*Y0 + tau_J + deltaW);

        # Call
        status = ccall((:update_attitude_FOF_SO3_continuous, AerialVehicleControl.CONT_LIB_PATH),
            Cint,
            (Ref{ref_state_qw_t}, Ref{dyn_state_qw_t}, Ref{con_state_qw_fof_t}, Cdouble,),
            R,
            S,
            C,
            dt)

        # test that the function is called correctly
        @test isapprox(status, 1);

        # Test that the memory that should be unchanged remains unchanged
        @test isapprox(C.status,  Cc.status );
        @test isapprox(C.gain_cw, Cc.gain_cw);
        @test isapprox(C.gain_cR, Cc.gain_cR);
        @test isapprox(C.param_a, Cc.param_a);
        @test isapprox(C.param_b, Cc.param_b);
        @test isapprox(C.param_c, Cc.param_c);
        @test isapprox(C.param_d, Cc.param_d);

        for i = 1:3
            # The torque should have changed in the controller, but
            # everything else should be the same
            @test !isapprox(C.torques[i], Cc.torques[i]);
            @test isapprox(C.gain_ki[i], Cc.gain_ki[i]);
            @test isapprox(C.measuredGyrorates[i], Cc.measuredGyrorates[i]);
            for j = 1:3
                @test isapprox(C.inertia[i][j], Cc.inertia[i][j]);
                @test isapprox(C.invinertia[i][j], Cc.invinertia[i][j]);
                @test isapprox(C.gain_Kw[i][j], Cc.gain_Kw[i][j]);
                @test isapprox(C.gain_Cw[i][j], Cc.gain_Cw[i][j]);
                @test isapprox(C.measuredDirections[i][j], Cc.measuredDirections[i][j]);
                @test isapprox(C.globalDirections[i][j], Cc.globalDirections[i][j]);
            end
        end

        # The state quaternion and attitude rates should have changed, but
        # everything in the reference struct should be the same
        for i = 1:4
            @test !isapprox(S.quaternion[i], Sc.quaternion[i]);
            @test  isapprox(R.quaternion[i], Rc.quaternion[i]);
        end
        for i = 1:3
            @test !isapprox(S.attrates[i], Sc.attrates[i]);
            @test  isapprox(R.attrates[i], Rc.attrates[i]);
            @test  isapprox(R.attaccel[i], Rc.attaccel[i]);
        end

        # Verify that the moments have been computed correctly
        for i = 1:3
            @test isapprox(tau_J[i,1], C.torques[i]);
        end

        # Verify that the estimator has performed a correct update
        qhatp1 = zeros(4,1);
        for i = 1:4
            qhatp1[i,1] = S.quaternion[i];
        end
        Rhatp1_C = quat_2_SO3(qhatp1);

        for i = 1:3
            for j = 1:3
                @test isapprox(Rhatp1_J[i,j], Rhatp1_C[i,j]);
            end
        end

        for i = 1:3
            @test isapprox(whatp1_J[i], S.attrates[i]);
        end

    end
end
