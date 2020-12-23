using LinearAlgebra

useVerbose  = [true,true,false]

macro verbose(level::Integer, string::String)
    @assert(0 < level < 4)
    if level > 0; padding = ""  ; end
    if level > 1; padding = "--"; end
    if level > 2; padding = "--"; end
    string = padding * string
    if useVerbose[level]
        return :( println($string) )
    end
end

function get_matrix_and_copy(in)
    return in, deepcopy(in), matrix_double_t(in)
end

function skew(u)
    S = zeros(3,3)
    S[1,2] = -u[3];
    S[1,3] = +u[2];
    S[2,1] = +u[3];
    S[2,3] = -u[1];
    S[3,1] = -u[2];
    S[3,2] = +u[1];
    return S;
end

function SO3_hat(u)
    return skew(u);
end

function SO3_vee(gfrak)
    out = zeros(3,1);
    out[1] = gfrak[3,2];
    out[2] = gfrak[1,3];
    out[3] = gfrak[2,1];
    return out;
end

function SU2_hat(u)
    L1 = [ 0.0   im;  im  0.0];
    L2 = [ 0.0 -1.0; 1.0  0.0];
    L3 = [  im  0.0; 0.0  -im];
    gfrak = L1*u[1] + L2*u[2] + L3*u[3];
    return gfrak
end

function SU2_vee(gfrak)
    out = [ imag(gfrak[1,2] + gfrak[2,1])/2.0;
            real(gfrak[2,1] - gfrak[1,2])/2.0;
            imag(gfrak[1,1] - gfrak[2,2])/2.0];
    return out
end

function quat_2_SU2(qX::Array{Float64, 2})
    X = zeros(2,2) + zeros(2,2)*im;
    X[1,1] = +qX[1]+qX[4]*im
    X[1,2] = -qX[3]+qX[2]*im
    X[2,1] = +qX[3]+qX[2]*im
    X[2,2] = +qX[1]-qX[4]*im
    return X
end

function SU2_2_quat(X::Array{Complex{Float64}, 2})
    q = zeros(4,1);
    q[1,1] = real(X[1,1])
    q[2,1] = imag(X[2,1])
    q[3,1] = real(X[2,1])
    q[4,1] = imag(X[1,1])
    return q
end

function quat_2_SO3(q::Array{Float64, 2})
    a, b, c, d = q[1,1], q[2,1], q[3,1], q[4,1];
    R = zeros(3,3);
    R[1,1] = a^2+b^2-c^2-d^2;
    R[1,2] = 2*(b*c-a*d);
    R[1,3] = 2*(b*d+a*c);
    R[2,1] = 2*(b*c+a*d);
    R[2,2] = a^2-b^2+c^2-d^2;
    R[2,3] = 2*(c*d-a*b);
    R[3,1] = 2*(b*d-a*c);
    R[3,2] = 2*(c*d+a*b);
    R[3,3] = a^2-b^2-c^2+d^2;
    return R
end

function quaternion_left_product(q)
    qim = q[2:4]
    qre = q[1]
    QL = zeros(4,4)
    QL[1,2:4]   = -qim;
    QL[2:4,1]   = +qim;
    QL[2:4,2:4] = skew(qim);
    QL = QL + I*qre;
    return QL
end


function quaternion_right_product(q)
    qim = q[2:4]
    qre = q[1]
    QR = zeros(4,4);
    QR[1,2:4]   = -qim;
    QR[2:4,1]   = +qim;
    QR[2:4,2:4] = -skew(qim);
    QR = QR + I*qre;
    return QR
end

function attitude_dynamics(q, w, J, tau)
    we      = zeros(4,1);
    we[2:4] = w;
    dq      = quaternion_left_product(q) * we / 2.0;
    dw      = (J) \ (-skew(w)*(J*w) + tau);
    return dq, dw
end

function ntuple_2_mat(tup, Nrow, Ncol)
    mat = zeros(Nrow, Ncol)
    for i = 1:Nrow
        for j = 1:Ncol
            mat[i,j] = tup[i][j];
        end
    end
    return mat;
end

function rand_PSD(N)
    J  = rand(N,N); J  = J + J'; J  = J + I*maximum(abs.(eigvals(J)));
    return J;
end
