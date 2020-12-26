using LinearAlgebra

# The matrix double struct
mutable struct matrix_double_t
    numRows::Cint
    numCols::Cint
    pData::Ptr{Cdouble}
    function matrix_double_t(data)
        @assert(0 < size(data,1))
        @assert(0 < size(data,2))
        return new(Cint(size(data,1)), Cint(size(data,2)), pointer(data))
    end
end

mutable struct ref_state_qw_t
    quaternion::NTuple{4, Cdouble}
    attrates::NTuple{3, Cdouble}
    attaccel::NTuple{3, Cdouble}
    thrust::Cdouble
    function ref_state_qw_t(
        q::Array{T,2},
        w::Array{T,2},
        a::Array{T,2},
        t::Cdouble
    ) where {T<:Number}
        @assert(size(q, 1) == 4)
        @assert(size(q, 2) == 1)
        @assert(size(w, 1) == 3)
        @assert(size(w, 2) == 1)
        @assert(size(a, 1) == 3)
        @assert(size(a, 2) == 1)
        @assert(t >= 0.0)
        qin = ntuple(i->Cdouble(q[i]), 4)
        win = ntuple(i->Cdouble(w[i]), 3)
        ain = ntuple(i->Cdouble(a[i]), 3)
        new(qin, win, ain, t)
    end
end

mutable struct dyn_state_qw_t
    quaternion::NTuple{4, Cdouble}
    attrates::NTuple{3, Cdouble}
    function dyn_state_qw_t(
        q::Array{T,2},
        w::Array{T,2}
    ) where {T<:Number}
        @assert(size(q, 1) == 4)
        @assert(size(q, 2) == 1)
        @assert(size(w, 1) == 3)
        @assert(size(w, 2) == 1)
        qin = ntuple(i->Cdouble(q[i]), 4)
        win = ntuple(i->Cdouble(w[i]), 3)
        new(qin, win)
    end
end

mutable struct con_state_qw_fsf_t
    status::Cint
    thrust::Cdouble
    torques::NTuple{3, Cdouble}
    inertia::NTuple{3,NTuple{3,Cdouble}}
    gain_kR::Cdouble
    gain_kc::Cdouble
    gain_kw::Cdouble
    gain_eps::Cdouble
    gain_L::Cdouble
    param_a::Cdouble
    param_b::Cdouble
    param_c::Cdouble
    param_d::Cdouble
    dist_Psi::Cdouble
    dist_Gamma::Cdouble
    dist_lyapunov::Cdouble
    function con_state_qw_fsf_t(
        J::Array{T,2},
        kR::Float64,
        kc::Float64,
        kw::Float64;
        eps::Float64 = 1e-2,
        L::Float64   = 1.0,
        a::Float64   = 1.0,
        b::Float64   = 1.0,
        c::Float64   = 1.0,
        d::Float64   = 1.0
    ) where {T<:Number}
        @assert(size(J,2) == 3)
        @assert(size(J,2) == 3)
        @assert(kR  > 0)
        @assert(kc  > 0)
        @assert(kw  > 0)
        @assert(eps > 0)
        @assert(L   > 0)
        @assert(a   > 0)
        @assert(b   > 0)
        @assert(c   > 0)
        @assert(d   > 0)
        Tin = ntuple(i->Cdouble(0.0), 3)
        Jin = ntuple(i->ntuple(j->Cdouble(J[i,j]), 3), 3)
        new(0, 0, Tin, Jin, kR, kc, kw, eps, L, a, b, c, d, 0, 0, 0)
    end
end

mutable struct con_state_qw_fof_t
    status::Cint
    torques::NTuple{3, Cdouble}
    inertia::NTuple{3,NTuple{3,Cdouble}}
    invinertia::NTuple{3,NTuple{3,Cdouble}}
    gain_Kw::NTuple{3,NTuple{3,Cdouble}}
    gain_ki::NTuple{3,Cdouble}
    gain_Cw::NTuple{3,NTuple{3,Cdouble}}
    gain_cw::Cdouble
    gain_cR::Cdouble
    param_a::Cdouble
    param_b::Cdouble
    param_c::Cdouble
    param_d::Cdouble
    measuredGyrorates::NTuple{3,Cdouble}
    measuredDirections::NTuple{3,NTuple{3,Cdouble}}
    globalDirections::NTuple{3,NTuple{3,Cdouble}}
    function con_state_qw_fof_t(
        J::Array{T,2},
        Kw::Array{T,2},
        ki::Array{T,2},
        Cw::Array{T,2},
        cw::Float64,
        cR::Float64,
        a::Float64,
        b::Float64,
        c::Float64,
        d::Float64,
        Y0::Array{T,2},
        Yi::Array{T,2},
        Vi::Array{T,2}
    ) where {T<:Number}
        # Check input dimensions
        @assert((3,3) == size(J))
        @assert((3,3) == size(Kw))
        @assert((1,3) == size(ki))
        @assert((3,3) == size(Cw))
        @assert((3,1) == size(Y0))
        @assert((3,3) == size(Yi))
        @assert((3,3) == size(Vi))

        # Check input characteristics
        @assert(cw > 0);
        @assert(cR > 0);
        @assert(a  > 0);
        @assert(b  > 0);
        @assert(c  > 0);
        @assert(d  > 0);

        # Form Ntuples
        invJ    = inv(J);
        T_in    = ntuple(i->Cdouble(0.0), 3);
        invJ_in = ntuple(i->ntuple(j->Cdouble(invJ[i,j]), 3), 3);
        J_in    = ntuple(i->ntuple(j->Cdouble(   J[i,j]), 3), 3);
        Kw_in   = ntuple(i->ntuple(j->Cdouble(Kw[i,j]), 3), 3);
        ki_in   = ntuple(i->Cdouble(ki[1,i]), 3);
        Cw_in   = ntuple(i->ntuple(j->Cdouble(Cw[i,j]), 3), 3);
        Y0_in   = ntuple(i->Cdouble(Y0[i,1]), 3);
        Yi_in   = ntuple(i->ntuple(j->Cdouble(Yi[j,i]), 3), 3);
        Vi_in   = ntuple(i->ntuple(j->Cdouble(Vi[j,i]), 3), 3);

        new(0, T_in, J_in, invJ_in, Kw_in, ki_in, Cw_in, cw, cR, a, b, c, d, Y0_in, Yi_in, Vi_in);
    end
end
