struct CurveCST{N, T, M} <: CurveFunction{N,T}
    bino::NTuple{M, T}
    hyp::Vector{T}
    function CurveCST(N::Int, N1::Real, N2::Real, ::Type{T}=Float64) where T
        M = N - 2
        bino = tuple((binomial(M, i) for i in 1:M)...)
        hyp = T[N1, N2]
        return new{N,T,M}(bino, hyp)
    end
end

num_hyper_param(::CurveCST) = 2
get_hyper_param(c::CurveCST) = copy(c.hyp)
set_hyper_param!(c::CurveCST, hyp::AbstractVector) = c.hyp .= hyp
set_hyper_param!(c::CurveCST, i::Int, val::Real) = c.hyp[i] = val

Base.broadcastable(cst::CurveCST) = Ref(cst)
Base.show(io::IO, cst::CurveCST{N,T,M}) where {N,T,M} = print(io, "CSTFunction($M, $T;N1=$(cst.hyp[1]), N2=$(cst.hyp[2])")


function _cst_shape_f(s::CurveCST{N,T,M}, i::Int, x::T) where {N,T,M}
    return if i == 0
        (1 - x)^M 
    elseif i == M
        x^M
    else
        s.bino[i] * x^i * (1 - x)^(M - i)
    end
end

function _cst_shape_d(s::CurveCST{N,T,M}, i::Int, x::T) where {N,T,M}
    return if i == 0 
        -M * (1 - x)^(M - 1)
    elseif i == M 
        M * x^(M - 1)
    else
        s.bino[i] * (i * x^(i - 1) * (1 - x)^(M - i) - (M - i) * x^i * (1 - x)^(M - i - 1))
    end
end

function _cst_class_f(c::CurveCST{N, T}, x::T) where {N,T}
    N1, N2 = c.hyp
    return x^N1 * (1 - x)^N2
end

function _cst_class_d(c::CurveCST{N, T}, x::T) where {N,T} 
    N1, N2 = c.hyp
    omx = 1-x
    return N1 * x^(N1 - 1) * omx^N2 - N2 * x^N1 * omx^(N2 - 1)
end


function curve_eval(cst::CurveCST{N,T,M}, param::P, x::X) where {N, T, M, X<:Real, P<:AbstractVector{T}}
    @assert length(param) == N
    a0 = param[N]
    yt = param[N-1]
    y = a0 * _cst_shape_f(cst, 0, x)
    @inbounds for i in 1:M
        y += param[i]*_cst_shape_f(cst, i, x)
    end
    y = y * _cst_class_f(cst, x) + yt*x
    return y
end

function curve_grad(cst::CurveCST{N,T,M}, param::P, x::X) where {N, T, M, X<:Real, P<:AbstractVector{T}}
    @assert length(param) == N
    a0 = param[N]
    yt = param[N-1]
    dy = a0 * _cst_shape_d(cst, 0, x)
    y = a0 * _cst_shape_f(cst, 0, x)
    @inbounds for i in 1:M
        dy += param[i] * _cst_shape_d(cst, i, x)
        y += param[i] * _cst_shape_f(cst, i, x)
    end
    dy = y * _cst_class_d(cst, x) + dy * _cst_class_f(cst, x) + yt
    return dy 
end

function curve_full(cst::CurveCST{N,T,M}, param::P, x::X) where {N, T, M, X<:Real, P<:AbstractVector{T}}
    @assert length(param) == N
    a0 = param[N]
    yt = param[N-1]
    dy = a0 * _cst_shape_d(cst, 0, x)
    y = a0 * _cst_shape_f(cst, 0, x)
    @inbounds for i in 1:M
        dy += param[i] * _cst_shape_d(cst, i, x)
        y += param[i] * _cst_shape_f(cst, i, x)
    end
    dy = y * _cst_class_d(cst, x) + dy * _cst_class_f(cst, x) + yt
    y = y * _cst_class_f(cst, x) + yt*x
    return y, dy
end


function curve_eval(cst::CurveCST{N,T,M}, param::P, x::X) where {N, T, M, X<:AbstractVector, P<:AbstractVector{T}}
    @assert length(param) == N
    a0 = param[N]
    yt = param[N-1]
    y = similar(x,T)
    @. y = a0 * _cst_shape_f(cst, 0, x)
    @inbounds for i in 1:M
        @. y += param[i]*_cst_shape_f(cst, i, x)
    end
    @. y = y * _cst_class_f(cst, x) + yt*x
    return y
end

function curve_grad(cst::CurveCST{N,T,M}, param::P, x::X) where {N, T, M, X<:AbstractVector, P<:AbstractVector{T}}
    @assert length(param) == N
    a0 = param[N]
    yt = param[N-1]
    y = similar(x,T)
    dy = similar(x,T)
    @. dy = a0 * _cst_shape_d(cst, 0, x)
    @. y = a0 * _cst_shape_f(cst, 0, x)
    @inbounds for i in 1:M
        @. dy += param[i] * _cst_shape_d(cst, i, x)
        @. y += param[i] * _cst_shape_f(cst, i, x)
    end
    @. dy = y * _cst_class_d(cst, x) + dy * _cst_class_f(cst, x) + yt
    return dy
end

function curve_full(cst::CurveCST{N,T,M}, param::P, x::X) where {N, T, M, X<:AbstractVector, P<:AbstractVector{T}}
    @assert length(param) == N
    a0 = param[N]
    yt = param[N-1]
    y = similar(x,T)
    dy = similar(x,T)
    @. dy = a0 * _cst_shape_d(cst, 0, x)
    @. y = a0 * _cst_shape_f(cst, 0, x)
    @inbounds for i in 1:M
        @. dy += param[i] * _cst_shape_d(cst, i, x)
        @. y += param[i] * _cst_shape_f(cst, i, x)
    end
    @. dy = y * _cst_class_d(cst, x) + dy * _cst_class_f(cst, x) + yt
    @. y = y * _cst_class_f(cst, x) + yt*x
    return y, dy
end

function matrix_LSQ(cst::CurveCST{N,T,M}, x::AbstractVector{T}) where {N,T,M}
    xs = convert(Vector{T}, x)
    A = Matrix{T}(undef, length(xs), N)
    @inbounds for ip in eachindex(xs)
        x = xs[ip]
        c = _cst_class_f(cst, x)
        for p in 1:M
            A[ip, p] = c * _cst_shape_f(cst, p, x) # cstShapeF(N, i, x)
        end
        A[ip, M + 1] = x # y_te
        A[ip, M + 2] = c * _cst_shape_f(cst, 0, x) # a0
    end
    return A
end