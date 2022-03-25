struct CurvePoly{N, T} <: CurveFunction{N, T}
    indices::NTuple{N, T}
    function CurvePoly(indices::NTuple{N, Real}, ::Type{T}=Float64) where {T, N}
        issorted(indices) || error("`indices` must be sortted")
        return new{N,T}(convert(NTuple{N, T}, indices))
    end
end

num_hyper_param(::CurvePoly)=0
get_hyper_param(::CurvePoly) = []
set_hyper_param!(::CurvePoly, ::AbstractVector) = nothing
set_hyper_param!(::CurvePoly, ::Int, ::Real) = nothing

function curve_eval(poly::CurvePoly{N,T}, param::P, x::X) where {N, T, X<:Real, P<:AbstractVector{T}}
    @assert length(param) == N
    y = zero(T)
    @inbounds for i in 1:N
        e = poly.indices[i]
        a = param[i]
        y += a*x^e
    end
    return y 
end

function curve_grad(poly::CurvePoly{N,T}, param::P, x::X) where {N, T, X<:Real, P<:AbstractVector{T}}
    @assert length(param) == N
    dy = zero(T)
    @inbounds for i in 1:N
        e = poly.indices[i]
        dy += e*x^(e-1)*param[i]
    end
    return dy 
end

function curve_full(poly::CurvePoly{N,T}, param::P, x::X) where {N, T, X<:Real, P<:AbstractVector{T}}
    @assert length(param) == N
    y = zero(T)
    dy = zero(T)
    @inbounds for i in 1:N
        e = poly.indices[i]
        a = param[i]
        y += a*x^e
        dy += e*x^(e-1)*a
    end
    return y, dy 
end

function curve_eval(poly::CurvePoly{N,T}, param::P, x::X) where {N, T, X<:AbstractVector, P<:AbstractVector{T}}
    @assert length(param) == N
    y = zeros(T, length(x))
    @inbounds for i in 1:N
        e = poly.indices[i]
        a = param[i]
        @. y += a*x^e
    end
    return y 
end

function curve_grad(poly::CurvePoly{N,T}, param::P, x::X) where {N, T, X<:AbstractVector, P<:AbstractVector{T}}
    @assert length(param) == N
    dy = zeros(T, length(x))
    @inbounds for i in 1:N
        e = poly.indices[i]
        a = param[i]
        @. dy += e*x^(e-1)*a
    end
    return dy 
end

function curve_full(poly::CurvePoly{N,T}, param::P, x::X) where {N, T, X<:AbstractVector, P<:AbstractVector{T}}
    @assert length(param) == N
    y = zeros(T, length(x))
    dy = zeros(T, length(x))
    @inbounds for i in 1:N
        e = poly.indices[i]
        a = param[i]
        @. y += a*x^e
        @. dy += e*x^(e-1)*a
    end
    return y, dy 
end

function curve_fit(poly::CurvePoly{N, T}, x::AbstractVector, y::AbstractVector) where {N, T}
    issorted(x) || error("`x` must be sorted")
    length(x) == length(y) || error("`x` and `y` must have the same length")
    xs = convert(Vector{T}, x)
    ys = convert(Vector{T}, y)
    A = Matrix{T}(undef, length(xs), N)
    @inbounds for i in 1:N
        e = poly.indices[i]
        for p in eachindex(xs)
            A[p, i] = xs[i]^e
        end
    end
    p = A\ys
    e = sqrt(maximum(abs2, A*p - ys))
    return p, e
end