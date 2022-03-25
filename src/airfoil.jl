"""
Parametric Airfoil Type `struct ParametricAirfoil{M, N, T}` 
where `M` is the parameterization method, `N` is the number of parameters and `T` is the floating-point type.
"""
struct ParametricAirfoil{N, T<:AbstractFloat, M<:Parameterization{N,T}}
    method::M
    params::Vector{T}
    function ParametricAirfoil(method::M) where {N, T, M<:Parameterization{N, T}}
        return new{N,T,M}(method, Vector{T}(undef,N))
    end
    function ParametricAirfoil(method::M, param) where {N, T, M<:Parameterization{N, T}}
        return new{N,T,M}(method, convert(Vector{T}, param))
    end
    function ParametricAirfoil(other::A, param) where {N, T, M, A<:ParametricAirfoil{N, T, M}}
        return new{N,T,M}(other.method, convert(Vector{T}, param))
    end
end

num_param(::ParametricAirfoil{N}) where {N} = N
set_param!(a::ParametricAirfoil, param) = a.params .= param
get_param(a::ParametricAirfoil) = copy(a.params)

function (a::ParametricAirfoil)(side::Symbol, x, with_grad::Bool=false)
    if with_grad
        airfoil_full(a.method, side, a.params, x)
    else
        airfoil_coord(a.method, side, a.params, x)
    end
end

function fit!(a::ParametricAirfoil, x, y)
    p, e = airfoil_fit(a.method, x, y)
    a.params .= p
    return a, e
end

function fit(m::Parameterization, x, y)
    p, e = airfoil_fit(m, x, y)
    a = ParametricAirfoil(m, p)
    return a, e
end

function linspace(start::T, stop::T, np::I) where {T <: AbstractFloat,I <: Integer}
    d = (stop - start) / (np - 1)
    p = Vector{T}(undef, np)
    p[1] = start
    p[end] = stop
    for i = 2:np - 1
        p[i] = start + (i - 1) * d
    end
    return p
end

function cosspace(start::T, stop::T, np::I) where {T <: AbstractFloat,I <: Integer}
    p = linspace(zero(T), convert(T, π), np)
    r = (stop - start) / 2
    for i = 1:np
        p[i] = start + r * (one(T) - cos(p[i]))
    end
    return p
end

function gen_airfoil(a::ParametricAirfoil{N, T}, np::Int=65; cos_space::Bool=true, orient::Symbol=:CW) where {N, T <: AbstractFloat}
    if cos_space
        xc = cosspace(zero(T), one(T), np)
    else
        xc = linspace(zero(T), one(T), np)
    end
    yu = a(:U, xc)
    yl = a(:L, xc)
    x = vcat(xc[end:-1:2],xc)
    y = if orient == :CW
        vcat(yl[end:-1:2],yu)
    elseif orient == :CCW
        vcat(yu[end:-1:2],yl)
    else
        throw(ArgumentError("invalid orient=$orient, must be :CW or :CCW"))
    end
    return x, y
end
