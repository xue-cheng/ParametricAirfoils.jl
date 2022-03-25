
abstract type Parameterization{N,T} end

num_param(::Parameterization{N}) where {N} = N
Base.eltype(::Parameterization{N,T}) where {N,T} = T

num_hyper_param(::Parameterization) = error("Not implemented")
get_hyper_param(::Parameterization) = error("Not implemented")
set_hyper_param!(::Parameterization, hyp) = error("Not implemented")
set_hyper_param!(::Parameterization, i, val) = error("Not implemented")

airfoil_fit(::Parameterization,x, y) = error("Not implemented")

airfoil_coord(::Parameterization,side,param,x) = error("Not implemented")
airfoil_grad(::Parameterization,side,param,x) = error("Not implemented")
airfoil_full(::Parameterization,side,param,x) = error("Not implemented")

function _split_points(::Type{T}, x::AbstractVector, y::AbstractVector) where T
    xs = convert(Vector{T}, x)
    ys = convert(Vector{T}, y)
    ile = argmin(xs)
    xu = xs[ile:end] .- xs[ile]
    yu = ys[ile:end] .- ys[ile]
    xl = xs[ile:-1:1] .- xs[ile]
    yl = ys[ile:-1:1] .- ys[ile]
    if yl[end-1] > yu[end-1]
        xu, xl = xl, xu
        yu, yl = yl, yu
    end
    issorted(xu) && issorted(xl) || error("surface points must be sortted")
    return xu,yu,xl,yl
end

include("cst.jl")
include("mct.jl")