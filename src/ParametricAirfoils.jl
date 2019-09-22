module ParametricAirfoils

using JuMP
using Printf

abstract type AbstractAirfoil{T<:AbstractFloat} end

dy_lower(::AbstractAirfoil{T}, xc::T) where {T} = error("no gradient information")
dy_upper(::AbstractAirfoil{T}, xc::T) where {T} = error("no gradient information")
fy_lower(::AbstractAirfoil{T}, xc::T) where {T} = error("no gradient information")
fy_upper(::AbstractAirfoil{T}, xc::T) where {T} = error("no gradient information")
x_lower(::AbstractAirfoil{T}, xc::T) where {T} = xc
x_upper(::AbstractAirfoil{T}, xc::T) where {T} = xc
@inline x_lower(a::AbstractAirfoil{T}, xc::R) where {T, R<:Real} = x_lower(a, convert(T, xc))
@inline x_upper(a::AbstractAirfoil{T}, xc::R) where {T, R<:Real} = x_upper(a, convert(T, xc))
@inline y_lower(a::AbstractAirfoil{T}, xc::R) where {T, R<:Real} = y_lower(a, convert(T,xc))
@inline y_upper(a::AbstractAirfoil{T}, xc::R) where {T, R<:Real} = y_upper(a, convert(T,xc))
@inline dy_lower(a::AbstractAirfoil{T}, xc::R) where {T, R<:Real} = dy_lower(a, convert(T,xc))
@inline dy_upper(a::AbstractAirfoil{T}, xc::R) where {T, R<:Real} = dy_upper(a, convert(T,xc))
@inline fy_lower(a::AbstractAirfoil{T}, xc::R) where {T, R<:Real} = fy_lower(a, convert(T,xc))
@inline fy_upper(a::AbstractAirfoil{T}, xc::R) where {T, R<:Real} = fy_upper(a, convert(T,xc))

include("CST/CST.jl")
include("NACA/NACA.jl")
include("utils.jl")

export AbstractAirfoil, CST, NACA, @NACA_str
export x_upper, x_lower, y_upper, y_lower, dy_upper, dy_lower, fy_upper, fy_lower
export gen_airfoil, n_upper, n_lower

end # module
