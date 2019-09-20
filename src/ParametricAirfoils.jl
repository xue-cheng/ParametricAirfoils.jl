module ParametricAirfoils

using JuMP
using Printf

abstract type AbstractAirfoil{T<:AbstractFloat} end

dy_lower(::AbstractAirfoil, xc) = error("no gradient information")
dy_upper(::AbstractAirfoil, xc) = error("no gradient information")
fy_lower(::AbstractAirfoil, xc) = error("no gradient information")
fy_upper(::AbstractAirfoil, xc) = error("no gradient information")
x_lower(::AbstractAirfoil{T}, xc::T) where {T} = xc
x_upper(::AbstractAirfoil{T}, xc::T) where {T} = xc

include("CST/CST.jl")
include("NACA/NACA.jl")
include("utils.jl")

export AbstractAirfoil, CST, NACA, @NACA_str
export x_upper, x_lower, y_upper, y_lower, dy_upper, dy_lower, fy_upper, fy_lower
export gen_airfoil

end # module
