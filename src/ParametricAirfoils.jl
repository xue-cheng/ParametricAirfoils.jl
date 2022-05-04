module ParametricAirfoils

using LinearAlgebra
using Printf

include("curve/_init.jl")
include("methods/_init.jl")
include("airfoil.jl")

export ParametricAirfoil
export CurveFunction, CurveCST, CurvePoly
export Parameterization, CST, MCT
export num_hyper_param, get_hyper_param, set_hyper_param!
export num_param, set_param!, get_param, fit, fit!, gen_airfoil

end # module
