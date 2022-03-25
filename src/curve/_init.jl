abstract type CurveFunction{N,T} end

num_param(::CurveFunction{N}) where {N} = N
num_hyper_param(::CurveFunction) = error("Not implemented")
get_hyper_param(::CurveFunction) = error("Not implemented")
set_hyper_param!(::CurveFunction, hyp) = error("Not implemented")
set_hyper_param!(::CurveFunction, i, val) = error("Not implemented")
curve_eval(::CurveFunction, param, x) = error("Not implemented")
curve_grad(::CurveFunction, param, x) = error("Not implemented")
curve_full(::CurveFunction, param, x) = error("Not implemented")



include("cst.jl")
include("poly.jl")