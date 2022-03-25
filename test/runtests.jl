using AirfoilDatabase
using Printf
using Test
using ParametricAirfoils


result = query_airfoil("N0012")
@assert length(result)==1
NACA0012 = result[1]
result = query_airfoil("RAE2822")
@assert length(result)==1
RAE2822 = result[1]

NCST = 10
NPoly = 10
indicies = tuple((max(0.5, i-1) for i in 1:NPoly)...)

include("testcst.jl")
include("testmct.jl")