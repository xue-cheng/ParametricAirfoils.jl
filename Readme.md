# ParametricAirfoils.jl #
A Parametric Airfoils package for julia.

![ci](https://github.com/xue-cheng/ParametricAirfoils.jl/actions/workflows/ci.yml/badge.svg)
[![Coverage Status](https://coveralls.io/repos/github/xue-cheng/ParametricAirfoils.jl/badge.svg)](https://coveralls.io/github/xue-cheng/ParametricAirfoils.jl)

## Examples ##

```julia
using ParametricAirfoils, AirfoilDatabase
########## STEP 1. Select a parameterization method ##########
# 1.1. CST Method
N = 9 # num of params for **EACH** curve
cst = CST(N) # == CST(N, Float64; N1=0.5, N2=1)
csts = CST(N; link_le=false) # Different `A0` for upper and lower surfaces, **NOT recommended**
# 1.2. Mean Camber + Thickness (MCT)
# 1.2.1 MCT(CST+CST)
mctcc = MCT(CurveCST(N, 1, 1, Float64), CurveCST(N, 0.5, 1.0, Float64))
# 1.2.2 MCT(CST+Poly)
indicies = tuple((max(0.5, i-1) for i in 1:N)...) # poly nomial indicies
mctcp = MCT(CurveCST(N, 1, 1, Float64), CurvePoly(indicies, Float64))
# 1.2.3 More curve types can be derived from `CurveFunction`, see `src/curve/_init.jl`
########## STEP 2. Create ParametricAirfoil ##########
# 2.1. Fitting airfoil data
# install AirfoilDatabase.jl
# run in REPL: `] add AirfoilDatabase`
result = query_airfoil("N0012")
@assert length(result)==1
NACA0012 = result[1]
airfoil, fiterr = fit(cst, NACA0012.x, NACA0012.y)
# `cst` can be replaced with any parameterization method
# 2.2. Create from parameters
param = get_param(airfoil) # any vector of length `num_param(cst)`
airfoil2 = ParametricAirfoil(cst, param)
########## STEP 3. Airfoil Coordinates ##########
xm = LinRange(0, 1, 33)
yu = airfoil(:U, xm) # Upper side
yl = airfoil(:L, xm) # Lower side
# OR with gradient information
yu, dyu =  airfoil(:U, xm, true)
# OR generate points at once for meshing 
x, y = gen_airfoil(airfoil)
```
