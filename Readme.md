# ParametricAirfoils.jl #
A Parametric Airfoils package for julia.

[![Build Status](https://app.travis-ci.com/xue-cheng/ParametricAirfoils.jl.svg?branch=master)](https://app.travis-ci.com/xue-cheng/ParametricAirfoils.jl) 
[![codecov](https://codecov.io/gh/xue-cheng/ParametricAirfoils.jl/branch/master/graph/badge.svg?token=8JqTeI0yto)](https://codecov.io/gh/xue-cheng/ParametricAirfoils.jl)

## Examples ##

```julia
using ParametricAirfoils, Plots

naca0012 = NACA"0012"s # `s` for sharp trailing edge
# NACA0015 = NACA"0015"

x1, y1 = gen_airfoil(naca0012, 129) 
@asser length(x1) == 2*129 - 1

cst = fit(:CST, x1, y1; N=9)
x2, y2 = gen_airfoil(cst, 129) 
x2r, y2r = gen_airfoil(cst, 129, orient=:CCW) # `orient` can be :CW or :CCW, default value is :CW
yu = y_upper(cst, 0.2)
dyu = dy_upper(cst, 0.2)
yu, dyu = fy_upper(cst, 0.2) # eval y & dy in one call
nx, ny = n_upper(cst, 0.3) # surface normal 
tx, ty = t_upper(cst, 0.3) # surface tangent
p = plot()
plot!(p, x1, y1, label="NACA 0012")
plot!(p, x2, y2, label="CST N=9")
```
