# ParametricAirfoils.jl #

A Parametric Airfoils package for julia.

## Examples ##

```julia
using ParametricAirfoils, JuMP, Plots

naca0012 = NACA"0012"s # `s` for sharp trailing edge

xᵤ, yᵤ, xₗ, yₗ = gen_airfoil(naca0012, 129)

gen_segment("NACA0012.dat", naca0012) # generate segment filefor mesh tools

cst = CST(9, xᵤ, yᵤ, xₗ, yₗ)

yᵤ₉ = map(xᵤ) do x
  y_upper(cst, x)
end

dyᵤ = map(xᵤ) do x
  dy_upper(cst, x)
end

y, dy = fy_upper(cst, 0.3)

plot(xᵤ, yᵤ₉ - yᵤ, label=:error)
```
