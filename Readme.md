# ParametricAirfoils.jl #

A Parametric Airfoils package for julia.

## Examples ##

```julia
using ParametricAirfoils, JuMP, Ipopt, Plots

naca0012 = NACA"0012"s # `s` for sharp trailing edge

xᵤ, yᵤ, xₗ, yₗ = gen_airfoil(naca0012, 129)

cst = CST(9, xᵤ, yᵤ, xₗ, yₗ, with_optimizer(Ipopt.Optimizer))

yᵤ₉ = map(xᵤ) do x
  y_upper(cst, x)
end

dyᵤ = map(xᵤ) do x
  dy_upper(cst, x)
end

y, dy = dy_upper(cst, 0.3)

plot(xᵤ, yᵤ₉ - yᵤ, label=:error)
```
