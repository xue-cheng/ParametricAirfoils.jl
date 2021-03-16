using NLopt
@inline cst_S(n, i, x) = binomial(n, i) * x^i * (1 - x)^(n - i)
@inline cst_C(n1, n2, x) = x^n1 * (1 - x)^n2
@inline function cst_dS(n, i, x)
    if i == 0 # (1-x)^n
        -n * (1 - x)^(n - 1)
    elseif i == n # x^n
        n * x^(n - 1)
    else
        binomial(n, i) * (i * x^(i - 1) * (1 - x)^(n - i) - (n - i) * x^i * (1 - x)^(n - i - 1))
    end
end

@inline function cst_dC(n1, n2, x) 
    n1 * x^(n1 - 1) * (1 - x)^n2 - n2 * x^n1 * (1 - x)^(n2 - 1)
end

function fit_cst(n::Int, dz::T,
    xu::AbstractVector{T},yu::AbstractVector{T},
    xl::AbstractVector{T},yl::AbstractVector{T},
    algorithm::Symbol; 
    c1::T=0.5,
    c2::T=1.0 ) where {T}
    m = Model(NLopt.Optimizer)
    set_optimizer_attribute(m, "algorithm", algorithm)
    lb = oneunit(T) * -2
    up = oneunit(T) * 2
    n_up = length(xu)
    n_lo = length(xl)
    @variable(m, lb <= A₀ <= up, start = 1.0)
    @variable(m, lb <= Aᵤ[1:n] <= up, start = 1.0)
    @variable(m, lb <= Aₗ[1:n] <= up, start = -1.0)
    @expression(m, yuc[j=1:n_up], cst_C(c1, c2, xu[j]) * ( A₀ * cst_S(n, 0, xu[j]) + sum(Aᵤ[i] * cst_S(n, i, xu[j]) for i = 1:n)) + xu[j] * dz)
    @expression(m, ylc[j=1:n_lo], cst_C(c1, c2, xl[j]) * (-A₀ * cst_S(n, 0, xl[j]) + sum(Aₗ[i] * cst_S(n, i, xl[j]) for i = 1:n)) - xl[j] * dz)
    @expression(m, fup, sum((yuc[j] - yu[j])^2 for j = 1:n_up))
    @expression(m, flo, sum((ylc[j] - yl[j])^2 for j = 1:n_lo))
    @objective(m, Min, fup + flo)
    JuMP.optimize!(m)
    status = termination_status(m)
    value(A₀), value.(Aᵤ), value.(Aₗ), status
end


