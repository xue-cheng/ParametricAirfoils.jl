include("utils.jl")

struct CST{T <: AbstractFloat,N,N1,N2} <: AbstractAirfoil{T}
    au::Vector{T} 
    al::Vector{T}
    a0::T
    yu::T
    yl::T
end

function CST(params::AbstractArray{T}; N::Int, N1::T=0.5, N2::T=1) where T
    @assert length(params) == 2N + 3
    au = params[1:N]
    al = params[N + 1:2N]
    a0 = params[2N + 1]
    yu = params[2N + 2]
    yl = params[2N + 3]
    CST{T,N,N1,N2}(au, al, a0, yu, yl)
end

"""
    fit(:CST, xu, yu, xl, yl; N=N, N1=N1, N2=N2)
fit CST model with order `N` and class function indices `N1` and `N2`. 
An optimizer is utilized to determin `N1` and `N2` if any of them is `:OPT`.
"""
function fit(v::Val{:CST}, xu, yu, xl, yl; N=5, N1=:OPT, N2=:OPT, kw...)
    @assert issorted(xu)
    @assert issorted(xl)
    @assert xu[1] == xl[1] == 0
    @assert yu[1] == yl[1] == 0
    xu1, yu1, xl1, yl1 = promote(xu, yu, xl, yl)
    T = eltype(xu1)
    if isa(N1, Number) && isa(N2, Number)
        cst_fit(N, convert(T, N1), convert(T, N2), xu1, yu1, xl1, yl1)
    else
        N1_ = isa(N1, Symbol) ? N1 : convert(T, N1)
        N2_ = isa(N2, Symbol) ? N2 : convert(T, N2)
        cst_fit_opt(N, N1_, N2_, xu1, yu1, xl1, yl1)
    end
end

    function lsfit(A, b)
    U, S, V = svd(A)
    V * ((transpose(U) * b) ./ S)
end

function cst_fit(N::Int, N1::T, N2::T,
    xu::AbstractVector{T}, yu::AbstractVector{T},
    xl::AbstractVector{T}, yl::AbstractVector{T} ) where {T}
    A, b = cstMatrix(xu, yu, xl, yl, N, N1, N2)
    params = lsfit(A, b)
    CST(params; N=N, N1=N1, N2=N2)
end

function cst_fit_opt(N::Int, N1::Symbol, N2::Symbol,
    xu::AbstractVector{T}, yu::AbstractVector{T},
    xl::AbstractVector{T}, yl::AbstractVector{T}
    ) where {T}
    function objective(x::Vector, ::Vector)
        N1_ = convert(T, x[1])
        N2_ = convert(T, x[2])
        A, b = cstMatrix(xu, yu, xl, yl, N, N1_, N2_)
        cst = CST(lsfit(A, b); N=N, N1=N1_, N2=N2_)
        yu1 = y_upper(cst, xu)
        yl1 = y_lower(cst, xl)
        @. yu1 = abs(yu1 - yu)
        @. yl1 = abs(yl1 - yl)
        convert(eltype(x), max(maximum(yu1), maximum(yl1)))
    end
    opt = Opt(:LN_COBYLA, 2)
    opt.lower_bounds = [0.1, 0.1]
    opt.upper_bounds = [1.0, 1.0]
    opt.xtol_rel = 1e-4
    opt.min_objective = objective
    (minf, minx, ret) = optimize(opt, [0.5, 1.0])
    numevals = opt.numevals # the number of function evaluations
    @info "CSTOPT: got [N1,N2]=$minx after $numevals iters"
    cst_fit(N, convert(T, minx[1]), convert(T, minx[2]), xu, yu, xl, yl)
end

function cst_fit_opt(N::Int, N1::Symbol, N2::T,
    xu::AbstractVector{T}, yu::AbstractVector{T},
    xl::AbstractVector{T}, yl::AbstractVector{T}) where {T}
    function objective(x::Vector, ::Vector)
        N1_ = convert(T, x[1])
        A, b = cstMatrix(xu, yu, xl, yl, N, N1_, N2)
        cst = CST(lsfit(A, b); N=N, N1=N1_, N2=N2)
        yu1 = y_upper(cst, xu)
        yl1 = y_lower(cst, xl)
        @. yu1 = abs(yu1 - yu)
        @. yl1 = abs(yl1 - yl)
        convert(eltype(x), max(maximum(yu1), maximum(yl1)))
    end
    opt = Opt(:LN_COBYLA, 1)
    opt.lower_bounds = [0.1]
    opt.upper_bounds = [1.0]
    opt.xtol_rel = 1e-4
    opt.min_objective = objective
    (minf, minx, ret) = optimize(opt, [0.5])
    numevals = opt.numevals # the number of function evaluations
    @info "CSTOPT: got N1=$(minx[1]) after $numevals iters"
    cst_fit(N, convert(T, minx[1]), N2, xu, yu, xl, yl)
end

function cst_fit_opt(N::Int, N1::T, N2::Symbol,
    xu::AbstractVector{T}, yu::AbstractVector{T},
    xl::AbstractVector{T}, yl::AbstractVector{T}) where {T}
    function objective(x::Vector, ::Vector)
        N2_ = convert(T, x[1])
        A, b = cstMatrix(xu, yu, xl, yl, N, N1, N2_)
        cst = CST(lsfit(A, b); N=N, N1=N1, N2=N2_)
        yu1 = y_upper(cst, xu)
        yl1 = y_lower(cst, xl)
        @. yu1 = abs(yu1 - yu)
        @. yl1 = abs(yl1 - yl)
        convert(eltype(x), max(maximum(yu1), maximum(yl1)))
    end
    opt = Opt(:LN_COBYLA, 1)
    opt.lower_bounds = [0.1]
    opt.upper_bounds = [1.0]
    opt.xtol_rel = 1e-4
    opt.min_objective = objective
    (minf, minx, ret) = optimize(opt, [1.0])
    numevals = opt.numevals # the number of function evaluations
    @info "CSTOPT: got N2=$(minx[1]) after $numevals iters"
    cst_fit(N, N1, convert(T, minx[1]), xu, yu, xl, yl)
end

function y_upper(cst::CST{T,N,N1,N2}, x::T)::T where {T <: AbstractFloat,N,N1,N2}
    y = cst.a0 * cstShapeF(N, 0, x)
    for i = 1:N
        y += cst.au[i] * cstShapeF(N, i, x)
    end
    y * cstClassF(N1, N2, x) + cst.yu * x
end

function y_lower(cst::CST{T,N,N1,N2}, x::T)::T where {T <: AbstractFloat,N,N1,N2}
    y = -cst.a0 * cstShapeF(N, 0, x)
    for i = 1:N
        y += cst.al[i] * cstShapeF(N, i, x)
    end
    y * cstClassF(N1, N2, x) + cst.yl * x
end

function dy_upper(cst::CST{T,N,N1,N2}, x::T)::T where {T <: AbstractFloat,N,N1,N2}
    dy = cst.a0 * cstShapeD(N, 0, x)
    y = cst.a0 * cstShapeF(N, 0, x)
    for i = 1:N
        dy += cst.au[i] * cstShapeD(N, i, x)
        y += cst.au[i] * cstShapeF(N, i, x)
    end
    y * cstClassD(N1, N2, x) + dy * cstClassF(N1, N2, x) + cst.yu
end

function dy_lower(cst::CST{T,N,N1,N2}, x::T)::T where {T <: AbstractFloat,N,N1,N2}
    dy = -cst.a0 * cstShapeD(N, 0, x)
    y = -cst.a0 * cstShapeF(N, 0, x)
    for i = 1:N
        dy += cst.al[i] * cstShapeD(N, i, x)
        y += cst.al[i] * cstShapeF(N, i, x)
    end
    y * cstClassD(N1, N2, x) + dy * cstClassF(N1, N2, x) + cst.yl
end


function fy_upper(cst::CST{T,N,N1,N2}, x::T)::Tuple{T,T} where {T <: AbstractFloat,N,N1,N2}
    dy = cst.a0 * cstShapeD(N, 0, x)
    y = cst.a0 * cstShapeF(N, 0, x)
    for i = 1:N
        dy += cst.au[i] * cstShapeD(N, i, x)
        y += cst.au[i] * cstShapeF(N, i, x)
    end
    C = cstClassF(N1, N2, x)
    dC = cstClassD(N1, N2, x)
    y * C + cst.yu * x, y * dC + dy * C + cst.yu
end

function fy_lower(cst::CST{T,N,N1,N2}, x::T)::Tuple{T,T} where {T <: AbstractFloat,N,N1,N2}
    dy = -cst.a0 * cstShapeD(N, 0, x)
    y = -cst.a0 * cstShapeF(N, 0, x)
    for i = 1:N
        dy += cst.al[i] * cstShapeD(N, i, x)
        y += cst.al[i] * cstShapeF(N, i, x)
    end
    C = cstClassF(N1, N2, x)
    dC = cstClassD(N1, N2, x)
    y * C + cst.yl * x, y * dC + dy * C + cst.yl
end
