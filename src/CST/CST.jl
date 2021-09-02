include("utils.jl")

mutable struct CST{T,N,N1,N2} <: AbstractAirfoil{T}
    class::CSTClassFunction{N1,N2}
    shape::CSTShapeFunction{N,T}
    au::Vector{T} 
    al::Vector{T}
    a0::T
    yu::T
    yl::T
end

function CST(::UndefInitializer, ::Type{T}, N::Int, N1=0.5, N2=1) where {T <: AbstractFloat}
    N1_ = convert(T, N1)
    N2_ = convert(T, N2)
    class = CSTClassFunction(N1_, N2_)
    shape = CSTShapeFunction(N, T)
    CST{T,N,N1_,N2_}(class, shape, Vector{T}(undef, N), Vector{T}(undef, N), 0, 0, 0)
end

function CST(a::AbstractVector{T}, N1=0.5, N2=1, T_TE=0, Y_TE=0) where {T <: AbstractFloat}
    N::Int, rem = divrem(length(a), 2)
    rem == 1 || error("CST params must be of length 2*N+1")
    N1_ = convert(T, N1)
    N2_ = convert(T, N2)
    class = CSTClassFunction(N1_, N2_)
    shape = CSTShapeFunction(N, T)
    yu = convert(T, Y_TE + T_TE / 2)
    yl = convert(T, Y_TE - T_TE / 2)
    a0 = a[1]
    au = a[2:N + 1]
    al = a[N + 2:2N + 1]
    CST{T,N,N1_,N2_}(class, shape, au, al, a0, yu, yl)
end

function set_param!(cst::CST{T,N,N1,N2}, a::AbstractVector{T}) where {T,N,N1,N2}
    length(a) == 2N + 1 || error("CST params must be of length 2*N+1")
    cst.a0 = a[1]
    cst.au = a[2:N + 1]
    cst.al = a[N + 2:2N + 1]
end

get_param(cst::CST{T}) where {T} = vcat(cst.a0, cst.au, cst.al)

function modify_te!(cst::CST{T}; T_TE=nothing, Y_TE=nothing) where T
    if T_TE === nothing
        T_TE = cst.yu - cst.yl
    end
    if Y_TE === nothing
        Y_TE = (cst.yu + cst.yl) / 2
    end
    cst.yu = convert(T, Y_TE + T_TE / 2)
    cst.yl = convert(T, Y_TE - T_TE / 2)
end


@inline shape_f(cst::CST, i, x) = shape_f(cst.shape, i, x)
@inline shape_d(cst::CST, i, x) = shape_d(cst.shape, i, x)
@inline class_f(cst::CST, x) = class_f(cst.class, x)
@inline class_d(cst::CST, x) = class_d(cst.class, x)

function Base.show(io::IO, cst::CST{T,N,N1,N2}) where {T,N,N1,N2}
    println(io, "CST{$T,$N,$N1,$N2}")
    println(io, "  Upper surface:", vcat(cst.a0, cst.au))
    println(io, "  Lower surface:", vcat(-cst.a0, cst.al))
    println(io, "  Trailing edge: yu=", cst.yu, " yl=", cst.yl)
end

"""
    fit(:CST, xu, yu, xl, yl; type=Float64, N=9, N1=0.5, N2=1.0)
"""
function fit(::Val{:CST}, xx, yy; T=Float64, N=9, N1=0.5, N2=1.0, kw...)
    if (!isempty(kw))
        @warn "Unused keyword arguments:" kw...
    end
    N1_ = convert(T, N1)
    N2_ = convert(T, N2)
    xu, yu, xl, yl = cst_preprocess(convert.(T, xx), convert.(T, yy))
    @assert issorted(xu) && issorted(xl)
    @assert xu[1] == xl[1] == yu[1] == yl[1] == 0
    @assert (xu[end] + xl[end])/2 ≈ 1
    class = CSTClassFunction(N1_, N2_)
    shape = CSTShapeFunction(N, T)
    au, al, a0, yu1, yl1 = cst_fit(class, shape, xu, yu, xl, yl)
    CST{T,N,N1,N2}(class, shape, au, al, a0, yu1, yl1)
end

"""
    fit!(cst, xu, yu, xl, yl)
    
"""
function fit!(cst::CST{T,N,N1,N2}, xx, yy) where {T,N,N1,N2}
    xu, yu, xl, yl = cst_preprocess(convert.(T, xx), convert.(T, yy))
    @assert issorted(xu) && issorted(xl)
    @assert xu[1] == xl[1] == yu[1] == yl[1] == 0
    @assert (xu[end] + xl[end])/2 ≈ 1
    cst.au, cst.al, cst.a0, cst.yu, cst.yl = cst_fit(cst.class, cst.shape, xu, yu, xl, yl)
    cst
end

function y_upper(cst::CST{T,N}, x::T)::T where {T <: AbstractFloat,N}
    y = cst.a0 * shape_f(cst, 0, x) 
    for i = 1:N
        y += cst.au[i] * shape_f(cst, i, x)
    end
    y * class_f(cst, x) + cst.yu * x
end

function y_lower(cst::CST{T,N}, x::T)::T where {T <: AbstractFloat,N}
    y = -cst.a0 * shape_f(cst, 0, x) 
    for i = 1:N
        y += cst.al[i] * shape_f(cst, i, x)
    end
    y * class_f(cst, x) + cst.yl * x
end

function dy_upper(cst::CST{T,N}, x::T)::T where {T <: AbstractFloat,N}
    dy = cst.a0 * shape_d(cst, 0, x)
    y = cst.a0 * shape_f(cst, 0, x)
    for i = 1:N
        dy += cst.au[i] * shape_d(cst, i, x)
        y += cst.au[i] * shape_f(cst, i, x)
    end
    y * class_d(cst, x) + dy * class_f(cst, x) + cst.yu
end

function dy_lower(cst::CST{T,N}, x::T)::T where {T <: AbstractFloat,N}
    dy = -cst.a0 * shape_d(cst, 0, x)
    y = -cst.a0 * shape_f(cst, 0, x)
    for i = 1:N
        dy += cst.al[i] * shape_d(cst, i, x)
        y += cst.al[i] * shape_f(cst, i, x)
    end
    y * class_d(cst, x) + dy * class_f(cst, x) + cst.yl
end


function fy_upper(cst::CST{T,N}, x::T)::Tuple{T,T} where {T <: AbstractFloat,N}
    dy = cst.a0 * shape_d(cst, 0, x)
    y = cst.a0 * shape_f(cst, 0, x)
    for i = 1:N
        dy += cst.au[i] * shape_d(cst, i, x)
        y += cst.au[i] * shape_f(cst, i, x)
    end
    C = class_f(cst, x)
    dC = class_d(cst, x)
    y * C + cst.yu * x, y * dC + dy * C + cst.yu
end

function fy_lower(cst::CST{T,N}, x::T)::Tuple{T,T} where {T <: AbstractFloat,N}
    dy = -cst.a0 * shape_d(cst, 0, x)
    y = -cst.a0 * shape_f(cst, 0, x)
    for i = 1:N
        dy += cst.al[i] * shape_d(cst, i, x)
        y += cst.al[i] * shape_f(cst, i, x)
    end
    C = class_f(cst, x)
    dC = class_d(cst, x)
    y * C + cst.yl * x, y * dC + dy * C + cst.yl
end
