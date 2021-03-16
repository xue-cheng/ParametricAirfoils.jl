include("utils.jl")

struct CST{T <: AbstractFloat,N,N1,N2} <: AbstractAirfoil{T}
    au::Vector{T} 
    al::Vector{T}
    a0::T
    yu::T
    yl::T
end

function CST(params::AbstractArray{T}; N::Int, N1::T=0.5, N2::T=0) where T
    @assert length(params) == 2N + 3
    au = params[1:N]
    al = params[N + 1:2N]
    a0 = params[2N + 1]
    yu = params[2N + 2]
    yl = params[2N + 3]
    CST{T,N,N1,N2}(au, al, a0, yu, yl)
end

"""
    fit(:CST, xu, yu, xl, yl; N=N, [N1=N1, N2=N2])
fit CST model with order `N`
"""
function fit(v::Val{:CST}, xu, yu, xl, yl;N1=0.5, N2=1.0, kw...)
    xu1, yu1, xl1, yl1 = promote(xu, yu, xl, yl)
    T = eltype(xu1)
    fit(
        v,
        xu1, yu1, xl1, yl1;
        N1=convert(T, N1),
        N2=convert(T, N2)
    )
end

function fit(
    ::Val{:CST}, 
    xu::AbstractVector{T},
    yu::AbstractVector{T},
    xl::AbstractVector{T},
    yl::AbstractVector{T};
    N::Int, N1::T=0.5, N2::T=1.0) where {T}
    @assert issorted(xu)
    @assert issorted(xl)
    @assert xu[1] == xl[1] == 0
    @assert yu[1] == yl[1] == 0
    A, b = cstMatrix(xu, yu, xl, yl, N, N1, N2)
    F = svd(A, alg=LinearAlgebra.QRIteration())
    b1 = transpose(F.U) * b
    y = b1 ./ F.S
    params = F.V * y
    CST(params; N=N, N1=N1, N2=N2)
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
