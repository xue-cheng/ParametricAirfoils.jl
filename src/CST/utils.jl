struct CSTShapeFunction{N,T <: AbstractFloat}
    bino::NTuple{N,Int}
    CSTShapeFunction(N::Int, ::Type{T}) where T = new{N,T}(tuple((binomial(N, i) for i in 1:N)...))
end

struct CSTClassFunction{N1,N2,T}
    CSTClassFunction(N1::T, N2::T) where T = new{N1,N2,T}()
end

shape_f(s::CSTShapeFunction{N,T}, i::Int, x::T) where {N,T} =
if i == 0
    (1 - x)^N 
elseif i == N
    x^N
else
    s.bino[i] * x^i * (1 - x)^(N - i)
end

shape_d(s::CSTShapeFunction{N,T}, i::Int, x::T) where {N,T} = 
if i == 0 
    -N * (1 - x)^(N - 1)
elseif i == N 
    N * x^(N - 1)
else
    s.bino[i] * (i * x^(i - 1) * (1 - x)^(N - i) - (N - i) * x^i * (1 - x)^(N - i - 1))
end

class_f(::CSTClassFunction{N1,N2,T},x::T) where {N1,N2,T} = x^N1 * (1 - x)^N2

class_d(::CSTClassFunction{N1,N2,T},x::T) where {N1,N2,T} = N1 * x^(N1 - 1) * (1 - x)^N2 - N2 * x^N1 * (1 - x)^(N2 - 1)


function cst_preprocess(xx::AbstractVector{T}, yy::AbstractVector{T}) where T
    xte = (xx[1] + xx[end]) / 2
    
    imin = argmin(xx)
    x0 = xx[imin]
    y0 = yy[imin]
    cc = xte - x0
    xu_ = @view(xx[imin:-1:1])
    yu_ = @view(yy[imin:-1:1])
    xl_ = @view(xx[imin:end] )
    yl_ = @view(yy[imin:end] )
    xu = (xu_ .- x0) ./ cc
    yu = (yu_ .- y0) ./ cc
    xl = (xl_ .- x0) ./ cc
    yl = (yl_ .- y0) ./ cc
    xu, yu, xl, yl
end

function cst_matrix(
    class::CSTClassFunction{N1,N2,T}, 
    shape::CSTShapeFunction{N,T},
    xu::V, yu::V, xl::V, yl::V)　where {T,V <: AbstractVector{T},N1,N2,N}
    NP = length(xu) + length(xl) - 2
    A = zeros(T, NP, 2N + 3)
    b = vcat(yu[2:end], yl[2:end])
    ip = 1
    @inbounds for x in @view xu[2:end]
        c = class_f(class, x)
        for i in 1:N
            A[ip, i] = c * shape_f(shape, i, x) # cstShapeF(N, i, x)
        end
        A[ip, 2N + 1] = c * shape_f(shape, 0, x) # a_le
        A[ip, 2N + 2] = x # y_tu
        ip += 1
    end
    @inbounds for x in @view xl[2:end]
        c = class_f(class, x)
        for i in 1:N
            A[ip, i + N] = c * shape_f(shape, i, x)
        end
        A[ip, 2N + 1] = -c * shape_f(shape, 0, x) # -a_le
        A[ip, 2N + 3] = x # y_tl
        ip += 1
    end
    return A, b
end

function cst_fit(
    class::CSTClassFunction{N1,N2}, 
    shape::CSTShapeFunction{N,T},
    xu::AbstractVector{T}, yu::AbstractVector{T},
    xl::AbstractVector{T}, yl::AbstractVector{T} ) where {T,N,N1,N2}
    A, b = cst_matrix(class, shape, xu, yu, xl, yl)
    U, S, V = svd(A)
    params = V * ((transpose(U) * b) ./ S)
    au = params[1:N]
    al = params[N + 1:2N]
    a0 = params[2N + 1]
    yu = params[2N + 2]
    yl = params[2N + 3]
    return au, al, a0, yu, yl
end