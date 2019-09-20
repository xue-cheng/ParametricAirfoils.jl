include("utils.jl")

struct CST{T<:AbstractFloat, N} <: AbstractAirfoil{T}
    A₀::T
    Aᵤ::Vector{T}
    Aₗ::Vector{T}
    dz::T
    C₁::T
    C₂::T
end

function CST(A::AbstractVector{T}; dz::T=0, c1::T=0.5, c2::T=1.0) where {T<:AbstractFloat}
    n, r = divrem(length(A), 2)
    if (r!=1)
        throw(ArgumentError("length of parameters must be `2N-1`"))
    end
    CST(A[1], A[2:n+1], A[n+2:end]; dz=dz, c1=c1, c2=c2)
end

function CST(A₀::T, Aᵤ::AbstractVector{T}, Aₗ::AbstractVector{T};
             dz::T=0, c1::T=0.5, c2::T=1.0) where {T<:AbstractFloat}
    if length(Aₗ) != length(Aᵤ)
        throw(ArgumentError("length(Aᵤ) != length(Aₗ)"))
    end
    n = length(Aᵤ)
    CST{T,n}(A₀, copy(Aᵤ), copy(Aₗ), dz, c1, c2)
end

function CST(n::Int,
             xᵤ::AbstractVector{T}, yᵤ::AbstractVector{T},
             xₗ::AbstractVector{T}, yₗ::AbstractVector{T},
             opt::OptimizerFactory; 
             c1::T=0.5, c2::T=1.0) where {T}
    if length(xᵤ)!=length(yᵤ)
        throw(ArgumentError("imbalanced upper-surface data: length(xᵤ)!=length(yᵤ)"))
    end
    if length(xₗ)!=length(yₗ)
        throw(ArgumentError("imbalanced lower-surface data: length(xₗ)!=length(yₗ)"))
    end
    if !(xᵤ[1] == xₗ[1] == 0)
        throw(ArgumentError("`xᵤ[1]` and `xₗ[1]` must be 0"))
    end
    if !(xᵤ[end] == xₗ[end] == 1)
        throw(ArgumentError("`xᵤ[end]` and `xₗ[end]` must be 0"))
    end
    if !isapprox(yᵤ[1], yₗ[1])
        throw(ArgumentError("`yᵤ[1]` must approximate `yₗ[1]`"))
    end
    # move LE to (0,0)
    oy = (yᵤ[1]+yₗ[1] )/2
    yₗ .-= oy
    yᵤ .-= oy
    #
    dz = (yᵤ[end] - yₗ[end])/2
    if dz < -eps(T)
        throw(ArgumentError("`yᵤ[end]` must be greater than `yₗ[end]`"))
    elseif dz < eps(T)
        dz = zero(T) # set to zero
    end
    A₀, Aᵤ, Aₗ, status = fit_cst(n, dz, xᵤ, yᵤ, xₗ, yₗ, opt; c1=c1, c2=c2)
    CST{T,n}(A₀, Aᵤ, Aₗ, dz, c1, c2)
end

function y_upper(cst::CST{T,N}, x::T)::T where{T,N}
    y = cst.A₀*cst_S(N,0,x)
    for i = 1:N
        y += cst.Aᵤ[i]*cst_S(N,i,x)
    end
    y*cst_C(cst.C₁, cst.C₂, x) + cst.dz*x
end

function y_lower(cst::CST{T,N}, x::T)::T where{T,N}
    y = -cst.A₀*cst_S(N,0,x)
    for i = 1:N
        y += cst.Aₗ[i]*cst_S(N,i,x)
    end
    y*cst_C(cst.C₁, cst.C₂, x)-cst.dz*x
end

function dy_upper(cst::CST{T,N}, x::T)::T where{T,N}
    dy = cst.A₀*cst_dS(N,0,x)
    y = cst.A₀*cst_S(N,0,x)
    for i = 1:N
        y += cst.Aᵤ[i]*cst_S(N,i,x)
        dy += cst.Aᵤ[i]*cst_dS(N,i,x)
    end
    y*cst_dC(cst.C₁, cst.C₂, x)+dy*cst_C(cst.C₁, cst.C₂,x)+cst.dz
end

function dy_lower(cst::CST{T,N}, x::T)::T where{T,N}
    dy = -cst.A₀*cst_dS(N,0,x)
    y = -cst.A₀*cst_S(N,0,x)
    for i = 1:N
        dy += cst.Aₗ[i]*cst_dS(N,i,x)
        y += cst.Aₗ[i]*cst_S(N,i,x)
    end
    y*cst_dC(cst.C₁, cst.C₂, x)+dy*cst_C(cst.C₁, cst.C₂,x)-cst.dz
end


function fy_upper(cst::CST{T,N}, x::T)::Tuple{T,T} where{T,N}
    dy = cst.A₀*cst_dS(N,0,x)
    y = cst.A₀*cst_S(N,0,x)
    for i = 1:N
        y += cst.Aᵤ[i]*cst_S(N,i,x)
        dy += cst.Aᵤ[i]*cst_dS(N,i,x)
    end
    C = cst_C(cst.C₁, cst.C₂, x)
    dC = cst_dC(cst.C₁, cst.C₂, x)
    y*C+cst.dz*x, y*dC+dy*C+cst.dz
end

function fy_lower(cst::CST{T,N}, x::T)::Tuple{T,T} where{T,N}
    dy = -cst.A₀*cst_dS(N,0,x)
    y = -cst.A₀*cst_S(N,0,x)
    for i = 1:N
        dy += cst.Aₗ[i]*cst_dS(N,i,x)
        y += cst.Aₗ[i]*cst_S(N,i,x)
    end
    C = cst_C(cst.C₁, cst.C₂, x)
    dC = cst_dC(cst.C₁, cst.C₂, x)
    y*C-cst.dz*x, y*dC+dy*C-cst.dz
end
