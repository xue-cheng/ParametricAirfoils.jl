function linspace(start::T, stop::T, np::I) where {T<:AbstractFloat, I<:Integer}
    d = (stop-start)/(np-1)
    p = Vector{T}(undef, np)
    p[1] = start
    p[end] = stop
    for i = 2:np-1
        p[i] = start + (i-1)*d
    end
    return p
end

function cosspace(start::T, stop::T, np::I) where {T<:AbstractFloat, I<:Integer}
    p = linspace(zero(T), convert(T,π), np)
    r = (stop-start)/2
    for i = 1:np
        p[i] = start + r*(one(T)-cos(p[i]))
    end
    return p
end

function gen_airfoil(a::AbstractAirfoil{T}, np::Int=65; cos_space::Bool = true) where {T<:AbstractFloat}
    if cos_space
        xc = cosspace(zero(T), one(T), np)
    else
        xc = linspace(zero(T), one(T), np)
    end
    xu = Vector{T}(undef, np)
    yu = Vector{T}(undef, np)
    xl = Vector{T}(undef, np)
    yl = Vector{T}(undef, np)
    @inbounds for i = 1:np
        xu[i] = x_upper(a, xc[i])
        yu[i] = y_upper(a, xc[i])
        xl[i] = x_lower(a, xc[i])
        yl[i] = y_lower(a, xc[i])
    end
    return xu,yu,xl,yl
end

function n_upper(a::AbstractAirfoil{T}, xc::T)::Tuple{T,T} where {T<:AbstractFloat}
    θ = atan(dy_upper(a, xc))
    return -sin(θ), cos(θ)
end

function n_lower(a::AbstractAirfoil{T}, xc::T)::Tuple{T,T} where {T<:AbstractFloat}
    θ = atan(dy_lower(a, xc))
    return sin(θ), -cos(θ)
end

@inline n_upper(a::AbstractAirfoil{T}, xc::R) where {T<:AbstractFloat, R<:Real} = n_upper(a, convert(T, xc))
@inline n_lower(a::AbstractAirfoil{T}, xc::R) where {T<:AbstractFloat, R<:Real} = n_lower(a, convert(T, xc))

function t_upper(a::AbstractAirfoil{T}, xc::T)::Tuple{T,T} where {T<:AbstractFloat}
    dy = dy_upper(a,xc)
    if isinf(dy)
        return zero(T), sign(dy)
    else
        n = sqrt(dy^2+oneunit(T))
        return 1/n, dy/n
    end
end

function t_lower(a::AbstractAirfoil{T}, xc::T)::Tuple{T,T} where {T<:AbstractFloat}
    dy = dy_lower(a,xc)
    if isinf(dy)
        return zero(T), sign(dy)
    else
        n = sqrt(dy^2+oneunit(T))
        return 1/n, dy/n
    end
end

@inline t_upper(a::AbstractAirfoil{T}, xc::R) where {T<:AbstractFloat, R<:Real} = t_upper(a, convert(T, xc))
@inline t_lower(a::AbstractAirfoil{T}, xc::R) where {T<:AbstractFloat, R<:Real} = t_lower(a, convert(T, xc))