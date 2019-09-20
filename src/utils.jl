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

function gen_airfoil(a::AbstractAirfoil{T}, np::Int=65; cos_space::Bool = true) where {T}
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