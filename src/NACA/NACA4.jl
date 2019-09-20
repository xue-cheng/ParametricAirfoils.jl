struct NACA4{T<:AbstractFloat, SYM, STE} <: NACABase{T}
    m::T
    p::T
    t::T
end

function NACA4(::Type{T}, num::AbstractString, sharp_TE::Bool) where {T<:AbstractFloat}
    if length(num) != 4
        throw(ArgumentError("invalid 4-dig NACA airfoil: NACA$num"))
    end
    t = parse(T, num[3:4]) / 100
    m = parse(T, num[1:1]) / 100
    p = parse(T, num[2:2]) / 10
    NACA4{T, m==0, sharp_TE}(m,p,t)
end

function y_t(a::NACA4{T, SYM, false}, x::T)::T where {T<:AbstractFloat,SYM}
    5*a.t*(.2969*√x-.1260*x-.3516*x^2+.2843*x^3-.1015*x^4)
end
function y_t(a::NACA4{T, SYM, true}, x::T)::T where {T<:AbstractFloat,SYM}
    5*a.t*(.2969*√x-.1260*x-.3516*x^2+.2843*x^3-.1036*x^4)
end

function dy_t(a::NACA4{T, SYM, false}, x::T)::T where {T<:AbstractFloat,SYM}
    5*a.t*(.2969/(2*√x)-.1260-.3516*2*x+.2843*3*x^2-.1015*4*x^3)
end

function dy_t(a::NACA4{T, SYM, true}, x::T)::T where {T<:AbstractFloat,SYM}
    5*a.t*(.2969/(2*√x)-.1260-.3516*2*x+.2843*3*x^2-.1036*4*x^3)
end

function y_c(a::NACA4{T}, x::T)::T where {T<:AbstractFloat}
    if x < a.p 
        a.m/a.p^2 * (2*a.p*x - x^2)
    else
        a.m/(1-a.p)^2 * (1 - 2*a.p + 2*a.p*x -x^2)
    end
end

function dy_c(a::NACA4{T}, x::T)::T where {T<:AbstractFloat}
    if x < a.p 
        2 * a.m/a.p^2 * (a.p - x)
    else
        2 * a.m/(1-a.p)^2 * (a.p*x -x)
    end
end

function y_upper(a::NACA4{T,true}, x::T)::T where {T<:AbstractFloat}
    y_t(a, x)
end

function y_lower(a::NACA4{T,true}, x::T)::T where {T<:AbstractFloat}
    -y_t(a, x)
end

function dy_upper(a::NACA4{T,true}, x::T)::T where {T<:AbstractFloat}
    dy_t(a, x)
end

function dy_lower(a::NACA4{T,true}, x::T)::T where {T<:AbstractFloat}
    -dy_t(a, x)
end

function fy_upper(a::NACA4{T,true}, x::T)::T where {T<:AbstractFloat}
    y_t(a, x), dy_t(a, x)
end

function fy_lower(a::NACA4{T,true}, x::T)::T where {T<:AbstractFloat}
    -y_t(a, x), -dy_t(a, x)
end


function y_upper(a::NACA4{T,false}, x::T)::T where {T<:AbstractFloat}
    y_c(a,x) + y_t(a, x) * cos(atan(dy_c(a,x)))
end

function y_lower(a::NACA4{T,false}, x::T)::T where {T<:AbstractFloat}
    y_c(a,x) - y_t(a, x) * cos(atan(dy_c(a,x)))
end

function x_upper(a::NACA4{T,false}, x::T)::T where {T<:AbstractFloat}
    x - y_t(a, x) * sin(atan(dy_c(a,x)))
end

function x_lower(a::NACA4{T,false}, x::T)::T where {T<:AbstractFloat}
    x + y_t(a, x) * sin(atan(dy_c(a,x)))
end
