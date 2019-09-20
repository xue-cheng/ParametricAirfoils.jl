abstract type NACABase{T<:AbstractFloat} <: AbstractAirfoil{T} end

include("NACA4.jl")

function NACA(::Type{T}, num::AbstractString, sharp_TE::Bool=false) where {T<:AbstractFloat}
    if length(num) == 4
        NACA4(T, num, sharp_TE)
    else
        throw(ArgumentError("unsupported airfoil NACA$num"))
    end
end

NACA(num, sharp_TE) = NACA(Float64, num, sharp_TE)

macro NACA_str(num)
    NACA(num, false)
end

macro NACA_str(num, s)
    if (s != "s")
        error("use NACA\"XXXX\"s for sharp trailing edge")
    end
    NACA(num, true)
end
