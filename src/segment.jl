using Printf

print2d(io::IO,x,y) = @printf(io, "%+.12f\t%+.12f\n"   , x, y)
printxy(io::IO,x,y) = @printf(io, "%+.12f\t%+.12f\t0\n", x, y)
printxz(io::IO,x,y) = @printf(io, "%+.12f\t0\t%+.12f\n", x, y)

function gen_segment(io::IO, a::AbstractAirfoil, np::Int=65; dim::Int=3, plane::Symbol=:XY)
    dim == 2 || dim == 3 || error("`dim` must be either 2 or 3")
    if dim == 3 && plane != :XY && plane != :XZ
        error("`plane` must be either :XY or :XZ")
    end
    prnt = if dim == 2
        print2d
    elseif plane == :XY
        printxy
    elseif plane == :XZ
        printxz
    end
    println(io, np * 2 - 1)
    xu, yu, xl, yl = gen_airfoil(a, np)
    for i in np:-1:2
        prnt(io, xl[i], yl[i])
    end
    for i in 1:np
        prnt(io, xu[i], yu[i])
    end
end

function gen_segment(fn::AbstractString, a::AbstractAirfoil, np::Int=65; dim::Int=3, plane::Symbol=:XY)
    open(fn, "w") do io
        gen_segment(io, a, np; dim=dim, plane=plane)
    end
end

gen_segment(a::AbstractAirfoil, np::Int=65; dim::Int=3, plane::Symbol=:XY) = gen_segment(stdout, a, np; dim=dim, plane=plane)
