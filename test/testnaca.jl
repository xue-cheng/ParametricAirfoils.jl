@testset "NACA 4" begin
function read_airfoil_data(fn::AbstractString)
    xs = Float64[]
    ys = Float64[]
    open(fn, "r") do io
        title = readline(io)
        while true
            line = strip(readline(io))
            if isempty(line)
                break
            end
            x, y = parse.(Float64, split(line))
            push!(xs, x)
            push!(ys, y)
        end
    end
    if (iseven(length(xs)))
        error("File Broken")
    end
    n = (length(xs) + 1) >> 1
    xu = xs[n:-1:1]
    yu = ys[n:-1:1]
    xl = xs[n:end]
    yl = ys[n:end]
    xu, yu, xl, yl
end
xu, yu, xl, yl = read_airfoil_data("NACA0015.dat")
naca0015 = NACA"0015"
yuaf = y_upper(naca0015, xu)
@test all(isapprox.(yu, yuaf, atol=0.001))
cst0015 = CST(9, xu, yu, xl, yl)
yuaf = y_upper(cst0015, xu)
@test all(isapprox.(yu, yuaf, atol=0.001))
end