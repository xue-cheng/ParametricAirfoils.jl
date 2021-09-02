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
ylaf = y_lower(naca0015, xl)
@test all(isapprox.(yu, yuaf, atol=1e-4))
@test all(isapprox.(yl, ylaf, atol=1e-4))
naca4312s = NACA"4312"s
xu, yu, xl, yl = read_airfoil_data("NACA4312S.dat")
xtest = collect(0:0.1:1)
xuaf = x_upper(naca4312s, xtest)
yuaf = y_upper(naca4312s, xtest)
xlaf = x_lower(naca4312s, xtest)
ylaf = y_lower(naca4312s, xtest)
@test all(isapprox.(xu, xuaf, atol=1e-4))
@test all(isapprox.(yu, yuaf, atol=1e-4))
@test all(isapprox.(xl, xlaf, atol=1e-4))
@test all(isapprox.(yl, ylaf, atol=1e-4))

end
