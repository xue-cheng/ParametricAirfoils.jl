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
end

@testset "CST" begin
naca = NACA"0015"
NP = 65
x0, y0 = gen_airfoil(naca, NP)
# original airfoil
cst = fit(:CST, x0, y0, N=6)
xtest = x0[NP:end]
yuaf = y_upper(cst, xtest)
ylaf = y_lower(cst, xtest)
@test all(isapprox.(y0[NP:-1:1], yuaf, atol=2e-4))
@test all(isapprox.(y0[NP:end], ylaf, atol=2e-4))
# fit cambered airfoil
naca = NACA"4412"s
x0, y0 = gen_airfoil(naca, NP)
fit!(cst, x0, y0)
xx, yy = gen_airfoil(cst, NP)
@test xx[end] == xx[1] == 1
@test issorted(xx[NP:-1:1]) && issorted(xx[NP:end])
end