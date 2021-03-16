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

using Printf

@testset "CST" begin
naca = NACA"0015"
xu, yu, xl, yl = gen_airfoil(naca)
dyu = dy_upper(naca, xu)
dyl = dy_lower(naca, xl)
# original airfoil
cst = fit(:CST, xu, yu, xl, yl; N=9)
yuaf = y_upper(cst, xu)
ylaf = y_lower(cst, xl)
@test all(isapprox.(yu, yuaf, atol=1e-4))
@test all(isapprox.(yl, ylaf, atol=1e-4))
dyuaf = dy_upper(cst, xu)
dylaf = dy_lower(cst, xl)
@test all(isapprox.(dyu, dyuaf, rtol=1e-2, atol=1e-2))
@test all(isapprox.(dyl, dylaf, rtol=1e-2, atol=1e-2))

# modify airfoil
@. yu += xu * 0.001
@. yl -= xl * 0.0015
cst = fit(:CST, xu, yu, xl, yl; N=9)
yuaf = y_upper(cst, xu)
ylaf = y_lower(cst, xu)
@test all(isapprox.(yu, yuaf, atol=1e-4))
@test all(isapprox.(yl, ylaf, atol=1e-4))

# shart TE
naca = NACA"0015"s
xu, yu, xl, yl = gen_airfoil(naca)
# original airfoil
cst = fit(:CST, xu, yu, xl, yl; N=9)
yuaf = y_upper(cst, xu)
ylaf = y_lower(cst, xl)
@test all(isapprox.(yu, yuaf, atol=1e-4))
@test all(isapprox.(yl, ylaf, atol=1e-4))

end