using AirfoilDatabase
using Printf
using Test
using ParametricAirfoils


result = query_airfoil("N0012")
@assert length(result)==1
NACA0012 = result[1]
result = query_airfoil("RAE2822")
@assert length(result)==1
RAE2822 = result[1]

NCST = 10
NPoly = 10
indicies = tuple((max(0.5, i-1) for i in 1:NPoly)...)

function test_method(method)
    for airfoil in (NACA0012, RAE2822)
        paf, err = fit(method, airfoil.x, airfoil.y)
        @info "$(airfoil.name) Fit error= $err" 
        ile = argmin(airfoil.x)
        x_up = airfoil.x[ile:-1:1]
        y_up = airfoil.y[ile:-1:1]
        x_lo = airfoil.x[ile:end]
        y_lo = airfoil.y[ile:end]

        yu = paf(:U, x_up)
        eu = maximum(abs, yu-y_up)
        yl = paf(:L, x_lo)
        el = maximum(abs, yl.-y_lo)
        @info "Maximum error= $eu(Upper) $el(Lower)"
        @test eu < 2e-4 
        @test el < 2e-4
        if airfoil === NACA0012
            y, dy = paf(:U, x_up, true)
            @test isinf(dy[1])
            @test yu == y
            t = 0.12
            x = x_up[2:end]
            dyt = similar(x)
            @. dyt = 5t*(.2969/2*x^-0.5-.126-.3516*2*x+.2843*3*x^2-.1015*4*x^3)
            @test all(isapprox.(dy[2:end], dyt, atol=5e-2, rtol=5e-2))
        end
        gen_airfoil(paf)
    end
end


@testset "CST($NCST+$NCST)" begin
    cst = CST(NCST, Float64)
    test_method(cst)
end



@testset "MCT(CST$NCST+CST$NCST)" begin
    mct = MCT(CurveCST(NCST, 1, 1, Float64), CurveCST(NCST, 0.5, 1.0, Float64))
    test_method(mct)
end


@testset "MCT(CST$NCST+Poly$NPoly)" begin
    mct = MCT(CurveCST(NCST, 1, 1, Float64), CurvePoly(indicies, Float64))
    test_method(mct)

end
