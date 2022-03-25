
@testset "CST($NCST+$NCST)" begin
    cst = CST(NCST, Float64)
    for airfoil in (NACA0012, RAE2822)
        paf, err = fit(cst, airfoil.x, airfoil.y)
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