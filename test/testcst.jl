
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
  
  dynaca = dy_upper(cst, xtest)
  dycst = dy_upper(cst, xtest)
  @test all(isapprox.(dynaca, dycst, rtol=.01, atol=.01))
  dynaca = dy_lower(cst, xtest)
  dycst = dy_lower(cst, xtest)
  @test all(isapprox.(dynaca, dycst, rtol=.01, atol=.01))
  
  # fit cambered airfoil
  cst = CST(undef, Float64, 5)
  a = get_param(cst)
  set_param!(cst, a)
  naca = NACA"4412" 
  x0, y0 = gen_airfoil(naca, NP)
  fit!(cst, x0, y0)
  xx, yy = gen_airfoil(cst, NP)
  @test xx[end] == xx[1] == 1
  @test issorted(xx[NP:-1:1]) && issorted(xx[NP:end])
  @test y_upper(cst, 1.0) > y_lower(cst, 1.0)
  modify_te!(cst, T_TE=0) # modify to sharp TE
  @test isapprox(y_upper(cst, 1.0), y_lower(cst, 1.0), atol=1e-6)
  end