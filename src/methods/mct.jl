struct MCT{N,T,A,B,C<:CurveFunction{A,T},D<:CurveFunction{B,T}} <: Parameterization{N,T}
    camber::C
    thickness::D
    function MCT(camber::C, thickness::D) where {T, A, B, C<:CurveFunction{A,T}, D<:CurveFunction{B,T}}
        N = A+B
        new{N,T,A,B,C,D}(camber, thickness)
    end
end

num_hyper_param(m::MCT) = num_hyper_param(m.camber)+num_hyper_param(m.thickness)
get_hyper_param(m::MCT) = vcat(get_hyper_param(m.camber), get_hyper_param(m.thickness))
function set_hyper_param!(m::MCT, hyp)
    nc = num_hyper_param(m.camber)
    nt = num_hyper_param(m.thickness)
    nc+nt == length(hyp) || error("length of `hyp` does not match requirement")
    set_hyper_param!(m.camber, view(hyp, 1:nc)) 
    set_hyper_param!(m.thickness, view(hyp, nc+1:nc+nt))
    return nothing
end

function set_hyper_param!(m::MCT, i::Int, val::Real)
    nc = num_hyper_param(m.camber)
    nt = num_hyper_param(m.thickness)
    if 0<i<=nc
        set_hyper_param!(m.camber, i, val)
    elseif nc<i<=nc+nt
        set_hyper_param!(m.thickness, i-nc, val)
    else
        error("invalid `i`")
    end
    return nothing
end

function side_sign(side::Symbol, ::Type{T}) where {T}
    if side == :U
        one(T)
    elseif side == :L
        -one(T)
    else
        error("Invalid `side=:$side`, must be `:U` or `:L`")
    end
end

for op in (:eval, :grad, :full)
    @eval function $(Symbol(:_mct_c_, op))(m::MCT{N,T,A,B}, param, x) where {N,T,A,B}
        return $(Symbol(:curve_, op))(m.camber, param[1:A], x)
    end
    @eval function $(Symbol(:_mct_t_, op))(m::MCT{N,T,A,B}, param, x) where {N,T,A,B}
        return $(Symbol(:curve_, op))(m.thickness,param[A+1:end], x)
    end
end

function airfoil_coord(m::MCT, side::Symbol, param, x)
    c = side_sign(side, eltype(m))
    yc = _mct_c_eval(m, param, x)
    yt = _mct_t_eval(m, param, x)
    return yc + c*yt
end

function airfoil_grad(m::MCT, side::Symbol, param, x)
    c = side_sign(side, eltype(m))
    yc = _mct_c_grad(m, param, x)
    yt = _mct_t_grad(m, param, x)
    return yc + c*yt
end

function airfoil_full(m::MCT, side::Symbol, param, x)
    c = side_sign(side, eltype(m))
    yc, dyc = _mct_c_full(m, param, x)
    yt, dyt = _mct_t_full(m, param, x)
    return yc + c*yt, dyc + c*dyt
end


function airfoil_fit(m::MCT{N, T}, x::AbstractVector, y::AbstractVector) where {N,T}
    xu,yu,xl,yl = _split_points(eltype(m), x, y)
    Acu = matrix_LSQ(m.camber, xu)
    Atu = side_sign(:U, T) * matrix_LSQ(m.thickness, xu)
    Acl = matrix_LSQ(m.camber, xl)
    Atl = side_sign(:L, T) * matrix_LSQ(m.thickness, xl)
    A = [Acu Atu; Acl Atl]
    b = vcat(yu, yl)
    p = A \ b
    yy = vcat(airfoil_coord(m, :U, p, xu), airfoil_coord(m, :L, p, xl))
    e = sqrt(maximum(abs2, yy - b))
    return p, e
end