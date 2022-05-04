struct CST{N,T,C,L} <: Parameterization{N,T}
    curve::CurveCST{C,T}
    function CST(n::Int, ::Type{T}=Float64;N1::Real=0.5, N2::Real=1, link_le::Bool=true) where {T}
        curve = CurveCST(n, N1, N2, T)
        N = link_le ? 2n-1 : 2n
        L = link_le
        return new{N,T,n,L}(curve)
    end
end

linked_le(::CST{N,T,C,L}) where{N,T,C,L} = L

num_hyper_param(c::CST) = num_hyper_param(c.curve)
get_hyper_param(c::CST) = get_hyper_param(c)
set_hyper_param!(c::CST, hyp) = set_hyper_param!(c, hyp)
set_hyper_param!(c::CST, i, val) = set_hyper_param!(c, i, val)

function _filt_param(cst::CST{N,T,C}, side::Symbol, param::AbstractVector) where {N,T,C} 
    return if side == :U
        view(param, 1:C)
    elseif side == :L
        if linked_le(cst)
            vcat(view(param, C+1:N), -param[C])
        else
            view(param, C+1:N)
        end
    else
        error("Invalid `side=:$side`, must be `:U` or `:L`")
    end
end


airfoil_coord(c::CST, side::Symbol, param, x) = curve_eval(c.curve, _filt_param(c, side, param), x)

airfoil_grad(c::CST, side::Symbol, param, x) = curve_grad(c.curve, _filt_param(c, side, param), x)

airfoil_full(c::CST, side::Symbol, param, x) = curve_full(c.curve, _filt_param(c, side, param), x)

function airfoil_fit(cst::CST{N,T,C},x::AbstractVector, y::AbstractVector) where {N,T,C}
    xu,yu,xl,yl = _split_points(eltype(cst), x, y)
    if linked_le(cst)
        Au = matrix_LSQ(cst.curve, xu)
        Al = matrix_LSQ(cst.curve, xl)
        A = [Au zeros(T, size(Au, 1), N-C);
             zeros(T, size(Al, 1), C-1) -Al[:, C] Al[:, 1:(C-1)]]
        b = vcat(yu, yl)
        @assert N == 2C-1
        p = A\b
        e = sqrt(maximum(abs2, A*p - b))
    else
        Au = matrix_LSQ(cst.curve, xu)
        Al = matrix_LSQ(cst.curve, xl)
        pu = Au\yu
        pl = Al\yl
        eu = sqrt(maximum(abs2, Au*pu - yu))
        el = sqrt(maximum(abs2, Al*pl - yu))
        p = vcat(pu, pl)
        e = max(eu, el)
    end
    return p, e
end
