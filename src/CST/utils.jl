cstShapeF(n, i, x) = binomial(n, i) * x^i * (1 - x)^(n - i)

cstClassF(N1, N2, x) = x^N1 * (1 - x)^N2

function cstShapeD(n, i, x)
    if i == 0 # (1-x)^n
        -n * (1 - x)^(n - 1)
    elseif i == n # x^n
        n * x^(n - 1)
    else
        binomial(n, i) * (
            i * x^(i - 1) * (1 - x)^(n - i) - 
            (n - i) * x^i * (1 - x)^(n - i - 1)
        )
    end
end

function cstClassD(N1, N2, x) 
    N1 * x^(N1 - 1) * (1 - x)^N2 - N2 * x^N1 * (1 - x)^(N2 - 1)
end

function cstMatrix(xu::V, yu::V, xl::V, yl::V, N, N1, N2)　where {T,V <: AbstractVector{T}}
    NP = length(xu) + length(xl) - 2
    A = zeros(T, NP, 2N + 3)
    b = vcat(yu[2:end], yl[2:end])
    ip = 1
    
    @inbounds for x in @view xu[2:end]
        c = cstClassF(N1, N2, x)
        for i in 1:N
            A[ip, i] = c * cstShapeF(N, i, x)
        end
        A[ip, 2N + 1] = c * cstShapeF(N, 0, x) # a_le
        A[ip, 2N + 2] = x # y_tu
        ip += 1
    end
    @inbounds for x in @view xl[2:end]
        c = cstClassF(N1, N2, x)
        for i in 1:N
            A[ip, i + N] = c * cstShapeF(N, i, x)
        end
        A[ip, 2N + 1] = -c * cstShapeF(N, 0, x) # -a_le
        A[ip, 2N + 3] = x # y_tl
        ip += 1
    end
    return A, b
end