function _sqsum(a)
    m, n = size(a)
    ans = Array(Float64, n, 1)
    for i = 1:n
        ans[i] = sum(abs2(a[:, i]))
    end
    return(ans)
end

# The function below is adapted from Dahua Lin's blog post here: http://julialang.org/blog/2013/09/fast-numeric/
function vecdist(a1::Array{Float64, 2}, a2::Array{Float64, 2})
    m, n = size(a1)
    sa = _sqsum(a1)
    sb = _sqsum(a2)
    r = sa .+ reshape(sb, 1, n)
    
    #Update C as alpha*A*B + beta*C
    # or the other three variants according to tA (transpose A) and tB. Returns the updated C.
    gemm!('T', 'N', -2.0, a1, a2, 1.0, r)
    for i = 1:length(r)
        if r[i] < 0
            if r[i] > -1e-9
                r[i] = 0
            else
                error("Negative squared distance in vecdist")
            end
        else
            r[i] = sqrt(r[i])
        end
    end
    return(r)
end


### functions for attractor reconstruction and distance calculation
function calc_dists(shadowmat::Array{Float64,2})
    return vecdist(shadowmat, shadowmat)
end


function precalc_manif_dists(Evals::AbstractVector{Int}, tau_vals::AbstractVector{Int}, vector::AbstractVector)
    shadowmats = (Int=>Dict)[tt => (Int=>Array{Float64,2})[E => construct_shadow(vector, E, tt) for E in Evals] for tt in tau_vals]
    distmats   = (Int=>Dict)[tt => (Int=>Array{Float64,2})[E => calc_dists(shadowmats[tt][E]) for E in Evals] for tt in tau_vals]
    return shadowmats, distmats
end


function construct_shadow(vector::AbstractVector, E::Int, tau_s::Int=1)
    if tau_s < 1
        throw(ArgumentError("tau_s must be greater than 0!"))
    end
    n::Int   = length(vector)
    lag::Int = 0
    shadowmat  = fill(NaN, E, n)
    
    for ii in 1:E
        shadowmat[ii, (lag+1):n] = vector[1:(n-lag)]
        lag += tau_s
    end
    return shadowmat
end
#### End attractor reconstruction and distance calculation functions
