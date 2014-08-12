using Base.LinAlg.BLAS 
using Distance

function _sqsum(a)
    m, n = size(a)
    ans = Array(Float64, n, 1)
    for i = 1:n
        ans[i] = sum(abs2(a[:, i]))
    end
    return(ans)
end

# The code below is adapted from Dahua Lin's blog post here: http://julialang.org/blog/2013/09/fast-numeric/
function vecdist(a1::Array{Float64, 2}, a2::Array{Float64, 2})
    m, n = size(a1)
    sa = _sqsum(a1)
    sb = _sqsum(a2)
    r = sa .+ reshape(sb, 1, n)
    
    #Update C as alpha*A*B + beta*C 
    # or the other three variants according to tA (transpose A) and tB. Returns the updated C.
    gemm!('T', 'N', -2.0, a1, a2, 1.0, r)
    for i = 1:length(r)
        r[i] = sqrt(r[i])
    end
    return(r)
end


x = [0.; 1.; 2.]
y = [1.; 0.; 5.]
z = [3.5; 2.2; 1.1]

mat = [x  y z]
println(mat)

ans1 = sqrt(pairwise(SqEuclidean(), mat))
println(ans1)

ans2 = vecdist(mat, mat)
println(ans2)

println(sqrt(2))