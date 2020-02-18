using Polynomials, Combinatorics, LinearAlgebra, Random, Statistics

"""
The moment (Prony) mapping
"""
𝔐(x,c) = k->sum([c[i]*x[i]^k for i=1:length(x)])

"""
The Fourier mapping
"""
𝔉(x,c) = k->sum([c[i]*exp(im*x[i]*k) for i=1:length(x)])

"""
Given measurement function and an approximation to the nodes,
    solve the linear least squares system.
"""
function solveLinear(f,Ω,x̃)
    V = [exp(im*x*k) for k=0:Ω,x=x̃]
    d = [f(k) for k=0:Ω]
    V\d
end

"""
The inverse moment mapping
"""
function 𝔓(mvec)
    d2 = length(mvec)
    if (isodd(d2))
        error("Length of moment vector should be even")
    end

    d = d2>>1

    # construct the Hankel
    H = [mvec[i+j+1] for i=0:d-1,j=0:d-1]
    rhs = -mvec[d+1:end]
    q = H\rhs
    p = Polynomials.Poly(vcat(q,1))
    xrec=roots(p)
    V = [x^k for k=0:d-1,x=xrec]
    crec = V\mvec[1:d]
    xrec,crec
end

"""
Worst case perturbation function
"""
𝔉ᵪ(x,c,ϵ,s::Array) = begin
    R = length(x)
    x₀ = x[s]
    c₀ = c[s]

    # move to origin
    μ = (minimum(x₀)+maximum(x₀))/2
    y₀ = x₀.-μ

    m₀=[𝔐(y₀,c₀)(j+0.) for j=0:(2*length(s)-1)]
    #m₁ = copy(m₀)
    m₀[end] = m₀[end]+ϵ
    y₁,c₁ = 𝔓(m₀)

    # move back
    x₁ = y₁ .+ μ

    xₙ = copy(vec(x))
    xₙ[s]=real.(x₁)

    cₙ = convert(typeof(c₁),copy(vec(c)))
    cₙ[s] = c₁
    𝔉(xₙ,cₙ)
end

𝔉ᵪ(x,c,ϵ,ℓ::Integer) = 𝔉ᵪ(x,c,ϵ,collect(1:ℓ))
𝔉ᵪ(x,c,ϵ) = 𝔉ᵪ(x,c,ϵ,length(x))
𝔉ᵪ(ℓ) = (x,c,ϵ)->𝔉ᵪ(x,c,ϵ,ℓ)

"""
Random perturbation function
"""
𝔉ᵣ(x,c,ϵ) = k->(𝔉(x,c)(k) + ϵ*((-1+2*rand(MersenneTwister(Int(floor(k)))))+im*(-1+2*rand(MersenneTwister(Int(floor(k)))))))

𝔉ᵣ(Ωₘ::Integer,seed::Integer) = begin
    rng = MersenneTwister(seed)
    seq = -1 .+ rand(rng,Ωₘ+1)

    (x,c,ϵ)->(k->(𝔉(x,c)(k)+ϵ*seq[Int(floor(k+1))]))
end

Hank(F,L,M) = [F(i+j) for i=1:L,j=1:M]

"""
The matrix pencil algorithm for recovery of complex exponentials from equispaced samples

F - function handle for sampling from the model
N - number of equispaced samples
n - the model order
L - the pencil parameter L (default = ceil(N/2))

Returns: array of n estimated frequencies, normalized to [0,2π]
"""
function MatrixPencil(F,N,n,L)
    # check we have enough samples
    if(N < 2*n)
        error("N should be at least 2n")
    end

    if (L < n+1)
        error("L should be at least n+1")
    end


    # construct the matrices H₁ and H₂
    M = N-L
    H = Hank(F,L,M)

    Rₗ = H[2:end,:] # lower part
    Rᵤ = H[1:end-1,:] # upper part

    # compute truncated SVDs
    S₁ = svd(Rᵤ;full=true)
    u₁ = S₁.U[:,1:n]
    Σ₁ = Diagonal(S₁.S[1:n])
    v₁ = S₁.V[:,1:n]

    S₂ = svd(Rₗ;full=true)
    u₂ = S₂.U[:,1:n]
    Σ₂ = Diagonal(S₂.S[1:n])
    v₂ = S₂.V[:,1:n]

    # compute the generalized eigenvalues of the corresponding *regular* pencil
    A = (u₂)'*u₁*Σ₁*(v₁)'*v₂
    B = Σ₂
    z = eigvals(B\A)
    mod.(2π.-angle.(z),2π)

end

MatrixPencil(F,N,n) = MatrixPencil(F,N,n,ceil(N/2))

"""
Amplitude function: all 1's
"""
all1s(ℓ,R) = ones(R)

"""
Amplitude function: i^j
"""
imag1s(ℓ,R) = [(im)^j for j=0:R-1]

"""
Amplitude function: plus minus 1
"""
opposite1s(ℓ,R) = [(-1)^j for j=0:R-1]

"""
Amplitude function: complex amplitudes
"""
mixed1s(ℓ,R) = [1+im*(-1)^j for j=0:R-1]

"""
Creates a clustered node configuration with a single cluster
and other nodes well-separated
(equidistributed in the remaining part of the circle)
"""
Ζ(ℓ,R,Δ,endpoint=π) = begin

    if (ℓ>R)
        error("ℓ should be at most R")
    end

    ρ = (ℓ-1)*Δ
    if (endpoint < ρ)
        error("Δ is too large")
    end

    δ = (endpoint-ρ)/(R-ℓ+1)

    r₁ = collect(range(0,step=Δ,length=ℓ))
    r₂ = collect(range(ρ+δ,step=δ,length=R-ℓ))
    [r₁;r₂].-π/2
end

distCirc(x₀::Number,x₁::Number) = min(mod(x₀-x₁,2π),mod(x₁-x₀,2π))

function housedorff_distance(test_points, candidates)
    n = length(test_points)
    m = length(candidates)
    @assert m==n "Err"
    t = reshape(test_points,(n,1))
    c = reshape(candidates,(n,1))
    return minimum([maximum(abs.([tt[i]-c[i] for i=1:n])) for tt in permutations(t)])
end

"""
Computes the minimal separation between the nodes in an array
"""
function min_sep(X::Array)
    # distance matrix
    D = [distCirc(x,y) for x in X, y in X]
    return minimum(D+diagm(Inf*ones(length(X))))
end

get_separation(X::Array) = [minimum(
            [distCirc(X[j],X[k]) for k in setdiff(1:length(X),j)])
        for j=1:length(X)]

get_separation(X::Array,Y::Array) = [minimum(
            [distCirc(X[j],Y[k]) for k in 1:length(Y)])
        for j=1:length(X)]

get_separation_ext(X::Array,Y::Array) = [findmin(
            [distCirc(X[j],Y[k]) for k in 1:length(Y)])
        for j=1:length(X)]

"""
Test the Matrix Pencil algo: maximal error
"""
testMP(ℓ,R,Δ,Ω,ϵ,clusterFun=Ζ,pertFun=𝔉ᵪ,amplitudeFun = all1s) = begin
    x = clusterFun(ℓ,R,Δ)
    c = amplitudeFun(ℓ,R)
    f = pertFun(x,c,ϵ)
    ϵ₀ = maximum(abs.([𝔉(x,c)(k+0.) - f(k+0.) for k=0:Ω]))
    x̃ = MatrixPencil(f,Ω,R)
    c̃ = solveLinear(f,Ω,x̃)
    Δx = housedorff_distance(exp.(im.*x),exp.(im.*x̃))/ϵ₀
    Δc = housedorff_distance(c,c̃)/ϵ₀


    # need to match the first two
    η = get_separation(x) # distance to closest neighbor
    σ = get_separation(x,x̃) # error for each node

    # success if all errors are smaller than corresponding separations
    succ = reduce((x,y)->x&&y,σ .< η/3)
    [Ω*Δx,Δc,Int(succ),ϵ₀]
end

"""
Test the Matrix Pencil algo: a specific node
"""
testMP_individual(ℓ,R,Δ,Ω,ϵ,node_index,clusterFun=Ζ,pertFun=𝔉ᵪ,amplitudeFun = all1s) = begin
    x = clusterFun(ℓ,R,Δ)
    c = amplitudeFun(ℓ,R)
    f = pertFun(x,c,ϵ)
    ϵ₀ = maximum(abs.([𝔉(x,c)(k+0.) - f(k+0.) for k=0:Ω]))
    x̃ = MatrixPencil(f,Ω,R)
    c̃ = solveLinear(f,Ω,x̃)


    # need to match the first two
    η = get_separation(x) # distance to closest neighbor
    σ̃ = get_separation_ext(x,x̃) # error for each node, including index
    Δx = σ̃[node_index][1]
    idx = σ̃[node_index][2]
    Δc = (abs(c̃[idx] - c[node_index]))

    succ = (σ̃[node_index][1] < η[node_index]/3)

    # success if all errors are smaller than corresponding separations
    [Ω*Δx/ϵ₀,Δc/ϵ₀,Int(succ),ϵ₀]
end

"""
Test the Matrix Pencil algorithm: accuracy for every node
"""
testMP_all(ℓ,R,Δ,Ω,ϵ,clusterFun=Ζ,pertFun=𝔉ᵪ,amplitudeFun = all1s) = begin
    x = clusterFun(ℓ,R,Δ)
    c = amplitudeFun(ℓ,R)
    f = pertFun(x,c,ϵ)
    ϵ₀ = maximum(abs.([𝔉(x,c)(k+0.) - f(k+0.) for k=0:Ω]))
    x̃ = MatrixPencil(f,Ω,R)
    c̃ = solveLinear(f,Ω,x̃)


    # need to match the first two
    η = get_separation(x) # distance to closest neighbor
    σ̃ = get_separation_ext(x,x̃) # error for each node, including index
    Δx = [sep[1] for sep=σ̃]
    Δc = [abs(c̃[σ̃[i][2]]-c[i]) for i=1:R]

    # success if all errors are smaller than corresponding separations
    succ = [Int(σ̃[i][1] < η[i]/3) for i=1:R]

    [Ω*Δx/ϵ₀,Δc/ϵ₀,succ,ϵ₀]
end

NrangeLog(low,high,n=100) = Int.(unique(floor.(10.0.^LinRange(low,high,n))))

"""
Run a single random experiment testing accuracy
"""
runSingle_Acc(ℓ,R,clusterFun, pertFun, amplitudeFun, ϵR,ΩR,ΔR) = begin
    ϵ = ϵR[rand(1:end)]
    Ω = ΩR[rand(1:end)]
    Δ = ΔR[rand(1:end)]

    res = testMP_all(ℓ,R,Δ,Ω,ϵ,clusterFun,pertFun,amplitudeFun)
    [inv(Ω*Δ),res...]
end

"""
Run a single random experiment testing accuracy
"""
runSingle(ℓ,R,clusterFun, pertFun, amplitudeFun, ϵR,ΩR,ΔR) = begin
    ϵ = ϵR[rand(1:end)]
    Ω = ΩR[rand(1:end)]
    Δ = ΔR[rand(1:end)]

    res = testMP(ℓ,R,Δ,Ω,ϵ,clusterFun,pertFun,amplitudeFun)

    [inv(Ω*Δ),res[3:4]...]
end

runSingle_ind(ℓ,R,clusterFun, pertFun, amplitudeFun, ni, ϵR,ΩR,ΔR) = begin
    ϵ = ϵR[rand(1:end)]
    Ω = ΩR[rand(1:end)]
    Δ = ΔR[rand(1:end)]

    res = testMP_individual(ℓ,R,Δ,Ω,ϵ,ni,clusterFun,pertFun,amplitudeFun)

    [inv(Ω*Δ),res[3:4]...]
end

cc(x) = (x==0) ? (:red) : (:blue)
ss(x) = (x==0) ? (:circle) : (:dtriangle)