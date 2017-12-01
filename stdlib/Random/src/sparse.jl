# This file is a part of Julia. License is MIT: https://julialang.org/license


## matrices

function sprand_IJ(r::AbstractRNG, m::Integer, n::Integer, density::AbstractFloat)
    ((m < 0) || (n < 0)) && throw(ArgumentError("invalid Array dimensions"))
    0 <= density <= 1 || throw(ArgumentError("$density not in [0,1]"))
    N = n*m

    I, J = Vector{Int}(), Vector{Int}() # indices of nonzero elements
    sizehint!(I, round(Int,N*density))
    sizehint!(J, round(Int,N*density))

    # density of nonzero columns:
    L = log1p(-density)
    coldensity = -expm1(m*L) # = 1 - (1-density)^m
    colsparsity = exp(m*L) # = 1 - coldensity
    iL = 1/L

    rows = Vector{Int}()
    for j in randsubseq(r, 1:n, coldensity)
        # To get the right statistics, we *must* have a nonempty column j
        # even if p*m << 1.   To do this, we use an approach similar to
        # the one in randsubseq to compute the expected first nonzero row k,
        # except given that at least one is nonzero (via Bayes' rule);
        # carefully rearranged to avoid excessive roundoff errors.
        k = ceil(log(colsparsity + rand(r)*coldensity) * iL)
        ik = k < 1 ? 1 : k > m ? m : Int(k) # roundoff-error/underflow paranoia
        randsubseq!(r, rows, 1:m-ik, density)
        push!(rows, m-ik+1)
        append!(I, rows)
        nrows = length(rows)
        Jlen = length(J)
        resize!(J, Jlen+nrows)
        @inbounds for i = Jlen+1:length(J)
            J[i] = j
        end
    end
    I, J
end


"""
    sprand([rng],[type],m,[n],p::AbstractFloat,[rfn])

Create a random length `m` sparse vector or `m` by `n` sparse matrix, in
which the probability of any element being nonzero is independently given by
`p` (and hence the mean density of nonzeros is also exactly `p`). Nonzero
values are sampled from the distribution specified by `rfn` and have the type `type`. The uniform
distribution is used in case `rfn` is not specified. The optional `rng`
argument specifies a random number generator, see [Random Numbers](@ref).

# Examples
```jldoctest
julia> rng = MersenneTwister(1234);

julia> sprand(rng, Bool, 2, 2, 0.5)
2×2 SparseMatrixCSC{Bool,Int64} with 2 stored entries:
  [1, 1]  =  true
  [2, 1]  =  true

julia> sprand(rng, Float64, 3, 0.75)
3-element SparseVector{Float64,Int64} with 1 stored entry:
  [3]  =  0.298614
```
"""
function sprand(r::AbstractRNG, m::Integer, n::Integer, density::AbstractFloat,
                rfn::Function, ::Type{T}=eltype(rfn(r,1))) where T
    N = m*n
    N == 0 && return spzeros(T,m,n)
    N == 1 && return rand(r) <= density ? sparse([1], [1], rfn(r,1)) : spzeros(T,1,1)

    I,J = sprand_IJ(r, m, n, density)
    Base.SparseArrays.sparse_IJ_sorted!(I, J, rfn(r,length(I)), m, n, +)  # it will never need to combine
end

function sprand(m::Integer, n::Integer, density::AbstractFloat,
                rfn::Function, ::Type{T}=eltype(rfn(1))) where T
    N = m*n
    N == 0 && return spzeros(T,m,n)
    N == 1 && return rand() <= density ? sparse([1], [1], rfn(1)) : spzeros(T,1,1)

    I,J = sprand_IJ(defaultRNG(), m, n, density)
    Base.SparseArrays.sparse_IJ_sorted!(I, J, rfn(length(I)), m, n, +)  # it will never need to combine
end

truebools(r::AbstractRNG, n::Integer) = ones(Bool, n)

sprand(m::Integer, n::Integer, density::AbstractFloat) = sprand(defaultRNG(),m,n,density)

sprand(r::AbstractRNG, m::Integer, n::Integer, density::AbstractFloat) = sprand(r,m,n,density,rand,Float64)
sprand(r::AbstractRNG, ::Type{T}, m::Integer, n::Integer, density::AbstractFloat) where {T} = sprand(r,m,n,density,(r, i) -> rand(r, T, i), T)
sprand(r::AbstractRNG, ::Type{Bool}, m::Integer, n::Integer, density::AbstractFloat) = sprand(r,m,n,density, truebools, Bool)
sprand(::Type{T}, m::Integer, n::Integer, density::AbstractFloat) where {T} = sprand(defaultRNG(), T, m, n, density)

"""
    sprandn([rng], m[,n],p::AbstractFloat)

Create a random sparse vector of length `m` or sparse matrix of size `m` by `n`
with the specified (independent) probability `p` of any entry being nonzero,
where nonzero values are sampled from the normal distribution. The optional `rng`
argument specifies a random number generator, see [Random Numbers](@ref).

# Examples
```jldoctest
julia> rng = MersenneTwister(1234);

julia> sprandn(rng, 2, 2, 0.75)
2×2 SparseMatrixCSC{Float64,Int64} with 3 stored entries:
  [1, 1]  =  0.532813
  [2, 1]  =  -0.271735
  [2, 2]  =  0.502334
```
"""
sprandn(r::AbstractRNG, m::Integer, n::Integer, density::AbstractFloat) = sprand(r,m,n,density,randn,Float64)
sprandn(m::Integer, n::Integer, density::AbstractFloat) = sprandn(defaultRNG(),m,n,density)


## vectors


### Rand Construction
sprand(n::Integer, p::AbstractFloat, rfn::Function, ::Type{T}) where {T} = sprand(defaultRNG(), n, p, rfn, T)
function sprand(r::AbstractRNG, n::Integer, p::AbstractFloat, rfn::Function, ::Type{T}) where T
    I = randsubseq(r, 1:convert(Int, n), p)
    V = rfn(r, T, length(I))
    SparseVector(n, I, V)
end

sprand(n::Integer, p::AbstractFloat, rfn::Function) = sprand(defaultRNG(), n, p, rfn)
function sprand(r::AbstractRNG, n::Integer, p::AbstractFloat, rfn::Function)
    I = randsubseq(r, 1:convert(Int, n), p)
    V = rfn(r, length(I))
    SparseVector(n, I, V)
end

sprand(n::Integer, p::AbstractFloat) = sprand(defaultRNG(), n, p, rand)

sprand(r::AbstractRNG, n::Integer, p::AbstractFloat) = sprand(r, n, p, rand)
sprand(r::AbstractRNG, ::Type{T}, n::Integer, p::AbstractFloat) where {T} = sprand(r, n, p, (r, i) -> rand(r, T, i))
sprand(r::AbstractRNG, ::Type{Bool}, n::Integer, p::AbstractFloat) = sprand(r, n, p, truebools)
sprand(::Type{T}, n::Integer, p::AbstractFloat) where {T} = sprand(defaultRNG(), T, n, p)

sprandn(n::Integer, p::AbstractFloat) = sprand(defaultRNG(), n, p, randn)
sprandn(r::AbstractRNG, n::Integer, p::AbstractFloat) = sprand(r, n, p, randn)


## _rand_pm1! (used in Base.LinAlg)

function Base._rand_pm1!(v)
    for i in eachindex(v)
        v[i] = rand(Bool) ? 1 : -1
    end
end
