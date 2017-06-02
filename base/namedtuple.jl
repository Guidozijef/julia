# This file is a part of Julia. License is MIT: https://julialang.org/license

@generated function namedtuple(::Type{NamedTuple{names}}, args...) where names
    N = length(names)
    if length(args) == N
        Expr(:new, :(NamedTuple{names,$(Tuple{args...})}), Any[ :(args[$i]) for i in 1:N ]...)
    else
        :(throw(ArgumentError("wrong number of arguments to named tuple constructor")))
    end
end

NamedTuple() = namedtuple(NamedTuple{()})

length(t::NamedTuple) = nfields(t)
start(t::NamedTuple) = 1
done(t::NamedTuple, iter) = iter > nfields(t)
next(t::NamedTuple, iter) = (getfield(t, iter), iter + 1)
endof(t::NamedTuple) = length(t)
getindex(t::NamedTuple, i::Int) = getfield(t, i)
getindex(t::NamedTuple, i::Symbol) = getfield(t, i)

function getindex(t::NamedTuple, I::AbstractVector)
    names = unique( Symbol[ isa(i,Symbol) ? i : fieldname(typeof(t),i) for i in I ] )
    namedtuple(NamedTuple{(names...)}, [ getfield( t, i ) for i in names ]...)
end

convert(::Type{NamedTuple{names,T}}, itr::NamedTuple{names,T}) where {names,T} = itr

function convert(::Type{NamedTuple{names,T}}, itr) where {names,T}
    namedtuple(NamedTuple{names}, T(itr)...)
end

function show(io::IO, t::NamedTuple)
    n = length(t)
    if n == 0
        print(io, "NamedTuple()")
    else
        print(io, "(")
        for i = 1:n
            print(io, fieldname(typeof(t),i), " = "); show(io, getfield(t,i))
            if n == 1
                print(io, ",")
            elseif i < n
                print(io, ", ")
            end
        end
        print(io, ")")
    end
end

eltype(::Type{NamedTuple{names,T}}) where {names,T} = eltype(T)

@generated function convert(::Type{Tuple}, n::NamedTuple)
    Expr(:tuple, Any[ Expr(:getfield, :n, f) for f = 1:nfields(n) ]...)
end

==(a::NamedTuple{n}, b::NamedTuple{n}) where {n} = Tuple(a) == Tuple(b)
==(a::NamedTuple, b::NamedTuple) = false

isequal(a::NamedTuple{n}, b::NamedTuple{n}) where {n} = isequal(Tuple(a), Tuple(b))
isequal(a::NamedTuple, b::NamedTuple) = false

_nt_names(::NamedTuple{names}) where {names} = names
_nt_names(::Type{T}) where {names,T<:NamedTuple{names}} = names

hash(x::NamedTuple, h::UInt) = xor(object_id(_nt_names(x)), hash(Tuple(x), h))

isless(a::NamedTuple{n}, b::NamedTuple{n}) where {n} = isless(Tuple(a), Tuple(b))
# TODO: case where one argument's names are a prefix of the other's

function map(f, nt::NamedTuple, nts::NamedTuple...)
    # this method makes sure we don't define a map(f) method
    _nt_map(f, nt, nts...)
end

@generated function _nt_map(f, nts::NamedTuple...)
    fields = _nt_names(nts[1])
    for x in nts[2:end]
        if _nt_names(x) != fields
            throw(ArgumentError("All NamedTuple arguments to map must have the same fields in the same order"))
        end
    end
    N = nfields(nts[1])
    M = length(nts)

    NT = NamedTuple{fields}
    args = Expr[:(f($(Expr[:(getfield(nts[$i], $j)) for i = 1:M]...))) for j = 1:N]
    quote
        namedtuple($NT, $(args...))
    end
end

# a version of `in` for the older world these generated functions run in
function sym_in(x, itr)
    for y in itr
        y === x && return true
    end
    return false
end

@generated function merge(a::NamedTuple{an}, b::NamedTuple{bn}) where {an, bn}
    names = Symbol[an...]
    for n in bn
        if !sym_in(n, an)
            push!(names, n)
        end
    end
    vals = map(names) do n
        if sym_in(n, bn)
            :(getfield(b, $(Expr(:quote, n))))
        else
            :(getfield(a, $(Expr(:quote, n))))
        end
    end
    names = (names...,)
    :(namedtuple(NamedTuple{$names}, $(vals...)))
end

@generated function structdiff(a::NamedTuple{an},
                               b::Union{NamedTuple{bn},Type{NamedTuple{bn}}}) where {an,bn}
    names = Symbol[]
    for n in an
        if !sym_in(n, bn)
            push!(names, n)
        end
    end
    vals = map(names) do n
        :(getfield(a, $(Expr(:quote, n))))
    end
    names = (names...,)
    :(namedtuple(NamedTuple{$names}, $(vals...)))
end
