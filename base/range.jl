# This file is a part of Julia. License is MIT: http://julialang.org/license

## 1-dimensional ranges ##

abstract Range{T} <: AbstractArray{T,1}

## ordinal ranges

abstract OrdinalRange{T,S} <: Range{T}
abstract AbstractUnitRange{T} <: OrdinalRange{T,Int}

# ordinal ranges rely on having an ordering
immutable HasOrder{b} end
HasOrder{T<:Real}(::Type{T}) = HasOrder{true}()
HasOrder{T}(::Type{T}) = HasOrder{false}()

immutable StepRange{T,S} <: OrdinalRange{T,S}
    start::T
    step::S
    stop::T

    function StepRange(start::T, step::S, stop::T)
        new(start, step, steprange_last(start,step,stop))
    end
end

# to make StepRange constructor inlineable, so optimizer can see `step` value
function steprange_last{T}(start::T, step, stop)
    if isa(start,AbstractFloat) || isa(step,AbstractFloat)
        throw(ArgumentError("StepRange should not be used with floating point"))
    end
    z = zero(step)
    step == z && throw(ArgumentError("step cannot be zero"))

    if stop == start
        last = stop
    else
        if (step > z) != (stop > start)
            last = steprange_last_empty(start, step, stop)
        else
            diff = stop - start
            if T<:Signed && (diff > zero(diff)) != (stop > start)
                # handle overflowed subtraction with unsigned rem
                if diff > zero(diff)
                    remain = -convert(T, unsigned(-diff) % step)
                else
                    remain = convert(T, unsigned(diff) % step)
                end
            else
                remain = steprem(start,stop,step)
            end
            last = stop - remain
        end
    end
    last
end

function steprange_last_empty{T<:Integer}(start::T, step, stop)
    # empty range has a special representation where stop = start-1
    # this is needed to avoid the wrap-around that can happen computing
    # start - step, which leads to a range that looks very large instead
    # of empty.
    if step > zero(step)
        last = start - one(stop-start)
    else
        last = start + one(stop-start)
    end
    last
end
# For types where x+one(x) may not be well-defined
steprange_last_empty(start, step, stop) = start - step

steprem(start,stop,step) = (stop-start) % step

StepRange{T,S}(start::T, step::S, stop::T) = StepRange{T,S}(start, step, stop)

immutable UnitRange{T<:Real} <: AbstractUnitRange{T}
    start::T
    stop::T
    UnitRange(start, stop) = new(start, unitrange_last(start,stop))
end
UnitRange{T<:Real}(start::T, stop::T) = UnitRange{T}(start, stop)

unitrange_last(::Bool, stop::Bool) = stop
unitrange_last{T<:Integer}(start::T, stop::T) =
    ifelse(stop >= start, stop, convert(T,start-one(stop-start)))
unitrange_last{T}(start::T, stop::T) =
    ifelse(stop >= start, convert(T,start+floor(stop-start)),
                          convert(T,start-one(stop-start)))

"""
    Base.OneTo(n)

Define an `AbstractUnitRange` that behaves like `1:n`, with the added
distinction that the lower limit is guaranteed (by the type system) to
be 1.
"""
immutable OneTo{T<:Integer} <: AbstractUnitRange{T}
    stop::T
    OneTo(stop) = new(max(zero(T), stop))
end
OneTo{T<:Integer}(stop::T) = OneTo{T}(stop)

colon(a::Real, b::Real) = colon(promote(a,b)...)

colon{T<:Real}(start::T, stop::T) = UnitRange{T}(start, stop)

range(a::Real, len::Integer) = UnitRange{typeof(a)}(a, oftype(a, a+len-1))

colon{T}(start::T, stop::T) = StepRange(start, one(stop-start), stop)

range{T}(a::T, len::Integer) = range(a, one(a-a), len)

# first promote start and stop, leaving step alone
# this is for non-numeric ranges where step can be quite different
colon{A<:Real,C<:Real}(a::A, b, c::C) = colon(convert(promote_type(A,C),a), b, convert(promote_type(A,C),c))

"""
    colon(start, [step], stop)

Called by `:` syntax for constructing ranges.
"""
colon{T<:Real}(start::T, step, stop::T) = StepRange(start, step, stop)

"""
    :(start, [step], stop)

Range operator. `a:b` constructs a range from `a` to `b` with a step size of 1, and `a:s:b`
is similar but uses a step size of `s`. These syntaxes call the function `colon`. The colon
is also used in indexing to select whole dimensions.
"""
colon{T<:Real}(start::T, step::T, stop::T) = StepRange(start, step, stop)
colon{T<:Real}(start::T, step::Real, stop::T) = StepRange(promote(start, step, stop)...)

colon{T}(start::T, step, stop::T) = StepRange(start, step, stop)

"""
    range(start, [step], length)

Construct a range by length, given a starting value and optional step (defaults to 1).
"""
range{T,S}(a::T, step::S, len::Integer) = _range(HasOrder(T), a, step, len)
_range{T,S}(::HasOrder{true}, a::T, step::S, len::Integer) = StepRange{T,S}(a, step, convert(T, a+step*(len-1)))
_range{T,S}(::HasOrder{false}, a::T, step::S, len::Integer) = StepRangeHiLo(promote(a, step)..., len)

## Roundoff-compensating ranges

# Use 2x arithmetic to achieve high precision.
# Necessary for ranges like 0.1:0.1:0.3, since 0.1+2*0.1 = 0.30000000000000004

immutable StepRangeHiLo{T} <: Range{T}
    # ref is the smallest-magnitude element in the range (not necessarily start or stop)
    ref_hi::T    # most significant bits
    ref_lo::T    # least significant bits

    # step is the separation between r[i] and r[i+1]
    step_hi::T   # ideally, has enough trailing zeros that (0:len-1)*step_hi is exact
    step_lo::T   # remaining bits of step

    offset::Int  # the index of ref
    len::Int     # length of the range

    function StepRangeHiLo(ref_hi, ref_lo, step_hi, step_lo, offset, len)
        len >= 0 || throw(ArgumentError("length cannot be negative, got $len"))
        1 <= offset <= max(len,1) || throw(ArgumentError("StepRangeHiLo: offset must be in [1,$len], got $offset"))
        len > 1 && step_hi == 0 && step_lo == 0 && throw(ArgumentError("range step cannot be 0"))
        new(ref_hi, ref_lo, step_hi, step_lo, offset, len)
    end
end

# promote types without risking a StackOverflowError
StepRangeHiLo(ref_hi, ref_lo, step_hi, step_lo, offset, len) =
    _srhl(promote(ref_hi, ref_lo, step_hi, step_lo)..., offset, len)
_srhl{T}(ref_hi::T, ref_lo::T, step_hi::T, step_lo::T, offset::Integer, len::Integer) =
    StepRangeHiLo{T}(ref_hi, ref_lo, step_hi, step_lo, Int(offset), Int(len))
_srhl(ref_hi, ref_lo, step_hi, step_lo, offset::Integer, len::Integer) =
    throw(ArgumentError("$ref_hi::$(typeof(ref_hi)), $ref_hi::$(typeof(ref_lo)), $ref_hi::$(typeof(step_hi)), and $ref_hi::$(typeof(step_lo)) cannot be promoted to a common type"))

function (::Type{StepRangeHiLo{T}}){T}(start::Integer, step::Integer, len::Integer, den::Integer)
    step == 0 && throw(ArgumentError("range step cannot be zero"))
    if len < 2
        return StepRangeHiLo{T}(div2(Int(start), 0, den)..., step/den, zero(T), 1, Int(len))
    end
    # index of smallest-magnitude value
    imin = clamp(round(Int, -start/step+1), 1, Int(len))
    # Compute smallest-magnitude element to 2x precision
    ref_n = start+(imin-1)*step  # this shouldn't overflow, so don't check
    ref_hi, ref_lo = div2(T(ref_n), 0, den)
    step_hi, step_lo = step2x(T, step, den, len)
    StepRangeHiLo{T}(ref_hi, ref_lo, step_hi, step_lo, imin, Int(len))
end

function (::Type{StepRangeHiLo{T}}){T}(start::T, step::T, len::Integer)
    step == 0 && throw(ArgumentError("range step cannot be zero"))
    StepRangeHiLo{T}(start, zero(T), step, zero(T), 1, Int(len))
end

function step2x{T}(::Type{T}, step_n::Integer, step_d::Integer, len::Integer)
    # This computes step=step_n/step_d to 2x precision...
    step_hi_pre, step_lo = div2(T(step_n), zero(T), step_d)
    # ...and then truncates enough low bits of step_hi to ensure that
    # multiplication by 0:(len-1) is exact
    nb = ceil(UInt, log2(len-1))
    step_hi = truncbits(step_hi_pre, nb)
    step_lo += step_hi_pre - step_hi
    step_hi, step_lo
end

function StepRangeHiLo(a::AbstractFloat, st::AbstractFloat, len::Real, divisor::AbstractFloat)
    T = promote_type(typeof(a), typeof(st), typeof(divisor))
    m = maxintfloat(T)
    if abs(a) <= m && abs(st) <= m && abs(divisor) <= m
        ia, ist, idivisor = round(Int, a), round(Int, st), round(Int, divisor)
        if ia == a && ist == st && idivisor == divisor
            # We can return the high-precision range
            return StepRangeHiLo{T}(ia, ist, Int(len), idivisor)
        end
    end
    # Fall back to 1x-precision range
    StepRangeHiLo{T}(T(a/divisor), zero(T), T(st/divisor), zero(T), 1, Int(len))
end

function colon{T<:Union{Float16,Float32,Float64}}(start::T, step::T, stop::T)
    step == 0 && throw(ArgumentError("range step cannot be zero"))
    len = max(0, floor(Int, (stop-start)/step) + 1)
    # Because len might be too small by 1 due to roundoff error, let's
    # see if the inputs have exact rational approximations (and if so,
    # perform all computations in terms of the rationals)
    step_n, step_d = rat(step)
    if T(step_n/step_d) == step
        start_n, start_d = rat(start)
        stop_n, stop_d = rat(stop)
        if T(start_n/start_d) == start && T(stop_n/stop_d) == stop
            den = lcm(start_d, step_d) # use same denominator for start and step
            m = maxintfloat(T)
            if abs(start*den) <= m && abs(step*den) <= m
                start_n = round(Int, start*den)
                step_n = round(Int, step*den)
                len = max(0, div(den*stop_n - stop_d*start_n + step_n*stop_d, step_n*stop_d))
                # Integer ops could overflow, so check that this makes sense
                if T(start_n/den) == start && T(step_n/den) == step &&
                    (isbetween(start, start + (len-1)*step, stop + step/2) &&
                     !isbetween(start, start + len*step, stop))
                    # Return a 2x precision range
                    return StepRangeHiLo{T}(start_n, step_n, len, den)
                end
            end
        end
    end
    # Fall back to the 1x precision range
    StepRangeHiLo{T}(start, zero(T), step, zero(T), 1, len)
end

colon{T<:AbstractFloat}(a::T, b::T) = colon(a, one(a), b)

colon{T<:Real}(a::T, b::AbstractFloat, c::T) = _colon(promote(a,b,c)...)
colon{T<:AbstractFloat}(a::T, b::AbstractFloat, c::T) = _colon(promote(a,b,c)...)
colon{T<:AbstractFloat}(a::T, b::Real, c::T) = _colon(promote(a,b,c)...)
_colon{T}(a::T, b::T, c::T) = colon(a, b, c)
_colon(a, b, c) = throw(ArgumentError("$a::$(typeof(a)), $b::$(typeof(b)), and $c::$(typeof(c)) cannot be promoted to a common type"))

range(a::AbstractFloat, len::Integer) = StepRangeHiLo(a,zero(a),one(a),zero(a),1,len)
range(a::AbstractFloat, st::AbstractFloat, len::Integer) = _range(promote(a, st)..., len)
range(a::Real, st::AbstractFloat, len::Integer) = range(float(a), st, len)
range(a::AbstractFloat, st::Real, len::Integer) = range(a, float(st), len)
function _range{T<:Union{Float16,Float32,Float64}}(a::T, st::T, len::Integer)
    start_n, start_d = rat(a)
    step_n, step_d = rat(st)
    if T(start_n/start_d) == a && T(step_n/step_d) == st
        den = lcm(start_d, step_d)
        m = maxintfloat(T)
        if abs(den*a) <= m && abs(den*st) <= m
            start_n = round(Int, den*a)
            step_n = round(Int, den*st)
            return StepRangeHiLo{T}(start_n, step_n, len, den)
        end
    end
    StepRangeHiLo{T}(a, zero(a), st, zero(st), 1, Int(len))
end
_range{T<:AbstractFloat}(a::T, st::T, len::Integer) = StepRangeHiLo{T}(a, zero(a), st, zero(st), 1, Int(len))
_range(a, st, len::Integer) = throw(ArgumentError("range: $a::$(typeof(a)) and $st::$(typeof(st)) cannot be promoted to a common type"))

## linspace and logspace

immutable LinSpace{T} <: Range{T}
    start::T
    stop::T
    len::Int
    lendiv::Int

    function LinSpace(start,stop,len)
        len >= 0 || throw(ArgumentError("linspace($start, $stop, $len): negative length"))
        if len == 1
            start == stop || throw(ArgumentError("linspace($start, $stop, $len): endpoints differ"))
            return new(start, stop, 1, 1)
        end
        new(start,stop,len,max(len-1,1))
    end
end

function LinSpace(start, stop, len::Integer)
    T = typeof((stop-start)/len)
    LinSpace{T}(start, stop, len)
end

"""
    linspace(start, stop, n=50)

Construct a range of `n` linearly spaced elements from `start` to `stop`.

```jldoctest
julia> linspace(1.3,2.9,9)
9-element LinSpace{Float64}:
 1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9
```
"""
linspace(start, stop, len::Real=50) = LinSpace(start, stop, Int(len))

# For Float16, Float32, and Float64, linspace returns a StepRangeHiLo
function linspace{T<:Union{Float16,Float32,Float64}}(start::T, stop::T, len::Integer)
    len < 2 && return _linspace1(T, start, stop, len)
    # Attempt to find exact rational approximations
    start_n, start_d = rat(start)
    stop_n, stop_d = rat(stop)
    den = lcm(start_d, stop_d)
    m = maxintfloat(T)
    if abs(den*start) <= m && abs(den*stop) <= m
        start_n = round(Int, den*start)
        stop_n = round(Int, den*stop)
        if T(start_n/den) == start && T(stop_n/den) == stop
            return linspace(T, start_n, stop_n, len, den)
        end
    end
    _linspace(start, stop, len)
end

function _linspace{T<:Union{Float16,Float32,Float64}}(start::T, stop::T, len)
    if !isfinite(start) || !isfinite(stop)
        len > 2 && throw(ArgumentError("start and stop must be finite, got $start and $stop"))
        return StepRangeHiLo{T}(start, zero(T), -start, stop, 1, len)
    end
    # Find the index that returns the smallest-magnitude element
    Δ, Δfac = stop-start, 1
    if ~isfinite(Δ)   # handle overflow
        Δ, Δfac = stop/len - start/len, Int(len)
    end
    tmin = -(start/Δ)/Δfac            # interpolation t such that return value is 0
    imin = round(Int, tmin*(len-1)+1)
    if 1 < imin < len
        # The smallest-magnitude element is in the interior
        t = (imin-1)/(len-1)
        ref = T((1-t)*start + t*stop)
        step = imin-1 < len-imin ? (ref-start)/(imin-1) : (stop-ref)/(len-imin)
    elseif imin <= 1
        imin = 1
        ref = start
        step = (Δ/(len-1))*Δfac
    else
        imin = Int(len)
        ref = stop
        step = (Δ/(len-1))*Δfac
    end
    if len == 2 && ~isfinite(step)
        # For very large endpoints where step overflows, exploit the
        # split-representation to handle the overflow
        return StepRangeHiLo{T}(start, zero(T), -start, stop, 1, 2)
    end
    # 2x calculations to get high precision endpoint matching while also
    # preventing overflow in ref_hi+(i-offset)*step_hi
    t, k = prevfloat(realmax(T)), max(imin-1, len-imin)
    step_hi_pre = clamp(step, max(-(t+ref)/k, (-t+ref)/k), min((t-ref)/k, (t+ref)/k))
    nb = ceil(UInt, log2(len-1))
    step_hi = truncbits(step_hi_pre, nb)
    x1_hi, x1_lo = add2((1-imin)*step_hi, ref)
    x2_hi, x2_lo = add2((len-imin)*step_hi, ref)
    a, b = (start - x1_hi) - x1_lo, (stop - x2_hi) - x2_lo
    step_lo = (b - a)/(len - 1)
    ref_lo = abs(ref) < eps(max(abs(start), abs(stop))) ? zero(T) : a - (1 - imin)*step_lo
    StepRangeHiLo(ref, ref_lo, step_hi, step_lo, imin, Int(len))
end

# linspace for rational numbers, start = start_n/den, stop = stop_n/den
# Note this returns a StepRangeHiLo
function linspace{T}(::Type{T}, start_n::Integer, stop_n::Integer, len::Integer, den::Integer)
    len < 2 && return _linspace1(T, start_n/den, stop_n/den, len)
    tmin = -start_n/(Float64(stop_n) - Float64(start_n))
    imin = round(Int, tmin*(len-1)+1)
    imin = clamp(imin, 1, Int(len))
    # Compute (1-t)*a and t*b separately in 2x precision (itp = interpolant)...
    dent = (den, len-1)  # represent products as a tuple to eliminate risk of overflow
    start_itp_hi, start_itp_lo = proddiv(T, (len-imin, start_n), dent)
    stop_itp_hi, stop_itp_lo = proddiv(T, (imin-1, stop_n), dent)
    # ...and then combine them to make ref
    ref_hi, ref_lo = add2(start_itp_hi, stop_itp_hi)
    ref_hi, ref_lo = add2(ref_hi, ref_lo + (start_itp_lo + stop_itp_lo))
    # Compute step to 2x precision without risking overflow...
    end_hi, end_lo = proddiv(T, (stop_n,), dent)
    beg_hi, beg_lo = proddiv(T, (-start_n,), dent)
    step_hi_pre, step_lo = add2(end_hi, beg_hi)
    step_hi_pre, step_lo = add2(step_hi_pre, step_lo + (end_lo + beg_lo))
    # ...and then truncate enough low-bits of step_hi to ensure that
    # multiplication by 0:len-1 is exact
    nb = ceil(UInt, log2(len-1))
    step_hi = truncbits(step_hi_pre, nb)
    step_lo += step_hi_pre - step_hi
    StepRangeHiLo(ref_hi, ref_lo, step_hi, step_lo, imin, Int(len))
end

# For len < 2
function _linspace1{T}(::Type{T}, start, stop, len)
    len >= 0 || throw(ArgumentError("linspace($start, $stop, $len): negative length"))
    if len <= 1
        len == 1 && (start == stop || throw(ArgumentError("linspace($start, $stop, $len): endpoints differ")))
        # Ensure that first(r)==start and last(r)==stop even for len==0
        return StepRangeHiLo(start, zero(T), start, -stop, 1, len)
    end
    throw(ArgumentError("should only be called for len < 2, got $len"))
end

linspace(start::Real, stop::Real, len::Integer) = _linspace(promote(start, stop)..., len)
_linspace{T<:AbstractFloat}(start::T, stop::T, len) = linspace(start, stop, len)
_linspace{T<:Integer}(start::T, stop::T, len) = linspace(Float64, start, stop, len, 1)
_linspace(start, stop, len) = throw(ArgumentError("$start::$(typeof(start)) and $stop::$(typeof(stop)) cannot be promoted to a common type"))

function show(io::IO, r::LinSpace)
    print(io, "linspace(")
    show(io, first(r))
    print(io, ',')
    show(io, last(r))
    print(io, ',')
    show(io, length(r))
    print(io, ')')
end

"""
`print_range(io, r)` prints out a nice looking range r in terms of its elements
as if it were `collect(r)`, dependent on the size of the
terminal, and taking into account whether compact numbers should be shown.
It figures out the width in characters of each element, and if they
end up too wide, it shows the first and last elements separated by a
horizontal elipsis. Typical output will look like `1.0,2.0,3.0,…,4.0,5.0,6.0`.

`print_range(io, r, pre, sep, post, hdots)` uses optional
parameters `pre` and `post` characters for each printed row,
`sep` separator string between printed elements,
`hdots` string for the horizontal ellipsis.
"""
function print_range(io::IO, r::Range,
                     pre::AbstractString = " ",
                     sep::AbstractString = ",",
                     post::AbstractString = "",
                     hdots::AbstractString = ",\u2026,") # horiz ellipsis
    # This function borrows from print_matrix() in show.jl
    # and should be called by show and display
    limit = get(io, :limit, false)
    sz = displaysize(io)
    if !haskey(io, :compact)
        io = IOContext(io, compact=true)
    end
    screenheight, screenwidth = sz[1] - 4, sz[2]
    screenwidth -= length(pre) + length(post)
    postsp = ""
    sepsize = length(sep)
    m = 1 # treat the range as a one-row matrix
    n = length(r)
    # Figure out spacing alignments for r, but only need to examine the
    # left and right edge columns, as many as could conceivably fit on the
    # screen, with the middle columns summarized by horz, vert, or diag ellipsis
    maxpossiblecols = div(screenwidth, 1+sepsize) # assume each element is at least 1 char + 1 separator
    colsr = n <= maxpossiblecols ? (1:n) : [1:div(maxpossiblecols,2)+1; (n-div(maxpossiblecols,2)):n]
    rowmatrix = r[colsr]' # treat the range as a one-row matrix for print_matrix_row
    A = alignment(io, rowmatrix, 1:m, 1:length(rowmatrix), screenwidth, screenwidth, sepsize) # how much space range takes
    if n <= length(A) # cols fit screen, so print out all elements
        print(io, pre) # put in pre chars
        print_matrix_row(io,rowmatrix,A,1,1:n,sep) # the entire range
        print(io, post) # add the post characters
    else # cols don't fit so put horiz ellipsis in the middle
        # how many chars left after dividing width of screen in half
        # and accounting for the horiz ellipsis
        c = div(screenwidth-length(hdots)+1,2)+1 # chars remaining for each side of rowmatrix
        alignR = reverse(alignment(io, rowmatrix, 1:m, length(rowmatrix):-1:1, c, c, sepsize)) # which cols of rowmatrix to put on the right
        c = screenwidth - sum(map(sum,alignR)) - (length(alignR)-1)*sepsize - length(hdots)
        alignL = alignment(io, rowmatrix, 1:m, 1:length(rowmatrix), c, c, sepsize) # which cols of rowmatrix to put on the left
        print(io, pre)   # put in pre chars
        print_matrix_row(io, rowmatrix,alignL,1,1:length(alignL),sep) # left part of range
        print(io, hdots) # horizontal ellipsis
        print_matrix_row(io, rowmatrix,alignR,1,length(rowmatrix)-length(alignR)+1:length(rowmatrix),sep) # right part of range
        print(io, post)  # post chars
    end
end

"""
    logspace(start::Real, stop::Real, n::Integer=50)

Construct a vector of `n` logarithmically spaced numbers from `10^start` to `10^stop`.

```jldoctest
julia> logspace(1.,10.,5)
5-element Array{Float64,1}:
   10.0
 1778.28
    3.16228e5
    5.62341e7
    1.0e10
```
"""
logspace(start::Real, stop::Real, n::Integer=50) = 10.^linspace(start, stop, n)

## interface implementations

size(r::Range) = (length(r),)

isempty(r::StepRange) =
    (r.start != r.stop) & ((r.step > zero(r.step)) != (r.stop > r.start))
isempty(r::AbstractUnitRange) = first(r) > last(r)
isempty(r::StepRangeHiLo) = length(r) == 0
isempty(r::LinSpace) = length(r) == 0

"""
    step(r)

Get the step size of a `Range` object.
```jldoctest
julia> step(1:10)
1

julia> step(1:2:10)
2

julia> step(2.5:0.3:10.9)
0.3

julia> step(linspace(2.5,10.9,85))
0.1
```
"""
step(r::StepRange) = r.step
step(r::AbstractUnitRange) = 1
step(r::StepRangeHiLo) = r.step_hi + r.step_lo
step(r::LinSpace) = (last(r)-first(r))/r.lendiv

unsafe_length(r::Range) = length(r)  # generic fallback

function unsafe_length(r::StepRange)
    n = Integer(div(r.stop+r.step - r.start, r.step))
    isempty(r) ? zero(n) : n
end
length(r::StepRange) = unsafe_length(r)
unsafe_length(r::AbstractUnitRange) = Integer(last(r) - first(r) + 1)
unsafe_length(r::OneTo) = r.stop
length(r::AbstractUnitRange) = unsafe_length(r)
length(r::OneTo) = unsafe_length(r)
length(r::StepRangeHiLo) = Integer(r.len)
length(r::LinSpace) = r.len

function length{T<:Union{Int,UInt,Int64,UInt64}}(r::StepRange{T})
    isempty(r) && return zero(T)
    if r.step > 1
        return checked_add(convert(T, div(unsigned(r.stop - r.start), r.step)), one(T))
    elseif r.step < -1
        return checked_add(convert(T, div(unsigned(r.start - r.stop), -r.step)), one(T))
    else
        checked_add(div(checked_sub(r.stop, r.start), r.step), one(T))
    end
end

function length{T<:Union{Int,Int64}}(r::AbstractUnitRange{T})
    @_inline_meta
    checked_add(checked_sub(last(r), first(r)), one(T))
end
length{T<:Union{Int,Int64}}(r::OneTo{T}) = T(r.stop)

length{T<:Union{UInt,UInt64}}(r::AbstractUnitRange{T}) =
    r.stop < r.start ? zero(T) : checked_add(last(r) - first(r), one(T))

# some special cases to favor default Int type
let smallint = (Int === Int64 ?
                Union{Int8,UInt8,Int16,UInt16,Int32,UInt32} :
                Union{Int8,UInt8,Int16,UInt16})
    global length

    function length{T <: smallint}(r::StepRange{T})
        isempty(r) && return Int(0)
        div(Int(r.stop)+Int(r.step) - Int(r.start), Int(r.step))
    end

    length{T <: smallint}(r::AbstractUnitRange{T}) = Int(last(r)) - Int(first(r)) + 1
    length{T <: smallint}(r::OneTo{T}) = Int(r.stop)
end

first{T}(r::OrdinalRange{T}) = convert(T, r.start)
first{T}(r::OneTo{T}) = one(T)
first(r::StepRangeHiLo) = unsafe_getindex(r, 1)
first(r::LinSpace) = r.start

last{T}(r::OrdinalRange{T}) = convert(T, r.stop)
last(r::StepRangeHiLo) = unsafe_getindex(r, length(r))
last(r::LinSpace) = r.stop

minimum(r::AbstractUnitRange) = isempty(r) ? throw(ArgumentError("range must be non-empty")) : first(r)
maximum(r::AbstractUnitRange) = isempty(r) ? throw(ArgumentError("range must be non-empty")) : last(r)
minimum(r::Range)  = isempty(r) ? throw(ArgumentError("range must be non-empty")) : min(first(r), last(r))
maximum(r::Range)  = isempty(r) ? throw(ArgumentError("range must be non-empty")) : max(first(r), last(r))

ctranspose(r::Range) = [x for _=1:1, x=r]
transpose(r::Range) = r'

# Ranges are immutable
copy(r::Range) = r


## iteration

start(r::Union{StepRangeHiLo,LinSpace}) = 1
done(r::Union{StepRangeHiLo,LinSpace}, i::Int) = length(r) < i
function next(r::Union{StepRangeHiLo,LinSpace}, i::Int)
    @_inline_meta
    unsafe_getindex(r, i), i+1
end

start(r::StepRange) = oftype(r.start + r.step, r.start)
next{T}(r::StepRange{T}, i) = (convert(T,i), i+r.step)
done{T,S}(r::StepRange{T,S}, i) = isempty(r) | (i < min(r.start, r.stop)) | (i > max(r.start, r.stop))
done{T,S}(r::StepRange{T,S}, i::Integer) =
    isempty(r) | (i == oftype(i, r.stop) + r.step)

start{T}(r::UnitRange{T}) = oftype(r.start + one(T), r.start)
next{T}(r::AbstractUnitRange{T}, i) = (convert(T, i), i + one(T))
done{T}(r::AbstractUnitRange{T}, i) = i == oftype(i, r.stop) + one(T)

start{T}(r::OneTo{T}) = one(T)

# some special cases to favor default Int type to avoid overflow
let smallint = (Int === Int64 ?
                Union{Int8,UInt8,Int16,UInt16,Int32,UInt32} :
                Union{Int8,UInt8,Int16,UInt16})
    global start
    global next
    start{T<:smallint}(r::StepRange{T}) = convert(Int, r.start)
    next{T<:smallint}(r::StepRange{T}, i) = (i % T, i + r.step)
    start{T<:smallint}(r::UnitRange{T}) = convert(Int, r.start)
    next{T<:smallint}(r::AbstractUnitRange{T}, i) = (i % T, i + 1)
    start{T<:smallint}(r::OneTo{T}) = 1
end

## indexing

function getindex{T}(v::UnitRange{T}, i::Integer)
    @_inline_meta
    ret = convert(T, first(v) + i - 1)
    @boundscheck ((i > 0) & (ret <= v.stop) & (ret >= v.start)) || throw_boundserror(v, i)
    ret
end

function getindex{T}(v::OneTo{T}, i::Integer)
    @_inline_meta
    @boundscheck ((i > 0) & (i <= v.stop)) || throw_boundserror(v, i)
    convert(T, i)
end

function getindex{T}(v::Range{T}, i::Integer)
    @_inline_meta
    ret = convert(T, first(v) + (i - 1)*step(v))
    ok = ifelse(step(v) > zero(step(v)),
                (ret <= v.stop) & (ret >= v.start),
                (ret <= v.start) & (ret >= v.stop))
    @boundscheck ((i > 0) & ok) || throw_boundserror(v, i)
    ret
end

function getindex(r::Union{StepRangeHiLo,LinSpace}, i::Integer)
    @_inline_meta
    @boundscheck checkbounds(r, i)
    unsafe_getindex(r, i)
end

# This is separate to make it useful even when running with --check-bounds=yes
function unsafe_getindex(r::StepRangeHiLo, i::Integer)
    # Use 2x arithmetic for high precision
    u = i - r.offset
    shift_hi, shift_lo = u*r.step_hi, u*r.step_lo
    x_hi, x_lo = add2(r.ref_hi, shift_hi)
    x_hi + (x_lo + (shift_lo + r.ref_lo))
end

function _getindex_hiprec(r::StepRangeHiLo, i::Integer)
    u = i - r.offset
    shift_hi, shift_lo = u*r.step_hi, u*r.step_lo
    x_hi, x_lo = add2(r.ref_hi, shift_hi)
    add2(x_hi, x_lo + (shift_lo + r.ref_lo))
end

function unsafe_getindex(r::LinSpace, i::Integer)
    d = r.lendiv
    j, a, b = ifelse(2i >= length(r), (i-1, r.start, r.stop), (length(r)-i, r.stop, r.start))
    lerpi(j, d, a, b)
end

# High-precision interpolation. Accurate for t ∈ [0.5,1], so that 1-t is exact.
function lerpi{T}(j::Integer, d::Integer, a::T, b::T)
    @_inline_meta
    t = j/d
    # computes (1-t)*a + t*b
    # T(fma(t, b, fma(-t, a, a)))
    T((1-t)*a + t*b)
end

getindex(r::Range, ::Colon) = copy(r)

function getindex{T<:Integer}(r::AbstractUnitRange, s::AbstractUnitRange{T})
    @_inline_meta
    @boundscheck checkbounds(r, s)
    f = first(r)
    st = oftype(f, f + first(s)-1)
    range(st, length(s))
end

function getindex{T}(r::OneTo{T}, s::OneTo)
    @_inline_meta
    @boundscheck checkbounds(r, s)
    OneTo(T(s.stop))
end

function getindex{T<:Integer}(r::AbstractUnitRange, s::StepRange{T})
    @_inline_meta
    @boundscheck checkbounds(r, s)
    st = oftype(first(r), first(r) + s.start-1)
    range(st, step(s), length(s))
end

function getindex{T<:Integer}(r::StepRange, s::Range{T})
    @_inline_meta
    @boundscheck checkbounds(r, s)
    st = oftype(r.start, r.start + (first(s)-1)*step(r))
    range(st, step(r)*step(s), length(s))
end

function getindex(r::Union{StepRangeHiLo,LinSpace}, s::OrdinalRange)
    @_inline_meta
    @boundscheck checkbounds(r, s)
    sl = length(s)
    ifirst = first(s)
    ilast = last(s)
    vfirst = unsafe_getindex(r, ifirst)
    vlast  = unsafe_getindex(r, ilast)
    return linspace(vfirst, vlast, sl)
end

show(io::IO, r::Range) = print(io, repr(first(r)), ':', repr(step(r)), ':', repr(last(r)))
show(io::IO, r::UnitRange) = print(io, repr(first(r)), ':', repr(last(r)))
show(io::IO, r::OneTo) = print(io, "Base.OneTo(", r.stop, ")")

=={T<:Range}(r::T, s::T) = (first(r) == first(s)) & (step(r) == step(s)) & (last(r) == last(s))
==(r::OrdinalRange, s::OrdinalRange) = (first(r) == first(s)) & (step(r) == step(s)) & (last(r) == last(s))
=={T<:Union{StepRangeHiLo,LinSpace}}(r::T, s::T) = (first(r) == first(s)) & (length(r) == length(s)) & (last(r) == last(s))

function ==(r::Range, s::Range)
    lr = length(r)
    if lr != length(s)
        return false
    end
    u, v = start(r), start(s)
    while !done(r, u)
        x, u = next(r, u)
        y, v = next(s, v)
        if x != y
            return false
        end
    end
    return true
end

intersect(r::OneTo, s::OneTo) = OneTo(min(r.stop,s.stop))

intersect{T1<:Integer, T2<:Integer}(r::AbstractUnitRange{T1}, s::AbstractUnitRange{T2}) = max(first(r),first(s)):min(last(r),last(s))

intersect{T<:Integer}(i::Integer, r::AbstractUnitRange{T}) =
    i < first(r) ? (first(r):i) :
    i > last(r)  ? (i:last(r))  : (i:i)

intersect{T<:Integer}(r::AbstractUnitRange{T}, i::Integer) = intersect(i, r)

function intersect{T1<:Integer, T2<:Integer}(r::AbstractUnitRange{T1}, s::StepRange{T2})
    if isempty(s)
        range(first(r), 0)
    elseif step(s) == 0
        intersect(first(s), r)
    elseif step(s) < 0
        intersect(r, reverse(s))
    else
        sta = first(s)
        ste = step(s)
        sto = last(s)
        lo = first(r)
        hi = last(r)
        i0 = max(sta, lo + mod(sta - lo, ste))
        i1 = min(sto, hi - mod(hi - sta, ste))
        i0:ste:i1
    end
end

function intersect{T1<:Integer, T2<:Integer}(r::StepRange{T1}, s::AbstractUnitRange{T2})
    if step(r) < 0
        reverse(intersect(s, reverse(r)))
    else
        intersect(s, r)
    end
end

function intersect(r::StepRange, s::StepRange)
    if isempty(r) || isempty(s)
        return range(first(r), step(r), 0)
    elseif step(s) < 0
        return intersect(r, reverse(s))
    elseif step(r) < 0
        return reverse(intersect(reverse(r), s))
    end

    start1 = first(r)
    step1 = step(r)
    stop1 = last(r)
    start2 = first(s)
    step2 = step(s)
    stop2 = last(s)
    a = lcm(step1, step2)

    # if a == 0
    #     # One or both ranges have step 0.
    #     if step1 == 0 && step2 == 0
    #         return start1 == start2 ? r : Range(start1, 0, 0)
    #     elseif step1 == 0
    #         return start2 <= start1 <= stop2 && rem(start1 - start2, step2) == 0 ? r : Range(start1, 0, 0)
    #     else
    #         return start1 <= start2 <= stop1 && rem(start2 - start1, step1) == 0 ? (start2:step1:start2) : Range(start1, step1, 0)
    #     end
    # end

    g, x, y = gcdx(step1, step2)

    if rem(start1 - start2, g) != 0
        # Unaligned, no overlap possible.
        return range(start1, a, 0)
    end

    z = div(start1 - start2, g)
    b = start1 - x * z * step1
    # Possible points of the intersection of r and s are
    # ..., b-2a, b-a, b, b+a, b+2a, ...
    # Determine where in the sequence to start and stop.
    m = max(start1 + mod(b - start1, a), start2 + mod(b - start2, a))
    n = min(stop1 - mod(stop1 - b, a), stop2 - mod(stop2 - b, a))
    m:a:n
end

function intersect(r1::Range, r2::Range, r3::Range, r::Range...)
    i = intersect(intersect(r1, r2), r3)
    for t in r
        i = intersect(i, t)
    end
    i
end

# findin (the index of intersection)
function _findin{T1<:Integer, T2<:Integer}(r::Range{T1}, span::AbstractUnitRange{T2})
    local ifirst
    local ilast
    fspan = first(span)
    lspan = last(span)
    fr = first(r)
    lr = last(r)
    sr = step(r)
    if sr > 0
        ifirst = fr >= fspan ? 1 : ceil(Integer,(fspan-fr)/sr)+1
        ilast = lr <= lspan ? length(r) : length(r) - ceil(Integer,(lr-lspan)/sr)
    elseif sr < 0
        ifirst = fr <= lspan ? 1 : ceil(Integer,(lspan-fr)/sr)+1
        ilast = lr >= fspan ? length(r) : length(r) - ceil(Integer,(lr-fspan)/sr)
    else
        ifirst = fr >= fspan ? 1 : length(r)+1
        ilast = fr <= lspan ? length(r) : 0
    end
    ifirst, ilast
end

function findin{T1<:Integer, T2<:Integer}(r::AbstractUnitRange{T1}, span::AbstractUnitRange{T2})
    ifirst, ilast = _findin(r, span)
    ifirst:ilast
end

function findin{T1<:Integer, T2<:Integer}(r::Range{T1}, span::AbstractUnitRange{T2})
    ifirst, ilast = _findin(r, span)
    ifirst:1:ilast
end

## linear operations on ranges ##

-(r::OrdinalRange) = range(-first(r), -step(r), length(r))
-(r::StepRangeHiLo) = StepRangeHiLo(-r.ref_hi, -r.ref_lo, -r.step_hi, -r.step_lo, r.offset, length(r))
-(r::LinSpace) = LinSpace(-r.start, -r.stop, r.len)

+(x::Number, r::AbstractUnitRange) = range(x + first(r), length(r))
+(x::Number, r::Range) = (x+first(r)):step(r):(x+last(r))
+(x::Number, r::StepRangeHiLo) = StepRangeHiLo(add2(r.ref_hi, r.ref_lo, x)..., r.step_hi, r.step_lo, r.offset, r.len)
function +(x::Number, r::LinSpace)
    LinSpace(x + r.start, x + r.stop, r.len)
end
+(r::Range, x::Number) = x + r

-(x::Number, r::Range)      = (x-first(r)):-step(r):(x-last(r))
-(x::Number, r::StepRangeHiLo) = +(x, -r)
function -(x::Number, r::LinSpace)
    LinSpace(x - r.start, x - r.stop, r.len)
end

-(r::Range, x::Number) = +(-x, r)

*(x::Number, r::OrdinalRange) = range(x*first(r), x*step(r), length(r))
*(x::Number, r::StepRangeHiLo)   = StepRangeHiLo(mul2(r.ref_hi, r.ref_lo, x)..., mul2(r.step_hi, r.step_lo, x)..., r.offset, r.len)
*(x::Number, r::LinSpace)     = LinSpace(x * r.start, x * r.stop, r.len)
*(r::Range, x::Number)        = x * r

/(r::OrdinalRange, x::Number) = range(first(r)/x, step(r)/x, length(r))
/(r::StepRangeHiLo, x::Number)   = StepRangeHiLo(div2(r.ref_hi, r.ref_lo, x)..., div2(r.step_hi, r.step_lo, x)..., r.offset, r.len)
/(r::LinSpace, x::Number)     = LinSpace(r.start / x, r.stop / x, r.len)

promote_rule{T1,T2}(::Type{UnitRange{T1}},::Type{UnitRange{T2}}) =
    UnitRange{promote_type(T1,T2)}
convert{T<:Real}(::Type{UnitRange{T}}, r::UnitRange{T}) = r
convert{T<:Real}(::Type{UnitRange{T}}, r::UnitRange) = UnitRange{T}(r.start, r.stop)

promote_rule{T1,T2}(::Type{OneTo{T1}},::Type{OneTo{T2}}) =
    OneTo{promote_type(T1,T2)}
convert{T<:Real}(::Type{OneTo{T}}, r::OneTo{T}) = r
convert{T<:Real}(::Type{OneTo{T}}, r::OneTo) = OneTo{T}(r.stop)

promote_rule{T1,UR<:AbstractUnitRange}(::Type{UnitRange{T1}}, ::Type{UR}) =
    UnitRange{promote_type(T1,eltype(UR))}
convert{T<:Real}(::Type{UnitRange{T}}, r::AbstractUnitRange) = UnitRange{T}(first(r), last(r))
convert(::Type{UnitRange}, r::AbstractUnitRange) = UnitRange(first(r), last(r))

promote_rule{T1a,T1b,T2a,T2b}(::Type{StepRange{T1a,T1b}},::Type{StepRange{T2a,T2b}}) =
    StepRange{promote_type(T1a,T2a),promote_type(T1b,T2b)}
convert{T1,T2}(::Type{StepRange{T1,T2}}, r::StepRange{T1,T2}) = r

promote_rule{T1a,T1b,UR<:AbstractUnitRange}(::Type{StepRange{T1a,T1b}},::Type{UR}) =
    StepRange{promote_type(T1a,eltype(UR)),promote_type(T1b,eltype(UR))}
convert{T1,T2}(::Type{StepRange{T1,T2}}, r::Range) =
    StepRange{T1,T2}(convert(T1, first(r)), convert(T2, step(r)), convert(T1, last(r)))
convert{T}(::Type{StepRange}, r::AbstractUnitRange{T}) =
    StepRange{T,T}(first(r), step(r), last(r))

promote_rule{T1,T2}(::Type{StepRangeHiLo{T1}},::Type{StepRangeHiLo{T2}}) =
    StepRangeHiLo{promote_type(T1,T2)}
convert{T}(::Type{StepRangeHiLo{T}}, r::StepRangeHiLo{T}) = r
function convert{T<:AbstractFloat,S}(::Type{StepRangeHiLo{T}}, r::StepRangeHiLo{S})
    # if start and step have a rational approximation in the old type,
    # then we transfer that rational approximation to the new type
    f, s = first(r), step(r)
    start_n, start_d = rat(f)
    step_n, step_d = rat(s)
    if S(start_n/start_d) == f && S(step_n/step_d) == s
        den = lcm(start_d, step_d)
        m = maxintfloat(T)
        if abs(f*den) <= m && abs(s*den) <= m
            start_n = round(Int, f*den)
            step_n = round(Int, s*den)
            if S(start_n/den) == f && S(step_n/den) == s
                return StepRangeHiLo{T}(start_n, step_n, length(r), den)
            end
        end
    end
    StepRangeHiLo{T}(r.ref_hi, r.ref_lo, r.step_hi, r.step_lo, r.offset, r.len)
end
convert{T}(::Type{StepRangeHiLo{T}}, r::StepRangeHiLo) =
    StepRangeHiLo{T}(r.ref_hi, r.ref_lo, r.step_hi, r.step_lo, r.offset, r.len)

promote_rule{F,OR<:OrdinalRange}(::Type{StepRangeHiLo{F}}, ::Type{OR}) =
    StepRangeHiLo{promote_type(F,eltype(OR))}
convert{T<:Union{Float16,Float32,Float64}}(::Type{StepRangeHiLo{T}}, r::OrdinalRange) =
    linspace(first(r), last(r), length(r))
convert{T<:AbstractFloat}(::Type{StepRangeHiLo{T}}, r::OrdinalRange) =
    StepRangeHiLo{T}(first(r), zero(T), step(r), zero(T), 1, length(r))
convert{T}(::Type{StepRangeHiLo}, r::OrdinalRange{T}) =
    convert(StepRangeHiLo{typeof(float(first(r)))}, r)

promote_rule{T1,T2}(::Type{LinSpace{T1}},::Type{LinSpace{T2}}) =
    LinSpace{promote_type(T1,T2)}
convert{T<:AbstractFloat}(::Type{LinSpace{T}}, r::LinSpace{T}) = r
convert{T<:AbstractFloat}(::Type{LinSpace{T}}, r::LinSpace) =
    LinSpace{T}(r.start, r.stop, r.len)

promote_rule{F,OR<:OrdinalRange}(::Type{LinSpace{F}}, ::Type{OR}) =
    LinSpace{promote_type(F,eltype(OR))}
convert{T<:AbstractFloat}(::Type{LinSpace{T}}, r::OrdinalRange) =
    linspace(convert(T, first(r)), convert(T, last(r)), length(r))
convert{T}(::Type{LinSpace}, r::OrdinalRange{T}) =
    convert(LinSpace{typeof(float(first(r)))}, r)

# Promote StepRangeHiLo to LinSpace
promote_rule{F,OR<:StepRangeHiLo}(::Type{LinSpace{F}}, ::Type{OR}) =
    LinSpace{promote_type(F,eltype(OR))}
convert{T<:AbstractFloat}(::Type{LinSpace{T}}, r::StepRangeHiLo) =
    LinSpace{T}(convert(T, first(r)), convert(T, last(r)), length(r))
convert{T<:AbstractFloat}(::Type{LinSpace}, r::StepRangeHiLo{T}) =
    convert(LinSpace{T}, r)


# +/- of ranges is defined in operators.jl (to be able to use @eval etc.)

## concatenation ##

function vcat{T}(rs::Range{T}...)
    n::Int = 0
    for ra in rs
        n += length(ra)
    end
    a = Array{T}(n)
    i = 1
    for ra in rs, x in ra
        @inbounds a[i] = x
        i += 1
    end
    return a
end

convert{T}(::Type{Array{T,1}}, r::Range{T}) = vcat(r)
collect(r::Range) = vcat(r)

reverse(r::OrdinalRange) = colon(last(r), -step(r), first(r))
reverse(r::StepRangeHiLo)   = StepRangeHiLo(r.ref_hi, r.ref_lo, -r.step_hi, -r.step_lo, length(r)-r.offset+1, length(r))
reverse(r::LinSpace)     = LinSpace(r.stop, r.start, length(r))

## sorting ##

issorted(r::AbstractUnitRange) = true
issorted(r::Range) = length(r) <= 1 || step(r) >= zero(step(r))

sort(r::AbstractUnitRange) = r
sort!(r::AbstractUnitRange) = r

sort(r::Range) = issorted(r) ? r : reverse(r)

sortperm(r::AbstractUnitRange) = 1:length(r)
sortperm(r::Range) = issorted(r) ? (1:1:length(r)) : (length(r):-1:1)

function sum{T<:Real}(r::Range{T})
    l = length(r)
    # note that a little care is required to avoid overflow in l*(l-1)/2
    return l * first(r) + (iseven(l) ? (step(r) * (l-1)) * (l>>1)
                                     : (step(r) * l) * ((l-1)>>1))
end

function sum(r::StepRangeHiLo)
    l = length(r)
    np, nn = l - r.offset, r.offset - 1
    sp, sn = sumints(np), sumints(nn)
    s_hi, s_lo = mul2(r.step_hi, r.step_lo, sp - sn)
    r_hi, r_lo = mul2(r.ref_hi, r.ref_lo, l)
    sm_hi, sm_lo = add2(s_hi, r_hi)
    add2(sm_hi, sm_lo + r_lo)[1]
end

# sum(1:n)
sumints(n::Integer) = iseven(n) ? (n+1)*(n>>1) : n*((n+1)>>1)


function mean{T<:Real}(r::Range{T})
    isempty(r) && throw(ArgumentError("mean of an empty range is undefined"))
    (first(r) + last(r)) / 2
end

median{T<:Real}(r::Range{T}) = mean(r)

function in(x, r::Range)
    n = step(r) == 0 ? 1 : round(Integer,(x-first(r))/step(r))+1
    n >= 1 && n <= length(r) && r[n] == x
end

in{T<:Integer}(x::Integer, r::AbstractUnitRange{T}) = (first(r) <= x) & (x <= last(r))
in{T<:Integer}(x, r::Range{T}) = isinteger(x) && !isempty(r) && x>=minimum(r) && x<=maximum(r) && (mod(convert(T,x),step(r))-mod(first(r),step(r)) == 0)
in(x::Char, r::Range{Char}) = !isempty(r) && x >= minimum(r) && x <= maximum(r) && (mod(Int(x) - Int(first(r)), step(r)) == 0)

### Numeric utilities

function rat(x)
    y = x
    a = d = 1
    b = c = 0
    m = maxintfloat(narrow(typeof(x)))
    while abs(y) <= m
        f = trunc(Int,y)
        y -= f
        a, c = f*a + c, a
        b, d = f*b + d, b
        max(abs(a), abs(b)) <= convert(Int,m) || return c, d
        oftype(x,a)/oftype(x,b) == x && break
        y = inv(y)
    end
    return a, b
end

narrow(::Type{Float64}) = Float32
narrow(::Type{Float32}) = Float16
narrow(::Type{Float16}) = Float16

truncbits(x::Float16, nb) = box(Float16, unbox(UInt16, _truncbits(box(UInt16, unbox(Float16, x)), nb)))
truncbits(x::Float32, nb) = box(Float32, unbox(UInt32, _truncbits(box(UInt32, unbox(Float32, x)), nb)))
truncbits(x::Float64, nb) = box(Float64, unbox(UInt64, _truncbits(box(UInt64, unbox(Float64, x)), nb)))
function _truncbits{U<:Unsigned}(xi::U, nb::Integer)
    @_inline_meta
    mask = typemax(U)
    xi & (mask << nb)
end
splitprec(x::Float16) = truncbits(x, UInt(5))
splitprec(x::Float32) = truncbits(x, UInt(11))
splitprec(x::Float64) = truncbits(x, UInt(26))

function add2{T<:Number}(u::T, v::T)
    if abs(v) > abs(u)
        return add2(v, u)
    end
    w = u + v
    w, (u-w) + v
end

add2(u, v) = _add2(promote(u, v)...)
_add2{T<:Number}(u::T, v::T) = add2(u, v)
_add2(u, v) = error("$u::$(typeof(u)) and $v::$(typeof(v)) cannot be promoted to a common type")

function add2(u, v, x)
    s_hi, s_lo = add2(u, x)
    s_hi, s_lo+v
end

function mul2(u_hi, u_lo, v::Integer)
    v == 0 && return zero(u_hi), zero(u_lo)
    nb = ceil(Int, log2(abs(v)))
    ut = truncbits(u_hi, nb)
    add2(ut*v, ((u_hi-ut) + u_lo)*v)
end

function _mul2{T<:Union{Float16,Float32,Float64}}(u_hi::T, u_lo::T, v::T)
    v == 0 && return u_hi*v, u_lo*v
    uhh, uhl = splitprec(u_hi)
    vh, vl = splitprec(v)
    z = u_hi*v
    add2(z, ((uhh*vh-z) + uhh*vl + uhl*vh) + uhl*vl)
end

_mul2(u_hi, u_lo, v) = u_hi*v, u_lo*v

mul2(u_hi, u_lo, v) = _mul2(promote(u_hi, u_lo, v)...)

function div2(u_hi, u_lo, v)
    hi = u_hi/v
    w_hi, w_lo = mul2(hi, zero(hi), v)
    lo = (((u_hi - w_hi) - w_lo) + u_lo)/v
    add2(hi, lo)
end

function proddiv(T, num, den)
    @_inline_meta
    v_hi = T(num[1])
    v_hi, v_lo = _prod(v_hi, zero(T), tail(num)...)
    _div(v_hi, v_lo, den...)
end
function _prod(v_hi, v_lo, x, y...)
    @_inline_meta
    _prod(mul2(v_hi, v_lo, x)..., y...)
end
_prod(v_hi, v_lo) = v_hi, v_lo
function _div(v_hi, v_lo, x, y...)
    @_inline_meta
    _div(div2(v_hi, v_lo, x)..., y...)
end
_div(v_hi, v_lo) = v_hi, v_lo

isbetween(a, x, b) = a <= x <= b || b <= x <= a
