# This file is a part of Julia. License is MIT: https://julialang.org/license

@test_throws TypeError NamedTuple{1,Tuple{}}
@test_throws TypeError NamedTuple{(),1}
@test_throws TypeError NamedTuple{(:a,1),Tuple{Int}}
@test_throws ErrorException NamedTuple{(:a,:b),Tuple{Int}}
@test_throws ErrorException NamedTuple{(:a,:b),Tuple{Int,Vararg{Int}}}
@test_throws ErrorException NamedTuple{(:a,),Union{Tuple{Int},Tuple{String}}}

@test (a=1,).a == 1
@test (a=2,)[1] == 2
@test (a=3,)[:a] == 3
@test (x=4, y=5, z=6).y == 5
@test (x=4, y=5, z=6).z == 6
@test_throws ErrorException (x=4, y=5, z=6).a
@test_throws BoundsError (a=2,)[0]
@test_throws BoundsError (a=2,)[2]

@test length(NamedTuple()) == 0
@test length((a=1,)) == 1
@test length((a=1, b=0)) == 2

@test (a=1,b=2) === (a=1,b=2)
@test (a=1,b=2) !== (b=1,a=2)

@test (a=1,b=2) == (a=1,b=2)
@test (a=1,b=2) != (b=1,a=2)

@test string((a=1,)) == "(a = 1,)"
@test string((name="", day=:today)) == "(name = \"\", day = :today)"
@test string(NamedTuple()) == "NamedTuple()"

@test hash((a = 1, b = "hello")) == hash(convert(NamedTuple{(:a,:b),Tuple{Int,String}}, (1, "hello")))
@test hash((a = 1, b = "hello")) != hash(convert(NamedTuple{(:a,:c),Tuple{Int,String}}, (1, "hello")))
@test hash((a = 1, b = "hello")) != hash(convert(NamedTuple{(:a,:b),Tuple{Int,String}}, (1, "helo")))

@test (x=4, y=5, z=6)[[1,3]] == (x=4, z=6)
@test (x=4, y=5, z=6)[[:z,:y]] == (z=6, y=5)
@test_throws BoundsError (x=4, y=5, z=6)[[1,5]]

@test convert(NamedTuple{(:a,:b),Tuple{Int8,Int16}}, (1,2)) === (a=Int8(1), b=Int16(2))
@test convert(NamedTuple{(:a,:b),Tuple{Int8,Int16}}, (x=3,y=4)) === (a=Int8(3), b=Int16(4))

@test eltype((a=[1,2], b=[3,4])) === Vector{Int}

@test Tuple((a=[1,2], b=[3,4])) == ([1,2], [3,4])
@test Tuple(NamedTuple()) === ()
@test Tuple((x=4, y=5, z=6)) == (4,5,6)
@test convert(Tuple, (x=4, y=5, z=6)) == (4,5,6)
@test collect((x=4, y=5, z=6)) == [4,5,6]

@test isless((a=1,b=2), (a=1,b=3))
@test_broken isless((a=1,), (a=1,b=2))
@test !isless((a=1,b=2), (a=1,b=2))
@test !isless((a=2,b=1), (a=1,b=2))
@test_throws MethodError isless((a=1,), (x=2,))

@test map(-, (x=1, y=2)) == (x=-1, y=-2)
@test map(+, (x=1, y=2), (x=10, y=20)) == (x=11, y=22)
@test map(string, (x=1, y=2)) == (x="1", y="2")
@test map(round, (x=1//3, y=Int), (x=3, y=2//3)) == (x=0.333, y=1)

@test merge((a=1, b=2), (a=10,)) == (a=10, b=2)
@test merge((a=1, b=2), (a=10, z=20)) == (a=10, b=2, z=20)
@test merge((a=1, b=2), (z=20,)) == (a=1, b=2, z=20)

@test Base.structdiff((a=1, b=2), (b=3,)) == (a=1,)
@test Base.structdiff((a=1, b=2, z=20), (b=3,)) == (a=1, z=20)
@test Base.structdiff((a=1, b=2, z=20), (b=3, q=20, z=1)) == (a=1,)
@test Base.structdiff((a=1, b=2, z=20), (b=3, q=20, z=1, a=0)) == NamedTuple()
@test Base.structdiff((a=1, b=2, z=20), NamedTuple{(:b,)}) == (a=1, z=20)

# splatting and implicit naming syntax

let d = [:a=>1, :b=>2, :c=>3]   # use an array to preserve order
    x = 10
    t = (x=1, y=20)
    @test (;d...) == (a=1, b=2, c=3)
    @test (;d..., a=10) == (a=10, b=2, c=3)
    @test (;d..., x, t.y) == (a=1, b=2, c=3, x=10, y=20)
    @test (;d..., :z=>20) == (a=1, b=2, c=3, z=20)
    @test (;a=10, d..., :c=>30) == (a=1, b=2, c=30)
    @test (;a=0, b=0, z=1, d..., x=4, y=5) == (a=1, b=2, z=1, c=3, x=4, y=5)
    @test (;t.x, t.y) == (x=1, y=20)
    @test (b=1, t.y) == (b=1, y=20)
    y = (w=30, z=40)
    @test (;t..., y...) == (x=1, y=20, w=30, z=40)
    @test (;t..., y=0, y...) == (x=1, y=0, w=30, z=40)
end

# syntax errors

@test expand(Main, parse("(a=1, 0)")) == Expr(:error, "invalid named tuple element \"0\"")
@test expand(Main, parse("(a=1, f(x))")) == Expr(:error, "invalid named tuple element \"f(x)\"")
@test expand(Main, parse("(; f(x))")) == Expr(:error, "invalid named tuple element \"f(x)\"")
@test expand(Main, parse("(;1=0)")) == Expr(:error, "invalid named tuple field name \"1\"")
@test expand(Main, parse("(a=1,a=2)")) == Expr(:error, "field name \"a\" repeated in named tuple")
@test expand(Main, parse("(a=1,b=0,a=2)")) == Expr(:error, "field name \"a\" repeated in named tuple")
@test expand(Main, parse("(c=1,a=1,b=0,a=2)")) == Expr(:error, "field name \"a\" repeated in named tuple")
@test expand(Main, parse("(c=1,a=1,b=0,d.a)")) == Expr(:error, "field name \"a\" repeated in named tuple")
@test expand(Main, parse("(c=1,a=1,b=0,c)")) == Expr(:error, "field name \"c\" repeated in named tuple")
@test expand(Main, parse("(;d.c,c)")) == Expr(:error, "field name \"c\" repeated in named tuple")

@test parse("(;)") == quote end
@test expand(Main, parse("(1,;2)")) == Expr(:error, "unexpected semicolon in tuple")

# inference tests

namedtuple_get_a(x) = x.a
@test Base.return_types(namedtuple_get_a, (NamedTuple,)) == Any[Any]
@test Base.return_types(namedtuple_get_a, (typeof((b=1,a="")),)) == Any[String]

namedtuple_fieldtype_a(x) = fieldtype(typeof(x), :a)
@test Base.return_types(namedtuple_fieldtype_a, (NamedTuple,)) == Any[Type]
@test Base.return_types(namedtuple_fieldtype_a, (typeof((b=1,a="")),)) == Any[Type{String}]

namedtuple_nfields(x) = nfields(x) === 0 ? 1 : ""
@test Union{Int,String} <: Base.return_types(namedtuple_nfields, (NamedTuple,))[1]

struct NamedTupleFieldType
    a::NamedTuple
end
namedtuple_field(x) = getfield(x.a, 1)
@test namedtuple_field((a=1, b=2, c=3)) == 1

# Associative tests
@test keys((a=1, b=2, c=3)) == (:a, :b, :c)
@test values((a=1, b=2, c=3)) == (1, 2, 3)
@test haskey((a=1, b=2, c=3), :a)
@test haskey((a=1, b=2, c=3), 1)
@test get((a=1, b=2, c=3), :a, 0) == 1
@test get((a=1, b=2, c=3), 1, 0) == 1
@test get(()->0, (a=1, b=2, c=3), :a) == 1
@test get(()->0, (a=1, b=2, c=3), :d) == 0