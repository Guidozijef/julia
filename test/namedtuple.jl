# This file is a part of Julia. License is MIT: https://julialang.org/license

@test_throws TypeError NamedTuple{1,Tuple{}}
@test_throws TypeError NamedTuple{(),1}
@test_throws TypeError NamedTuple{(:a,1),Tuple{Int}}
@test_throws ErrorException NamedTuple{(:a,:b),Tuple{Int}}
@test_throws ErrorException NamedTuple{(:a,:b),Tuple{Int,Vararg{Int}}}

@test (a=1,).a == 1
@test (a=2,)[1] == 2
@test (a=3,)[:a] == 3
@test (x=4, y=5, z=6).y == 5
@test (x=4, y=5, z=6).z == 6

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

@test convert(NamedTuple{(:a,:b),Tuple{Int8,Int16}}, (1,2)) === (a=Int8(1), b=Int16(2))
@test convert(NamedTuple{(:a,:b),Tuple{Int8,Int16}}, (x=3,y=4)) === (a=Int8(3), b=Int16(4))

@test eltype((a=[1,2], b=[3,4])) === Vector{Int}

@test Tuple((a=[1,2], b=[3,4])) == ([1,2], [3,4])
@test Tuple(NamedTuple()) === ()
@test Tuple((x=4, y=5, z=6)) == (4,5,6)

@test isless((a=1,b=2), (a=1,b=3))
@test_broken isless((a=1,), (a=1,b=2))
@test !isless((a=1,b=2), (a=1,b=2))
@test !isless((a=2,b=1), (a=1,b=2))
@test_throws MethodError isless((a=1,), (x=2,))

@test map(-, (x=1, y=2)) == (x=-1, y=-2)
@test map(+, (x=1, y=2), (x=10, y=20)) == (x=11, y=22)
@test map(string, (x=1, y=2)) == (x="1", y="2")
@test map(round, (x=1//3, y=Int), (x=3, y=2//3)) == (x=0.333, y=1)
