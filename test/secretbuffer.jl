using Base: SecretBuffer, SecretBuffer!, shred!, isshredded
using Test, Random

@testset "SecretBuffer" begin
    @testset "original unmodified" begin
        str = "foobar"
        secret = SecretBuffer(str)

        @test read(secret, String) == str
        seekstart(secret)

        @test shred!(secret) === secret
        @test read(secret, String) == ""
        @test str == "foobar"
    end

    @testset "finalizer" begin
        v = UInt8[1, 2]
        secret_a = SecretBuffer!(v)
        secret_b = secret_a

        secret_a = nothing
        GC.gc()

        @test all(iszero, v)
        @test !isshredded(secret_b)

        # TODO: ideally we'd test that the finalizer warns from GC.gc(), but that is harder
        @test_logs (:warn, r".*SecretBuffer was `shred!`ed by the GC.*") finalize(secret_b)
        @test isshredded(secret_b)
        secret_b = nothing
        GC.gc()
    end

    @testset "initializers" begin
        s1 = SecretBuffer("setec astronomy")
        data2 = [0x73, 0x65, 0x74, 0x65, 0x63, 0x20, 0x61, 0x73, 0x74, 0x72, 0x6f, 0x6e, 0x6f, 0x6d, 0x79]
        s2 = SecretBuffer!(data2)
        @test all(==(0x0), data2)
        @test s1 == s2

        ptr3 = Base.unsafe_convert(Cstring, "setec astronomy")
        s3 = Base.unsafe_SecretBuffer!(ptr3)
        @test Base.unsafe_string(ptr3) == ""
        @test s1 == s2 == s3

        shred!(s1); shred!(s2); shred!(s3)
    end

    @testset "copiers" begin
        s1 = SecretBuffer()
        write(s1, "hello world")
        seekstart(s1)

        s2 = copy(s1)
        write(s2, 'c')
        seekstart(s2)

        @test read(s1) == codeunits("hello world")
        @test read(s2) == codeunits("cello world")

        shred!(s1)
        @test isshredded(s1)
        @test !isshredded(s2)
        shred!(s2)

        # Copying into a bigger destination
        s3 = SecretBuffer()
        s4 = SecretBuffer()
        write(s3, "original")
        seekstart(s3)
        write(s4, randstring(1234))
        s4data = s4.data
        copy!(s4, s3)
        @test s3.data == s4.data
        @test read(s3) == read(s4) == codeunits("original")
        @test all(iszero, s4data)
        shred!(s3); shred!(s4)

        # Copying into a smaller destination
        s5 = SecretBuffer()
        s6 = SecretBuffer("sekrit")
        str = randstring(321)
        write(s5, str)
        seekstart(s5)
        copy!(s6, s5)
        @test read(s5) == read(s6) == codeunits(str)
        shred!(s5); shred!(s6)
    end
end
