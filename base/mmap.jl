# This file is a part of Julia. License is MIT: http://julialang.org/license

module Mmap

export mmap_array, mmap_bitarray

type Stream # <: IO
    viewhandle::Ptr{Void}
    mmaphandle::Ptr{Void} # only needed on windows
    len::Int
    offset::Int
    name

    Stream(v::Ptr{Void},m::Ptr{Void},len::Int,o::Int,n) = new(v,m,len,o,n)

    @unix_only begin
    function Stream(len::Integer, prot::Integer, flags::Integer, fd, offset::Integer, name::AbstractString="")
        pagesize::Int = ccall(:jl_getpagesize, Clong, ())
        # Check that none of the computations will overflow
        if len < 0
            throw(ArgumentError("requested size must be ≥ 0, got $len"))
        end
        if len > typemax(Int)-pagesize
            throw(ArgumentError("requested size must be ≤ $(typemax(Int)-pagesize), got $len"))
        end
        # Set the offset to a page boundary
        offset_page::FileOffset = div(offset,pagesize)*pagesize
        len_page::Int = (offset-offset_page) + len
        # Mmap the file
        p = ccall(:jl_mmap, Ptr{Void}, (Ptr{Void}, Csize_t, Cint, Cint, Cint, FileOffset), C_NULL, len_page, prot, flags, fd, offset_page)
        systemerror("memory mapping failed", reinterpret(Int,p) == -1)
        # Also return a pointer that compensates for any adjustment in the offset
        stream = new(p,C_NULL,len_page,Int(offset-offset_page),name)
        finalizer(stream,close)
        return stream
    end
    end # @unix_only

    @windows_only begin
        function Stream(len::Integer, hdl, readonly::Bool, create::Bool, offset::Integer, name)
            granularity::Int = ccall(:jl_getallocationgranularity, Clong, ())
            if len < 0
                throw(ArgumentError("requested size must be ≥ 0, got $len"))
            end
            if len > typemax(Int)-granularity
                throw(ArgumentError("requested size must be ≤ $(typemax(Int)-granularity), got $len"))
            end
            # Set the offset to a page boundary
            offset_page::FileOffset = div(offset, granularity)*granularity
            szfile = convert(Csize_t, len + offset)
            szarray = szfile - convert(Csize_t, offset_page)
            access = readonly ? 4 : 2
            if create
                flprotect = readonly ? 0x02 : 0x04
                mmaphandle = ccall(:CreateFileMappingW, stdcall, Ptr{Void}, (Cptrdiff_t, Ptr{Void}, Cint, Cint, Cint, Cwstring),
                    hdl, C_NULL, flprotect, szfile>>32, szfile&typemax(UInt32), name)
            else
                mmaphandle = ccall(:OpenFileMappingW, stdcall, Ptr{Void}, (Cint, Cint, Cwstring),
                    access, true, name)
            end
            if mmaphandle == C_NULL
                error("could not create file mapping: $(Base.FormatMessage())")
            end
            viewhandle = ccall(:MapViewOfFile, stdcall, Ptr{Void}, (Ptr{Void}, Cint, Cint, Cint, Csize_t),
                mmaphandle, access, offset_page>>32, offset_page&typemax(UInt32), szarray)
            if viewhandle == C_NULL
                error("could not create mapping view: $(Base.FormatMessage())")
            end
            stream = new(viewhandle,mmaphandle,szarray,Int(offset-offset_page),name)
            finalizer(stream,close)
            return stream
        end
    end # @windows_only
end

Base.show(io::IO, s::Stream) = print(io, "Mmap.Stream(", s.name, ")")

function Base.close(m::Stream)
    @unix_only systemerror("munmap", ccall(:munmap,Cint,(Ptr{Void},Int),m.viewhandle,m.len) != 0)
    @windows_only begin
        status = ccall(:UnmapViewOfFile, stdcall, Cint, (Ptr{Void},), m.viewhandle)!=0
        status |= ccall(:CloseHandle, stdcall, Cint, (Ptr{Void},), m.mmaphandle)!=0
        if !status
            error("could not unmap view: $(Base.FormatMessage())")
        end
    end
end

# msync flags
const MS_ASYNC = 1
const MS_INVALIDATE = 2
const MS_SYNC = 4

sync!(m::Stream, flags::Integer=MS_SYNC) = sync!(m.viewhandle, m.len, flags)
sync!{T}(A::Array{T}) = sync!(pointer(A), length(A)*sizeof(T))
sync!(B::BitArray)    = sync!(pointer(B.chunks), length(B.chunks)*sizeof(UInt64))
function sync!(p::Ptr, len::Integer, flags::Integer=MS_SYNC)
    @unix_only systemerror("msync", ccall(:msync, Cint, (Ptr{Void}, Csize_t, Cint), p, len, flags) != 0)
    @windows_only begin
        status = ccall(:FlushViewOfFile, stdcall, Cint, (Ptr{Void}, Csize_t), p, len)!=0
        status || error("could not msync: $(Base.FormatMessage())")
    end
end

@unix_only begin

const SEEK_SET = Cint(0)
const SEEK_CUR = Cint(1)
const SEEK_END = Cint(2)
# Before mapping, grow the file to sufficient size
# (Required if you're going to write to a new memory-mapped file)
#
# Note: a few mappable streams do not support lseek. When Julia
# supports structures in ccall, switch to fstat.
function grow!(len::Integer, fd::Integer, offset::FileOffset)
    # Save current file position so we can restore it later
    cpos = ccall(:jl_lseek, FileOffset, (Cint, FileOffset, Cint), fd, 0, SEEK_CUR)
    systemerror("lseek", cpos < 0)
    filelen = ccall(:jl_lseek, FileOffset, (Cint, FileOffset, Cint), fd, 0, SEEK_END)
    systemerror("lseek", filelen < 0)
    if (filelen < offset + len)
        systemerror("pwrite", ccall(:jl_pwrite, Cssize_t, (Cint, Ptr{Void}, UInt, FileOffset), fd, Int8[0], 1, offset + len - 1) < 1)
    end
    cpos = ccall(:jl_lseek, FileOffset, (Cint, FileOffset, Cint), fd, cpos, SEEK_SET)
    systemerror("lseek", cpos < 0)
    return
end

const PROT_READ  = Cint(1)
const PROT_WRITE = Cint(2)
const MAP_SHARED = Cint(1)
const F_GETFL    = Cint(3)
# Determine a stream's read/write mode, and return prot & flags
# appropriate for mmap
# We could use isreadonly here, but it's worth checking that it's readable too
function settings(s::IO)
    mode = ccall(:fcntl,Cint,(Cint,Cint),fd(s),F_GETFL)
    systemerror("fcntl F_GETFL", mode == -1)
    mode = mode & 3
    if mode == 0
        prot = PROT_READ
    elseif mode == 1
        prot = PROT_WRITE
    else
        prot = PROT_READ | PROT_WRITE
    end
    if prot & PROT_READ == 0
        throw(ArgumentError("mmap requires read permissions on the file (choose r+)"))
    end
    flags = MAP_SHARED
    return prot, flags, (prot & PROT_WRITE) > 0
end

end # @unix_only

# @windows_only
type SharedMemSpec
    name :: AbstractString
    readonly :: Bool
    create :: Bool
end

# Mmapped-array constructor
function mmap_array{T,N}(::Type{T}, dims::NTuple{N,Integer}, s::Union(IO,SharedMemSpec), offset::FileOffset=position(s); grow::Bool=true)
    applicable(fd,s) || throw(ArgumentError("fd() is not defined for $s"))
    n = 1
    for (i, d) in enumerate(dims)
        if d < 0
            throw(ArgumentError("dimension size must be ≥ 0, got $d size for dimension $i"))
        end
        n *= d
    end
    len = prod(dims)*sizeof(T)
    @unix_only begin
        isa(s,SharedMemSpec) && throw(ArgumentError("$s is not valid on unix"))
        prot, flags, iswrite = settings(s)
        if iswrite && grow
            grow!(len, fd(s), offset)
        end
        mm = Stream(len, prot, flags, fd(s), offset, isdefined(s,:name) ? s.name : "")
    end # @unix_only

    @windows_only begin
        if isa(s,IO)
            hdl::Int = Base._get_osfhandle(RawFD(fd(s))).handle
            hdl == -1 && throw(ArgumentError("could not get handle for file to map: $(Base.FormatMessage())"))
            name = Ptr{Cwchar_t}(C_NULL)
            ro = isreadonly(s)
            create = true
        else
            # shared memory
            hdl = -1
            name = utf16(s.name)
            ro = s.readonly
            create = s.create
        end
        mm = Stream(len, hdl, ro, create, offset, name)
    end # @windows_only

    A = pointer_to_array(convert(Ptr{T}, UInt(mm.viewhandle)+mm.offset), dims)
    return A
end

# Mmapped-bitarray constructor
mmap_bitarray{N}(::Type{Bool}, dims::NTuple{N,Integer}, s::IOStream, offset::FileOffset) =
    mmap_bitarray(dims, s, offset)
mmap_bitarray{N}(::Type{Bool}, dims::NTuple{N,Integer}, s::IOStream) = mmap_bitarray(dims, s, position(s))

function mmap_bitarray{N}(dims::NTuple{N,Integer}, s::IOStream, offset::FileOffset=position(s))
    iswrite = !isreadonly(s)
    n = prod(dims)
    nc = Base.num_bit_chunks(n)
    chunks = mmap_array(UInt64, (nc,), s, offset)
    if iswrite
        chunks[end] &= Base._msk_end(n)
    else
        if chunks[end] != chunks[end] & Base._msk_end(n)
            throw(ArgumentError("the given file does not contain a valid BitArray of size $(join(dims, 'x')) (open with \"r+\" mode to override)"))
        end
    end
    B = BitArray{N}(ntuple(N,i->0)...)
    B.chunks = chunks
    B.len = n
    if N != 1
        B.dims = dims
    end
    return B
end

end # module