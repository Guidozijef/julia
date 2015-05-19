# This file is a part of Julia. License is MIT: http://julialang.org/license

module Mmap

# platform-specific mmap utilities
pagesize() = Int(@unix ? ccall(:jl_getpagesize, Clong, ()) : ccall(:jl_getallocationgranularity, Clong, ()))

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
function settings(s)
    mode = ccall(:fcntl,Cint,(Cint,Cint),s,F_GETFL)
    systemerror("fcntl F_GETFL", mode == -1)
    mode = mode & 3
    prot = mode == 0 ? PROT_READ : mode == 1 ? PROT_WRITE : PROT_READ | PROT_WRITE
    if prot & PROT_READ == 0
        throw(ArgumentError("mmap requires read permissions on the file (choose r+)"))
    end
    flags = MAP_SHARED
    return prot, flags, (prot & PROT_WRITE) > 0
end
end # @unix_only

@windows_only begin
type SharedMemSpec <: IO
    name::AbstractString
    readonly::Bool
    create::Bool
end

Base.fd(sh::SharedMemSpec) = -2 # -1 == INVALID_HANDLE_VALUE

const INVALID_HANDLE_VALUE = -1
gethandle(io::SharedMemSpec) = INVALID_HANDLE_VALUE
function gethandle(io::IO)
    handle = Base._get_osfhandle(RawFD(fd(io))).handle
    systemerror("could not get handle for file to map: $(Base.FormatMessage())", handle == -1)
    return Int(handle)
end

settings(sh::SharedMemSpec) = utf16(sh.name), sh.readonly, sh.create
settings(io::IO) = Ptr{Cwchar_t}(C_NULL), isreadonly(io), true

# Memory mapped file constants
const PAGE_READONLY          = UInt32(0x02)
const PAGE_READWRITE         = UInt32(0x04)
const PAGE_WRITECOPY         = UInt32(0x08)

const PAGE_EXECUTE_READ      = UInt32(0x20)
const PAGE_EXECUTE_READWRITE = UInt32(0x40)
const PAGE_EXECUTE_WRITECOPY = UInt32(0x80)

const FILE_MAP_COPY          = UInt32(0x01)
const FILE_MAP_WRITE         = UInt32(0x02)
const FILE_MAP_READ          = UInt32(0x04)
const FILE_MAP_EXECUTE       = UInt32(0x20)
end # @windows_only

# core impelementation of mmap
type Stream <: IO
    ptr::Ptr{Void}    # pointer to mmapped-memory
    handle::Ptr{Void} # only needed on windows for file mapping object
    len::Int          # amount of memory mapped
    offset::FileOffset
    pos::Int64
    name

    function Stream{T<:IO}(io::T, len::Integer=filesize(io), offset::Integer=0; finalize::Bool=true, grow::Bool=true)
        # check inputs
        isopen(io) || throw(ArgumentError("$io must be open to mmap"))
        applicable(fd,io) || throw(ArgumentError("method `fd(::$T)` doesn't exist, unable to mmap $io"))
        # Check that none of the computations will overflow
        len >= 0 || throw(ArgumentError("requested size must be ≥ 0, got $len"))
        ps = pagesize()
        len < typemax(Int)-ps || throw(ArgumentError("requested size must be < $(typemax(Int)-ps), got $len"))

        offset >= 0 || throw(ArgumentError("requested offset must be ≥ 0, got $offset"))
        offset > filesize(io) && throw(ArgumentError("requested offset is beyond size of file: $offset > $(filesize(io))"))

        # Set the offset to a page boundary
        offset_page::FileOffset = div(offset,ps)*ps
        len_page::Int = (offset-offset_page) + len

        # platform-specific internals
         @unix_only begin
            file_desc = fd(io)
            prot, flags, iswrite = settings(file_desc)
            iswrite && grow && grow!(len, file_desc, offset)
            # mmap the file
            ptr = ccall(:jl_mmap, Ptr{Void}, (Ptr{Void}, Csize_t, Cint, Cint, Cint, FileOffset), C_NULL, len_page, prot, flags, file_desc, offset_page)
            systemerror("memory mapping failed", reinterpret(Int,ptr) == -1)
            handle = C_NULL
        end # @unix_only

        @windows_only begin
            hdl::Int = gethandle(io)
            name, readonly, create = settings(io)
            szfile = convert(Csize_t, len + offset)
            handle = create ? ccall(:CreateFileMappingW, stdcall, Ptr{Void}, (Cptrdiff_t, Ptr{Void}, Cint, Cint, Cint, Cwstring),
                                    hdl, C_NULL, readonly ? PAGE_READONLY : PAGE_READWRITE, szfile>>32, szfile&typemax(UInt32), name) :
                                  ccall(:OpenFileMappingW, stdcall, Ptr{Void}, (Cint, Cint, Cwstring),
                                    readonly ? FILE_MAP_READ : FILE_MAP_WRITE, true, name)
            handle == C_NULL && error("could not create file mapping: $(Base.FormatMessage())")
            ptr = ccall(:MapViewOfFile, stdcall, Ptr{Void}, (Ptr{Void}, Cint, Cint, Cint, Csize_t),
                            handle, readonly ? FILE_MAP_READ : FILE_MAP_WRITE, offset_page>>32, offset_page&typemax(UInt32), szfile - convert(Csize_t, offset_page))
            ptr == C_NULL && error("could not create mapping view: $(Base.FormatMessage())")
        end # @windows_only
        offset = FileOffset(offset-offset_page)
        stream = new(ptr,handle,len_page,offset,offset+1,isdefined(io,:name) ? io.name : "")
        finalize && finalizer(stream,close)
        return stream
    end
end

Stream(file::AbstractString, len::Integer=filesize(file), offset::Integer=0; finalize::Bool=true, grow::Bool=true) =
    open(io->Stream(io, len, offset; finalize=finalize, grow=grow), file, "r+")
Base.show(io::IO, s::Stream) = print(io, "Mmap.Stream(", s.name, ")")

function Base.close(m::Stream)
    @unix_only systemerror("munmap", ccall(:munmap,Cint,(Ptr{Void},Int),m.ptr,m.len) != 0)
    @windows_only begin
        status = ccall(:UnmapViewOfFile, stdcall, Cint, (Ptr{Void},), m.ptr)!=0
        status |= ccall(:CloseHandle, stdcall, Cint, (Ptr{Void},), m.handle)!=0
        status || error("could not unmap view: $(Base.FormatMessage())")
    end
    return
end

function Base.read(from::Stream, ::Type{UInt8})
    from.pos > from.len && throw(EOFError())
    byte = unsafe_load(convert(Ptr{UInt8},from.ptr), from.pos)
    from.pos += 1
    return byte
end

function Base.peek(from::Stream)
    from.pos > from.len && throw(EOFError())
    return unsafe_load(from.ptr, from.pos)
end

# msync flags for unix
const MS_ASYNC = 1
const MS_INVALIDATE = 2
const MS_SYNC = 4

sync!(m::Stream, flags::Integer=MS_SYNC) = sync!(m.ptr, m.len, flags)
sync!{T}(A::Array{T}) = sync!(pointer(A), length(A)*sizeof(T))
sync!(B::BitArray)    = sync!(pointer(B.chunks), length(B.chunks)*sizeof(UInt64))
function sync!(p::Ptr, len::Integer, flags::Integer=MS_SYNC)
    @unix_only systemerror("msync", ccall(:msync, Cint, (Ptr{Void}, Csize_t, Cint), p, len, flags) != 0)
    @windows_only systemerror("could not FlushViewOfFile: $(Base.FormatMessage())",
                    ccall(:FlushViewOfFile, stdcall, Cint, (Ptr{Void}, Csize_t), p, len) == 0)
end

# Mmapped-array constructor
function Array{T,N}(::Type{T}, dims::NTuple{N,Integer}, io::IO, offset=position(io); grow::Bool=true)
    n = 1
    for (i, d) in enumerate(dims)
        d < 0 && throw(ArgumentError("dimension size must be ≥ 0, got $d size for dimension $i"))
        n *= d
    end
    mm = Stream(io, n*sizeof(T), offset; finalize=false, grow=grow)
    A = pointer_to_array(convert(Ptr{T}, UInt(mm.ptr)+mm.offset), dims)
    finalizer(A, x->close(mm))
    return A
end
Array{T,N}(::Type{T}, dims::NTuple{N,Integer}, file::AbstractString, offset::Integer=0; grow::Bool=true) = open(io->Array(T,dims,io,offset;grow=grow),file, "r+")
Array(io::IO, len::Integer=filesize(io), offset::Integer=0; grow::Bool=true) = Array(UInt8, (len,), io, offset; grow=grow)
Array(file::AbstractString, len::Integer=filesize(file), offset::Integer=0; grow::Bool=true) = open(io->Array(UInt8,(len,),io,offset;grow=grow),file, "r+")

# Mmapped-bitarray constructor
function BitArray{N}(dims::NTuple{N,Integer}, io::IO, offset::Integer=position(io); grow::Bool=true)
    n = prod(dims)
    nc = Base.num_bit_chunks(n)
    chunks = Mmap.Array(UInt64, (nc,), io, offset)
    if !isreadonly(io)
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
BitArray{N}(dims::NTuple{N,Integer}, file::AbstractString, offset::Integer=0; grow::Bool=true) = open(io->BitArray(dims,io,offset;grow=grow),file, "r+")

end # module