module aHDPs

using BioSequences
using CircularArrays

export findcores


const FORMULA = CircularVector(collect("XPHHPPHHPXHPPHHPPH"))
const ALPHABET = Set(aa for aa in aa"ARNDCQEGHILKMFPSTWYV")
const AA_IDS = Dict(a => i for (i, a) in enumerate(ALPHABET))

const POLAR = Set(aa for aa in aa"EDKRHQNTSAG")
const HYDROPHOBIC = Set(aa for aa in aa"VMCILFWYAG")


"""
    findcores(sequence; window_size = $(length(FORMULA)))

Returns ranges where the given `sequences` matches the `aHDPs.FORMULA`. By
default, we use the entire formula to match. The default `window_size` of
18, the length of the formula.

# Example
```julia-repl
julia> ll37 = aa"LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTES"
37aa Amino Acid Sequence:
LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTES

julia> window_size = 12
12

julia> core_idx = findcores(ll37; window_size = window_size)
1-element Vector{UnitRange{Int64}}:
 11:32

julia> core_seqs = map(core_idx) do idx
           ll37[idx]
       end
1-element Vector{LongAA}:
 EKIGKEFKRIVQRIKDFLRNLV
```
"""
function findcores(target::AASeq; window_size = length(FORMULA))
    pattern = generate_pattern(window_size)
    scores = conv(pattern, onehot(target))
    
    match_indices = findall(≈(1), scores)
    ranges = map(Tuple.(match_indices)) do (i,)
        i:i+window_size-1
    end

    merge_adjacent(unique(ranges))
end


"""
    onehot([T], sequence)

Computes a "one-hot" encoding representation of a given amino acid `sequence`.
Optionally, a target value type `T` can be supplied.
"""
function onehot(::Type{T}, seq::AASeq) where T
    enc = zeros(T, length(ALPHABET), length(seq))
    for (i, a) in enumerate(seq)
        enc[:, i] = encode(T, a)
    end
    enc
end

onehot(seq::AASeq) = onehot(Float64, seq)


"""
    encode([T], sequence)

Encodes an amino acid or a collection of amino acids as a one-hot vector.
"""
function encode(::Type{T}, aa::AminoAcid) where T
    enc = zeros(T, length(ALPHABET))
    encode!(enc, aa)
    enc
end

function encode(::Type{T}, aas::Union{AbstractSet{AminoAcid},AASeq}) where T
    enc = zeros(T, length(ALPHABET))
    for a in aas
        encode!(enc, a)
    end
    enc
end

encode(x) = encode(Float64, x)
encode!(dst::AbstractArray{T}, aa::AminoAcid) where T = dst[AA_IDS[aa]] = one(T)

"""
    generate_pattern(window_size)

Computes a one-hot representation of the `aHDPs.FORMULA`. Returns a 3D array
where the third axis stores each pattern window.
"""
function generate_pattern(window_size)
    enc = ones(length(ALPHABET), window_size, length(FORMULA))

    ps = encode(POLAR)
    hs = encode(HYDROPHOBIC)

    for j in axes(enc, 3)
        for (i, x) in enumerate(FORMULA[j:j+window_size-1])
            if x == 'P'
                enc[:, i, j] = ps
            elseif x == 'H'
                enc[:, i, j] = hs
            end
        end
    end

    enc ./ window_size
end

"""
    conv(pattern, sequence)

Performs a 1D convolution of one-hot encoded `pattern` and `sequence` along all
windows of the pattern.
"""
function conv(pattern::AbstractArray{T}, seq::AbstractMatrix{T}) where T
    n, w = size(seq, 2), size(pattern, 2)
    res = similar(pattern, n - w, size(pattern, 3))

    for i in axes(res, 1)
        @inbounds res[i, :] = sum(pattern .* view(seq, :, i:i+w-1); dims = (1, 2))
    end

    res
end

"""
    merge_adjacent(ranges)

Merges adjacent `UnitRange`s.
"""
function merge_adjacent(ranges)
    adjacent_segments, start = [], 1
    sort!(ranges)

    for gap in findall(!=(1), diff(first.(ranges)))
        push!(adjacent_segments, start:gap)
        start = gap + 1
    end

    push!(adjacent_segments, start:length(ranges))

    map(adjacent_segments) do idx
        (; start), (; stop) = ranges[idx][[1, end]]
        start:stop
    end
end

end # module aHDPs
