module aHDPs

using BioSequences
using CircularArrays
using SnoopPrecompile

export findcores


const FORMULA = CircularVector(collect("XPHHPPHHPXHPPHHPPH"))
const ALPHABET = Set(aa for aa in aa"ARNDCQEGHILKMFPSTWYV")
const AA_IDS = Dict(a => i for (i, a) in enumerate(ALPHABET))

const POLAR = Set(aa for aa in aa"EDKRHQNTSAG")
const HYDROPHOBIC = Set(aa for aa in aa"VMCILFWYAG")


"""
    findcores(sequence; window_size = min(length(sequence), $(length(FORMULA))))

Returns ranges where the given `sequences` matches the `aHDPs.FORMULA`. By
default, we use the entire formula to match. The default `window_size` is the
smaller of the length of the target and the length of the formula.

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
function findcores(target::AASeq; window_size = min(length(target), length(FORMULA)))
    pattern = generate_pattern(window_size)
    scores = conv(pattern, onehot(target))

    match_indices = findall(≈(window_size), scores)
    ranges = map(Tuple.(match_indices)) do (i,)
        i:i+window_size-1
    end

    isempty(ranges) ? ranges : merge_adjacent(unique(ranges))
end


"""
    onehot(sequence)
    onehot(aas)
    onehot(aa)

Computes a "one-hot" encoding representation of a given amino acid `sequence`.
Optionally, a target value type `T` can be supplied.
"""
function onehot(seq::AASeq)
    enc = falses(length(ALPHABET), length(seq))
    for (i, aa) in enumerate(seq)
        encode!(view(enc, :, i), aa)
    end
    enc
end

function onehot(aas::Union{AbstractArray{AminoAcid}, AbstractSet{AminoAcid}})
    enc = falses(length(ALPHABET))
    for aa in aas
        encode!(enc, aa)
    end
    enc
end

function onehot(aa::AminoAcid)
    enc = falses(length(ALPHABET))
    encode!(enc, aa)
    enc
end

function encode!(dst::AbstractVector{Bool}, aas::Union{AbstractSet{AminoAcid},AASeq})
    for a in aas
        encode!(dst, a)
    end
end

function encode!(dst::AbstractVector{Bool}, aa::AminoAcid)
    if aa == AA_X
        dst[:] .= true
    else
        dst[AA_IDS[aa]] = true
    end
end


"""
    generate_pattern(window_size)

Computes a one-hot representation of the `aHDPs.FORMULA`. Returns a 3D array
where the third axis stores each pattern window.
"""
function generate_pattern(window_size)
    enc = trues(length(ALPHABET), window_size, length(FORMULA))

    ps = onehot(POLAR)
    hs = onehot(HYDROPHOBIC)

    for j in axes(enc, 3)
        for (i, x) in enumerate(FORMULA[j:j+window_size-1])
            if x == 'P'
                enc[:, i, j] = ps
            elseif x == 'H'
                enc[:, i, j] = hs
            end
        end
    end

    enc
end

"""
    conv(pattern, sequence)

Performs a 1D boolean, convolution of one-hot encoded `pattern` and `sequence`
along all windows of the pattern.
"""
function conv(pattern::BitArray, seq::BitArray) where T
    n, w = size(seq, 2), size(pattern, 2)
    res = zeros(Int, max(1, n - w + 1), size(pattern, 3))

    for i in axes(res, 1)
        @inbounds res[i, :] = count(any(pattern .& view(seq, :, i:i+w-1); dims = 1); dims = 2)
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


@precompile_setup begin
    ll37 = aa"LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTES"

    @precompile_all_calls begin
        core_idx = findcores(ll37)
        core_idx = findcores(ll37; window_size = 12)
    end
end

end # module aHDPs
