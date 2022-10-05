# aHDPs

Yount et al., PNAS, 2019

## Example
```julia
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
