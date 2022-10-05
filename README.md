# aHDPs

Provides a basic tool for identifying αHDPs as described in [Yount et al., PNAS,
2019](https://www.pnas.org/doi/10.1073/pnas.1819250116) by applying the
alpha-core formula described in that paper.

## Example
```julia
using aHDPs
using BioSequences


ll37 = aa"LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTES"
core_ranges = findcores(ll37; window_size = 12)

core_seqs = map(core_ranges) do idx
    ll37[idx]
end
```
