using AlphaCore
using BioSequences
using CairoMakie


lp_007321 = aa"MYIRKASFILLITLVILSQTMASQLFDELADDDSNAREESWGALTKLRKKVAHRTCEQTKKCQTTGGHETCANYICCHIVVFGKNTTPSCTVICYFKLDTKIACTCISK"
ll37 = aa"LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTES"

window_size = 12
target = lp_007321

core_idx = findcores(target; window_size = window_size)
core_seqs = map(core_idx) do idx
    target[idx]
end

