include("iola.jl")
using DSP

fs = 44100
T = 10
t = 0:1/fs:T
f = 440
sig = sin.(2π*f*t)

len = 1152

diffs = iola.analyze(sig, (s) -> 20*log10.(abs.(iola.transforms.STMDCT(s, Windows.cosine(len)))),
                    fs, div(len,2))
