include("iola.jl")
using DSP
using WAV

len = 1152

(y, fs, nbits, opt) = wavread("audio.wav")

diffs = IOLA.analyze(y[:,1], (s) -> 20*log10.(abs.(IOLA.transforms.STMDCT(s, Windows.cosine(len)))),
                      convert(Int, fs), div(len, 2))

(y, fs, nbits, opt) = wavread("audio-128.wav")

diffs128 = IOLA.analyze(y[:,1], (s) -> 20*log10.(abs.(IOLA.transforms.STMDCT(s, Windows.cosine(len)))),
                      convert(Int, fs), div(len, 2))

(y, fs, nbits, opt) = wavread("audio-320.wav")

diffs320 = IOLA.analyze(y[:,1], (s) -> 20*log10.(abs.(IOLA.transforms.STMDCT(s, Windows.cosine(len)))),
                      convert(Int, fs), div(len, 2))
