include("iola.jl")
using DSP
using WAV
using Gnuplot

len = 1152
hop = div(len,2)
period = sqrt(hop)/2

(y, fs, nbits, opt) = wavread("audio.wav")

diffs = IOLA.analyze(y[:,1], (s) -> 20*log10.(abs.(IOLA.transforms.STMDCT(s, Windows.cosine(len)))),
                      convert(Int, fs), hop)

(y, fs, nbits, opt) = wavread("audio-128.wav")

diffs128 = IOLA.analyze(y[:,1], (s) -> 20*log10.(abs.(IOLA.transforms.STMDCT(s, Windows.cosine(len), hop))),
                      convert(Int, fs), hop)

(y, fs, nbits, opt) = wavread("audio-320.wav")

diffs320 = IOLA.analyze(y[:,1], (s) -> 20*log10.(abs.(IOLA.transforms.STMDCT(s, Windows.cosine(len), hop))),
                      convert(Int, fs), hop)
@gp("set angles radians",
    "set polar",
    "set grid polar 15. lt -1 dt 0 lw 0.5",
    "unset border",
    "unset param",
    "unset xtics",
    2π*mod.(diffs128[:,2], period)/period, diffs128[:,1], "pt 7 t '128'",
    2π*mod.(diffs320[:,2], period)/period, diffs320[:,1], "pt 7 t '320'",
    2π*mod.(diffs[:,2], period)/period, diffs[:,1], "pt 7 t 'none'"
   )

