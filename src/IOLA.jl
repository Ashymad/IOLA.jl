module IOLA
using Statistics

export analyze, Codec

module Codec
import MDCT, DSP, ...IOLA

module Transform
@enum TransformType begin
    MDCT = 1
end
end

struct CodecParams
    transform::Transform.TransformType
    window::Function
    length::Integer
    hop::Integer
end

@enum CodecType begin
    MP3 = 1
    AAC = 2
    AC3 = 3
    OGG = 4
    WMA = 5
end

function getparams(type::CodecType)
    return if type == MP3
        CodecParams(Transform.MDCT,
                    DSP.Windows.cosine,
                    1152,
                    div(1152, 2))
    elseif type == AAC
        CodecParams(Transform.MDCT,
                    (n) -> IOLA.Windows.KBD(n, 4),
                    2048,
                    div(2048, 2))
    elseif type == OGG
        CodecParams(Transform.MDCT,
                    IOLA.Windows.slope,
                    2048,
                    div(2048, 2))
    elseif type == AC3
        CodecParams(Transform.MDCT,
                    (n) -> IOLA.Windows.KBD(n, 5),
                    512,
                    div(512, 2))
    elseif type == WMA
        CodecParams(Transform.MDCT,
                    DSP.Windows.cosine,
                    4096,
                    div(4096, 2))
                    
    else
        error("Codec $type not yet implemented");
    end
end

function gettransform(params::CodecParams)
    if params.transform == Transform.MDCT
        win = params.window(params.length)
        mdct_fn = let mdct_plan = MDCT.plan_mdct(win)
            function(X::AbstractArray{T}) where T<:Number
                mdct_plan*X
            end
        end
        return IOLA.Utils.shorttermify(mdct_fn, win, params.hop)
    else
        error("Transform $(params.transform) not yet implemented!")
    end
end

end

module Windows
import DSP

function KBD(n::Integer, α::Real)
    if isodd(n)
        throw(ArgumentError("KBD window length must be even"))
    end
    N = div(n, 2)
    kais = DSP.Windows.kaiser(N+1, α)
    out = zeros(n)
    out[n] = (out[1] = kais[1])
    @inbounds for i = 2:N
        out[n+1-i] = (out[i] = kais[i] + out[i-1])
    end
    return sqrt.(out./out[N])
end

function slope(n::Integer)
    return sin.(0.5π.*sin.(π/n.*((1:n).-0.5)).^2)
end

end

module Utils

using Statistics
export findperiod, radiusofmean

function fast_sum_abs_log10_abs!(array::AbstractArray{SignalType}) where SignalType <: Number
    N = length(array)
    while N > 1
        for it = 1:2:N
            a1 = abs(array[it])
            if a1 == 0; return Inf; end
            it2 = div(it+1,2)
            if it+1 <= N
                a2 = abs(array[it+1])
                if a2 == 0; return Inf; end
                if (a1 >= 1 && a2 >= 1) || (a1 < 1 && a2 < 1)
                    a1 *= a2
                else
                    a1 /= a2
                end
                if a1 == 0 || a1 == Inf
                    return sum(abs.(log10.(view(array,1:it2-1)))) .+
                        sum(abs.(log10.(view(array,it:N))))
                end
            end
            array[it2] = a1
        end
        N = div(N+1,2)
    end
    return abs(log10(array[1]))
end

function shorttermify(fun::Function, window::AbstractVector{SignalType},
                      hop::Integer) where SignalType <: Number
    return let winlen = length(window), outlen = length(fun(window))
        function(signal::AbstractVector{SignalType})
            siglen = length(signal)
            N = div(siglen, hop)-div(winlen, hop)+1
            pad = div(siglen - (N-1)*hop - winlen, 2)
            signal_padded = zeros(SignalType, 2*pad + siglen + 1)
            signal_padded[pad+1:siglen+pad] = signal
            out = zeros(outlen, N)

            for i = 1:N
                sind = (i-1)*hop+1
                @views out[:,i] = fun(signal_padded[sind:(sind+winlen-1)].*window)
            end
            return out
        end
    end
end

function findperiod(vector::AbstractVector{T},
                    domain::AbstractVector) where T<:Number
    mod_val = typemax(T)
    p_min = 0
    for p in domain
        per = mod.(vector, p)
        val = std(per; corrected=false)/mean(per)
        if val <= mod_val
            mod_val = val
            p_min = p
        end
    end
    return p_min
end

function radiusofmean(radiuses::AbstractVector,
                      angles::AbstractVector)
    M = length(radiuses)
    sqrt(sum(radiuses.*sin.(angles))^2 + sum(radiuses.*cos.(angles))^2)/M
end

end

function energydiff(signal::AbstractVector{SignalType}, transform::Function,
                    segment_length::Integer, start_index::Integer,
                    end_index::Integer) where SignalType <: Number
    N = end_index-start_index+1
    energies = zeros(N)
    out = zeros(N-1)
    for i = start_index:end_index
        arind = i-start_index+1
        segment = transform(view(signal, i:(i+segment_length-1)))
        energies[arind] = Utils.fast_sum_abs_log10_abs!(segment)./segment_length
        if i > start_index
            out[arind-1] = abs(energies[arind-1] - energies[arind])
        end
    end
    return out
end

function analyze(signal::AbstractVector{SignalType}, transform::Function,
                 segment_length::Integer, overlap::Integer,
                 normalize::Bool = false) where SignalType <: Number
    siglen = length(signal)
    return let N = div(siglen-overlap, segment_length), out = zeros(N, 2)
        Threads.@threads for i = 1:N
            sind = (i-1)*segment_length+1
            segment_diffs = energydiff(signal, transform, segment_length, sind, sind+overlap)
            if normalize 
                segment_diffs = (segment_diffs.-mean(segment_diffs))./std(segment_diffs)
            end
            out[i, :] = [findmax(segment_diffs)...]
        end
        out
    end
end

end
