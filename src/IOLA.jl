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
end

function getParams(type::CodecType)
    return Dict{CodecType, CodecParams}(
        MP3 => CodecParams(Transform.MDCT,
                   DSP.Windows.cosine,
                   1152,
                   576)
    )[type]
end

function getTransformFunction(params::CodecParams)
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

module Utils

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

end

function energyDiff(signal::AbstractVector{SignalType}, transform::Function,
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
                 segment_length::Integer, overlap::Integer) where SignalType <: Number
    siglen = length(signal)
    N = div(siglen-overlap, segment_length)
    out = zeros(N, 2)
    Threads.@threads for i = 1:N
        sind = (i-1)*segment_length+1
        segment_diffs = energyDiff(signal, transform, segment_length, sind, sind+overlap)
        segment_diffs = (segment_diffs.-mean(segment_diffs))./std(segment_diffs)
        out[i, :] = [findmax(segment_diffs)...]
    end
    return out
end

end
