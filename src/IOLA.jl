module IOLA
using Statistics

export energyDiff, analyze, transforms

struct Codec
    transform::Function
    window::Function
    length::Integer
    hop::Integer
end

module utils

function fast_sum_abs_log10_abs!(vector::AbstractVector{SignalType}) where SignalType <: Number
    N = length(vector)
    while N > 1
        for it = 1:2:N
            if it+1 <= N
                a1 = vector[it]
                a2 = vector[it+1]
                if ((a1 >= 1 || a1 <= -1) && (a2 >= 1 || a2 <= -1)) ||
                    (a1 < 1  && a1 > -1 && a2 > -1 && a2 < 1)
                    vector[div(it+1,2)] = a1*a2
                else
                    vector[div(it+1,2)] = a1/a2
                end
            else
                vector[div(it+1,2)] = vector[it]
            end
        end
        N = div(N+1,2)
    end
    return abs(log10(abs(vector[1])))
end

end

module transforms
import MDCT

function STMDCT(signal::AbstractVector{SignalType},
                window::AbstractVector{SignalType},
                hop::Integer=div(length(window),2),
                mdct_fn::Function=MDCT.mdct) where SignalType <: Number
    winlen = length(window)
    siglen = length(signal)
    N = div(siglen, hop)-div(winlen, hop)+1
    out = zeros(div(winlen,2), N)
    pad = div(siglen - (N-1)*hop - winlen, 2)

    signal_padded = [zeros(SignalType, pad); signal; zeros(SignalType, pad+1)]
    
    for i = 1:N
        sind = (i-1)*hop+1
        out[:,i] = mdct_fn(signal_padded[sind:(sind+winlen-1)].*window)
    end
    return out
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
        segment = transform(signal[i:(i+segment_length-1)])
        energies[arind] = mean(abs.(log10.(abs.(segment))))
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
