module IOLA
using Statistics

export energyDiff, analyze

struct Codec
    transform::Function
    window::Function
    length::Integer
    hop::Integer
end

module transforms
import MDCT

export STMDCT

function STMDCT(signal::AbstractVector{SignalType},
                window::AbstractVector{SignalType},
                hop::Integer=div(length(window),2)) where SignalType <: Number
    winlen = length(window)
    siglen = length(signal)
    N = div(siglen, hop)-div(winlen, hop)+1
    out = zeros(div(winlen,2), N)
    pad = div(siglen - (N-1)*hop - winlen, 2)

    signal_padded = [zeros(SignalType, pad); signal; zeros(SignalType, pad+1)]
    
    for i = 1:N
        sind = (i-1)*hop+1
        out[:,i] = MDCT.mdct(signal_padded[sind:(sind+winlen-1)].*window)
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
        energies[arind] = sqrt(mean(segment.^2))
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
    for i = 1:N
        sind = (i-1)*segment_length+1
        segment_diffs = energyDiff(signal, transform, segment_length, sind, sind+overlap)
        segment_diffs = (segment_diffs.-mean(segment_diffs))./std(segment_diffs)
        out[i, :] = [findmax(segment_diffs)...]
    end
    return out
end
    
end
