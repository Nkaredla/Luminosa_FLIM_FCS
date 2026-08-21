function [countsPerFrame, meanBrightnessPhotons] = ...
        bnpBrightEmitterCounts(state, opts)
%BNPBRIGHTEMITTERCOUNTS Count active emitters above the focus threshold.
%
%   Brightness is the expected detector-array signal after focused
%   illumination and finite-array capture. Intrinsic molecular brightness
%   remains shared through state.brightnessScale.

    amplitudes = reshape(sum(state.contributions, 3), ...
        state.nCandidates, state.nFrames);
    bright = amplitudes >= opts.bnpInFocusPhotonThreshold & state.active;
    countsPerFrame = sum(bright, 1).';
    selectedBrightness = amplitudes(bright);
    if isempty(selectedBrightness)
        meanBrightnessPhotons = NaN;
    else
        meanBrightnessPhotons = mean(selectedBrightness);
    end
end
