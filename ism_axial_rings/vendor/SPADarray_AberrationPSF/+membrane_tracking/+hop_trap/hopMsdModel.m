function msdUm2 = hopMsdModel(parameters, lagTimeS)
%HOPMSDMODEL Confined short-range motion plus long-range diffusion.
%
%   parameters = [Dmacro, plateauAmplitude, confinementTime].

    macroDiffusion = parameters(1);
    amplitude = parameters(2);
    confinementTime = parameters(3);
    msdUm2 = 4 * macroDiffusion * lagTimeS + ...
        amplitude * (1 - exp(-lagTimeS / confinementTime));
end
