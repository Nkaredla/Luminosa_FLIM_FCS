function [hEff, hExc, eta] = spadEffectivePSFArrayInterface(sim, coeffs, ...
        stageZ, emitterHeights)
%SPADEFFECTIVEPSFARRAYINTERFACE Per-channel effective PSF across an interface.
%
% [hEff, hExc, eta] = spadEffectivePSFArrayInterface(sim, coeffs, stageZ, ...
%                                                    emitterHeights)
%
% The interface counterpart of spadEffectivePSFArray, which is hard-wired to
% the homogeneous psfBessel and therefore cannot represent an oil objective
% looking into water.
%
% WHY A SEPARATE FUNCTION IS NEEDED
%
% In the homogeneous model the pupil is parameterised as
%   sinTheta = (NA / nMedium) * rho
% which requires NA <= nMedium and so cannot express NA 1.45 in n = 1.33. The
% stratified model instead uses
%   q = NA * rho,  cosSample = sqrt(1 - (q/nSample)^2)
% so the pupil spans the FULL aperture and components with q > nSample get a
% complex cosSample. Those are the supercritical-angle components: a dipole
% near the interface has near-field content with in-plane momentum above
% n_sample*k0, and at the water-glass boundary it couples into propagating
% waves in the higher-index glass. This is the reciprocal of TIRF and is
% collected, not lost.
%
% That channel matters here. Supercritical-angle fluorescence decays over
% roughly 100-200 nm from the surface, so its share of the signal is a steep
% axial reporter exactly where the membrane sits, and is nearly absent for
% internalised dye several hundred nm up. Capping the aperture at n_sample
% would discard the most distance-sensitive part of the measurement, not a
% marginal one.
%
% The excitation is passed through the same stratified model rather than a
% homogeneous one, because with a full-aperture NA 1.45 focus the
% supercritical components also produce an evanescent excitation contribution
% near the surface. That is a second distance-dependent term and it belongs in
% the product.
%
% INPUTS
%   sim             defaultParams-style struct. Requires nSample, nGlass,
%                   nImmersion, coverslipThicknessUm and the detector fields.
%                   sim.NA is the FULL objective NA and is not capped.
%   coeffs          aberration coefficients, e.g. coeffStruct(sim, zeros(...))
%   stageZ          objective/stage focus position(s) relative to the
%                   interface, um. Scalar is broadcast.
%   emitterHeights  emitter height(s) above the interface, um.
%
% OUTPUTS
%   hEff(y,x,c,k)   effective excitation-times-detection PSF for channel k
%   hExc(y,x,c)     excitation PSF
%   eta(y,x,c,k)    per-channel collection efficiency
%
% The third dimension indexes the (stageZ, emitterHeight) condition pairs in
% the order psfBesselAirInterface returns them.

    if nargin < 3 || isempty(stageZ); stageZ = 0; end
    if nargin < 4 || isempty(emitterHeights); emitterHeights = 0; end
    for required = {'psfBesselAirInterface', ...
            'detectorCollectionEfficiencyExplicit'}
        if exist(required{1}, 'file') ~= 2
            error('spadEffectivePSFArrayInterface:MissingDependency', ...
                '%s.m is not on the MATLAB path.', required{1});
        end
    end
    if ~isfield(sim, 'nSample') || isempty(sim.nSample)
        error('spadEffectivePSFArrayInterface:NoSampleIndex', ...
            'sim.nSample must be set (1.33 for a water-like medium).');
    end

    coeffs = coeffStruct(sim, coeffs);

    hExc = psfBesselAirInterface(sim, coeffs, sim.lamExc, stageZ, ...
        emitterHeights);
    hDet = psfBesselAirInterface(sim, coeffs, sim.lamEm, stageZ, ...
        emitterHeights);

    % detectorCollectionEfficiencyExplicit integrates the detection PSF over
    % each finite detector element, so it is geometry-agnostic and is reused
    % unchanged.
    eta = detectorCollectionEfficiencyExplicit(sim, hDet);

    if ndims(eta) ~= 4
        error('spadEffectivePSFArrayInterface:EtaShape', ...
            'Expected eta to be [ny nx nCond nCh], got %s.', ...
            mat2str(size(eta)));
    end
    if size(hExc, 3) ~= size(eta, 3)
        error('spadEffectivePSFArrayInterface:CondMismatch', ...
            ['Excitation has %d conditions but the efficiency stack has ' ...
             '%d.'], size(hExc, 3), size(eta, 3));
    end

    % Confocal/ISM detection: the signal on a channel is the excitation at the
    % emitter times that channel's collection efficiency.
    hEff = double(hExc) .* double(eta);
end
