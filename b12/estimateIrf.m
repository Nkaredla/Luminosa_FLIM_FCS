function [irf, source, irfParameters, irfTau, irfFit, irfChi] = ...
        estimateIrf(globalCurve, dtNs, suppliedIrf, minimumPhotons, ...
        useGpu, cascadeCount, substartCount, tauSeed)
%ESTIMATEIRF Return a normalised IRF on the selected B12 gate.
% A supplied calibration IRF is preferred. Otherwise the existing global
% shifted-gamma estimator is applied to the pooled B12 decay.

    if ~isempty(suppliedIrf)
        irf = max(double(suppliedIrf(:)), 0);
        if numel(irf) ~= numel(globalCurve)
            error('The supplied IRF has %d bins; the B12 gate has %d bins.', ...
                numel(irf), numel(globalCurve));
        end
        irf = irf ./ max(sum(irf), eps);
        source = 'supplied calibration IRF';
        irfParameters = nan(1, 3);
        irfTau = nan(1, numel(tauSeed));
        irfFit = nan(size(globalCurve));
        irfChi = NaN;
        return;
    end

    if sum(globalCurve) < minimumPhotons
        error('Only %d B12 photons are available; %d are required for IRF estimation.', ...
            round(sum(globalCurve)), minimumPhotons);
    end

    % These structures belong to the external legacy API and are unpacked
    % immediately; they are not passed between B12 helper functions.
    fitHeader = struct('MeasDesc_Resolution', dtNs * 1e-9);
    options = struct('useGPU', useGpu, 'nCasc', cascadeCount, 'nSub', substartCount);
    legacyOutput = Calc_mIRF_Global_GammaShifted_fast( ...
        fitHeader, globalCurve, tauSeed, options);

    irf = max(double(legacyOutput.IRF(:)), 0);
    irf = irf ./ max(sum(irf), eps);
    source = 'global B12 shifted-gamma estimate';
    irfParameters = double(legacyOutput.pIRF(:)).';
    irfTau = double(legacyOutput.tauFit(:)).';
    irfFit = double(legacyOutput.fit(:));
    irfChi = double(legacyOutput.chi);
end
