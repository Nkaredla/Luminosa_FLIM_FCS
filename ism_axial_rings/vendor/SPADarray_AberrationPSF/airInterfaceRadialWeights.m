function [radialWeight, diagnostics] = airInterfaceRadialWeights(sim, rho, quadrature)
%AIRINTERFACERADIALWEIGHTS Radial pupil weights for layered air-interface PSF.
%
%   radialWeight = airInterfaceRadialWeights(sim, rho, quadrature)
%
%   The historical layered model used the immersion-side Debye solid-angle
%   measure, (NA/nImm)^2 rho/cos(thetaImm).  That has the correct homogeneous
%   limit but is not the right default for an air-side emitter because the
%   emitting angular spectrum is referenced to the sample-side medium.  The
%   default mode therefore uses the sample-side solid-angle/Weyl projection:
%
%       (NA/nSample)^2 rho/cos(thetaSample)
%
%   with complex cos(thetaSample) retained above the critical angle.  The
%   legacy immersion weighting remains available for controlled comparisons:
%
%       sim.interfaceRadialWeightMode = 'immersionDebye'

    if nargin < 3 || isempty(quadrature)
        quadrature = ones(size(rho));
    end

    rho = double(rho(:)).';
    quadrature = double(quadrature(:)).';
    q = sim.NA * rho;
    nImm = sim.nImmersion;
    nSample = sim.nSample;
    cosImm = positiveRootLocal(1 - (q/nImm).^2);
    cosSample = positiveRootLocal(1 - (q/nSample).^2);

    mode = radialWeightMode(sim);
    legacyMeasure = (sim.NA/nImm)^2 * rho ./ safeDenominator(cosImm);
    sampleMeasure = (sim.NA/nSample)^2 * rho ./ safeDenominator(cosSample);

    switch lower(mode)
        case {'samplesolidangle','sample','sampleemission','weyl'}
            measure = sampleMeasure;
            description = 'sample-side solid-angle / scalar Weyl projection';
        case {'immersiondebye','immersion','legacy'}
            measure = legacyMeasure;
            description = 'legacy immersion-side Debye solid-angle projection';
        otherwise
            error('airInterfaceRadialWeights:BadMode', ...
                ['Unknown interfaceRadialWeightMode "%s". Use ' ...
                '''sampleSolidAngle'' or ''immersionDebye''.'], mode);
    end

    radialWeight = measure .* quadrature;
    legacyWeight = legacyMeasure .* quadrature;
    sampleWeight = sampleMeasure .* quadrature;

    supercritical = q > nSample;
    areaDen = sum(max(rho,0));
    if areaDen > 0
        supercriticalAreaFraction = sum(rho(supercritical)) / areaDen;
    else
        supercriticalAreaFraction = 0;
    end
    contrastDiff = norm(abs(sampleWeight(:)) - abs(legacyWeight(:))) / ...
        max(norm(abs(legacyWeight(:))), eps);

    diagnostics = struct();
    diagnostics.mode = mode;
    diagnostics.description = description;
    diagnostics.q = q;
    diagnostics.cosImmersion = cosImm;
    diagnostics.cosSample = cosSample;
    diagnostics.supercriticalMask = supercritical;
    diagnostics.supercriticalAreaFraction = supercriticalAreaFraction;
    diagnostics.legacyWeight = legacyWeight;
    diagnostics.sampleSolidAngleWeight = sampleWeight;
    diagnostics.selectedWeight = radialWeight;
    diagnostics.sampleVsLegacyRelativeAbsDifference = contrastDiff;
    diagnostics.finiteSelectedWeight = all(isfinite(real(radialWeight(:)))) && ...
        all(isfinite(imag(radialWeight(:))));
end

function mode = radialWeightMode(sim)
    if isfield(sim, 'interfaceRadialWeightMode') && ...
            ~isempty(sim.interfaceRadialWeightMode)
        mode = char(sim.interfaceRadialWeightMode);
    else
        mode = 'sampleSolidAngle';
    end
end

function den = safeDenominator(den)
    den = complex(den);
    tiny = abs(den) < 1e-9;
    if any(tiny)
        phase = angle(den(tiny));
        phase(~isfinite(phase)) = 0;
        den(tiny) = 1e-9 * exp(1i*phase);
    end
end

function root = positiveRootLocal(value)
    root = sqrt(complex(value,0));
    root(imag(root)<0) = conj(root(imag(root)<0));
end
