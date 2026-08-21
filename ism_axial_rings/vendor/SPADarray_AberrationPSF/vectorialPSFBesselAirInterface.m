function intensity = vectorialPSFBesselAirInterface( ...
        sim, coeffs, wavelength, stageZ, emitterHeight)
%VECTORIALPSFBESSELAIRINTERFACE Vectorial air-interface PSF.
%
%   This is the vectorial counterpart to scalarPSFBesselAirInterface.
%
%   sim.airInterfaceDirection controls the physical direction:
%       'excitation' : objective/immersion -> glass -> air focus.
%                      Uses immersion-side Debye measure.
%       'collection' : air emitter -> glass -> immersion/objective.
%                      Uses sample-side solid-angle measure.
%
%   The glass is treated as a single-pass thick layer: phase is retained,
%   Fabry-Perot multiple reflections are not. The recorded axial stage/focus
%   displacement is referenced to the immersion by default (consistent with
%   scalarPSFBesselAirInterface, since the piezo moves the coverslip through
%   the immersion gap); set sim.airInterfaceStageMedium = 'glass' for the
%   legacy convention. Intensity defaults to the x/y-averaged (unpolarized /
%   isotropic) vector field; set sim.vectorialPolarizationMode = 'x' or 'y'
%   for a linearly polarized input. This is an effective isotropic-bead model,
%   not a full oriented-dipole Green tensor.

    validateScalarOptics(sim);
    coeffs = coeffStruct(sim, coeffs);
    if ~isscalar(wavelength) || ~isfinite(wavelength) || wavelength <= 0
        error('vectorialPSFBesselAirInterface:InvalidWavelength', ...
            'wavelength must be a finite positive scalar in micrometers.');
    end

    stageZ = double(stageZ(:)).';
    emitterHeight = double(emitterHeight(:)).';
    if isscalar(stageZ) && numel(emitterHeight) > 1
        stageZ = repmat(stageZ, size(emitterHeight));
    elseif isscalar(emitterHeight) && numel(stageZ) > 1
        emitterHeight = repmat(emitterHeight, size(stageZ));
    end
    if numel(stageZ) ~= numel(emitterHeight) || ...
            any(~isfinite(stageZ)) || any(~isfinite(emitterHeight))
        error('vectorialPSFBesselAirInterface:BadConditions', ...
            'stageZ and emitterHeight must be finite and scalar-expandable.');
    end
    if any(emitterHeight < 0)
        error('vectorialPSFBesselAirInterface:EmitterBelowInterface', ...
            'Air-side emitter heights must be nonnegative.');
    end

    psfKey = buildVectorialAirPsfKey(sim, coeffs, wavelength, stageZ, emitterHeight);
    cachedIntensity = vectorialAirPsfCache(psfKey);
    if ~isempty(cachedIntensity)
        intensity = cachedIntensity;
        return;
    end

    nImm = sim.nImmersion;
    nGlass = sim.nGlass;
    nSample = sim.nSample;
    nDesignGlass = sim.nDesignGlass;
    k0 = 2*pi/wavelength;

    rho = linspace(0, 1, sim.Nr);
    phi = linspace(0, 2*pi, sim.Nphi+1); phi(end) = [];
    [R, P] = ndgrid(rho, phi);

    W = wavefrontWaves(R, P, coeffs);
    aberrationPhase = exp(1i*2*pi*(sim.lamRef/wavelength).*W);

    q = sim.NA * rho;
    cosImm = positiveRoot(1 - (q/nImm).^2);
    cosGlass = positiveRoot(1 - (q/nGlass).^2);
    cosSample = positiveRoot(1 - (q/nSample).^2);
    cosDesignGlass = positiveRoot(1 - (q/nDesignGlass).^2);
    sinSample = q / nSample;

    direction = airInterfaceDirection(sim);
    switch direction
        case 'excitation'
            [tpIG, tsIG] = transmissionCoefficients(nImm, cosImm, nGlass, cosGlass);
            [tpGS, tsGS] = transmissionCoefficients(nGlass, cosGlass, nSample, cosSample);
            tp = tpIG .* tpGS;
            ts = tsIG .* tsGS;
        case 'collection'
            [tpSG, tsSG] = transmissionCoefficients(nSample, cosSample, nGlass, cosGlass);
            [tpGI, tsGI] = transmissionCoefficients(nGlass, cosGlass, nImm, cosImm);
            tp = tpSG .* tpGI;
            ts = tsSG .* tsGI;
        otherwise
            error('vectorialPSFBesselAirInterface:BadDirection', ...
                'airInterfaceDirection must be excitation or collection.');
    end

    coverPhase = k0 * ( ...
        nGlass*cosGlass*sim.coverslipThicknessUm - ...
        nDesignGlass*cosDesignGlass*sim.designCoverslipThicknessUm);
    common = sqrt(cosImm(:)) .* exp(1i*coverPhase(:)) .* aberrationPhase;

    c = cos(P); s = sin(P);
    pupil = vectorialLayeredPupilComponents( ...
        common, c, s, tp(:), ts(:), cosSample(:), sinSample(:));

    ms = -sim.M:sim.M;
    pm = componentSpectra(pupil, ms, sim.Nphi);

    [X, Y] = meshgrid(sim.x, sim.y);
    radius = hypot(X, Y);
    psi = atan2(Y, X);
    rv = radius(:);
    roundedRadius = round(rv, 12);
    uniqueRadius = unique(roundedRadius);
    [~, radiusIndex] = ismember(roundedRadius, uniqueRadius);

    argument = k0 * sim.NA * (uniqueRadius * rho);
    besselTables = besselTablesForVectorialAir( ...
        wavelength, sim.NA, sim.Nr, uniqueRadius, sim.M, argument);

    quadrature = trapezoidWeights(rho);
    radialWeight = vectorialAirRadialWeights(sim, rho, quadrature, cosImm, direction);
    [nStage, cosStage] = stageMedium(sim, nImm, cosImm, nGlass, cosGlass);
    axial = exp(1i*k0*( ...
        (nStage*cosStage(:))*stageZ + ...
        (nSample*cosSample(:))*emitterHeight));

    ix = synthesizeIntensity(pm.x, ms, radialWeight, axial, ...
        besselTables, radiusIndex, psi);
    iy = synthesizeIntensity(pm.y, ms, radialWeight, axial, ...
        besselTables, radiusIndex, psi);
    combined = combinePolarizations(ix, iy, polarizationModeString(sim));
    intensity = reshape(combined, numel(sim.y), numel(sim.x), []);

    vectorialAirPsfCache(psfKey, intensity);
end

function W = wavefrontWaves(R, P, coeffs)
    Z = zernikeModes(R, P);
    W = zeros(size(R));
    fn = fieldnames(coeffs);
    for k = 1:numel(fn)
        if ~isfield(Z, fn{k})
            error('vectorialPSFBesselAirInterface:UnknownMode', ...
                'Unknown Zernike mode "%s".', fn{k});
        end
        W = W + coeffs.(fn{k}) * Z.(fn{k});
    end
end

function pupil = vectorialLayeredPupilComponents(common, c, s, tp, ts, cs, ss)
    % x-polarized entrance field: cos(phi) p - sin(phi) s
    pupil.x.ex = common .* (tp .* cs .* c.^2 + ts .* s.^2);
    pupil.x.ey = common .* ((tp .* cs - ts) .* s .* c);
    pupil.x.ez = common .* (-tp .* ss .* c);

    % y-polarized entrance field: sin(phi) p + cos(phi) s
    pupil.y.ex = common .* ((tp .* cs - ts) .* s .* c);
    pupil.y.ey = common .* (tp .* cs .* s.^2 + ts .* c.^2);
    pupil.y.ez = common .* (-tp .* ss .* s);
end

function pm = componentSpectra(pupil, ms, nPhi)
    pm.x = oneInputSpectra(pupil.x, ms, nPhi);
    pm.y = oneInputSpectra(pupil.y, ms, nPhi);
end

function out = oneInputSpectra(pupil, ms, nPhi)
    out.ex = azimuthalCoefficients(pupil.ex, ms, nPhi);
    out.ey = azimuthalCoefficients(pupil.ey, ms, nPhi);
    out.ez = azimuthalCoefficients(pupil.ez, ms, nPhi);
end

function pm = azimuthalCoefficients(pupil, ms, nPhi)
    spectrum = fft(pupil, [], 2) / nPhi;
    pm = complex(zeros(size(pupil, 1), numel(ms)));
    for k = 1:numel(ms)
        pm(:,k) = spectrum(:, mod(ms(k), nPhi) + 1);
    end
end

function intensity = synthesizeIntensity(pm, ms, radialWeight, axial, ...
        besselTables, radiusIndex, psi)
    ex = synthesizeFields(pm.ex, ms, radialWeight, axial, ...
        besselTables, radiusIndex, psi);
    ey = synthesizeFields(pm.ey, ms, radialWeight, axial, ...
        besselTables, radiusIndex, psi);
    ez = synthesizeFields(pm.ez, ms, radialWeight, axial, ...
        besselTables, radiusIndex, psi);
    intensity = abs(ex).^2 + abs(ey).^2 + abs(ez).^2;
end

function fields = synthesizeFields(pm, ms, radialWeight, axial, ...
        besselTables, radiusIndex, psi)
    nCondition = size(axial, 2);
    nPixel = numel(psi);
    fields = complex(zeros(nPixel, nCondition));
    psiVector = psi(:);
    for k = 1:numel(ms)
        m = ms(k);
        J = besselTables{abs(m)+1};
        if m < 0
            J = ((-1)^(-m)) * J;
        end
        basis = (pm(:,k).*radialWeight(:)) .* axial;
        radial = J * basis;
        angular = (2*pi) * (1i^m) * exp(1i*m*psiVector);
        fields = fields + angular .* radial(radiusIndex,:);
    end
end

function w = trapezoidWeights(rho)
    w = ones(size(rho));
    if numel(rho) > 1
        w([1 end]) = 0.5;
        w = w * (rho(2) - rho(1));
    end
end

function root = positiveRoot(value)
    root = sqrt(complex(value, 0));
    root(imag(root) < 0) = conj(root(imag(root) < 0));
end

function [tp, ts] = transmissionCoefficients(n1, c1, n2, c2)
    tp = 2*n1*c1 ./ (n2*c1 + n1*c2);
    ts = 2*n1*c1 ./ (n1*c1 + n2*c2);
end

function J = besselTablesForVectorialAir(wavelength, NA, Nr, uniqueRadius, M, argument)
    persistent keys vals
    if isempty(keys), keys = {}; vals = {}; end
    for i = 1:numel(keys)
        ki = keys{i};
        if ki.wavelength == wavelength && ki.NA == NA && ki.Nr == Nr && ...
                ki.M == M && isequal(ki.uniqueRadius, uniqueRadius)
            J = vals{i};
            return;
        end
    end
    J = cell(M+1, 1);
    for a = 0:M
        J{a+1} = besselj(a, argument);
    end
    keys{end+1} = struct('wavelength', wavelength, 'NA', NA, 'Nr', Nr, ...
        'M', M, 'uniqueRadius', uniqueRadius); %#ok<AGROW>
    vals{end+1} = J; %#ok<AGROW>
    if numel(keys) > 6
        keys(1) = []; vals(1) = [];
    end
end

function key = buildVectorialAirPsfKey(sim, coeffs, wavelength, stageZ, emitterHeight)
    key = struct( ...
        'model', 'vectorialAirInterface', ...
        'wavelength', wavelength, ...
        'NA', sim.NA, ...
        'nImmersion', sim.nImmersion, ...
        'nGlass', sim.nGlass, ...
        'nSample', sim.nSample, ...
        'nDesignGlass', sim.nDesignGlass, ...
        'coverslipThicknessUm', sim.coverslipThicknessUm, ...
        'designCoverslipThicknessUm', sim.designCoverslipThicknessUm, ...
        'lamRef', sim.lamRef, ...
        'Nr', sim.Nr, 'Nphi', sim.Nphi, 'M', sim.M, ...
        'x', sim.x(:).', 'y', sim.y(:).', ...
        'radialMode', radialWeightModeString(sim), ...
        'direction', airInterfaceDirection(sim), ...
        'stageMedium', stageMediumString(sim), ...
        'polarizationMode', polarizationModeString(sim), ...
        'coeffs', coeffs, ...
        'stageZ', stageZ, 'emitterHeight', emitterHeight);
end

function [nStage, cosStage] = stageMedium(sim, nImm, cosImm, nGlass, cosGlass)
    switch stageMediumString(sim)
        case 'glass'
            nStage = nGlass;
            cosStage = cosGlass;
        case 'immersion'
            nStage = nImm;
            cosStage = cosImm;
        otherwise
            error('vectorialPSFBesselAirInterface:BadStageMedium', ...
                'airInterfaceStageMedium must be glass or immersion.');
    end
end

function s = stageMediumString(sim)
    % Default 'immersion': the piezo moves the coverslip through the immersion
    % gap, so the stage defocus accrues in the immersion. Matches the scalar
    % model. Set 'glass' only for the legacy convention.
    if isfield(sim, 'airInterfaceStageMedium') && ...
            ~isempty(sim.airInterfaceStageMedium)
        s = lower(char(sim.airInterfaceStageMedium));
    else
        s = 'immersion';
    end
    switch s
        case {'glass','coverslip'}
            s = 'glass';
        case {'immersion','water','objective'}
            s = 'immersion';
    end
end

function direction = airInterfaceDirection(sim)
    if isfield(sim, 'airInterfaceDirection') && ~isempty(sim.airInterfaceDirection)
        direction = lower(char(sim.airInterfaceDirection));
    else
        direction = 'collection';
    end
    switch direction
        case {'excitation','illumination','objective_to_sample'}
            direction = 'excitation';
        case {'collection','emission','sample_to_objective'}
            direction = 'collection';
    end
end

function radialWeight = vectorialAirRadialWeights(sim, rho, quadrature, cosImm, direction)
    switch direction
        case 'excitation'
            nImm = sim.nImmersion;
            denom = cosImm;
            tiny = abs(denom) < 1e-12;
            denom(tiny) = 1e-12;
            radialWeight = (sim.NA/nImm)^2 * rho ./ denom .* quadrature;
            radialWeight = radialWeight(:);
        case 'collection'
            radialWeight = airInterfaceRadialWeights(sim, rho, quadrature);
        otherwise
            error('vectorialPSFBesselAirInterface:BadDirection', ...
                'Unsupported air-interface direction "%s".', direction);
    end
end

function s = radialWeightModeString(sim)
    if isfield(sim, 'interfaceRadialWeightMode') && ...
            ~isempty(sim.interfaceRadialWeightMode)
        s = char(sim.interfaceRadialWeightMode);
    else
        s = 'sampleSolidAngle';
    end
end

function s = polarizationModeString(sim)
    if isfield(sim, 'vectorialPolarizationMode') && ...
            ~isempty(sim.vectorialPolarizationMode)
        s = char(sim.vectorialPolarizationMode);
    else
        s = 'xyAverage';
    end
end

function I = combinePolarizations(ix, iy, mode)
    % How the two entrance polarizations combine into the recorded intensity.
    switch lower(mode)
        case {'xyaverage','unpolarized','isotropic','circular'}
            I = 0.5*(ix + iy);     % unpolarized / isotropic-bead detection
        case 'x'
            I = ix;                % linearly x-polarized input
        case 'y'
            I = iy;                % linearly y-polarized input
        otherwise
            error('vectorialPSFBesselAirInterface:BadPolarizationMode', ...
                'vectorialPolarizationMode must be xyAverage, x, or y.');
    end
end

function val = vectorialAirPsfCache(key, putVal)
    persistent K V
    capacity = 24;
    if isempty(K), K = {}; V = {}; end
    if nargin == 2
        K{end+1} = key; V{end+1} = putVal; %#ok<AGROW>
        if numel(K) > capacity, K(1) = []; V(1) = []; end
        val = putVal;
        return;
    end
    val = [];
    for i = 1:numel(K)
        if isequaln(K{i}, key)
            val = V{i};
            return;
        end
    end
end
