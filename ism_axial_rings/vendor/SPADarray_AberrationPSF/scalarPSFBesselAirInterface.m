function intensity = scalarPSFBesselAirInterface( ...
        sim, coeffs, wavelength, stageZ, emitterHeight)
%SCALARPSFBESSELAIRINTERFACE Layered scalar field for an air-side emitter.
%
%   intensity(y,x,condition) is evaluated for a water-immersion objective,
%   glass substrate, and sample-side air layer. stageZ moves the complete
%   sample/interface relative to the objective. emitterHeight is measured
%   from the glass/air interface into air.
%
%   Separate s and p scalar fields are propagated and their intensities are
%   averaged. Complex sample-side kz retains evanescent decay above the
%   critical angle. This is a Gibson-Lanni-style scalar interface model,
%   not a vectorial dipole model.  The radial weighting defaults to the
%   sample-side angular spectrum / solid-angle projection.  Set
%   sim.interfaceRadialWeightMode = 'immersionDebye' to reproduce the legacy
%   immersion-side Debye weighting.

    validateScalarOptics(sim);
    coeffs = coeffStruct(sim,coeffs);
    if ~isscalar(wavelength) || ~isfinite(wavelength) || wavelength <= 0
        error('scalarPSFBesselAirInterface:InvalidWavelength', ...
            'wavelength must be a finite positive scalar in micrometers.');
    end

    stageZ = double(stageZ(:)).';
    emitterHeight = double(emitterHeight(:)).';
    if isscalar(stageZ) && numel(emitterHeight) > 1
        stageZ = repmat(stageZ,size(emitterHeight));
    elseif isscalar(emitterHeight) && numel(stageZ) > 1
        emitterHeight = repmat(emitterHeight,size(stageZ));
    end
    if numel(stageZ) ~= numel(emitterHeight) || ...
            any(~isfinite(stageZ)) || any(~isfinite(emitterHeight))
        error('scalarPSFBesselAirInterface:BadConditions', ...
            'stageZ and emitterHeight must be finite and scalar-expandable.');
    end
    if any(emitterHeight < 0)
        error('scalarPSFBesselAirInterface:EmitterBelowInterface', ...
            'Air-side emitter heights must be nonnegative.');
    end

    % Full-PSF memoization. The output depends only on the optics indices,
    % NA, pupil discretization, coverslip terms, lateral grid, radial-weight
    % mode, aberration coeffs, wavelength, and the (stageZ, emitterHeight)
    % conditions -- NOT on the detector layout. So across the detector-pitch
    % scan (same coeffs/optics/grid, only detXY changes) identical PSF calls
    % return instantly. KEY MAINTENANCE: if you add any new sim dependency to
    % the computation below, add it to buildAirPsfKey or the cache goes stale.
    psfKey = buildAirPsfKey(sim, coeffs, wavelength, stageZ, emitterHeight);
    cachedIntensity = airPsfCache(psfKey);
    if ~isempty(cachedIntensity)
        intensity = cachedIntensity;
        return;
    end

    nImm = sim.nImmersion;
    nGlass = sim.nGlass;
    nSample = sim.nSample;
    nDesignGlass = sim.nDesignGlass;
    k0 = 2*pi/wavelength;

    rho = linspace(0,1,sim.Nr);
    phi = linspace(0,2*pi,sim.Nphi+1);
    phi(end) = [];
    [R,P] = ndgrid(rho,phi);

    W = zeros(size(R));
    Z = zernikeModes(R,P);
    names = fieldnames(coeffs);
    for k = 1:numel(names)
        W = W+coeffs.(names{k})*Z.(names{k});
    end
    aberrationPhase = exp(1i*2*pi*(sim.lamRef/wavelength).*W);

    q = sim.NA*rho;
    cosImm = positiveRoot(1-(q/nImm).^2);
    cosGlass = positiveRoot(1-(q/nGlass).^2);
    cosSample = positiveRoot(1-(q/nSample).^2);
    cosDesignGlass = positiveRoot(1-(q/nDesignGlass).^2);

    [tpIG,tsIG] = transmissionCoefficients( ...
        nImm,cosImm,nGlass,cosGlass);
    [tpGS,tsGS] = transmissionCoefficients( ...
        nGlass,cosGlass,nSample,cosSample);
    tp = tpIG.*tpGS;
    ts = tsIG.*tsGS;

    coverPhase = k0*( ...
        nGlass*cosGlass*sim.coverslipThicknessUm- ...
        nDesignGlass*cosDesignGlass*sim.designCoverslipThicknessUm);
    common = sqrt(cosImm(:)).*exp(1i*coverPhase(:));
    pupilP = common.*tp(:).*aberrationPhase;
    pupilS = common.*ts(:).*aberrationPhase;

    ms = -sim.M:sim.M;
    pmP = azimuthalCoefficients(pupilP,ms,sim.Nphi);
    pmS = azimuthalCoefficients(pupilS,ms,sim.Nphi);

    [X,Y] = meshgrid(sim.x,sim.y);
    radius = sqrt(X.^2+Y.^2);
    psi = atan2(Y,X);
    rv = radius(:);
    roundedRadius = round(rv,12);
    uniqueRadius = unique(roundedRadius);
    [~,radiusIndex] = ismember(roundedRadius,uniqueRadius);

    argument = k0*sim.NA*(uniqueRadius*rho);
    % Bessel tables depend only on {wavelength, NA, Nr, uniqueRadius} -- never
    % on the aberration coeffs -- so they are identical on every forward
    % evaluation of the whole calibration. Cache them (the dominant per-eval
    % cost) keyed exactly on those invariants.
    besselTables = besselTablesForAirInterface( ...
        wavelength, sim.NA, sim.Nr, uniqueRadius, sim.M, argument);
    quadrature = ones(size(rho));
    if numel(rho) > 1
        quadrature([1 end]) = 0.5;
        quadrature = quadrature*(rho(2)-rho(1));
    end
    radialWeight = airInterfaceRadialWeights(sim, rho, quadrature);
    % Stage/focus defocus is referenced to the immersion by default (the piezo
    % moves the coverslip through the immersion gap). Configurable, and shared
    % with the vectorial model so the two stay consistent.
    [nStage, cosStage] = stageMediumScalar(sim, nImm, cosImm, nGlass, cosGlass);
    axial = exp(1i*k0*( ...
        (nStage*cosStage(:))*stageZ+ ...
        (nSample*cosSample(:))*emitterHeight));

    fieldP = synthesizeFields(pmP,ms,radialWeight,axial,besselTables, ...
        radiusIndex,psi);
    fieldS = synthesizeFields(pmS,ms,radialWeight,axial,besselTables, ...
        radiusIndex,psi);
    intensity = reshape(0.5*(abs(fieldP).^2+abs(fieldS).^2), ...
        numel(sim.y),numel(sim.x),[]);
    airPsfCache(psfKey, intensity);
end

function root = positiveRoot(value)
    root = sqrt(complex(value,0));
    root(imag(root)<0) = conj(root(imag(root)<0));
end

function [tp,ts] = transmissionCoefficients(n1,c1,n2,c2)
    tp = 2*n1*c1./(n2*c1+n1*c2);
    ts = 2*n1*c1./(n1*c1+n2*c2);
end

function pm = azimuthalCoefficients(pupil,ms,nPhi)
    spectrum = fft(pupil,[],2)/nPhi;
    pm = complex(zeros(size(pupil,1),numel(ms)));
    for k = 1:numel(ms)
        pm(:,k) = spectrum(:,mod(ms(k),nPhi)+1);
    end
end

function fields = synthesizeFields(pm,ms,radialWeight,axial,besselTables, ...
        radiusIndex,psi)
    nCondition = size(axial,2);
    nPixel = numel(psi);
    fields = complex(zeros(nPixel,nCondition));
    psiVector = psi(:);
    for k = 1:numel(ms)
        m = ms(k);
        J = besselTables{abs(m)+1};         % cached besselj(|m|, argument)
        if m < 0
            J = ((-1)^(-m))*J;
        end
        basis = (pm(:,k).*radialWeight(:)).*axial;
        radial = J*basis;
        angular = (2*pi)*(1i^m)*exp(1i*m*psiVector);
        fields = fields+angular.*radial(radiusIndex,:);
    end
end

function J = besselTablesForAirInterface(wavelength,NA,Nr,uniqueRadius,M,argument)
%   besselj(|m|, argument) tables for |m| = 0..M, memoized on the exact
%   invariants of the argument grid. besselj(|m|) covers both +m and -m.
    persistent keys vals
    if isempty(keys), keys = {}; vals = {}; end
    for i = 1:numel(keys)
        ki = keys{i};
        if ki.wavelength == wavelength && ki.NA == NA && ki.Nr == Nr && ...
                isequal(ki.uniqueRadius, uniqueRadius)
            J = vals{i};
            return;
        end
    end
    J = cell(M+1,1);
    for a = 0:M
        J{a+1} = besselj(a, argument);
    end
    keys{end+1} = struct('wavelength',wavelength,'NA',NA,'Nr',Nr, ...
        'uniqueRadius',uniqueRadius); %#ok<AGROW>
    vals{end+1} = J; %#ok<AGROW>
    if numel(keys) > 6        % keep a few wavelengths/grids; drop oldest
        keys(1) = []; vals(1) = [];
    end
end

function key = buildAirPsfKey(sim,coeffs,wavelength,stageZ,emitterHeight)
%   Every input the PSF output depends on. Keep this list complete: a
%   missing field would let a stale PSF be returned after that field changes.
    key = struct( ...
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
        'stageMedium', stageMediumStringScalar(sim), ...
        'coeffs', coeffs, ...
        'stageZ', stageZ, 'emitterHeight', emitterHeight);
end

function s = radialWeightModeString(sim)
    if isfield(sim,'interfaceRadialWeightMode') && ...
            ~isempty(sim.interfaceRadialWeightMode)
        s = char(sim.interfaceRadialWeightMode);
    else
        s = 'sampleSolidAngle';
    end
end

function [nStage, cosStage] = stageMediumScalar(sim, nImm, cosImm, nGlass, cosGlass)
    switch stageMediumStringScalar(sim)
        case 'glass'
            nStage = nGlass; cosStage = cosGlass;
        otherwise
            nStage = nImm;   cosStage = cosImm;     % immersion (default)
    end
end

function s = stageMediumStringScalar(sim)
    s = 'immersion';
    if isfield(sim, 'airInterfaceStageMedium') && ~isempty(sim.airInterfaceStageMedium)
        s = lower(char(sim.airInterfaceStageMedium));
    end
    switch s
        case {'glass','coverslip'},               s = 'glass';
        case {'immersion','water','objective'},   s = 'immersion';
        otherwise,                                s = 'immersion';
    end
end

function val = airPsfCache(key,putVal)
%   LRU memo of whole PSF volumes. Get with one arg, put with two.
%   The cap must hold one fit iteration's PSFs (base + finite-difference
%   columns, x2 wavelengths ~ a dozen-plus) so the detector-pitch scan can
%   reuse the identical iteration-1 PSFs across pitch candidates -- the (b)
%   win. Each entry is a full [ny x nx x nConditions] volume, so this trades
%   memory for time; lower it if RAM is tight (the Bessel cache below is the
%   low-memory primary win and is unaffected). Flush with
%   `clear scalarPSFBesselAirInterface`.
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
