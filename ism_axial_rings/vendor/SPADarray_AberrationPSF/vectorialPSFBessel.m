function psf = vectorialPSFBessel(sim, coeffs, wavelength)
%VECTORIALPSFBESSEL Vectorial Richards-Wolf/Bessel homogeneous PSF.
%
%   psf(y,x,z) is the orientation-averaged intensity from x- and y-polarized
%   entrance-pupil fields. Coefficients are RMS waves at sim.lamRef.

    validateScalarOptics(sim);
    coeffs = coeffStruct(sim, coeffs);
    if ~isscalar(wavelength) || ~isfinite(wavelength) || wavelength <= 0
        error('vectorialPSFBessel:InvalidWavelength', ...
            'wavelength must be a finite positive scalar in micrometers.');
    end

    rho = linspace(0, 1, sim.Nr);
    phi = linspace(0, 2*pi, sim.Nphi+1); phi(end) = [];
    [R, P] = ndgrid(rho, phi);

    W = wavefrontWaves(R, P, coeffs);
    aberrationPhase = exp(1i*2*pi*(sim.lamRef/wavelength).*W);

    sinT = (sim.NA / sim.nMedium) * rho;
    cosT = positiveRoot(1 - sinT.^2);
    apod = sqrt(cosT(:));
    common = apod .* aberrationPhase;

    c = cos(P); s = sin(P);
    ct = cosT(:); st = sinT(:);
    pupil = vectorialPupilComponents(common, c, s, ct, st);

    ms = -sim.M:sim.M;
    pm = componentSpectra(pupil, ms, sim.Nphi);

    [X, Y] = meshgrid(sim.x, sim.y);
    radius = hypot(X, Y);
    psi = atan2(Y, X);
    rv = radius(:);
    roundedRadius = round(rv, 12);
    uniqueRadius = unique(roundedRadius);
    [~, radiusIndex] = ismember(roundedRadius, uniqueRadius);

    k0 = 2*pi / wavelength;
    argument = k0 * sim.NA * (uniqueRadius * rho);
    besselTables = besselTablesForVectorial( ...
        wavelength, sim.NA, sim.Nr, uniqueRadius, sim.M, argument);

    quadrature = trapezoidWeights(rho);
    smax = sim.NA / sim.nMedium;
    radialWeight = (smax^2) * rho ./ max(abs(cosT), 1e-12) .* quadrature;

    polMode = polarizationModeString(sim);
    psf = zeros(numel(sim.y), numel(sim.x), numel(sim.z));
    for iz = 1:numel(sim.z)
        axial = exp(1i * k0 * sim.nMedium * cosT(:) * sim.z(iz));
        ix = synthesizeIntensity(pm.x, ms, radialWeight, axial, ...
            besselTables, radiusIndex, psi);
        iy = synthesizeIntensity(pm.y, ms, radialWeight, axial, ...
            besselTables, radiusIndex, psi);
        psf(:,:,iz) = reshape(combinePolarizations(ix, iy, polMode), ...
            numel(sim.y), numel(sim.x));
    end

    total = sum(psf(:));
    if total > 0, psf = psf / total; end
end

function W = wavefrontWaves(R, P, coeffs)
    Z = zernikeModes(R, P);
    W = zeros(size(R));
    fn = fieldnames(coeffs);
    for k = 1:numel(fn)
        if ~isfield(Z, fn{k})
            error('vectorialPSFBessel:UnknownMode', ...
                'Unknown Zernike mode "%s".', fn{k});
        end
        W = W + coeffs.(fn{k}) * Z.(fn{k});
    end
end

function pupil = vectorialPupilComponents(common, c, s, ct, st)
    pupil.x.ex = common .* (ct .* c.^2 + s.^2);
    pupil.x.ey = common .* ((ct - 1) .* s .* c);
    pupil.x.ez = common .* (-st .* c);

    pupil.y.ex = common .* ((ct - 1) .* s .* c);
    pupil.y.ey = common .* (ct .* s.^2 + c.^2);
    pupil.y.ez = common .* (-st .* s);
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

function pm = azimuthalCoefficients(pupil, ms, nPhi)
    spectrum = fft(pupil, [], 2) / nPhi;
    pm = complex(zeros(size(pupil,1), numel(ms)));
    for k = 1:numel(ms)
        pm(:,k) = spectrum(:, mod(ms(k), nPhi) + 1);
    end
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

function s = polarizationModeString(sim)
    if isfield(sim, 'vectorialPolarizationMode') && ...
            ~isempty(sim.vectorialPolarizationMode)
        s = char(sim.vectorialPolarizationMode);
    else
        s = 'xyAverage';
    end
end

function I = combinePolarizations(ix, iy, mode)
    switch lower(mode)
        case {'xyaverage','unpolarized','isotropic','circular'}
            I = 0.5*(ix + iy);
        case 'x'
            I = ix;
        case 'y'
            I = iy;
        otherwise
            error('vectorialPSFBessel:BadPolarizationMode', ...
                'vectorialPolarizationMode must be xyAverage, x, or y.');
    end
end

function root = positiveRoot(value)
    root = sqrt(complex(value, 0));
    root(imag(root) < 0) = conj(root(imag(root) < 0));
end

function J = besselTablesForVectorial(wavelength, NA, Nr, uniqueRadius, M, argument)
    persistent keys vals
    if isempty(keys), keys = {}; vals = {}; end
    for i = 1:numel(keys)
        ki = keys{i};
        if ki.wavelength == wavelength && ki.NA == NA && ki.Nr == Nr && ...
                isequal(ki.uniqueRadius, uniqueRadius) && ki.M == M
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
