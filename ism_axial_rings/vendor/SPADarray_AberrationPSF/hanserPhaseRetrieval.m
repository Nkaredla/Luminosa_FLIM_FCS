function res = hanserPhaseRetrieval(stack, planeZ, params)
%HANSERPHASERETRIEVAL Retrieve a pupil function from a 3-D PSF stack (Hanser).
%
%   res = hanserPhaseRetrieval(stack, planeZ, params)
%
%   Iterative (Gerchberg-Saxton / Hanser-Gustafsson) phase retrieval of the
%   complex pupil from a measured 3-D intensity PSF stack:
%
%     A_z = F^{-1}{ P . exp(i 2*pi k_z z) }          pupil -> image amplitude
%     A_z'= |I_z|^{1/2} . exp(i arg A_z)             enforce measured magnitude
%     P_z = F{ A_z' } . exp(-i 2*pi k_z z)           image -> pupil, de-focus
%     P   = aperture .* mean_z(P_z)                   average + aperture support
%
%   The converged pupil phase is the effective wavefront; it is decomposed
%   onto Zernike modes (same basis/units as the forward model). Piston, tilt,
%   and defocus are fitted (to absorb centering/focus) but only the requested
%   modes are reported.
%
%   params.propagatorMode:
%     'homogeneous'  (default) free-amplitude retrieval with a single-medium
%                    defocus propagator k_z = sqrt((nMedium/lambda)^2 - k_r^2).
%     'airInterface' constrain the pupil AMPLITUDE to the known water/glass/air
%                    apodization (sqrt(cosImm) * Fresnel(s/p) / cosSample Weyl
%                    Jacobian) -- which gives the broad effective-NA~1 PSF -- and
%                    retrieve only the phase. The reported coefficients are then
%                    the instrument aberration RELATIVE to the known interface
%                    phase (cover/height/Fresnel), so the retrieval no longer
%                    underestimates by absorbing the interface into a generic
%                    pupil. Needs nImmersion, nGlass, nSample, nDesignGlass,
%                    coverslipUm, designCoverslipUm, emitterHeightUm.
%
%   Effective NA: the geometric NA (1.2) is rarely the effective collection
%   aperture at a glass/air interface -- total internal reflection caps the
%   propagating cone near nSample (~1.0) and only thin-layer SAF reaches the
%   full NA. params.effectiveNA sets the pupil-support radius directly; with
%   params.fitEffectiveNA the retrieval is run over params.effectiveNACandidates
%   (or an auto 0.8..NA scan) and the aperture giving the lowest PSF residual is
%   kept. A smaller effective NA broadens the model PSF to match the data.
%
%   INPUTS
%     stack   : [ny nx nPlane] measured PSF intensity (background subtracted,
%               registered). Planes should already be the kept (non-shifted)
%               ones; the bead is re-centred internally.
%     planeZ  : [1 nPlane] plane z relative to focus [um] (z=0 at focus).
%     params  : struct with fields
%                 NA, nMedium, wavelengthUm, dxUm        (optics / sampling)
%                 lamRefUm    reference wavelength for coeff units [um]
%                 effectiveNA          fixed pupil-support NA ([] = NA)
%                 fitEffectiveNA       scan the aperture for the best NA (false)
%                 effectiveNACandidates  NA values to scan ([] = auto)
%                 fitModes    cellstr of modes to report
%                 padSize     FFT grid size (default 96)
%                 nIter       iterations (default 40)
%                 verbose     (default false)
%
%   OUTPUT res: pupil, pupilPhase, aperture, effectiveNA, naScan (NA vs
%               relResidual), coeffs (struct), coeffVector, modelStack
%               (|A|^2 in the input frame), relResidual, errorHistory.

    if nargin < 3, params = struct(); end
    params = fillDefaults(params);
    [ny0, nx0, nP] = size(stack);
    planeZ = planeZ(:).';
    if numel(planeZ) ~= nP
        error('hanserPhaseRetrieval:PlaneZ', 'planeZ must have one z per plane.');
    end

    % --- centre the bead on the focal-plane centroid (removes a tilt term) --
    [stackC, shiftXY] = centreStack(stack);

    % --- pad to the FFT grid -----------------------------------------------
    N = params.padSize;
    I = padStack(stackC, N);                       % [N N nP]
    planeEnergy = zeros(1, nP);
    for p = 1:nP
        e = sum(sum(I(:,:,p)));
        planeEnergy(p) = e;
        if e > 0, I(:,:,p) = I(:,:,p) / e; end      % equal energy per plane (model)
    end
    measAmp = sqrt(max(I, 0));
    % Use every plane, but weight the pupil average by photon signal so the
    % dim, noisy far-from-focus planes do not get equal say as the bright ones.
    planeWeight = planeEnergy / max(sum(planeEnergy), eps);

    % --- pupil grid and defocus propagator (NA-independent) ----------------
    dk  = 1 / (N * params.dxUm);                    % cycles/um per pixel
    ax  = ((0:N-1) - floor(N/2)) * dk;
    [KX, KY] = meshgrid(ax, ax);
    KR  = hypot(KX, KY);
    kMed = params.nMedium / params.wavelengthUm;
    kz = sqrt(max(kMed^2 - KR.^2, 0));              % axial frequency (propagating)
    defocus = zeros(N, N, nP);
    for p = 1:nP
        defocus(:,:,p) = exp(1i * 2*pi * kz * planeZ(p));
    end

    % --- retrieve at one or several effective NAs; keep the best ------------
    % The air sample limits the effective collection aperture: TIR at the
    % glass/air interface caps the propagating cone near nSample (~1.0) and only
    % thin-layer SAF reaches the geometric NA. Scanning the aperture radius lets
    % the data pick the effective NA that matches the measured PSF width instead
    % of forcing NA=1.2 (which makes the model too compact).
    ctx = struct('KR',KR,'KX',KX,'KY',KY,'N',N,'nP',nP,'defocus',defocus, ...
        'measAmp',measAmp,'planeWeight',planeWeight,'ny0',ny0,'nx0',nx0, ...
        'shiftXY',shiftXY,'stack',stack);
    naList = resolveEffectiveNAList(params);
    best = []; scanNA = nan(1, numel(naList)); scanRes = nan(1, numel(naList));
    for i = 1:numel(naList)
        ri = retrieveAtNA(naList(i), ctx, params);
        scanNA(i) = naList(i); scanRes(i) = ri.relResidual;
        if params.verbose && numel(naList) > 1
            fprintf('[hanser] NA_eff %.3f -> rel.residual %.4f\n', ...
                naList(i), ri.relResidual);
        end
        if isempty(best) || ri.relResidual < best.relResidual
            best = ri;
        end
    end

    res = best;
    res.naScan = struct('NA', scanNA, 'relResidual', scanRes);
    res.centreShiftXY = shiftXY;
    res.params = params;
    if params.verbose
        fprintf('[hanser] best NA_eff %.3f, %d iters, final amp.err %.4f, rel.residual %.4f\n', ...
            res.effectiveNA, params.nIter, res.errorHistory(end), res.relResidual);
    end
end

% ========================================================================
function res = retrieveAtNA(naEff, ctx, params)
    % One Gerchberg-Saxton/Hanser retrieval at a fixed effective NA (aperture).
    KR = ctx.KR; N = ctx.N; nP = ctx.nP;
    kCut = naEff / params.wavelengthUm;             % aperture radius [cyc/um]
    aperture = KR <= kCut;

    % air-interface base pupil: the known water/glass/air apodization (Fresnel
    % x sqrt(cosImm) x Weyl Jacobian) gives the broad PSF, and the
    % cover/height/Fresnel phase is the known interface wavefront. In this mode
    % the retrieval constrains the amplitude to this base and retrieves only the
    % phase; the reported coefficients are the residual instrument aberration.
    isAir = strcmpi(params.propagatorMode, 'airInterface');
    if isAir
        [basePupil, baseAmp, basePhase] = airInterfaceBasePupil(KR, aperture, params);
    else
        basePupil = double(aperture); baseAmp = double(aperture); basePhase = zeros(N);
    end

    P = basePupil;                                  % init at the base pupil
    errorHistory = nan(1, params.nIter);
    for it = 1:params.nIter
        accum = zeros(N, N);
        sse = 0;
        for p = 1:nP
            amp    = pup2img(P .* ctx.defocus(:,:,p));
            ampNew = ctx.measAmp(:,:,p) .* exp(1i * angle(amp));
            back   = img2pup(ampNew) .* conj(ctx.defocus(:,:,p));
            accum  = accum + ctx.planeWeight(p) * back;     % signal-weighted average
            mhat = abs(amp);
            mhat = mhat / max(sqrt(sum(mhat(:).^2)), eps);  % unit L2, like measAmp
            d = mhat - ctx.measAmp(:,:,p);
            sse = sse + ctx.planeWeight(p) * sum(d(:).^2);
        end
        if isAir
            P = baseAmp .* exp(1i*angle(accum)) .* aperture;  % fix amplitude, retrieve phase
        else
            P = accum .* aperture;                            % free amplitude + phase
        end
        errorHistory(it) = sqrt(sse);
    end

    % --- decompose pupil phase onto Zernikes --------------------------------
    phase = angle(P);
    if isAir
        aberrationPhase = wrapToPiLocal(phase - basePhase);   % residual instrument phase
        aberrationPhase(~aperture) = 0;
    else
        aberrationPhase = phase;
    end
    [coeffs, coeffVec, basisInfo] = decomposeZernikePhase( ...
        aberrationPhase, ctx.KX, ctx.KY, kCut, aperture, params);

    % --- forward model in the input frame, for comparison -------------------
    modelStack = forwardModelStack(P, ctx.defocus, N, ctx.ny0, ctx.nx0, ctx.shiftXY);
    rel = relResidual(ctx.stack, modelStack);

    res = struct();
    res.pupil = P;
    res.pupilPhase = phase;
    res.propagatorMode = params.propagatorMode;
    res.effectiveNA = naEff;
    res.aperture = aperture;
    if isAir
        res.aberrationPhase = aberrationPhase;
        res.basePhase = basePhase;
        res.basePupil = basePupil;
    end
    res.coeffs = coeffs;
    res.coeffVector = coeffVec;
    res.fitModes = params.fitModes;
    res.coeffNmRms = coeffVec * params.lamRefUm * 1000;
    res.modelStack = modelStack;
    res.relResidual = rel;
    res.errorHistory = errorHistory;
    res.basisInfo = basisInfo;
end

function naList = resolveEffectiveNAList(params)
    if params.fitEffectiveNA
        if ~isempty(params.effectiveNACandidates)
            naList = params.effectiveNACandidates(:).';
        else
            naList = linspace(0.8, params.NA, 9);   % auto scan up to geometric NA
        end
    elseif ~isempty(params.effectiveNA)
        naList = params.effectiveNA(:).';
    else
        naList = params.NA;
    end
    naList = naList(isfinite(naList) & naList > 0 & naList <= params.NA);
    if isempty(naList), naList = params.NA; end
end

% ========================================================================
function params = fillDefaults(params)
    d = struct('NA',1.2,'nMedium',1.33,'wavelengthUm',0.488,'dxUm',0.05, ...
        'lamRefUm',0.520,'padSize',96,'nIter',40,'verbose',false, ...
        'propagatorMode','homogeneous', ...   % 'homogeneous' | 'airInterface'
        'effectiveNA',[],'fitEffectiveNA',false,'effectiveNACandidates',[], ...
        'nImmersion',1.33,'nGlass',1.518,'nSample',1.0003,'nDesignGlass',1.518, ...
        'coverslipUm',190,'designCoverslipUm',190,'emitterHeightUm',0, ...
        'fitModes',{{'astig_x','astig_y','coma_x','coma_y', ...
        'trefoil_x','trefoil_y','spherical', ...
        'secondary_astig_x','secondary_astig_y', ...
        'quadrafoil_x','quadrafoil_y'}});
    fn = fieldnames(d);
    for k = 1:numel(fn)
        if ~isfield(params, fn{k}) || isempty(params.(fn{k}))
            params.(fn{k}) = d.(fn{k});
        end
    end
    if ischar(params.fitModes) || isstring(params.fitModes)
        params.fitModes = cellstr(params.fitModes);
    end
end

function [basePupil, baseAmp, basePhase] = airInterfaceBasePupil(KR, aperture, params)
    % Known water/glass/air Cartesian effective pupil: aplanatic apodization
    % sqrt(cosImm), s/p-averaged Fresnel transmission (sample->glass->immersion),
    % the cover-glass and emitter-height phases, and the 1/cosSample Weyl
    % Jacobian. Its magnitude is the broad effective-NA~1 apodization; its phase
    % is the known interface wavefront.
    lambda = params.wavelengthUm; k0 = 2*pi/lambda;
    q = KR * lambda;                                  % = k_parallel / k0
    nImm = params.nImmersion; nGlass = params.nGlass;
    nSample = params.nSample; nDesign = params.nDesignGlass;
    cosImm = posRootLocal(1 - (q/nImm).^2);
    cosGlass = posRootLocal(1 - (q/nGlass).^2);
    cosSample = posRootLocal(1 - (q/nSample).^2);
    cosDesign = posRootLocal(1 - (q/nDesign).^2);
    [tpSG, tsSG] = fresnelTLocal(nSample, cosSample, nGlass, cosGlass);
    [tpGI, tsGI] = fresnelTLocal(nGlass, cosGlass, nImm, cosImm);
    t = 0.5*(tpSG.*tpGI + tsSG.*tsGI);                % scalar s/p average
    coverPhase = k0*(nGlass.*cosGlass*params.coverslipUm - ...
        nDesign.*cosDesign*params.designCoverslipUm);
    heightFactor = exp(1i*k0*nSample.*cosSample*params.emitterHeightUm);
    cosSampleSafe = cosSample;
    cosSampleSafe(abs(cosSampleSafe) < 1e-6) = 1e-6;
    basePupil = sqrt(cosImm) .* t .* exp(1i*coverPhase) .* heightFactor ./ cosSampleSafe;
    basePupil(~aperture) = 0;
    mx = max(abs(basePupil(:)));
    if mx > 0, basePupil = basePupil / mx; end
    baseAmp = abs(basePupil);
    basePhase = angle(basePupil);
end

function [tp, ts] = fresnelTLocal(n1, c1, n2, c2)
    tp = 2*n1.*c1 ./ (n2.*c1 + n1.*c2);
    ts = 2*n1.*c1 ./ (n1.*c1 + n2.*c2);
end

function r = posRootLocal(v)
    r = sqrt(complex(v, 0));
    r(imag(r) < 0) = conj(r(imag(r) < 0));
end

function w = wrapToPiLocal(p)
    w = mod(p + pi, 2*pi) - pi;
end

function A = pup2img(P)
    A = fftshift(ifft2(ifftshift(P)));
end

function P = img2pup(A)
    P = fftshift(fft2(ifftshift(A)));
end

function [out, shiftXY] = centreStack(stack)
    [ny, nx, nP] = size(stack);
    [~, fp] = max(squeeze(sum(sum(stack,1),2)));    % focal plane
    w = stack(:,:,fp); w = max(w - 0.1*max(w(:)), 0);
    [X, Y] = meshgrid(1:nx, 1:ny);
    s = sum(w(:));
    if s > 0
        cx = sum(sum(w.*X))/s; cy = sum(sum(w.*Y))/s;
    else
        cx = (nx+1)/2; cy = (ny+1)/2;
    end
    shiftX = (nx+1)/2 - cx; shiftY = (ny+1)/2 - cy;
    shiftXY = [shiftX, shiftY];
    out = zeros(size(stack));
    for p = 1:nP
        out(:,:,p) = max(fourierShift2D(stack(:,:,p), shiftX, shiftY), 0);
    end
end

function I = padStack(stack, N)
    [ny, nx, nP] = size(stack);
    I = zeros(N, N, nP);
    r0 = floor((N - ny)/2); c0 = floor((N - nx)/2);
    if r0 < 0 || c0 < 0
        % crop instead if the data is larger than the pad grid
        rr = floor((ny - N)/2) + (1:N);
        cc = floor((nx - N)/2) + (1:N);
        I = stack(rr, cc, :);
        return;
    end
    I(r0+(1:ny), c0+(1:nx), :) = stack;
end

function model = forwardModelStack(P, defocus, N, ny0, nx0, shiftXY)
    nP = size(defocus, 3);
    model = zeros(ny0, nx0, nP);
    r0 = floor((N - ny0)/2); c0 = floor((N - nx0)/2);
    for p = 1:nP
        A = pup2img(P .* defocus(:,:,p));
        Ipad = abs(A).^2;
        if r0 >= 0 && c0 >= 0
            crop = Ipad(r0+(1:ny0), c0+(1:nx0));
        else
            rr = floor((N - ny0)/2); % unlikely branch
            crop = Ipad(max(rr,0)+(1:ny0), max(c0,0)+(1:nx0));
        end
        % undo the centring shift so the model aligns with the input frame
        model(:,:,p) = max(fourierShift2D(crop, -shiftXY(1), -shiftXY(2)), 0);
    end
end

function r = relResidual(meas, model)
    num = 0; den = 0;
    for p = 1:size(meas,3)
        m = meas(:,:,p); f = model(:,:,p);
        sm = sum(m(:)); sf = sum(f(:));
        if sm > 0, m = m/sm; end
        if sf > 0, f = f/sf; end
        num = num + sum((m(:)-f(:)).^2);
        den = den + sum(m(:).^2);
    end
    r = sqrt(num / max(den, eps));
end

function [coeffs, coeffVec, info] = decomposeZernikePhase(phase, KX, KY, ...
        kCut, aperture, params)
    rho = hypot(KX, KY) / kCut;
    phi = atan2(KY, KX);
    Z = zernikeModes(rho, phi);                     % repo basis (no piston)
    modeNames = zernikeModeNames(6);
    idx = aperture(:);
    % phase (radians at wavelength) -> waves at lamRef
    W = phase * params.wavelengthUm / (2*pi*params.lamRefUm);
    B = ones(sum(idx), numel(modeNames)+1);          % column 1 = piston
    for k = 1:numel(modeNames)
        zk = Z.(modeNames{k});
        B(:, k+1) = zk(idx);
    end
    c = B \ W(idx);                                  % least-squares fit
    full = struct();
    for k = 1:numel(modeNames)
        full.(modeNames{k}) = c(k+1);
    end
    coeffs = struct();
    coeffVec = zeros(1, numel(params.fitModes));
    for k = 1:numel(params.fitModes)
        m = params.fitModes{k};
        v = 0; if isfield(full, m), v = full.(m); end
        coeffs.(m) = v; coeffVec(k) = v;
    end
    info = struct('allCoeffs', full, 'pistonWaves', c(1), ...
        'modeNames', {modeNames});
end
