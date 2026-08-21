function fisher = multiEmitterLocalizationFisher(model, rows, cols, ...
        detector, opts)
    import membrane_tracking.focused_ism.*

    nEmitters = model.nEmitters;
    nPixels = numel(rows) * numel(cols);
    nParameters = 3 * nEmitters + 1;
    J = zeros(nPixels, nParameters);
    mu = model.background * ones(numel(rows), numel(cols));

    for e = 1:nEmitters
        [probability, dPdx, dPdy] = gaussianPixelModel( ...
            model.xUm(e), model.yUm(e), opts.psfSigmaUm, ...
            detector.xEdges(cols(1):(cols(end)+1)), ...
            detector.yEdges(rows(1):(rows(end)+1)));
        mu = mu + model.amplitude(e) * probability;
        positionBase = 2 * (e - 1);
        J(:,positionBase + 1) = model.amplitude(e) * dPdx(:);
        J(:,positionBase + 2) = model.amplitude(e) * dPdy(:);
        J(:,2*nEmitters + e) = probability(:);
    end
    J(:,end) = 1;

    mu = max(mu(:), opts.minExpectedCount);
    information = J.' * bsxfun(@rdivide, J, mu);
    information = 0.5 * (information + information.');
    [covariance, ok] = invertPositiveDefinite(information);

    covarianceXY = nan(2, 2, nEmitters);
    if ok
        for e = 1:nEmitters
            indices = 2 * (e - 1) + (1:2);
            covarianceXY(:,:,e) = covariance(indices, indices);
        end
    end
    fisher = struct('informationFull', information, ...
        'covarianceFull', covariance, 'covarianceXY', covarianceXY, ...
        'isValid', ok);
end
