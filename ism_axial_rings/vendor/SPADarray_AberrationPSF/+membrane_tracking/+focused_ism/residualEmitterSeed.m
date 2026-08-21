function seedXY = residualEmitterSeed(data, mu, rows, cols, existingXY, ...
        detector, opts)
    import membrane_tracking.focused_ism.*

    residual = data - mu;
    [~, order] = sort(residual(:), 'descend');
    minSeparationUm = opts.minResolvedEmitterSeparationPx * ...
        detector.pixelSizeUm;
    seedXY = [];
    for k = 1:numel(order)
        if residual(order(k)) <= max(1, ...
                opts.detectionThresholdSigma * sqrt(max(mu(order(k)), 1)))
            break;
        end
        [localRow, localCol] = ind2sub(size(data), order(k));
        candidateXY = [detector.xCenters(cols(localCol)), ...
            detector.yCenters(rows(localRow))];
        if isempty(existingXY) || all(hypot( ...
                existingXY(:,1) - candidateXY(1), ...
                existingXY(:,2) - candidateXY(2)) >= minSeparationUm)
            seedXY = candidateXY;
            return;
        end
    end
end
