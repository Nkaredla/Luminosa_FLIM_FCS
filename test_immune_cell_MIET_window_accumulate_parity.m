function test_immune_cell_MIET_window_accumulate_parity()
%TEST_IMMUNE_CELL_MIET_WINDOW_ACCUMULATE_PARITY MEX vs MATLAB equivalence.
%
% The MEX scatter-add must produce bit-identical output to the MATLAB
% fallback. This runs the same synthetic acquisition both ways by moving the
% MEX aside, and compares the compact cube, anchor indices, intensity and the
% photon-contribution counters.
%
% The MEX is restored by onCleanup, so a failure cannot leave it disabled.

    mexName = 'immune_cell_MIET_window_accumulate';
    here = fileparts(mfilename('fullpath'));
    mexFile = fullfile(here, [mexName '.' mexext]);
    if ~isfile(mexFile)
        fprintf(['%s: MEX not built, nothing to compare. Build it with\n' ...
                 '  mex -O %s.c\n'], mfilename, mexName);
        return;
    end

    [ptu, ismResult, imageSize] = syntheticAcquisition();
    mask = true(imageSize);
    mask(1, 1) = false;          % exercise a de-selected anchor
    mask(end, end) = false;

    windows = {[2 2], [4 4]};
    blockSizes = [7 1e6];        % force multi-block and single-block paths

    for w = 1:numel(windows)
        for b = 1:numel(blockSizes)
            windowSize = windows{w};
            blockSize = blockSizes(b);

            % --- with the MEX ---
            clear(mexName);
            assert(exist(mexName, 'file') == 3, ...
                'MEX should be visible before it is moved aside.');
            [cubeMex, idxMex, intMex, metaMex] = ...
                immune_cell_MIET_reassigned_sliding_tcspc( ...
                ptu, ismResult, mask, windowSize, blockSize);

            % --- without the MEX ---
            hidden = [mexFile '.hidden'];
            movefile(mexFile, hidden);
            restore = onCleanup(@() restoreMex(hidden, mexFile));
            clear(mexName);
            rehash;
            assert(exist(mexName, 'file') ~= 3, ...
                'MEX still visible after being moved aside.');
            [cubeM, idxM, intM, metaM] = ...
                immune_cell_MIET_reassigned_sliding_tcspc( ...
                ptu, ismResult, mask, windowSize, blockSize);
            clear restore;   % triggers restoreMex
            clear(mexName);
            rehash;

            label = sprintf('window [%d %d], block %g', ...
                windowSize(1), windowSize(2), blockSize);
            assert(isequal(cubeMex, cubeM), ...
                'Compact cube differs between MEX and MATLAB (%s).', label);
            assert(isequal(idxMex, idxM), ...
                'Anchor indices differ (%s).', label);
            assert(isequal(intMex, intM), ...
                'Intensity differs (%s).', label);
            assert(metaMex.windowPhotonContributionCount == ...
                   metaM.windowPhotonContributionCount, ...
                'windowPhotonContributionCount differs (%s): %g vs %g.', ...
                label, metaMex.windowPhotonContributionCount, ...
                metaM.windowPhotonContributionCount);
            assert(metaMex.uniqueContributingPhotonCount == ...
                   metaM.uniqueContributingPhotonCount, ...
                'uniqueContributingPhotonCount differs (%s): %g vs %g.', ...
                label, metaMex.uniqueContributingPhotonCount, ...
                metaM.uniqueContributingPhotonCount);

            fprintf('  %-28s identical (%d anchors, %g contributions)\n', ...
                label, numel(idxMex), metaMex.windowPhotonContributionCount);
        end
    end

    fprintf('test_immune_cell_MIET_window_accumulate_parity: PASS\n');
end

function restoreMex(hidden, mexFile)
    if isfile(hidden)
        if isfile(mexFile)
            delete(hidden);
        else
            movefile(hidden, mexFile);
        end
    end
end

function [ptu, ismResult, imageSize] = syntheticAcquisition()
    rng(11);
    imageSize = [12 14];
    gateLength = 8;
    photonCount = 4000;

    % Two detectors with distinct sub-pixel shifts, so the reassignment and
    % the detector lookup are both exercised.
    channelIDs = [0 3];
    shiftsToCenter = [0 0; 0.5 -0.5];

    rowList = randi(imageSize(1), photonCount, 1);
    columnList = randi(imageSize(2), photonCount, 1);
    timeList = randi(gateLength, photonCount, 1);
    channelList = channelIDs(randi(2, photonCount, 1)).';

    ptu = struct();
    ptu.head = struct('ImgHdr_PixX', imageSize(1), ...
        'ImgHdr_PixY', imageSize(2), 'MeasDesc_Resolution', 0.1e-9);
    ptu.Ngate = gateLength;
    ptu.Resolution_ns = 0.1;
    ptu.im_col = uint16(rowList);
    ptu.im_line = uint16(columnList);
    ptu.im_tcspc = uint16(timeList);
    ptu.im_chan = uint8(channelList);

    ismResult = struct('aprImage', zeros(imageSize), ...
        'channelIDs', channelIDs, 'shiftsToCenter', shiftsToCenter);
end
