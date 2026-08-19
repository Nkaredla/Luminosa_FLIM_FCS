function test_immune_cell_MIET_sliding_synthetic
%TEST_IMMUNE_CELL_MIET_SLIDING_SYNTHETIC Verify 2x2/step-1 accumulation.

    imageSize = [4 5];
    gateLength = 3;
    nativeCube = zeros([imageSize gateLength], 'uint16');
    rowList = [];
    columnList = [];
    timeList = [];
    channelList = [];
    for column = 1:imageSize(2)
        for row = 1:imageSize(1)
            for timeBin = 1:gateLength
                count = mod(row + 2 * column + timeBin, 3) + 1;
                nativeCube(row, column, timeBin) = uint16(count);
                rowList = [rowList; repmat(row, count, 1)]; %#ok<AGROW>
                columnList = [columnList; repmat(column, count, 1)]; %#ok<AGROW>
                timeList = [timeList; repmat(timeBin, count, 1)]; %#ok<AGROW>
                channelList = [channelList; zeros(count, 1)]; %#ok<AGROW>
            end
        end
    end

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
        'channelIDs', 0, 'shiftsToCenter', [0 0]);

    [compact, indices, intensity, meta] = ...
        immune_cell_MIET_reassigned_sliding_tcspc( ...
        ptu, ismResult, true(imageSize), [2 2], 17);
    expectedIndices = find([true(3, 4) false(3, 1); false(1, 5)]);
    assert(isequal(indices, expectedIndices));
    for compactIndex = 1:numel(indices)
        [row, column] = ind2sub(imageSize, indices(compactIndex));
        expected = squeeze(sum(sum(uint32(nativeCube( ...
            row:row+1, column:column+1, :)), 1), 2));
        observed = squeeze(uint32(compact(compactIndex, 1, :)));
        assert(isequal(observed, expected));
        assert(intensity(row, column) == sum(expected));
    end
    assert(all(intensity(end, :) == 0));
    assert(all(intensity(:, end) == 0));
    assert(isequal(meta.windowCentreOffsetPixels, [0.5 0.5]));
    assert(meta.overlappingWindowsSharePhotons);

    [compact4, indices4, intensity4, meta4] = ...
        immune_cell_MIET_reassigned_sliding_tcspc( ...
        ptu, ismResult, true(imageSize), [4 4], 13);
    expectedIndices4 = find([true(1, 2) false(1, 3); false(3, 5)]);
    assert(isequal(indices4, expectedIndices4));
    for compactIndex = 1:numel(indices4)
        [row, column] = ind2sub(imageSize, indices4(compactIndex));
        expected = squeeze(sum(sum(uint32(nativeCube( ...
            row:row+3, column:column+3, :)), 1), 2));
        observed = squeeze(uint32(compact4(compactIndex, 1, :)));
        assert(isequal(observed, expected));
        assert(intensity4(row, column) == sum(expected));
    end
    assert(isequal(meta4.windowCentreOffsetPixels, [1.5 1.5]));
    assert(meta4.windowArea == 16);
    fprintf('test_immune_cell_MIET_sliding_synthetic: PASS\n');
end
