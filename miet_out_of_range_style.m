function style = miet_out_of_range_style()
%MIET_OUT_OF_RANGE_STYLE Shared colours and wording for pixels without a MIET height.
%
% Used by MIET_HEIGHT_STATUS_RGB and every figure legend, so the GUI and the
% batch figures describe out-of-range pixels identically.

    style = struct();
    style.aboveColour = [1 1 1];
    style.belowColour = [0.5 0.5 0.5];
    style.aboveLabel = 'white: above calibration ceiling (beyond MIET range)';
    style.belowLabel = 'grey: below calibration floor';
end
