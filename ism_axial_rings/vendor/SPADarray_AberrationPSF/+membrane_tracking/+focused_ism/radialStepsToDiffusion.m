function D = radialStepsToDiffusion(radialSquaredDisplacement, dtS)
    import membrane_tracking.focused_ism.*

    if isempty(radialSquaredDisplacement)
        D = NaN;
    else
        D = mean(radialSquaredDisplacement) / (4 * dtS);
    end
end
