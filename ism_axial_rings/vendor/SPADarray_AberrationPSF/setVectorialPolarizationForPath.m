function simOut = setVectorialPolarizationForPath(simIn, pathName)
%SETVECTORIALPOLARIZATIONFORPATH Select vectorial polarization for one path.
%
%   sim.excitationPolarizationMode controls the illumination PSF.
%   sim.collectionPolarizationMode controls the detection/collection PSF.
%
%   If the path-specific field is absent, sim.vectorialPolarizationMode is
%   used as a legacy/global fallback. The final value is written back into
%   sim.vectorialPolarizationMode because vectorialPSFBessel* reads that
%   field internally.

    simOut = simIn;
    pathName = lower(char(pathName));
    switch pathName
        case {'excitation','illumination'}
            fieldName = 'excitationPolarizationMode';
            fallback = 'x';
        case {'collection','emission','detection'}
            fieldName = 'collectionPolarizationMode';
            fallback = 'xyAverage';
        otherwise
            error('setVectorialPolarizationForPath:BadPath', ...
                'pathName must be excitation or collection.');
    end

    if isfield(simIn, fieldName) && ~isempty(simIn.(fieldName))
        mode = simIn.(fieldName);
    elseif isfield(simIn, 'vectorialPolarizationMode') && ...
            ~isempty(simIn.vectorialPolarizationMode)
        mode = simIn.vectorialPolarizationMode;
    else
        mode = fallback;
    end
    simOut.vectorialPolarizationMode = char(mode);
end
