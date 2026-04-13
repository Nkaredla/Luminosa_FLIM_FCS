function mode = luminosa_miet_normalize_orientation(mode)
% Canonicalize dipole-orientation labels used by the local MIET wrappers.

    if nargin == 0 || isempty(mode)
        mode = 'fast_rotating';
        return;
    end

    if isnumeric(mode) || islogical(mode)
        mode = luminosa_miet_orientation_from_value(double(mode));
        return;
    end

    mode = lower(strtrim(char(mode)));
    mode = strrep(mode, '-', ' ');
    mode = strrep(mode, '_', ' ');
    mode = regexprep(mode, '\s+', ' ');

    switch mode
        case {'fast rotating', 'fastly rotating', 'freely rotating', 'randomly oriented', 'random orientation', 'random'}
            mode = 'fast_rotating';
        case {'random fixed', 'fixed random', 'randomly fixed'}
            mode = 'random_fixed';
        case 'parallel'
            mode = 'parallel';
        case 'vertical'
            mode = 'vertical';
        otherwise
            mode = 'fast_rotating';
    end
end
