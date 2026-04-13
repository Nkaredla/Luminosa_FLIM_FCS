function mode = luminosa_miet_orientation_from_value(value)
% Map popup-menu index to the canonical local MIET dipole-orientation mode.

    switch round(double(value))
        case 2
            mode = 'random_fixed';
        case 3
            mode = 'parallel';
        case 4
            mode = 'vertical';
        otherwise
            mode = 'fast_rotating';
    end
end
