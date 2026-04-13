function value = luminosa_miet_orientation_to_value(mode)
% Map a dipole-orientation label to the popup-menu index used by Luminosa.

    switch luminosa_miet_normalize_orientation(mode)
        case 'random_fixed'
            value = 2;
        case 'parallel'
            value = 3;
        case 'vertical'
            value = 4;
        otherwise
            value = 1;
    end
end
