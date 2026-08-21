function subset = subsampleForGp(n, maximum)
    import membrane_tracking.fluctuating_miet.*

% A dense GP is O(n^3). Thin uniformly in index, which preserves the
% temporal span and therefore the ability to identify the field
% correlation time.
    if n <= maximum
        subset = (1:n).';
    else
        subset = unique(round(linspace(1, n, maximum))).';
    end
end
