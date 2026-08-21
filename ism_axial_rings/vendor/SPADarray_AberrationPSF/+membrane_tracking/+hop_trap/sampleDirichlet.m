function sample = sampleDirichlet(concentration)
%SAMPLEDIRICHLET Sample a Dirichlet vector using core MATLAB randg.

    concentration = max(double(concentration(:).'), 1e-12);
    sample = randg(concentration);
    total = sum(sample);
    if ~(total > 0) || ~isfinite(total)
        sample = concentration / sum(concentration);
    else
        sample = sample / total;
    end
end
