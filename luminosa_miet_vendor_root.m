function vendorRoot = luminosa_miet_vendor_root()
% Resolve the preferred MIET-GUI vendor checkout shipped alongside Luminosa.

    baseDir = fileparts(mfilename('fullpath'));
    candidates = { ...
        fullfile(baseDir, 'external', 'miet-gui'), ...
        fullfile(baseDir, 'vendor', 'miet-gui')};

    vendorRoot = '';
    for idx = 1:numel(candidates)
        if isfolder(candidates{idx})
            vendorRoot = candidates{idx};
            addpath(genpath(vendorRoot), '-end');
            return;
        end
    end
end
