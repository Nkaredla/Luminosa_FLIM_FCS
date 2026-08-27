cd('D:\MATLAB\server\Luminosa\Luminosa_FLIM_FCS');
r = run_immune_cell_MIET_slb_anchored([], [1 2 4], ...
        struct('slbOpts', struct('doPerPixelNative', true)));
disp('ANCH_DONE');
