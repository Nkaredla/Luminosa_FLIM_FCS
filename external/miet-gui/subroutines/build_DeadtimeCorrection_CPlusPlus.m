% Build the mex files for DeadTimeCorrection using Microsoft Visual Studio compiler
%
% Note: The /MT option statically links in the runtime, so it does not need
% to be deployed with the application

mex OPTIMFLAGS="/O2 /Oy- /DNDEBUG /MT" ...
    circConv_BoxcarSet.cpp

% Note: Here we use openmp parallelization
mex COMPFLAGS="$COMPFLAGS /openmp /MT" ... ...
    OPTIMFLAGS="/O2 /Oy- /DNDEBUG" ...
    computeWW.cpp