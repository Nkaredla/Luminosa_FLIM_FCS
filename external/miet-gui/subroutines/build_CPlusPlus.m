% Build the mex files for DeadTimeCorrection

% Windows with Microsoft Visual Studio compiler
if ispc
% Note: The /MT option statically links in the runtime, so it does not need to be deployed with the application.    
mex OPTIMFLAGS="/O2 /Oy- /DNDEBUG /MT" ...
    computeWeightingFunction.cpp

% Unix (assuming gcc > 4.5)
elseif isunix 
    mex OPTIMFLAGS="-O2 -static-libgcc -static-libstdc++" ...
    computeWeightingFunction.cpp
end