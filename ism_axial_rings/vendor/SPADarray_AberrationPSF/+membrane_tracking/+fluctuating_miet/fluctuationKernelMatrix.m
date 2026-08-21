function K = fluctuationKernelMatrix(gpData, logParameters, opts)
    import membrane_tracking.fluctuating_miet.*

%FLUCTUATIONKERNELMATRIX Prior covariance of the membrane height field.
%
%   'sqexp'    k = A^2 * exp(-d^2/(2*L^2)) * exp(-|dt|/T)
%              This is a useful empirical alternative.
%   'helfrich' uses the exact discrete covariance of the simulation:
%       sum_m v_m cos(q_m.(x-x')) exp(-|t-t'|/tau_m) / sum_m v_m.
%              Only its overall rms amplitude is fitted.
    K = fluctuationKernelCross(gpData.position, gpData.timeS, ...
        gpData.position, gpData.timeS, logParameters, opts, gpData.modes);
    K = 0.5 * (K + K.');
end
