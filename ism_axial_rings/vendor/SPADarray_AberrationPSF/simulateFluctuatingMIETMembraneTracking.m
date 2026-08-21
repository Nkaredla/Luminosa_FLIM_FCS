function out = simulateFluctuatingMIETMembraneTracking(varargin)
%SIMULATEFLUCTUATINGMIETMEMBRANETRACKING MIET/ISM SPT on a fluctuating membrane.
%
%   out = simulateFluctuatingMIETMembraneTracking()
%   out = simulateFluctuatingMIETMembraneTracking('bendingModulusKT',20)
%
%   The membrane is a mean paraboloid plus a thermally fluctuating field:
%
%       z(r,t) = h0 + 0.5*kappa*|r|^2 + dh(r,t)
%       dh(r,t) = sum_m [ a_m(t) cos(q_m.r) + b_m(t) sin(q_m.r) ]
%
%   Mode variances follow the Helfrich spectrum. For the retained
%   half-plane real sine/cosine basis,
%       S(q) = 2 / ( A_patch * ( kappaBend*q^4 + tension*q^2 ) )
%   in units where kB*T = 1, so kappaBend is in kT and tension in kT/um^2.
%   Each mode relaxes as an independent Ornstein-Uhlenbeck process with
%   tau(q) proportional to 1/(kappaBend*q^3 + tension*q), normalised so the
%   fundamental mode relaxes in fluctuationRelaxationTimeS. Amplitudes are
%   advanced with the exact OU transition, not an Euler approximation.
%
%   The molecule performs intrinsic Brownian motion on the INSTANTANEOUS
%   surface. The Ito drift of the Laplace-Beltrami operator on a height
%   graph z=h(r) is
%
%       b = -p * ( trace(H)*s - p'*H*p ) / s^2 ,  p = grad h, H = hess h,
%       s = 1 + |p|^2
%
%   which reduces to the static-paraboloid drift when h is quadratic. This
%   identity was verified against a numerical divergence of
%   (1/sqrt(det g)) d_j ( sqrt(det g) g^{ij} ) to a relative error of 1e-7.
%
%   WHAT IS AND IS NOT IDENTIFIABLE
%   A fluctuation mode whose wavelength is comparable to the observation
%   DIAMETER is quadratic across the observed patch and is therefore
%   confounded with kappa. Shorter modes average out; longer modes are
%   absorbed into h0. Monte Carlo with a 10 nm rms field: the curvature
%   uncertainty peaks when the field correlation length equals the
%   observation radius, and an analysis that ignores field correlation
%   reports error bars 1.2x to 2.3x too small. This code therefore reports
%   an amplitude-profile interval for kappa in addition to the conditional
%   Fisher error bar.
%
%   Estimation differs from the static-membrane code in four ways:
%     1. Curvature uses a moment-corrected (regression-dilution) estimator.
%        Subtracting trace(C) from rhat^2 removes the MEAN bias but NOT the
%        attenuation. Monte Carlo: at 35 nm lateral precision the
%        mean-corrected estimator returns kappa 28% too small; the moment
%        correction recovers it to 0.3%.
%     2. The fluctuation field is marginalised as a Gaussian process rather
%        than ignored. The default kernel uses the same discrete Helfrich
%        modes and mode-specific relaxation times as the simulator.
%     3. Input noise enters the GP diagonal (NIGP), so heteroscedastic
%        localisation precision is handled.
%     4. The diffusion step likelihood is conditioned on the actual
%        anisotropic Mahalanobis tracking gate.
%
%   NOT MODELLED: instrument response function, intra-frame motion blur,
%   MIET calibration nonlinearity, membrane material-flow advection, or
%   more than one bright molecule per frame. Gate conditioning corrects the
%   geometric tracking gate, not focus/localisation selection.
%
%   See also SIMULATEMIETCURVEDMEMBRANETRACKING.
%
%   Implementation files: +membrane_tracking/+fluctuating_miet

    out = membrane_tracking.fluctuating_miet.runSimulation(varargin{:});
end
