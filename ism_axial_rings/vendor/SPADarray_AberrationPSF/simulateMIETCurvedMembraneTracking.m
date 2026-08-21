function out = simulateMIETCurvedMembraneTracking(varargin)
%SIMULATEMIETCURVEDMEMBRANETRACKING Curved-membrane MIET/ISM particle tracking.
%
%   out = simulateMIETCurvedMembraneTracking()
%   out = simulateMIETCurvedMembraneTracking('curvaturePerUm',-0.8)
%   out = simulateMIETCurvedMembraneTracking('writeOutputs',true)
%
%   This function simulates and analyzes membrane-bound molecules on the
%   rotationally symmetric local surface
%
%       z(x,y) = zTip + 0.5 * kappa * (x^2 + y^2).
%
%   The symmetry axis, membrane tip, focused excitation, and center of the
%   ISM detector all coincide at x=y=0. Signed kappa is the common principal
%   curvature at the tip: kappa<0 makes the centered tip a local maximum.
%
%   Brownian motion is generated intrinsically on the surface using the
%   inverse surface metric and the Ito drift of the Laplace-Beltrami
%   operator. It is not a planar random walk followed by axial projection.
%
%   The requested local linear MIET calibration is
%
%       lifetimeNs = lifetimeAtSubstrateNs + lifetimeSlopeNsPerUm * zUm.
%
%   Signal TCSPC microtimes follow a repetition-period-truncated exponential.
%   Background microtimes are uniform. Lateral position is fitted from the
%   finite-channel ISM microimage, lifetime is fitted from the channel and
%   microtime mixture likelihood, and height follows from the calibration.
%   The curvature fit pools all accepted photon microtimes. Diffusion is
%   estimated in local tangent coordinates using the fitted surface metric.
%
%   The low-concentration analysis fits at most one emitter per frame. Frames
%   with a poor single-emitter Poisson fit are rejected. Keep nMolecules low
%   enough that simultaneous occupancy of the focus remains uncommon.
%
%   Important outputs:
%     out.photonEvents       channel-resolved TCSPC photon table
%     out.localizationTable per-frame x, y, lifetime, height, and Fisher CRBs
%     out.curvature          fitted tip height and signed apex curvature
%     out.trackTable         nearest-neighbor low-density trajectories
%     out.diffusion          surface MSD and Fisher/MLE diffusion estimates
%
%   The image model is stroboscopic: one position is used for each frame.
%   Intra-frame motion blur and a finite instrument response function are not
%   included. All distances are sample-equivalent micrometers.
%
%   Implementation files: +membrane_tracking/+curved_miet

    out = membrane_tracking.curved_miet.runSimulation(varargin{:});
end
