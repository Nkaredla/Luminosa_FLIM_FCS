function out = simulateMembraneDiffusionParticleTracking(varargin)
%SIMULATEMEMBRANEDIFFUSIONPARTICLETRACKING Focused-excitation ISM membrane SPT.
%
%   out = simulateMembraneDiffusionParticleTracking()
%   out = simulateMembraneDiffusionParticleTracking('diffusionUm2PerS',0.15)
%   out = simulateMembraneDiffusionParticleTracking('writeOutputs',true)
%
%   This is an end-to-end focused-ISM particle-tracking simulation:
%     1. Brownian diffusion of membrane-bound molecules in a 2-D box.
%     2. A stationary Gaussian excitation focus with saturating photophysics.
%     3. Photon assignment to finite hexagonal channels of a Luminosa ISM
%        detector using an inverted Gaussian detection PSF.
%     4. Poisson detector-channel microimages with background and dark counts.
%     5. BIC/Laplace model-order scoring plus joint multi-emitter ISM fits.
%     6. Joint Fisher-information / Cramer-Rao localization uncertainty.
%     7. JPDA-style probabilistic tracking, ambiguity-filtered MSD estimates,
%        and a hard-gate-conditioned Brownian step-likelihood estimate.
%     8. Optional experimental finite beta-Bernoulli joint trajectory sampler
%        evaluated directly on raw detector counts. A fixed or calibrated
%        shared monomer brightness constrains the emitter-count posterior.
%
%   All detector coordinates are sample-equivalent micrometers. The excitation
%   focus is fixed at the origin; this is not a raster-scanned ISM image.
%
%   The image model is stroboscopic: all photons in a frame are emitted from
%   that frame's stored molecule position. Continuous-exposure motion blur is
%   therefore not included. Use a frame interval short enough that
%   sqrt(4*D*dtS) is small compared with psfSigmaUm and laserWaistUm when
%   interpreting the CRB as a static-emitter localization bound.
%
%   IMPORTANT DEFAULTS
%     diffusionUm2PerS           1.00 um^2/s
%     nMolecules                 40
%     moleculeConcentration...   [] (set to override nMolecules)
%     boxSizeUm                  [2 2]
%     boundaryCondition          periodic
%     dtS                        0.001 s
%     nFrames                    1000
%     detectorLayout             honeycomb23
%     detectorPitchUm            0.18 um, sample-equivalent
%     laserWaistUm               0.30 um, excitation 1/e^2 radius
%     psfSigmaUm                 0.14 um, detection intensity sigma
%     laserPowerW                1e-5 W
%     excitationRatePerW         5e11 s^-1/W at the beam center
%     illuminationMode           gaussian
%     quantumYield               0.7
%     collectionEfficiency       0.8
%     detectorQuantumEfficiency  0.6
%     blinkOffRateS              0 (constant-brightness default)
%     localizationMethod         jointISMPoissonBayes
%     trackingMethod             jpda
%     gateCensoringCorrection    true
%     nTrackingRefinement...     3
%
%   OUTPUT
%     out.trajectories.positionsUm          [nMolecules x 2 x nFrames]
%     out.frames                            [nDetector x nFrames] microimages
%     out.photonEvents                      photon table (empty if storage disabled)
%     out.localizationTable                 fitted positions and CRB per spot
%     out.trackTable                        linked localization table
%     out.diffusion                         MSD, Fisher, and diffusion results
%       .nIndependentStepsForFisher         non-overlapping steps used by MLE
%     out.densityDiagnostics                overlap/association diagnostics
%     out.experimentalBnp                   raw-count joint sampler result
%       .frameEmitterCountPosterior         per-frame count posterior summary
%       .posteriorMapActiveEmitterCount     MAP total active-emitter count
%       .countDiagnosticsPassed             count-specific internal checks
%     out.figureSummary                     diagnostic figure handle
%     out.bnpCountFigure                    count-posterior figure handle
%
%   Implementation files: +membrane_tracking/+focused_ism
%   Concentration sweep: verifyBrightnessConstrainedBnpCounting

    out = membrane_tracking.focused_ism.runSimulation(varargin{:});
end
