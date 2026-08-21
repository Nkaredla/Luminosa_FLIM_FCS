function out = simulateHopTrapMembraneTracking(varargin)
%SIMULATEHOPTRAPMEMBRANETRACKING Infer free, hop, or trap membrane diffusion.
%
%   out = simulateHopTrapMembraneTracking()
%   out = simulateHopTrapMembraneTracking('diffusionMode','trap')
%   out = simulateHopTrapMembraneTracking( ...
%       'referencePreset','stedFcsNrkHop','nFrames',3000)
%
%   The simulation follows low-concentration molecules on a flat periodic
%   membrane. A molecule can diffuse freely, cross periodic Voronoi
%   compartment boundaries with a finite probability (hop diffusion), enter
%   a slowly diffusing state (trap diffusion), or experience both effects.
%   Focused excitation, finite-channel ISM detection, Poisson photons,
%   maximum-likelihood localization, Fisher information, and trajectory
%   linking reuse the focused-ISM implementation in +membrane_tracking.
%
%   Two evidence calculations are deliberately kept separate:
%     1. A step-likelihood BIC compares free diffusion with a two-state
%        diffusion HMM and is the primary test for transient trapping.
%     2. An MSD quasi-BIC compares a linear free model with a confined-plus-
%        long-range model and is the primary test for hop diffusion.
%
%   A finite sticky HDP-HMM approximation provides posterior mobility-state
%   counts and state diffusion coefficients. It detects persistent mobility
%   states; it does not by itself identify spatial membrane fences.
%
%   Named parameter presets are based on Schneider et al. ACS Nano 2018 and
%   related STED-FCS/STED-FLCS work. See the package README for references,
%   parameter mappings, model assumptions, and interpretation.
%
%   Important outputs:
%     out.trajectories      true positions, trap states, and compartments
%     out.localizationTable focused-ISM localizations and Fisher CRBs
%     out.trackTable        linked low-density trajectories
%     out.modelComparison   BIC evidence and final diffusion classification
%     out.stickyBnp         posterior mobility-state summary
%
%   Implementation files: +membrane_tracking/+hop_trap

    out = membrane_tracking.hop_trap.runSimulation(varargin{:});
end
