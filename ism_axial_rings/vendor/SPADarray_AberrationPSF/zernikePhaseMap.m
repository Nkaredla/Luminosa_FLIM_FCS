function phase = zernikePhaseMap(sim, coeffs, wavelength)
%--------------------------------------------------------------------------
% zernikePhaseMap
%
% PURPOSE
%   Convert recovered Zernike coefficients into a pupil phase map.
%
% INPUTS
%   sim        : simulation structure from defaultParams()
%   coeffs     : aberration coefficient struct or numeric vector
%   wavelength : wavelength for phase in radians [um], default sim.lamRef
%
% OUTPUT
%   phase : structure with fields
%       rho, phi     : pupil coordinates
%       wavesRef     : wavefront in waves at sim.lamRef
%       radians      : phase in radians at wavelength
%       coeffs       : coefficient structure used for the map
%--------------------------------------------------------------------------

    if nargin < 3 || isempty(wavelength)
        wavelength = sim.lamRef;
    end

    coeffs = coeffStruct(sim, coeffs);

    rho = linspace(0, 1, sim.Nr);
    phi = linspace(0, 2*pi, sim.Nphi+1);
    phi(end) = [];
    [R, P] = ndgrid(rho, phi);

    Z = zernikeModes(R, P);
    W = zeros(size(R));
    fn = fieldnames(coeffs);

    for k = 1:numel(fn)
        W = W + coeffs.(fn{k}) * Z.(fn{k});
    end

    phase.rho = R;
    phase.phi = P;
    phase.wavesRef = W;
    phase.radians = 2*pi*(sim.lamRef/wavelength).*W;
    phase.coeffs = coeffs;
    phase.wavelength = wavelength;
end
