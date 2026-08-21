function sim = defaultParams()
%--------------------------------------------------------------------------
% defaultParams
%
% PURPOSE
%   Construct a structure containing all default simulation settings.
%
% FIELDS
%   modeOrder   : list of aberration modes used in the fit / simulation
%   lamExc      : excitation wavelength [um]
%   lamEm       : emission wavelength [um]
%   lamRef      : reference wavelength used to scale OPD coefficients [um]
%   NA          : objective numerical aperture
%   objectiveMagnification : nominal objective magnification
%   nMedium     : refractive index of immersion/sample medium
%   fovXY       : lateral field of view [um]
%   nzRange     : axial extent of simulation [um]
%   nx          : number of x,y pixels
%   nz          : number of z planes
%   beadRadius  : fluorescent bead radius [um]
%   detPitch    : detector pitch in sample units [um]
%   detSize     : detector pixel size projected to sample units [um]
%   detectorLayout : detector-layout flag ('honeycomb23' default, or 'square25')
%   nDet        : number of detector channels
%   Nr          : number of radial pupil samples
%   Nphi        : number of azimuthal pupil samples
%   M           : highest azimuthal harmonic retained
%   x,y,z       : coordinate grids
%   dx          : lateral sampling pitch [um]
%   detXY       : detector center coordinates [nDet x 2]
%   obj         : 3D bead object
%--------------------------------------------------------------------------

    % Ordered list of aberration modes
    sim.modeOrder = zernikeModeNames(6);

    % Wavelengths in microns
    sim.lamExc = 0.488;
    sim.lamEm = 0.520;
    sim.lamRef = 0.520;

    % Objective / medium parameters
    sim.NA = 1.2;
    sim.nMedium = 1.33;
    sim.sampleGeometry = 'homogeneous';
    sim.nImmersion = 1.33;
    sim.nGlass = 1.518;
    sim.nSample = 1.0003;
    sim.nDesignGlass = 1.518;
    sim.coverslipThicknessUm = 190;
    sim.designCoverslipThicknessUm = 190;
    sim.beadBottomHeightUm = 0;
    sim.airBeadAxialSamples = 3;
    sim.objectiveMagnification = 60;
    sim.objectiveDescription = '60x 1.2 NA water';
    sim.diffractionModel = 'vectorial Richards-Wolf/Bessel';
    sim.includesVectorialPolarization = true;
    sim.vectorialPolarizationMode = 'xyAverage';     % legacy/global fallback
    sim.excitationPolarizationMode = 'x';            % linearly polarized scan focus
    sim.collectionPolarizationMode = 'xyAverage';    % isotropic bead detection
    sim.modelAccuracyNote = ['Default forward model includes vectorial ' ...
        'Richards-Wolf polarization terms. Air-interface emission remains an ' ...
        'isotropic-orientation effective model unless a dipole orientation is supplied.'];
    sim.interfaceRadialWeightMode = 'sampleSolidAngle';
    sim.airInterfaceStageMedium = 'immersion';   % piezo moves coverslip through immersion

    % Spatial simulation size
    sim.fovXY = 1.8;     % lateral field of view [um]
    sim.nzRange = 3.5;   % axial range [um]

    % Grid resolution
    sim.nx = 27;         % number of pixels in x and y
    sim.nz = 5;          % number of z slices

    % Bead size
    sim.beadRadius = 0.08; % [um]
    sim.beadSubsamples = [3 3 3];

    % Detector geometry
    sim.detPitch = 0.18; % detector pitch [um]
    sim.detSize = sim.detPitch;  % detector flat-to-flat width [um], touching honeycomb cells
    sim.detectorLayout = 'honeycomb23'; % detector-layout flag; use 'square25' for legacy 5x5
    sim.detectorPixelShape = 'hex';
    sim.detectorHexRadius = sim.detSize / sqrt(3);

    % Pupil discretization
    sim.Nr = 36;         % radial samples in pupil
    sim.Nphi = 72;       % azimuthal samples in pupil
    sim.M = 6;           % azimuthal harmonics kept: m = -M ... +M

    % Lateral and axial coordinate grids
    sim.x = linspace(-sim.fovXY/2, sim.fovXY/2, sim.nx);
    sim.y = sim.x;
    sim.z = linspace(-sim.nzRange/2, sim.nzRange/2, sim.nz);

    % Lateral pixel size
    sim.dx = abs(sim.x(2)-sim.x(1));

    % Detector center coordinates. The default is the 23-channel Luminosa
    % honeycomb layout; callers can replace sim.detXY with data-derived
    % detector positions before running the forward model.
    [sim.detXY, sim.detectorIndexGrid, sim.detectorLayoutInfo] = ...
        detectorLayout(sim.detectorLayout, sim.detPitch);
    sim.nDet = size(sim.detXY, 1);
    sim.detectorGridSize = size(sim.detectorIndexGrid);
    sim.arrayN = sim.detectorGridSize(1); % legacy display row count

    sim.detectorSubsamples = 7;
    sim.detectorImageInverted = true;

    % Build the 3D bead object
    sim.obj = beadObject3D(sim);
end
