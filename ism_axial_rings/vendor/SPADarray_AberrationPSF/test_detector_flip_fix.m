function test_detector_flip_fix()
% Test script to verify the detector image flip fix
% Compares the simple and explicit calculation methods for no aberration case

% Set up simulation parameters
sim = defaultParams();

% Test with no aberrations (empty coefficient structure)
coeffs_none = struct();

fprintf('Testing detector image consistency between methods...\n');

% Generate stacks using both methods
fprintf('Generating stack with normalizedStack (simple method)...\n');
stack_simple = normalizedStack(sim, coeffs_none);

fprintf('Generating stack with normalizedStackExplicitDetector (explicit method)...\n'); 
stack_explicit = normalizedStackExplicitDetector(sim, coeffs_none);

% Compare the two methods
fprintf('Comparing results...\n');

% Calculate correlation between corresponding channels
correlations = zeros(size(sim.detXY,1), 1);
for k = 1:size(sim.detXY,1)
    img_simple = stack_simple(:,:,k);
    img_explicit = stack_explicit(:,:,k);
    
    % Normalize for correlation calculation
    img_simple_norm = (img_simple - mean(img_simple(:))) / std(img_simple(:));
    img_explicit_norm = (img_explicit - mean(img_explicit(:))) / std(img_explicit(:));
    
    correlations(k) = corr2(img_simple_norm, img_explicit_norm);
end

% Display results
fprintf('\nCorrelation coefficients between simple and explicit methods:\n');
fprintf('Channel\tCorrelation\n');
fprintf('-------\t-----------\n');
for k = 1:length(correlations)
    fprintf('%d\t\t%.4f\n', k, correlations(k));
end

fprintf('\nMean correlation: %.4f\n', mean(correlations));
fprintf('Min correlation:  %.4f\n', min(correlations));

% Test detector shift sign consistency
fprintf('\nTesting detector shift sign consistency...\n');
checkDetectorShiftSign(sim);

% Additional check: compare center channel intensities
center_idx = ceil(size(stack_simple,3)/2);
center_simple = stack_simple(:,:,center_idx);
center_explicit = stack_explicit(:,:,center_idx);

fprintf('Center channel intensity sum - Simple: %.6f, Explicit: %.6f\n', ...
    sum(center_simple(:)), sum(center_explicit(:)));

% Check if the fix resolved the flip issue
if mean(correlations) > 0.95
    fprintf('\n✓ SUCCESS: High correlation between methods indicates flip issue is resolved!\n');
else
    fprintf('\n✗ WARNING: Low correlation suggests remaining inconsistencies.\n');
end

end

function sim = defaultParams()
sim.modeOrder = {'tilt_x','tilt_y','defocus','astig_x','astig_y','coma_x','coma_y','spherical'};
sim.lamExc = 0.488;
sim.lamEm = 0.520;
sim.lamRef = 0.520;
sim.NA = 1.2;
sim.nMedium = 1.33;
sim.fovXY = 1.8;
sim.nzRange = 0.8;
sim.nx = 27;
sim.nz = 5;
sim.beadRadius = 0.08;
sim.detPitch = 0.18;
sim.detSize = 0.10;
sim.arrayN = 5;
sim.Nr = 36;
sim.Nphi = 72;
sim.M = 5;
sim.x = linspace(-sim.fovXY/2, sim.fovXY/2, sim.nx);
sim.y = sim.x;
sim.z = linspace(-sim.nzRange/2, sim.nzRange/2, sim.nz);
sim.dx = abs(sim.x(2)-sim.x(1));
coords = ((0:sim.arrayN-1) - (sim.arrayN-1)/2) * sim.detPitch;
[gx, gy] = meshgrid(coords, coords);
sim.detXY = [gx(:) gy(:)];
sim.obj = beadObject3D(sim);
end

function obj = beadObject3D(sim)
[X, Y] = meshgrid(sim.x, sim.y);
obj = zeros(numel(sim.y), numel(sim.x), numel(sim.z));
for iz = 1:numel(sim.z)
    obj(:,:,iz) = double(X.^2 + Y.^2 + sim.z(iz).^2 <= sim.beadRadius^2);
end
obj = obj / sum(obj(:));
end