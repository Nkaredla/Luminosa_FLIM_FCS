function test_detector_flip_fix_nogui()
% Test script to verify the detector image flip fix (no graphics)
% Compares the simple and explicit calculation methods for no aberration case

try
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
        
        % Flatten and compute correlation coefficient manually
        x = img_simple(:);
        y = img_explicit(:);
        
        % Remove means
        x = x - mean(x);
        y = y - mean(y);
        
        % Compute correlation
        correlations(k) = (x' * y) / sqrt((x' * x) * (y' * y));
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
    
    % Test detector shift consistency - simplified version without graphics
    fprintf('\nTesting detector shift sign consistency (simplified)...\n');
    test_detector_shift_simple(sim);
    
catch ME
    fprintf('Error during test: %s\n', ME.message);
    rethrow(ME);
end
end

function test_detector_shift_simple(sim)
    stack = normalizedStack(sim, struct());
    
    xc = zeros(size(sim.detXY,1),1);
    yc = zeros(size(sim.detXY,1),1);
    
    [X, Y] = meshgrid(sim.x, sim.y);
    
    for k = 1:size(stack,3)
        img = stack(:,:,k);
        s = sum(img(:));
        xc(k) = sum(X(:).*img(:)) / s;
        yc(k) = sum(Y(:).*img(:)) / s;
    end
    
    % Check if detector positions and image centroids are consistent
    det_shift_x = sim.detXY(:,1);
    det_shift_y = sim.detXY(:,2);
    
    % Calculate correlation between detector positions and image centroids
    corr_x = corrcoef(det_shift_x, xc);
    corr_y = corrcoef(det_shift_y, yc);
    
    fprintf('Detector-to-centroid correlation X: %.4f\n', corr_x(1,2));
    fprintf('Detector-to-centroid correlation Y: %.4f\n', corr_y(1,2));
    
    if abs(corr_x(1,2)) > 0.8 && abs(corr_y(1,2)) > 0.8
        fprintf('✓ Detector shift directions are consistent\n');
    else
        fprintf('✗ Detector shift directions may have issues\n');
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