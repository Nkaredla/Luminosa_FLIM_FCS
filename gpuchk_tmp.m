cd('D:\MATLAB\server\Luminosa\Luminosa_FLIM_FCS');
fprintf('MATLAB %s\n', version);
fprintf('Parallel Computing Toolbox: %d\n', license('test','Distrib_Computing_Toolbox'));
try
  n = gpuDeviceCount;
  fprintf('gpuDeviceCount = %d\n', n);
  if n > 0
    g = gpuDevice(1);
    fprintf('  name %s\n  compute %s\n  total %.1f GB, available %.1f GB\n', ...
      g.Name, g.ComputeCapability, g.TotalMemory/1e9, g.AvailableMemory/1e9);
    fprintf('  double support: %d\n', g.SupportsDouble);
  end
catch e
  fprintf('no GPU: %s\n', e.message);
end
fprintf('maxNumCompThreads = %d\n', maxNumCompThreads);
disp('GPU_OK');
