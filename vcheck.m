cd('D:\MATLAB\server\Luminosa\Luminosa_FLIM_FCS');
d = dir('*.m'); bad = 0;
for k = 1:numel(d)
  try
    nargin(erase(d(k).name,'.m'));
  catch e
    if contains(e.identifier,'m_unterminated') || contains(e.message,'Error:')
      fprintf('CANNOT PARSE %s : %s\n', d(k).name, e.message); bad = bad+1;
    end
  end
end
fprintf('unparseable: %d of %d\n', bad, numel(d));
fs = {'biexp_slb_basis.m','biexp_slb_pattern_batch.m','immune_cell_MIET_biexp_slb.m', ...
      'immune_cell_MIET_biexp_vp_run.m','benchmark_biexp_speed.m'};
n=0; for k=1:numel(fs)
  m=checkcode(fs{k},'-struct');
  for j=1:numel(m); fprintf('%s L%d: %s\n',fs{k},m(j).line,m(j).message); n=n+1; end
end
fprintf('LINT=%d\n',n); disp('VC_OK');
