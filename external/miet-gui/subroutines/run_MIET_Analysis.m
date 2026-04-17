%clear all
close all

load metals

%fn
dye.lamex = 0.488;
dye.lamem    = 0.542;
dye.qy       = 0.78;
dye.tau_free = 2.9;
dye.CurveType = 'maximum';
%
layers.n0    = [1.52 titan(wavelength == dye.lamem.*1e3) gold(wavelength == dye.lamem.*1e3) 1.46];
layers.n0    = 1.52;
layers.n1    = 1.33;
layers.n     = 1.33;  
layers.d0    = [2 10 30].*1e-3;
% layers.d0    = [];
layers.d1    = [];
layers.d     = 0.02;

mic.NA       = 1.49;
% mic.pattern  = 'radial';
mic.pattern  = 'linear';
% mic.pattern  = 'azimuthal';
mic.focpos   = 0;

% dname = 'U:\Arindam\170808_graphenescan\results_10\';
dname = 'U:\Arindam\181122_Hairpinstalk.sptw\analysis';
% dname = 'U:\Narain\140114\';
% dname = 'U:\Narain\from Anna\28-01-2014\';
fnames = dir([dname '*.ptu']);
pointidx = zeros(numel(fnames,1));
for i = 1:numel(fnames)
    pointidx(i)=~isempty(strfind(fnames(i).name,'hairpin'));  
end
fnames = fnames(logical(pointidx));
% flag='MIET';
% flag = '3D_MIET';
flag = 'lifetime';


for j=1:numel(fnames)
    name = fnames(j).name
%     IRF='U:\Narain\140305\Point_011.ht3';
    IRF=[];
   MIET_Analysis(name, IRF, dye, layers, mic, flag,1,1); 
  % MIET_Analysis(fnames, IRF, dye, layers, mic, flag,1,1); 
end

for j = 1:numel(fnames)
     name = fnames(j).name;
     MIET_viewer([dname name])
end