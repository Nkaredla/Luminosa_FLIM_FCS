function MIET_plotter(dirname)

fnames = dir([dirname '*.ht3']);
pointidx = zeros(numel(fnames,1));
for i = 1:numel(fnames)
    pointidx(i)=~isempty(strfind(fnames(i).name,'Image'));
end
fnames = fnames(logical(pointidx));

theta = []; ax_pos = []; tau = []; photons = []; tau2 = [];
for j = 27:32%numel(fnames)
    name = [dirname fnames(j)name];
%     if exist([name(1:end-4),'_miet_analysis.mat'],'file')
%         load(([name(1:end-4),'_miet_analysis.mat']));
%     else
%         return
%     end
    if exist([name(1:end-4),'_sm_flim.mat'],'file')
        load(([name(1:end-4),'_sm_flim.mat']));
    else
        return
    end
   
    theta = [theta; SM.theta(:)];
    lifetimes = [];
    
    for num = 1:numel(FLIM.life_mat)
        
        lifetimes(num) = FLIM.life_mat{num}(end);
    end
%         tau = [tau; FLIM.life_matav(:)];
    tau = [tau; lifetimes(:)];
%     tau2 = [tau2; FLIM.life_matav2(:)];
%     tmppos = [];
%     for i = 1:numel(SM.theta)
%         tmppos(i) = max(max(posimm(:,:,i)));
%     end
%     ax_pos = [ax_pos; tmppos(:)];
    photons = [photons; FLIM.photon_sm(:)];
end


colors = distinguishable_colors(size(lifecurve,2),'w');
% colors = hsv(size(lifecurve,2));
% al_res = pi/2/(size(lifecurve,2));
al_meas = unique(theta);

[~,theta_pos] = fuzzylocate(theta,al_meas);

for i = 1:size(lifecurve,2)
    plot(z,lifecurve(:,i),'Color',colors(i,:))
    ind = theta_pos==i;
    plot(ax_pos(ind),tau(ind),'o','Color',colors(i,:))
    hold on
end
hold off

zmax = max(ax_pos);
zmin = min(ax_pos);
% zmin = 4;
% zmax = 16;

mean_orient = []; cnt = []; std_orient = []; 
height = zmin:(zmax-zmin)/9:zmax;
% [r_axpos] = fuzzylocate(ax_pos,height);

[~,r_axpos] = histc(ax_pos,height);

for i = 1:max(r_axpos)
    mean_orient(i) = mean(theta(r_axpos==i));
    std_orient(i) = std(theta(r_axpos==i));
    cnt(i) =sum(r_axpos==i);    
end
figure
errorbar(unique(r_axpos),mean_orient*180/pi,std_orient*180/pi);

for i = 1:numel(al)
    mean_z(i) = pmean(ax_pos(theta==al(i)));
end

n = histc(ax_pos,0.5:8.5);
[AX,H1,H2] = plotyy(height(2:end-1),cnt(2:end-1),height(2:end-1),mean_orient(2:end-1)*180/pi)
hold(AX(2),'on')
errorbar(AX(2),height(2:end-1),mean_orient(2:end-1)*180/pi,std_orient(2:end-1)./sqrt(cnt(2:end-1))*180/pi,'bo');
hold off
xlim(AX(2),[0 zmax])
xlim([0 zmax])
ylim(AX(2),[50 80])
ylabel(AX(2),'\theta (degrees)')
ylabel(AX(1),'# molecules')
xlabel('height (nm)')



return

tau_mean = []; tau_std = []; cnt = [];
for i = 1:numel(al_meas)
    photon_mean(i) = mean(photons(theta_pos==i));
    photon_std(i)  = std(photons(theta_pos==i));
    tau_mean(i) = mean(tau(theta_pos==i));
    tau_std(i)  = std(tau(theta_pos==i));
    cnt(i) = sum(theta_pos==i);
end


