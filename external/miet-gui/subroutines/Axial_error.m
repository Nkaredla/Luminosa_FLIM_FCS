function [z_center, dz] = Axial_error(z,alpha,lifecurve,tau,theta)

theta_mean = mean(theta); del_theta = std(theta);
tau_mean = mean(tau); del_tau = std(tau);

[~, theta_idx] = fuzzylocate(theta_mean,alpha);
[~,tau_idx] = fuzzylocate(tau_mean,lifecurve(:,theta_idx));
z_center = z(tau_idx);

tmp = gradient(lifecurve(:,theta_idx)); % dtau/dh
dh_dtau = 1/tmp(tau_idx);

tmp = gradient(lifecurve(tau_idx,:)); % dtau/dtheta;
dh_dtheta = tmp(theta_idx).*dh_dtau;

dz = sqrt(dh_dtau^2*del_tau^2+dh_dtheta^2*del_theta^2);

dtau_dtheta = (lifecurve(tau_idx,theta_idx+1)-lifecurve(tau_idx,theta_idx-1))./2/(alpha(theta_idx+1)-alpha(theta_idx-1));
dtau_dh = (lifecurve(tau_idx+1,theta_idx)-lifecurve(tau_idx-1,theta_idx))./2/(z(tau_idx+1)-z(tau_idx-1));

dz = (1./dtau_dh)*sqrt(del_tau^2+dtau_dtheta^2*del_theta^2);