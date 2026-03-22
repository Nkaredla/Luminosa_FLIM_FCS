function [auto, autotime] = tttr2xfcs(y, num, Ncasc, Nsub)

% Function [auto, autotime] = tttr2xfcs(y, num, Ncasc, Nsub) calculates
% the correlation and crosscorrelation between photon streams with arrival 
% times y, where num is a matrix indicating which photon corresponds to
% which photon stream
%[autocorr, autotime] = tttr2xfcs(y, num, ncasc, nsub)

%y is a vector that contains the increasing arrival times of the
%detected photons, every number in the vector corresponds to one detected photon

%num is typically a n x m - matrix; here n is the number of detected photons and
%m is the number of different detection channels. num(j,k) is equal to
%1 if the jth photon was detected in the kth channel, otherwise the
%value is zero. In the case of FLCS, num(j,k) should be the weight of the jth photon in
%the kth detection channel, the weight as calculated by the FLCS-TCSPC filter
%as described in our paper.

%ncasc and nsub define the time values at which the autocorrelation
%function is calculated: ncasc logarithmically placed cascades of nsub
%linearly placed time values (good numbers are ncasc=20 and nsub=10)

%auto is then automatically a (Ncasc*Nsub) x m x m Array, whereby
%auto(:,j,k) gives the correlation function of the jth channel against
%the kth channel 

dt = max(y)-min(y);
y = round(y(:));
if size(num,1)<size(num,2)
    num = num';
end
autotime = zeros(Ncasc*Nsub,1);
auto = zeros(Ncasc*Nsub,size(num,2),size(num,2));
shift = 0;
delta = 1;
for j=1:Ncasc
    [y, k] = unique(y);
	tmp = cumsum(num);
    num = diff([zeros(1,size(num,2)); tmp(k,:)]);
    for k=1:Nsub
        shift = shift + delta;
        lag = round(shift/delta);
        [tmp,i1,i2] = intersect(y,y+lag);
        auto(k+(j-1)*Nsub,:,:) = num(i1,:)'*num(i2,:)/delta;
        autotime(k+(j-1)*Nsub) = shift;
    end
    y = round(0.5*y);
    delta = 2*delta;
end
for j=1:size(auto,1)
    auto(j,:,:) = auto(j,:,:)*dt./(dt-autotime(j));
end

