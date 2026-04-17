% Reading the TCSPC histogram
 
 [bin,tcspcdata,head] = Harp_tcspc(name); % name represents the raw data file (.ht3 or .ptu)
 % plotting the TCSPC curve
 
 resolution = head.Resolution; % the resolution of each time channel
 tau        = (1:numel(bin))*resolution; % tcspc time bins
 figure
 title('TCSPC decay')
 semilogy(tau,tcspcdata);
 xlabel('\tau (ns)')
 ylabel('Counts')
 
 % --- TCSPC fitting ---
 
 % What we perform below is a tail fitting. 
 [~, pos] = max(tcspcdata); % time bin with maximum counts
 cutoff   = 1; % 1 ns cutoff
 shift = round(cutoff./resolution); % time bin from which we start fitting
 tcspc = sum(tcspcdata(shift+max(pos):end,:),2); % 
 t   = [0:numel(tcspc)-1]*resolution;
 p = [2]; % initial guesses in ns. One can give multiple values for multiexponential fitting.
 
 figure
 p = Simplex('ExpFun',p,[],[],[],[],t,tcspc,1,[],1); % Fitting exponential functions
 
 disp(['The fitted values are  ' num2str(p.') ' ns'])
 [~, c, ~, z] = ExpFun(p, t,tcspc,1);
 mean_tau = sum(c(2:end).*p)/sum(c(2:end));
 clf
 plot(t,tcspc,'o',t,z,'r','LineWidth',2)
 axis tight
 text('Position',[numel(tcspc)*resolution*0.7,z(1)*0.9],'String',['\tau = ',mnum2str(mean_tau,[],2),' ns'],'FontSize',16)
 xlabel('\tau (ns)','FontSize',16)
 ylabel('Counts','FontSize',16)
 set(gca,'FontSize',14)
 
 % Average arrival time of all photons
 
 counts_end = mean(tcspc(end-11:end-1)); % Offset
 av_arriv_time = sum(tcspc.*t.'-counts_end*mean(t))./(sum(tcspc)-counts_end*mean(t)); % Average photon arrival time
 
 
