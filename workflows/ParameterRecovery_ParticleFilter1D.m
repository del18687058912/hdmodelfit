%% init
addpath('D:\My_Documents\GitHub\hdmodelfit\shared');

%% generate the behavioral data (AHV profile)
cfg = []; cfg.t = [0 100];
ahv = GenerateAHVTrajectory(cfg); % (simulated) experimentally observed AHV

cfg_hd = []; 
cfg_hd.hd0 = 90; cfg_hd.gain = [1 1]; % these are the parameters we want to recover, used here to generate frates
hd = AHVtoHD(cfg_hd,ahv); % simulated **internal** HD

plot(hd,'.');

%% make some tuning curves of cells
cfg_tc = [];
cfg_tc.pfd = [100 280]; % preferred firing direction of cells (length of vector --> nCells)
cfg_tc.maxfr = [45 35]; % peak firing rate
cfg_tc.sd = [30 40]; % standard deviation (in deg)
tc = GenerateHDTuningCurves(cfg_tc);

plot(tc.xbin,tc.tc);

%% simulate mean firing rates
sdf = GenerateSDFfromTC([],hd,tc);
plot(sdf);

%% now recover parameters -- particle filtering approach -- HD0 only
global ahv tc sdf;
cfg_pf.in = [180 1 1]; % initial conditions
cfg_pf.lb = [0 0.1 0.1];
cfg_pf.ub = [359 5 5];
cfg_pf.minDelta = [1 0.01 0.01]; % don't resample finer than this
cfg_pf.nP = 1000; % number of particles
cfg_pf.nT = 100; % number of time steps
cfg_pf.nCells = 2; % number of cells
cfg_pf.dt = 1; % not currently used, but will be important when running on real spikes

% initialize "live" particle matrix
particle_vals = nan(cfg_pf.nT+1,cfg_pf.nP);
particle_vals(1,:) = linspace(cfg_pf.lb(1),cfg_pf.ub(1),cfg_pf.nP); % initial distro of param values

% also some other matrices for housekeeping
particle_hd0 = nan(cfg_pf.nT,cfg_pf.nP);
particle_p = nan(cfg_pf.nT,cfg_pf.nP);

%% run filter

for iT = 1:cfg_pf.nT
    
    % AHV up to this time point
    this_AHV = ahv;
    this_AHV.tvec = this_AHV.tvec(1:iT);
    this_AHV.data = this_AHV.data(1:iT);
    
    clear HD_est;
    for iP = cfg_pf.nP:-1:1
        
        % for each particle, obtain HD_est from known AHV
        cfg_hd = []; cfg_hd.hd0 = particle_vals(iT,iP);
        this_hd = AHVtoHD(cfg_hd,this_AHV);
        HD_est(iP) = this_hd.data(iT);

        % for each particle, compute **expected** FR given HD_est
        for iC = 1:cfg_pf.nCells
            exp_fr(iP,iC) = cfg_pf.dt*interp1(tc.xbin,tc.tc(iC,:),HD_est(iP),'linear','extrap');
        end
        
        % get p(spikes|HD_est) by comparing to **observed** FR
        obs_fr = round(sdf.data(:,iT));
        % note, rounding is to simulate that we will be getting integer
        % spike counts when run on real data
        
        for iC = 1:cfg_pf.nCells
            % assume Poisson distributed spike counts
            this_p(iP,iC) = (((cfg_pf.dt * exp_fr(iP,iC)).^obs_fr(iC))./factorial(obs_fr(iC))).*(exp(-cfg_pf.dt.*exp_fr(iP,iC)));
        end
        
    end % of particle loop
    this_p = prod(this_p,2); % assume spike counts independent
    
   % plot the result
   subplot(221)
   plot(particle_vals(iT,:),this_p,'.k');
   title(sprintf('t %d',iT));
   
   % resample - simple method that doesn't generate new parameter values
   keep_idx = randsample(1:length(this_p),cfg_pf.nP,true,this_p);
   particle_vals(iT+1,:) = particle_vals(iT,keep_idx); 
   
   % plot comparison
   subplot(222)
   h1 = hist(particle_vals(iT,:),particle_vals(1,:));
   h2 = hist(particle_vals(iT+1,:),particle_vals(1,:));
   h(1) = plot(particle_vals(1,:),h1,'k'); hold on;
   h(2) = plot(particle_vals(1,:),h2,'r');
   legend(h,{'previous','resampled'});
   hold off;
    
   pause;
   
end % of time steps
    
    
    



