%% generate the behavioral data (AHV profile)
cfg = []; cfg.t = [0 100];
ahv = GenerateAHVTrajectory(cfg); % (simulated) experimentally observed AHV

cfg_hd = []; 
cfg_hd.hd0 = 90; cfg_hd.gain = [0.5 1]; % these are the parameters we want to recover, used here to generate spikes
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

%% now recover parameters -- optimization approach
global ahv tc sdf;
in = [180 1 1]; % initial conditions
lb = [0 0.1 0.1];
ub = [359 5 5];

opts = optimoptions('fmincon','Display','off');

hd0_init = 0:2:359; % try to avoid local minima by using many starting points

clear x_min fval;
for iI = length(hd0_init):-1:1
    
    fprintf('Starting point %d/%d\n',iI,length(hd0_init));
    in0 = [hd0_init(iI) in(2) in(3)];
    [x_min{iI},fval(iI),exitflag,output] = fmincon(@HDerrfun,in0,[],[],[],[],lb,ub,[],opts);

end
[~,min_idx] = min(fval);
fprintf('Minimum at hd0 = %.2f, gamma_l = %.2f, gamma_r = %.2f\n',x_min{min_idx}(1),x_min{min_idx}(2),x_min{min_idx}(3));

