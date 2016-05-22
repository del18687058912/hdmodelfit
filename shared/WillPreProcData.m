function data = WillPreProcData(cfg_in,data)
% function data = WillPreProcData(cfg_in,data)
%
% preprocess HD data

cfg_def = [];
cfg_def.dt = 1/60;
cfg_def.smoothwin = 11;
cfg_def.subsample_factor = 7;
cfg_def.debug = 0;

cfg = ProcessConfig(cfg_def,cfg_in);


fnames = fieldnames(data);
for iF = 1:length(fnames)
   
    this_data = data.(fnames{iF});
    
    %% need to unwrap HD to interpolate missing data
    [~,hd_unwrapped] = HDtoAHV(this_data.obs_hd);
    
    tvec_interp = (0:max(this_data.bin_idx))*cfg.dt;
    hd_unwrapped_interp = interp1(hd_unwrapped.tvec,hd_unwrapped.data,tvec_interp,'linear');
    hd_unwrapped_interp = medfilt1(hd_unwrapped_interp,cfg.smoothwin);
    hd_unwrapped_interp = smooth(hd_unwrapped_interp,cfg.smoothwin); % why does this take so long?
    
    hd_unwrapped_interp = tsd(tvec_interp,hd_unwrapped_interp');
    
    if cfg.debug
        %% could consider kalman filter instead of smoothing + linear interpolation... but this looks ok for now
        figure;
        plot(hd_unwrapped_interp,'.r'); hold on;
        plot(hd_unwrapped,'.k');
    end
    
    %% put it back in range
    hd_wrapped_interp = hd_unwrapped_interp;
    hd = wrapHD(hd_wrapped_interp.data);
    hd = tsd(hd_unwrapped_interp.tvec,hd);
    
    data.(fnames{iF}).hd = hd;
 
    %% subsample
    dt_ss = cfg.dt*cfg.subsample_factor;
    hd_ss = hd;
    hd_ss.tvec = hd_ss.tvec(1:cfg.subsample_factor:end);
    hd_ss.data = hd_ss.data(1:cfg.subsample_factor:end);
    
    if cfg.debug
        figure;
        plot(hd_ss,'.r'); hold on;
        plot(this_data.obs_hd,'.k');
    end
    
    %% get AHV
    ahv_ss = HDtoAHV(hd_ss);
    
    data.(fnames{iF}).ahv_ss = ahv_ss;
    data.(fnames{iF}).hd_ss = hd_ss;
    
end