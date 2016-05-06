function hd = AHVtoHD(cfg_in,ahv)
% function hd = AHVtoHD(cfg_in,ahv)
%
% ahv: input AHV tsd in deg/s

cfg_def = [];
cfg_def.gain = [1 1]; % L, R gain
cfg_def.noisetype = 'none'; % 'none','gauss'
cfg_def.noisesd = 1;
cfg_def.hd0 = 0; % starting HD (deg)

cfg = ProcessConfig(cfg_def,cfg_in);

%
if ~CheckTSD(ahv)
   error('Input not a well-formed TSD.'); 
end

dt = diff(cat(2,0,ahv.tvec));
dahv = zeros(size(ahv.data));

left_idx = ahv.data < 0;
right_idx = ahv.data >= 0;

dahv(left_idx) = dt(left_idx).*ahv.data(left_idx).*cfg.gain(1);
dahv(right_idx) = dt(right_idx).*ahv.data(right_idx).*cfg.gain(2);

hd_unwrapped = cfg.hd0 + cumsum(dahv);

hd_wrapped = nan(size(hd_unwrapped));
plus_idx = hd_unwrapped >= 0;
hd_wrapped(plus_idx) = rem(hd_unwrapped(plus_idx),360);
minus_idx = hd_unwrapped < 0;
hd_wrapped(minus_idx) = rem(hd_unwrapped(minus_idx),-360)+360;

hd = tsd(ahv.tvec,hd_wrapped);

