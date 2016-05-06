function err = HDerrfun(in)
%
% in(1): hd0
% in(2): gain_left
% in(3): gain_right

global ahv tc sdf;

cfg_hd = []; cfg_hd.hd0 = in(1); cfg_hd.gain = [in(2) in(3)];
hd = AHVtoHD(cfg_hd,ahv);

predicted_sdf = GenerateSDFfromTC([],hd,tc);

err = (predicted_sdf.data-sdf.data).^2;
err = nansum(err(:));


