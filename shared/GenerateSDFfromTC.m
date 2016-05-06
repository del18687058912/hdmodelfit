function sdf_out = GenerateSDFfromTC(cfg_in,hd,tc)
% function sdf_out = GenerateSDFfromTC(cfg_in,hd,tc)
%

cfg_def = [];

cfg = ProcessConfig(cfg_def,cfg_in);

nCells = size(tc.tc,1);

for iC = nCells:-1:1
    
    sdf_out(iC,:) = interp1(tc.xbin,tc.tc(iC,:),hd.data,'linear');
    
end

sdf_out = tsd(hd.tvec,sdf_out);