function data_out = WillDataLoader(cfg_in,sno)
% function data_out = WillDataLoader(cfg,sno)
%
% data inventory and loader

cfg_def.dt = 1/60;
cfg_def.columns = 'FGHI';
cfg_def.fd = 'D:\My_Documents\Dropbox\projects\HDfit\data';

cfg = ProcessConfig(cfg_def,cfg_in);

p = pwd;
cd(cfg.fd);

switch sno
    case 1 % WB84 4-15
        sess = {'std','dark','laser'};
        fn = {'WB84 4-15 s1 ST7 c1 ST8 c1 std.xls','WB84 4-15 s2 ST7 c1 ST8 c1 dark.xlsx','WB84 4-15 s3 ST7 c1 ST8 c1 darklaseron.xlsx'};
        sn = {'WB84 4-15 s1 ST7 c1 ST8 c1 std.','Sheet1','Sheet1'};
        nCells = 2;
        
    case 2 % WB84 3-23
        
        
    case 3 % WB84 3-26
        
        
    case 4 % WB89 7-1
        
    case 5 % WB89 7-2
        
    case 6 % WB89 7-7
        
        
    case 7 % WB95 8-28
        
    case 8 % WB95 9-1
        sess = {'dark','laser'};
        fn = {'WB95 9-1 s2 ST1c1 ST8c1 dark.xls','WB95 9-1 s3 ST1c1 ST8c1 darklaseron.xls'};
        sn = {'WB95 9-1 s2 ST1c1 ST8c1 dark.tx','WB95 9-1 s3 ST1c1 ST8c1 darklas'};
        nCells = 2;
    case 9 % WB95 9-7
        
end

for iF = 1:length(fn)
   
    data_out.(sess{iF}).fn = fn{iF};
    
    bin_idx = xlsread(fn{iF},sn{iF},'A1:A65536')';
    data_out.(sess{iF}).bin_idx = bin_idx;
     
    for iC = 1:nCells
    
        this_column = cfg.columns(iC);
        this_range = cat(2,this_column,'1:',this_column,'65536');
        data_out.(sess{iF}).obs_fr(iC,:) = xlsread(fn{iF},sn{iF},this_range);
        
    end
    
    this_column = cfg.columns(nCells+1);
    this_range = cat(2,this_column,'1:',this_column,'65536');
    obs_hd = xlsread(fn{iF},sn{iF},this_range);

    orig_tvec = (bin_idx-1)*cfg.dt;
    data_out.(sess{iF}).obs_hd = tsd(orig_tvec,obs_hd);

end

cd(pwd);