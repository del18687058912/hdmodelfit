%%
clear all; clear global; pack
addpath('D:\My_Documents\GitHub\hdmodelfit\shared');
%addpath('C:\Users\mvdm\Documents\GitHub\hdmodelfit\shared'); % isidro

cd('D:\My_Documents\Dropbox\projects\HDfit\data');

%% test
global param; % to access internals of crossval();
global param_count;
global param_hist;
param_count = 0;

%out = hdfit_crossval_func([]);

%% multirun

clear ALL_out;

sessions_to_run = 1:12;
targets_to_run = {'laser','std'};
%filters_to_run = {'smooth','kalman','kalmanwrapped'};
filters_to_run = {'smooth','kalman','kalmanwrapped'};

for iS = length(sessions_to_run):-1:1
    
    fprintf('\n\n*** SESSION %d/%d ***\n\n',iS,length(sessions_to_run));
    
    for iT = 1:length(targets_to_run)
        
        for iF = 1:length(filters_to_run)
            
            fprintf('\n\n--> %s\n\n',filters_to_run{iF});

            this_cfg = [];
            this_cfg.session = sessions_to_run(iS);
            this_cfg.target_session = targets_to_run{iT};
            this_cfg.mode = filters_to_run{iF};
            
            this_out = hdfit_crossval_func(this_cfg);
            
            ALL_out(iS).(targets_to_run{iT}).(filters_to_run{iF}) = this_out;
        end
    end
    
    save ALL_out;
end

   