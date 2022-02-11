% script to analyse results from the transport model and compare them to
% the traditional hydrograph separation by 2 end member mixing.

clear variables
addpath(genpath(fullfile('..','Source')))

%% SETTINGS

% select a simulation folder
mainfolder='saved_results'; %can be temp_results or saved_results
foldertoopen = '20220122230538';

% some flags here
use_saved_query=1; %use a saved query if it exists
plot_single=0; %plot a single simulation?
plot_multiple=1; %plot a results for all simulation?

% additional settings
event_selection=[3]; %select which events to display in the single plot


%% query the results
% open the file with the simulation summary and generate a table
load(fullfile(mainfolder,foldertoopen,'parameter_summary'))
T=array2table(M,'VariableNames',['ID',fieldvar']);
N=size(T,1);  
q = true(N,1);

% prepare the query
if use_saved_query == 1 && exist(fullfile(mainfolder,foldertoopen,'q_types.m'),'file')
    run(fullfile(mainfolder,foldertoopen,'q_types.m'))
else
    %build a query here
    q = q & T.z_Qmin==0.5; %condition on transport
    q = q & T.z_Qmax==0.8; %condition on transport
end

%% make plots
if plot_single==1
    % make an illustrative plot for a single simulation
    if ~exist('simul','var'), simul=1; end
    load(sprintf('%s%s/output_%05d',mainfolder,foldertoopen,simul))
    flag_plot.Cout=0;         %show plot with residuals
    flag_plot.TTDs=0;         %show plot with the selection of TTDs
    flag_plot.agemore=0;      %show plot with additional age timeseries
    flag_plot.events=1;       %show plot with individual events
    flag_plot.reactions=0;    %show the kinetics with selected TTDs
    run('plot_results_isoweathering')
end
if plot_multiple==1
    run('plot_results_multiple')
end

% eof   