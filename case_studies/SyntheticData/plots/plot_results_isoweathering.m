% script to plot the results of a model single-run. This script is
% standalone but it can be automatically called at the end of the model
% starter if any plot_flag is set to 1

%--------------------------------------------------------------------------
% load the data
%--------------------------------------------------------------------------
% reload the results if not already done in the model starter
if ~exist('C_Qm','var')
    
    % edit here if you want to modify the selected output file or the plot flags
    case_study = fscanf(fopen(fullfile('..','..','case_study_name.txt')),'%s');
    load(fullfile('..','..','case_studies',case_study,'results','all_output')) %loading the default output file 'all_output'
    flag_plot.Cout = 1;    %main plot with the full concentration timeseries
    flag_plot.TTDs = 1;    %plot TTDs for selected days
    flag_plot.agestats = 1;      %plot with additional age timeseries
    flag_plot.events = 1;       %plot with timeseries of the individual events
    flag_plot.eventstats = 1;   %plot with other stats during the events (only if flag_plot.events == 1)
    flag_plot.deltavsC = 1;     %plot with deltavsC in batch reaction and simulations (only if flag_plot.events == 1)
end


% some useful stuff:
% - check if spinup is to be shown or not 
% data.show_spinup = 0; 
if data.show_spinup==0
    t_1=data.dates(1)+data.ini_shift*data.dt/24; %avoid showing the spinup
    C_Qm(1:data.ini_shift) = NaN; %discard results during spinup
    d_f(1:data.ini_shift) = NaN; %discard results during spinup
    Fyw(1:data.ini_shift,:) = NaN; %discard results during spinup
else
    t_1=data.dates(1);
end
t_2=data.dates(end);

%--------------------------------------------------------------------------
% CALL THE FIGURE SCRIPTS
%--------------------------------------------------------------------------

% FIGURE: SHOW OUTPUT C_Q
if flag_plot.Cout==1
    run('plot_Cout.m')
end

% FIGURE: SHOW THE SELECTED TTDs
if flag_plot.TTDs==1 && ~isempty(data.index_datesel)==1
    run('plot_TTDs.m') 
end
    
% select the Fyw you want to display
%ii_sel = [2,3,4]; %index of the ywt that will be shown
ii_sel = [1:5]; %index of the ywt that will be shown


% FIGURE: ADDITIONAL AGE STATISTICS
if flag_plot.agestats==1
    run('plot_agestats.m')  
end

% FIGURES WITH EVENTS

    % select one or more periods to show:
    events = data.events;
%     events={'22-Jan-2020','03-Feb-2020';...
%         '22-Jun-2020','01-Jul-2020';...
%         '03-Dec-2020','19-Dec-2020';...
%         '16-Jan-2021','30-Jan-2021';...
%         '31-Oct-2021','21-Nov-2021';...
%         };
%     
    % select just some of those events
    event_sel = [3,4]; %just a selection of all the events
    event_sel = [1:5]; %just a selection of all the events
    

    % some plot settings
    use_addaxis = 0; %option to use the addaxis function
    show_TTDs = 0;
    cols = lines(2);
    cols4ev = lines(size(events,1)+1);
    cols4ev(1:end-1,:) = cols4ev(2:end,:); %nicer

    % INDIVIDUAL EVENTS
    if flag_plot.events==1
        run('plot_events.m')
    end

    % ALL C-Q, delta-Q and age-Q relationships; and all C-age, delta-age relationships
    if flag_plot.eventstats==1
        run('plot_eventstats.m')
    end

    % SHOW C-DELTA RELATIONSHIPS
    if flag_plot.deltavsC == 1
        show_batch = 1; % show results on top of a batch-reactor simulation?
        run('plot_deltavsC.m')
    end


%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%  END OF FILE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%