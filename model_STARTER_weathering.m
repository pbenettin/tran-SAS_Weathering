% This script imports data from a case study and runs the SAS weathering
% model to compute reactive solute and isotope transport and transit time
% distributions

% TO RUN THE MODEL, you need to prepare a case study first:
% - create a new folder inside the 'case_studies' folder
% - prepare a config file for this case study
% - copy the new folder name into the text file case_study_name
% After the case study folder is created, you can run this STARTER

% The STARTER is all automatic and you do not need to modify it
% It is organized as follows:
% 1) initial settings (e.g. choose model, case study, etc... )
% 2) load data from the case study and modify the timestep
% 3) set parameters for the specific case study and model
% 4) (if desired) select dates in which to save age computations
% 5) (if required) set/generate initial conditions
% 6) run the model
% 7) call some simple results visualization
% 8) export results to csv file

clear variables
addpath(genpath('Source'))


%--------------------------------------------------------------------------
% 1-GENERAL SIMULATION SETTINGS
%--------------------------------------------------------------------------

% identify the case study from textfile case_study_name
data.case_study = fscanf(fopen('case_study_name.txt'),'%s');
data.case_study_path = fullfile('case_studies',data.case_study);

% create a 'results' directory in the case study folder (if it does not exist)
[status,msg] = mkdir(data.case_study_path,'results');

% run the configuration file to define all settings
run(fullfile(data.case_study_path,'config_file'))


%--------------------------------------------------------------------------
% 2-DATA IMPORT and AGGREGATION
%--------------------------------------------------------------------------

% all model data (i.e. hydrologic variables, settings etc) is passed to the
% model through a structure named 'data'
data.f_thresh=f_thresh;       
data.save_output=save_output; 
data.extract_agedistrib=extract_agedistrib;
data.outfilename = outfilename;

% hourly data are imported from csv file 'datasetName' and aggregated to
% timestep dt_aggr
data = fdata_aggregate(data,fullfile(data.case_study_path,datasetName),dt_aggr);


%--------------------------------------------------------------------------
% 3-PARAMETERS
%--------------------------------------------------------------------------

% parameters are defined within the config_file and displayed here

% (automatic) unify the parameters into one series and print some info to screen
Pars = [Param.SAS,Param.chem];
data.SASQl=length(SASQparamnames);
data.SASETl=length(SASETparamnames);
data.paramnames=[SASQparamnames,SASETparamnames,otherparamnames,chemparamnames];  
fprintf('\n')
fprintf('Model:   %s\n',ModelName)
fprintf('dt:      %d h\n',dt_aggr)
fprintf('SAS Q:   %s\n',data.SASQName)
fprintf('SAS ET:  %s\n',data.SASETName)
fprintf('parameters:\n')
for i=1:length(Pars)
    fprintf(1,'%9s = %7.4f \n',char(data.paramnames(i)),Pars(i));
end



%--------------------------------------------------------------------------
% 4-SELECT THE AGE DISTRIBUTIONS TO EXTRACT
%--------------------------------------------------------------------------
% define a vector of dates in matlab format. Discharge age distributions 
% will be evaluated on the closest available day and will be stored into 
% the matrix age_matr. If not interested in evaluating age distributions,
% simply leave data.datesel={}

% dates are defined in the config file

% % otherwise, select automatically a set of dates
% data.datesel = cellstr(datestr(data.dates(end-365*(24/data.dt):(24/data.dt)*190:end),dateform));

% (automatic) check that desired dates actually exist and find the index
% of the closest time in the dataset
dateform='dd-mmm-yyyy'; %date format
tmp=zeros(size(data.datesel));
for i=1:size(data.datesel,1)
    if datenum(data.datesel(i))<data.dates(1) ||...
            datenum(data.datesel(i))>data.dates(end)
        warning(['the date ',data.datesel{i},...
            ' selected for streamflow age evaluation ',...
            'does not exist in the dataset'])
    else
        [~,tmp(i)]=min(abs(data.dates-...
            datenum(data.datesel(i),dateform))); %closest date in the dataset
    end
end
data.index_datesel=tmp(tmp>0);


%--------------------------------------------------------------------------
% 5-INITIAL CONDITIONS
%--------------------------------------------------------------------------
 
% if needed, create the spinup through an external function
if create_spinup==1
    data=fgenerate_spinup(data,period_rep,n_rep); %use an external function
else
    data.ini_shift=0; %no shift in the dataset if no spinup is used
end



%--------------------------------------------------------------------------
% 6-RUN THE MODEL
%--------------------------------------------------------------------------

% let's go for the model run
fprintf('\ncalculating model output... ')
tic
out = feval(ModelName,Pars,data); %run the model 'ModelName' using parameters 'Pars' and data 'data'
toc



%--------------------------------------------------------------------------
% 7-SHOW SOME OUTPUT
%--------------------------------------------------------------------------

if display_output==1
    
    % load the results that were just created
    load(fullfile(data.case_study_path,'results',data.outfilename))
    
    % plot the output
    if flag_plot.active == 1
        % check if a script 'plot_results_isoweathering' exists within the
        % case study folder, and in case execute it
        if exist(fullfile(data.case_study_path,'plot_results.m'),'file') ~= 0
            fprintf('\nplotting case-study-specific results...\n')
            run(fullfile(data.case_study_path,'plot_results'));
        else % otherwise run the default plotting script
            fprintf('\nplotting default results...\n')
            run(fullfile('Source','plotting_scripts','plot_results_isoweathering')); %run a script that plots some results
        end
    end
    
end

%--------------------------------------------------------------------------
% 8-PRINT SOME OUTPUT TO CSVFILE
%--------------------------------------------------------------------------

if save_output == 1 && export_output == 1
    
    % load the results that were just created
    load(fullfile(data.case_study_path,'results',data.outfilename))
    
    % make a table with various timeseries (removing the spinup)
    per = data.ini_shift+1:length(data.dates);
    T = array2table(datetime(data.dates(per),'ConvertFrom','datenum'),...
        'VariableNames',{'timestamp'}); %time stamp
    T.J = data.J(per);   %precipitation [mm/h]
    T.Q = data.Q(per);   %streamflow [mm/h]
    T.ET = data.ET(per); %ET [mm/h]
    T.wi = data.wi(per); %index of catchment wetness [mm/h] 
    T.C_J = data.C_J(per); %concentration in precipitation
    T.C_Qm = C_Qm(per);  %concentration for the major isotope [umol/l]
    T.C_Qm_NP = C_Qm_NP(per); %concentration for the major isotope in the absence of precipitation [umol/l] 
    T.d_f = d_f(per); %isotope ratio [permil]
    for i=1:size(Fyw,2) %the young water fractions (name is F_(threshold in days))
        tmp = ['F_',num2str(ywt(i),'%.0f')];
        T.(tmp) = Fyw(per,i);
    end
    
    % print the table to csv file
    writetable(T,fullfile(data.case_study_path,'results',[data.outfilename,'_timeseries.csv']));
    
end

%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%  END OF FILE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%