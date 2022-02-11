% this script cannot be launched alone. It is called from other scripts. 
% It takes settings and simply combines them and launches the model

% check that this is not run standalone
if ~exist('datasetName','var')
    error('Script ''model_launch.m'' cannot be run standalone. It needs to be called from other scripts')
end

% store some variables into the structure 'data'
data.f_thresh=f_thresh;       %pass storage threshold information to the model
data.save_output=save_output;  %pass saving information to the model
data.simul_ID=simul_ID; %pass ID information to the model
data.foldername=foldername; %pass the folder where to store the data
data.outfilename = 'sprintf(''%s/output_%05d'',data.foldername,data.simul_ID)';

% HYDROLOGIC DATA IMPORT
% hourly data are imported from csv file 'datasetName' and aggregated to
% timestep dt_aggr
data = fdata_aggregate(data,fullfile(data.case_study_path,datasetName),dt_aggr);

% % define here the precipitation event that we focus on
% % example of multiple events for synthetic data with J interarr 10 days ('rare' datasets)
% events={...
%     '22-Jan-2020','03-Feb-2020';...
%     '03-Dec-2020','18-Dec-2020';...
%     };
% data.event_list = events;

% select dates to compute TTDs (must give this, but does not matter)
dateform='dd-mmm-yyyy'; %date format
data.datesel = cellstr(datestr(data.dates(end-365*(24/data.dt):(24/data.dt)*190:end),dateform));
data.index_datesel=length(data.J);

% SET THE transport PARAMETERS

% select the SAS functions for Q and ET
data.SASQName='fSAS_pltv'; %time variant power-function SAS
data.SASETName='fSAS_pl'; %power-function SAS

% transport parameter names and values
SASQparamnames={'z_Qmin','z_Qmax'};   Param.SAS(1:2)=[z_Qmin,z_Qmax];  %[-]
SASETparamnames={'k_ET'};    Param.SAS(3)=1;  %[-]
otherparamnames(1)={'S0'};   Param.SAS(4)=S0; %[mm] 

% set chemical parameters and values (2 reactions, for major and rare isotopes)
chemparamnames(1) = {'Clim1'};  Param.chem(1) = chem1; %[mg/l] limiting C for major isotope reaction 1
chemparamnames(2) = {'Clim2'};  Param.chem(2) = chem2; %[mg/l] limiting C for major isotope reaction 2
chemparamnames(3) = {'lambda1'}; Param.chem(3) = chem3; %[1/h] kinetic constant for major isotope reaction 1
chemparamnames(4) = {'lambda2'}; Param.chem(4) = chem4; %[1/h] kinetic constant for major isotope reaction 2
chemparamnames(5) = {'alph_f'}; Param.chem(5) = chem5; %[-] fractionation factor rare isotope reaction 2
chemparamnames(6) = {'rlim1'};  Param.chem(6) = chem6; %[-] factor to compute limiting C for rare isotope reaction 1
chemparamnames(7) = {'rlim2'};  Param.chem(7) = chem7; %[-] factor to compute limiting C for rare isotope reaction 2

% (automatic) unify the parameters into one series Pars
Pars = [Param.SAS,Param.chem];
data.SASQl=length(SASQparamnames);
data.SASETl=length(SASETparamnames);
data.paramnames=[SASQparamnames,SASETparamnames,otherparamnames,chemparamnames];  

% 5-INITIAL CONDITIONS

% assign the concentration of the initial storage
data.C_S0=0; %this has no relevance right now   

% spinup settings (only used if create_spinup=1)
period_rep=[1,365*2]; %starting and ending day of the datasets which will be repeated
n_rep=3; % (integer) number of times that period_rep will be repeated; if n_rep=0, then no spinup is created
data.show_spinup = false; %show spinup along with the results?

% (automatic) if needed, create the spinup
if create_spinup==1
    data=fgenerate_spinup(data,period_rep,n_rep); %use an external function
else
    data.ini_shift=0; %no shift in the dataset if no spinup is used
end


% RUN THE MODEL with the current settings

% let's go for the model run
fprintf('model run: %d/%d (%.1f%%) \n',simul_ID,N,simul_ID/N*100)
%tic
feval(ModelName,Pars,data); %run the model 'ModelName' using parameters 'Pars' and data 'data'
%toc





%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%  END OF FILE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%