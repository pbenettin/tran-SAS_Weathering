% this file stores all the main model configurations

% GENERAL SIMULATION SETTINGS
% select the CASE STUDY and the MODEL
datasetName = 'SyntheticData.csv';  %this is the fast-responding hydrology
ModelName = 'SAS_EFs_isoweathering';  %weathering model with isotope fractionation
outfilename = 'all_output';     %output filename (both for .mat and .csv outputs)

% select the AGGREGATION timestep and the STORAGE threshold (numerical parameters)
dt_aggr = 4;          %dt for the computations [h] (must be integer)
f_thresh = 1;          %fraction of rank storage [-] after which the storage is sampled uniformly
                       %(just leave f_thresh=1 if not interested in this option)

% some further options
create_spinup = true;        %generate a spinup period
extract_agedistrib = false;  %save TTDs over selected days
save_output = true;          %save model variables
display_output = true;       %display the output at the end
export_output = false;        %print output to csvfile


% MODEL PARAMETERS
% select the SAS
data.SASQName='fSAS_pltv'; %time-variant power-law SAS
data.SASETName='fSAS_step'; %step SAS

% set parameter names and values
SASQparamnames={'zmin_Q','zmax_Q'}; Param.SAS(1:2)=[.3, 1]; %[-]
SASETparamnames={'k_ET'};    Param.SAS(3)=1;  %[-]
otherparamnames={'S0'};      Param.SAS(4)=1000; %[mm]

% set chemical parameters and values (2 reactions, for major and rare isotopes)
data.frac_reaction = @c_n_FF; %rare isotope reaction model (Fixed Fractionation model)
chemparamnames(1) = {'Clim1'};  Param.chem(1) = 800; %[umol/l] limiting C for major isotope reaction 1
chemparamnames(2) = {'Clim2'};  Param.chem(2) = 150; %[umol/l] limiting C for major isotope reaction 2
chemparamnames(3) = {'lambda1'}; Param.chem(3) = 1/(800*24); %[1/h] kinetic constant for major isotope reaction 1 (k/Clim in the paper)
chemparamnames(4) = {'lambda2'}; Param.chem(4) = 1/(300*24); %[1/h] kinetic constant for major isotope reaction 2 (k/Clim in the paper)
chemparamnames(5) = {'alph_f'}; Param.chem(5) = 0.998; %[-] fractionation factor rare isotope reaction 2
chemparamnames(6) = {'rlim1'};  Param.chem(6) = 1E-02; %[-] factor to compute limiting C for rare isotope reaction 1
chemparamnames(7) = {'rlim2'};  Param.chem(7) = 0.9995E-02; %[-] factor to compute limiting C for rare isotope reaction 2

% SELECT SOME AGE DISTRIBUTIONS TO EXTRACT
% choose dates for which TTDs are saved (only used if extract_agedistrib=1)
data.datesel={...
    '24-Jun-2020';...
    '29-Jun-2020';...
    };

% spinup settings (only used if create_spinup=1)
data.C_S0 = 0; %concentration of the initial storage (no relevance for weathering model)
period_rep = [1,365*2]; %starting and ending day of the datasets which will be repeated
n_rep = 3; % (integer) number of times that period_rep will be repeated; if n_rep=0, then no spinup is created
data.show_spinup = false; %show spinup along with the results?

% display some stuff (optional flags)
flag_plot.active = 1;  %to activate/disactivate all plotting
flag_plot.Cout = 1;       %main plot with the full concentration timeseries
flag_plot.TTDs = 0;       %plot TTDs for selected days
flag_plot.agestats = 0;   %plot with additional age timeseries
flag_plot.events = 1;     %plot with timeseries of the individual events
flag_plot.eventstats = 0; %plot with other stats during the events
flag_plot.deltavsC = 0;   %plot with deltavsC in batch reaction and simulations
