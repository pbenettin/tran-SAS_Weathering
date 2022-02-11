%**************************************************************************
% this is a model starter to prepare and call the model launcher. It can
% generate multiple parameter combinations from lists of parameter values
%**************************************************************************

% some initialization
clear variables
home_folder = pwd;
cd('..')
addpath(genpath('Source'))

%% ALL PARAMETERS AND SETTINGS
% LISTS WITH PARAMETER VALUES

% transport parameters (SASQ is pltv, SASET does not matter)
list.z_Qmin = [0.2,0.3,0.5]; %min k_Q value [-]
list.z_Qmax = [0.8,1,1.2]; %max k_Q value [-]
list.S0 = [300,1000,3000]; %[mm]

% chemical parameters
data.frac_reaction = @c_n_FF; %rare isotope reaction
list.chem1 = 800; %[umol/l] Clim1 limiting C for major isotope reaction 1
list.chem2 = 150; %[umol/l] Clim2 limiting C for major isotope reaction 2
list.chem3 = 1/(800*24); %[1/h] lambda1 kinetic constant for major isotope reaction 1 (k/Clim in the paper)
list.chem4 = 1/(300*24); %[1/h] lambda2 kinetic constant for major isotope reaction 2 (k/Clim in the paper)
list.chem5 = 0.998; %[-] alph_f fractionation factor rare isotope reaction 2
list.chem6 = 10^-2; %[-] rlim1 factor to compute limiting C for rare isotope reaction 1
list.chem7 = 0.9995*10^-2; %[-] rlim2 factor to compute limiting C for rare isotope reaction 2

% test a single run
%{
data.frac_reaction = @c_n_FF; %rare isotope reaction
list.z_Qmin = 0.3; %min k_Q value [-]
list.z_Qmax = 0.8; %max k_Q value [-]
list.S0 = 1000; %[mm]
list.chem1 = 300; %[umol/l] Clim1 limiting C for major isotope reaction 1
list.chem2 = 100; %[umol/l] Clim2 limiting C for major isotope reaction 2
list.chem3 = 1/(100*24); %[1/h] lambda1 kinetic constant for major isotope reaction 1 (k/Clim in the paper)
list.chem4 = 1/(200*24); %[1/h] lambda2 kinetic constant for major isotope reaction 2 (k/Clim in the paper)
list.chem5 = 0.98; %[-] alph_f fractionation factor rare isotope reaction 2
list.chem6 = 10^-2; %[-] rlim1 factor to compute limiting C for rare isotope reaction 1
list.chem7 = 10^-2; %[-] rlim2 factor to compute limiting C for rare isotope reaction 2
%}

% GENERAL SIMULATION SETTINGS

% select the CASE STUDY and the MODEL 
data.case_study = 'Syntheticdata_fast2';
data.case_study_path = fullfile('..','case_studies',data.case_study);
datasetName='Syntheticdata_fast2.csv';  %text file with data (check formatting specs)
ModelName='SAS_EFs_isoweathering';  %more advanced weathering model

% select the AGGREGATION timestep and the STORAGE threshold (numerical parameters)
dt_aggr=4;          %dt for the computations [h] (must be integer)
f_thresh=1;          %fraction of rank storage [-] after which the storage is sampled uniformly
                     %(just leave f_thresh=1 if not interested in this option)                                        
                     
% some further options
create_spinup=1;     %generate a spinup period: 1=yes
save_output=1;       %save model variables: 1=yes
load_output=0;       %display the output at the end: 1=yes 

%% (automatic) Setup the folders and nested loops and run them

% create a new temporary folder where to store all the simulations
tmp=datevec(now);
foldername=fullfile(home_folder,'temp_results',sprintf('%d%02d%02d%02d%02d%02.0f',tmp));
mkdir(foldername)
copyfile(fullfile(home_folder,'q_types_model.m'),fullfile(foldername,'q_types.m'));

% store fieldname into a variable
fieldvar=fieldnames(list);

% count how many combinations
N=1; %just an initialization
for i=1:numel(fieldvar)
    N = N*length(list.(fieldvar{i}));
end
fprintf('\nnumber of parameter combinations is %.f \n',N)

% generate a temp file with nested loops automatically generated from structure list
tempfilename='temp_loops_script.m';
fid=fopen(tempfilename,'w');
fprintf(fid,'%% temporary file to launch the model within a nested loop of variables depth. Regenerates everytime\n');
fprintf(fid,'%% last generated: %s\n',datestr(now));
fprintf(fid,'\n%% preallocate a table with the parameters of each simulation\n');
fprintf(fid,'M=zeros(N,length(fieldvar)+1);\n');
fprintf(fid,'\nsimul_ID=0;\n');
for i=1:numel(fieldvar) 
    fprintf(fid,'for %s = list.%s\n',fieldvar{i},fieldvar{i});
end
fprintf(fid,'simul_ID=simul_ID+%d;\n',1);
fprintf(fid,'M(simul_ID,1)=simul_ID;\n');
for i=1:numel(fieldvar) 
    fprintf(fid,'M(simul_ID,%d)=%s;\n',i+1,fieldvar{i});
end
fprintf(fid,'run(''run_multiple/model_launch.m'')\n'); 
for i=1:numel(fieldvar) 
    fprintf(fid,'end\n');
end
fprintf(fid,'save([foldername,''/parameter_summary''],''M'',''fieldvar'');\n');
fclose(fid);


% run the tempfile with the nested loops and delete it afterwards
run(tempfilename)
fprintf('saved to file %s\n',foldername);
delete(tempfilename)

cd(home_folder)


%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%  END OF FILE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%