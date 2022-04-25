# _tran_-SAS_Weathering v1.0

This code combines the time-variable transit time distributions provided by the basic _tran_-SAS model with chemical kinetics equations for weathering reactions.

The code is used to generate realistic solute and isotope dynamics in streams during storm events.

The analysis of these results are included in a paper submitted for publication on a scientific journal.

The code is entirely written in Matlab. Any Matlab version from 2016a should work, on any operative system.


## Running the model

The current model version can be run by executing the script `model_STARTER_weathering.m`. A test synthetic dataset (`SyntheticData`) is provided as case study. For testing, you can explore the synthetic case study and modify the model parameters and simulation settings through the configuration file (case_studies > SyntheticData > `config_file.m`).

To run the model on your own data, you need to prepare a new case study:

- create a new directory inside the 'case_studies' folder
- prepare a .csv data table formatted as SyntheticData.csv (see also the associated README file)
- fill-in a config_file for your case study
- copy the new folder name into the text file case_study_name, such that the model starter selects the right case study
