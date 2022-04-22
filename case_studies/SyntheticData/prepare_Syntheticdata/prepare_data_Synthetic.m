% generate synthetic hydrologic data for running reactive transport models

clear variables
clc
addpath('functions')

% -------------------------------------------------------------------------
% ASSIGN MAIN PARAMETERS AND GENERAL SETTINGS
% -------------------------------------------------------------------------

flag_dataexport = 0; %export data as textfile (quite slow file creation... )
filename='SyntheticData'; %export filename

% timeseries duration and dt
date0 = '01-Jan-2020';   %starting date of simulations
param.Nday = 365*2;      %days of simulation [d]

% important parametrers
param.Pyear = 1200;   %mean annual precip [mm/y]
param.interarr = 10;   %mean interarrival between events [d]
choice_parset = 1;    %1 is fast-responding system, 2 is slow responding
flag_transient = 1;   % 1 is transient fluxes, 0 is constant fluxes

% some other settings
%rng('default'); %set the random number generator
rng(1); %option to generate always the same random number sequence
plot_timeseries = 0;
plot_partitioning = 0;
plot_cumtimeseries = 0;

% -------------------------------------------------------------------------
% SET OTHER STUFF (can keep fixed)
% -------------------------------------------------------------------------
if choice_parset == 1 %fast-responding system
    % fast responding system (with ~70% soil water contributions)
    param.frech = 0.3;    %fraction of leakage that goes to gw [-]
    param.Z = 500;        %soil depth [mm]
    param.Ksat = 250;     %saturated hydraulic conductivity [mm/d]
    param.c = 8;          %leakage exponent [-]
    param.Kgw = 5;        %coefficient to produce gw flow [mm/d]
    param.cgw = 5;        %exponent to produce gw flow [-]
    param.Sgwmax = 200;   %max gw value (to normalize Sgw/Sgwmax)
end

if choice_parset == 2 %slow system
    % slow responding system (with ~40% soil water contributions)
    param.frech = 0.6;    %fraction of leakage that goes to gw [-]
    param.Z = 1200;       %soil depth [mm]
    param.Ksat = 50;      %saturated hydraulic conductivity [mm/d]
    param.c = 3;          %leakage exponent [-]
    param.Kgw = 50;       %coefficient to produce gw flow [mm/d]
    param.cgw = 3;        %exponent to produce gw flow [-]
    param.Sgwmax = 50;    %max gw value (to normalize Sgw/Sgwmax)
end

if choice_parset == 3  %other testing
    param.frech = 0.3;    %fraction of leakage that goes to gw [-]
    param.Z = 500;        %soil depth [mm]
    param.Ksat = 100;     %saturated hydraulic conductivity [mm/d]
    param.c = 4;          %leakage exponent [-]
    param.Kgw = 5;        %coefficient to produce gw flow [mm/d]
    param.cgw = 5;        %exponent to produce gw flow [-]
    param.Sgwmax = 200;   %max gw value (to normalize Sgw/Sgwmax)
end

% rainfall generation parameters
param.ndown = 6;      %number of downscalings of daily precip (5 gives a dt of 45 min)

% soil properties and groundwater balance parameters
param.n = 0.42;       %porosity [-]
param.sfc = 0.4;      %soil field capacity [-]

% ET parameters [IRRELEVANT IF ETMAX = 0]
param.ETmax = (0*param.Pyear/365);   %max yearly ET [mm/d]
param.ETampl = 0;               %amplitude of ET seasonality [-]
param.sstar = 0.35;             %water stress threshold on soil moisture [-]
param.sw = 0.2;                 %wilting point [-]
param.sh = 0.1;                 %hygroscopic point [-]
param.ETpart = 0.2;             %value of E/ET

% initial conditions
s0 = 0.7;                 %initial soil moisture [-]

% 18O isotope generation
param.m18O = -8;        %mean d18O value [permil]
param.amp18O = 2;       %amplitude of the d18O sinusoid [permil]
param.phase18O = 90;    %phase of the sinusoid [d]
param.lmwl1 = 8;        %slope LMWL [-]
param.lmwl2 = 10;       %intercept LMWL [permil]
param.noiseAmplitude = 0.5; %mean noise in d18O [permil]
param.noisecorr = 0.95; %lag-1 correlation in the d18O noise [-]



% -------------------------------------------------------------------------
% RAINFALL SERIES GENERATION
% -------------------------------------------------------------------------

% compute parameters for poisson rainfall generation (at daily time step)
lambda_P=1/param.interarr;  %poisson parameter for precipitation [d^-1]
gamma_P=(365/param.Pyear/param.interarr); %rain parameter [mm^-1] (inverse of the daily mean)

% generate rainfall at dtgen timesteps
[J,dtgen] = generate_rainfall(param.Nday,lambda_P,gamma_P,param.ndown); %J in [mm/d], dtgen in [d]
Nt = param.Nday/dtgen;            %number of timesteps
time = ((0:Nt-dtgen)*dtgen)';     %time in days from 0 to end-dt

% now interpolate to get hourly values
dtgen_h = 1/24; %now hourly time steps, measured in days
time_h = (time(1):dtgen_h:time(end))';
J_h = diff([0;interp1(time,cumsum(J)*dtgen,time_h)])/dtgen_h;
dtgen = dtgen_h;
time = time_h;
J = J_h;

% -------------------------------------------------------------------------
% POTENTIAL E and T GENERATION
% -------------------------------------------------------------------------
% ET0 = param.ETmax*ones(Nt,1); % mm/d constant
ET0 = (param.ETmax+param.ETampl*sin(2*pi*(time-90)/365)); %mm/d seasonal
E0 = param.ETpart*ET0;
T0 = (1-param.ETpart)*ET0;


% -------------------------------------------------------------------------
% SOIL MOISTURE BALANCE
% -------------------------------------------------------------------------
% (the equation is n*Z*ds/dt = J-E-T-L-R)
% simple Euler Forward at timesteps dtgen
data = sbalance_EF(s0,dtgen,J,E0,T0,param); %data table with moisture [-] and fluxes [mm/d]
data.time = time; % time in [d]

% -------------------------------------------------------------------------
% HYDROLOGIC MODEL
% -------------------------------------------------------------------------
% simple 2-bucket model to produce Q
data = bucket_model(data,param,dtgen);


% -------------------------------------------------------------------------
% ISOTOPES IN RAINFALL
% -------------------------------------------------------------------------
% add new columns d18O and d2H to the data table
data = generate_inputisotopetracer(data,param);



% -------------------------------------------------------------------------
% MAKE SOME PLOTS
% -------------------------------------------------------------------------

% some settings
xlm=[time(end)-365,time(end)]; % x-axis limits: for simplicity show only a subportion

% plotting the fluxes in a time series
if plot_timeseries == 1
    
    figure
    
    s1=subplot(3,1,1);
    hold on
    p1 = plot(time,data.J); %show rainfall
    ylabel('rainfall [mm/d]')
    set(gca,'TickDir','out','XLim',xlm,'XAxisLocation','top','YDir','reverse')
    s2=subplot(3,1,2);
    hold on
    plot(time,data.s)
    ylabel('s [-]') %show moisture
    set(gca,'TickDir','out','XLim',xlm)
    xlabel('time [d]')
    s3=subplot(3,1,3);
    hold on
    plot(time,data.R) %overland flow
%     plot(time,data.R+(1-param.frech)*data.L) %soil water flow
    plot(time,data.Q-(data.R+(1-param.frech)*data.L)) %gw flow
    plot(time,data.Q) %total runoff
    ylabel('Q [mm/d]')
    set(gca,'TickDir','out','XLim',xlm,'YScale','linear')
    xlabel('time [d]')
    linkaxes([s1,s2,s3],'x')
end

% plot with flow and partitioning
if plot_partitioning == 1
    figure(100)
    subplot(2,1,1)
    hold on
    title('Q and Qs/Q')
    plot(time,data.Q) %Q
    ylabel('Q [mm/d]') %
    set(gca,'TickDir','out','XLim',xlm)
    xlabel('time [d]')    
    subplot(2,1,2)
    hold on
    plot(time,data.f) %ratio Qs/Q
    plot(time,(data.s-min(data.s))./(max(data.s)-min(data.s))) %soil moisture
    tmp = cumsum(data.J-data.E-data.T-data.R-data.Q)*dtgen;
    plot(time,tmp/max(tmp)) %normalized totalstorage variations
    ylabel('f [-]') %
    set(gca,'TickDir','out','XLim',xlm)
    xlabel('time [d]')
    ylim([0,1])
end

% plot with cumulative balance
if plot_cumtimeseries == 1
    figure
    hold on
    plot(time,cumsum(data.J)*dtgen,'DisplayName','J')
    plot(time,cumsum(data.Q)*dtgen,'DisplayName','Q')
    plot(time,cumsum(data.R+(1-param.frech)*data.L)*dtgen,'DisplayName','Qs')
    plot(time,cumsum(data.Q-(data.R+(1-param.frech)*data.L))*dtgen,'DisplayName','Qgw')
    plot(time,cumsum(data.E+data.T)*dtgen,'DisplayName','ET')
    set(gca,'TickDir','out','XLim',xlm)
    xlabel('time [d]')
    ylabel('flux [mm]')
    box on
    legend(gca,'Location','NorthWest')
end




%--------------------------------------------------------------------------
% some final stuff before exporting the data
%--------------------------------------------------------------------------
Cout = -999*ones(size(data.J));
dates = datestr(datenum(date0)+time,'yyyy-mm-dd HH:MM');
if flag_transient ~= 1
    data.J = repmat(mean(data.J),size(data.J)); %mm/d
    data.Q = repmat(mean(data.Q),size(data.J)); %mm/d
    data.E = repmat(mean(data.E),size(data.J)); %mm/d
    data.T = repmat(mean(data.T),size(data.J)); %mm/d
    data.f = repmat(mean(data.f),size(data.J)); %-
end


%--------------------------------------------------------------------------
% export to textfile for tranSAS
%--------------------------------------------------------------------------

if flag_dataexport == 1
    fprintf('Exporting data to file %s...',[filename,'.csv'])
   
    % export as fixed-delimited text (former tran-SAS standard)
%     fid=fopen(['../',filename,'.dat'],'w');
%     
%     fprintf(fid,'%25s %10s %10s %11s %10s %10s %14s\n',... %headers
%         'date [yyyy-mm-dd HH:MM]','J [mm/h]','Q [mm/h]','ET [mm/h]',...
%         'Cin [-]','wi [-]','measC_Q [-]');
%     for i=1:length(J)
%         sprintf('line %.0f, %.1f%%',i,i/length(J)*100)
%         fprintf(fid,'%25s %10.5f %10.5f %11.5f %10.5f %10.5f %14.5f\n',...
%             datestr(dates(i,:),'yyyy-mm-dd HH:MM'),data.J(i)/24,data.Q(i)/24,...
%             (data.E(i)+data.T(i))/24,...
%             data.d18O(i),data.f(i),Cout(i));
%     end
%     fclose(fid);
    
    % print a csv file with data in mm/h
    data.date = datetime(datestr(dates));
    data.J = data.J/24;
    data.Q = data.Q/24;
    data.ET = (data.E+data.T)/24;
    data.measC_Q = NaN(size(data.J));
    data.Cin = data.d18O;
    data.wi = data.f;
    data4print = data(:,{'date','J','Q','ET','Cin','wi','measC_Q'});
    writetable(data4print,fullfile('..',[filename,'.csv']));
    
    % print a readme file with some info
    fid = fopen(fullfile('..',[filename,'_README.txt']),'w');
    fprintf(fid,'# Automatic README file generated for ''%s''.csv\n\n',filename);
    fprintf(fid,'The file ''%s''.csv includes tabular data at hourly timesteps:\n\n',filename);
    fprintf(fid,'- date, formatted as datetime dd-mmm-yyyy HH:MM:SS (e.g. 01-Jan-2020 00:00:00)\n');
    fprintf(fid,'- J, precipitation in mm/h\n');
    fprintf(fid,'- Q, streamflow in mm/h\n');
    fprintf(fid,'- ET, evapotranspiration in mm/h\n');
    fprintf(fid,'- Cin, input tracer concentration (here, d18O in permil)\n');
    fprintf(fid,'- wi, indicator of wetness index, non-dimensional and confined within [0-1]. Here, it is the ratio between soil runoff versus total runoff\n');
    fprintf(fid,'- measC_Q, output sampled concentration (here it''s all NaN)\n\n');
    fprintf(fid,'Missing values are indicated as NaN');
    fclose(fid);

    fprintf(' done\n')
end


