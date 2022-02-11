% space to test and implementing chemical kinetics
addpath(genpath(fullfile('..','Source')))


%--------------------------------------------------------------------------
% CASE OF MAJOR AND FRACTIONATING ISOTOPE
%--------------------------------------------------------------------------
    
% here we simplify to two reactions where Clim1>Clim2 (so it will be
% reaction 2 to cause precipitation and isotopic fractionation)

% create an age vector
dt=24; %time step [hours]
N=500*24/dt; %number of timesteps (at dt=1 hour)
T=(0:dt:(N-1)*dt)'; %age vector in dt timesteps

% PARAMETERS -------------------------

% select the rare isotope reaction
frac_reaction = @c_n_FF; %fixed-fractionation model

% parameters for the major species
par.Clim=[800,150]; %limit concentration for each reaction
par.lambda=[1/(800*24), 1/(300*24)]; %[1/h] %kinetic constants for each reaction
par.C0=0; %initial condition 
%
% coefficients for the rare isotope
par.alph_f=[1, .998];
par.rlim=[10^-2, 0.9995*10^-2];
par.r0=10^-2;

% from hereon, all is automatic
%-------------------------------------

% compute useful variables
K_m=sum(par.lambda); %overall kinetic constant for isotope m
CLIM_m=(par.lambda*par.Clim')/K_m; %overall limiting concentration for isotope m
CMAX_n = max(frac_reaction(linspace(0,10/par.lambda(2),1/par.lambda(2)),par)); %maximum concentration value for rare isotope
dMAX = (CMAX_n/CLIM_m-par.rlim(1))/par.rlim(1)*1000; %estimate of maximum delta value

% compute the solution for each species (using external functions c_m.m and <frac_reaction>.m)
C_m = c_m(T,par);
C_n = frac_reaction(T,par);

% compute the isotope ratios and delta values
r_f = C_n./C_m; 
d_f = (r_f-par.r0)/par.r0*1000; 

% -------------------------------------------------------------------------
% FIGURES
% -------------------------------------------------------------------------

% SHOW the REACTIONS
figure(1)

s1=subplot(3,1,1,'NextPlot','add');
title('reaction for major isotope (m)')
% plot([1/K_nf/24,1/K_nf/24],[0 max(par.Clim)],...
%     'Color',[.9 .9 .9],'HandleVisibility','off') %this is the timescale of the overall kinetic
plot([0 T(end)/24],[par.Clim(1) par.Clim(1)],'--k','DisplayName','Clim 1')
plot([0 T(end)/24],[par.Clim(2) par.Clim(2)],'-.k','DisplayName','Clim 2')
plot([0 T(end)/24],[CLIM_m CLIM_m],'-k','DisplayName','CLIM')
plot(T./24,C_m,'b-','DisplayName','major isotope');
ylim([0 1.05*par.Clim(1)])
ylim([0 +Inf])
xlim([0 max(T)/24])
xlabel('age [d]')
ylabel('Concentration')
%legend(gca,'Location','SE')

s2=subplot(3,1,2,'NextPlot','add');
axes(s2)
title('reaction for fractionating isotope (n)')
%plot([0 T(end)/24],[CLIM_af CLIM_af],'-k','DisplayName','CLIM')
% plot([p_start/24 p_start/24],[0 rlim(1)*par.Clim(1)],...
%     'Color',[.9 .9 .9],'HandleVisibility','off') %this is the time when precipitation starts    
plot(T./24,C_n,'r','DisplayName','fractionating isotope');
ylim([0 1.05*max(C_n)])
ylim([0 +Inf])
xlim([0 max(T)/24])
xlabel('age [d]')
ylabel('Concentration')
%legend(gca,'Location','SE')

s3=subplot(3,1,3,'NextPlot','add');
axes(s3)
title('isotope ratio')
% plot([0 T(end)/24],[dLIM_af dLIM_af],'-k','DisplayName','\deltaLIM')
plot(T./24,d_f,'g','LineWidth',1,'DisplayName','\delta');
xlabel('age [d]')
ylabel(['\delta [',char(8240),']'])
ylim([-Inf +Inf])
xlim([0 max(T)/24])
%legend(gca,'Location','SE')

linkaxes([s1,s2,s3],'x')


% SHOW the delta-C RELATIONSHIP
figure(2)
set(gca,'NextPlot','add','TickDir','out','box','on');
title('\delta VS C')
plot(C_m/CLIM_m,d_f,'r-','LineWidth',1);
% xlim([0 1])
xlabel('C/CLIM')
ylabel(['\delta [',char(8240),']'])




