% space to test and implementing chemical kinetics
clear variables
addpath(genpath(fullfile('..','Source')))

% set the flags for what to show
plot_TTDshapes = 1;
plot_allTTDs = 0;
plot_reactions = 1;
plot_Da = 0;

%--------------------------------------------------------------------------
% MAJOR AND FRACTIONATING ISOTOPE WITH AGE DISTRIBUTION
%--------------------------------------------------------------------------

% here we have reactions and mixing of different parcels with different ages

% PARAMETERS -------------------------
% select the rare isotope reaction
frac_reaction = @c_n_FF; %fixed-fractionation model

% parameters for the major species (m)
par.Clim=[800,150]; %limit concentration for each reaction
par.lambda=[1/(800*24), 1/(300*24)]; %[1/h] %kinetic constants for each reaction
par.C0=0; %initial condition on parcel's concentration

% parameters for the rare isotope (n)
par.alph_f=[1, .998];
par.rlim=[10^-2, 0.9995*10^-2];
par.r0=10^-2;

% create an age vector
dT=24*1; %time step [hours]
N=200/par.lambda(2)/dT; %number of timesteps
T=(0:dT:(N-1)*dT)'; %age vector in dt timesteps

% define parameters of the gamma distributions to be used
a = [0.25, 0.5, 1, 5]; %'shape' parameter
Ng1 = length(a); %number of shape parameter values to be used
Ng2 = 100; %number of scale parameter values to be used
Ng22 = 40;
gm = logspace(log10(24*1),log10(24*5000),Ng2);
% from hereon, all is automatic
%-------------------------------------


% chemical equations
% compute useful variables
Lambda_m=sum(par.lambda); %overall kinetic constant for isotope m
CLIM_m=(par.lambda*par.Clim')/Lambda_m; %overall limiting concentration for isotope m
Lambda_af=sum(par.alph_f.*par.lambda); %overall kinetic constant for isotope n with fractionation
CLIM_af=((par.alph_f.*par.lambda)*(par.rlim.*par.Clim)')/Lambda_af; %overall limiting concentration for isotope n
dLIM_af = ((CLIM_af./CLIM_m)-par.r0)/par.r0*1000; %limiting delta values (useful for plots)

% compute the solution for each species (using external functions c_m.m and <frac_reaction>.m)
C_m = c_m(T,par);
C_n = frac_reaction(T,par);

% compute the isotope ratios and delta values
r_f=C_n./C_m; %fractionating only after the flag is activated
d_f = (r_f-par.r0)/par.r0*1000; %fractionating only upon precipitation

% show an example of the selected gamma distributions (with same mean)
if plot_TTDshapes == 1
    figure
    ax_T = axes('NextPlot','add','TickDir','out');
    sel_gm = 1; %mean [years]
    for i= 1:Ng1
        cumTTD=cdf('gamma',T/24/365,a(i),sel_gm/a(i)); %matlab built-in function 'cdf'
        TTD = diff([cumTTD;cumTTD(end)])/(dT/24/365);
        plot(T/24/365,TTD,'Parent',ax_T,...
            'LineWidth',1,...
            'DisplayName',sprintf('alpha = %.2f',a(i)))
    end
    legend(ax_T)
    set(ax_T,...
        'XLim',[0 sel_gm*3],...
        'YLim',[0 4])
    xlabel('Age [y]')
    ylabel('pdf [1/y]')
end

% start a loop to compute C for various TTDs
fprintf('\n Computing %d gamma distributions\n',Ng1*Ng2)
if plot_allTTDs == 1
    figure(101)
    ax = axes('Nextplot','add');
    xlabel('age [d]')
end
C_Qmstat = zeros(Ng2,Ng1,1);
C_Qnstat = zeros(Ng2,Ng1,1);
Da = zeros(Ng2,Ng1,1);
mage = zeros(Ng2,Ng1,1);
for i = 1:Ng1
    for j = 1:Ng2
        cumTTD=cdf('gamma',T,a(i),gm(j)/a(i)); %matlab built-in function 'cdf'
        TTD = diff([cumTTD;cumTTD(end)]);
        mage(j,i) = gm(j);
        Da(j,i) = gm(j)*Lambda_m; %Damkohler number = flow time scale / reaction time scale
        C_Qmstat(j,i) = C_m'*TTD;
        C_Qnstat(j,i) = C_n'*TTD;
        if plot_allTTDs == 1
            plot(T/24,cumTTD,'Color',[.6 .6 .6],'Parent',ax)
        end
    end
end
if plot_allTTDs == 1
    plot(T/24,C_m/CLIM_m,'b','LineWidth',2,'Parent',ax)
end
R_Q = C_Qnstat./C_Qmstat;
d_Q = (R_Q-par.r0)./par.r0*1000;


% FIGURES

% some common settings
dotstyle = '-k'; %symbol and color for the marker in all plots
col = [.6 .6 .6];
agemax = 365*5; %days

% PLOT with C_m and delta in batch VS mixing experiments
if plot_reactions == 1
    
    figure(1)
    set(gcf,'color','white');
    
    s1=subplot(3,1,1,'TickDir','out','box','on','NextPlot','add');
    title('Reaction for major isotope')
    % plot([1/Lambda_nf/24,1/Lambda_nf/24],[0 max(Clim)],...
    %     'Color',[.9 .9 .9],'HandleVisibility','off') %this is the timescale of the overall kinetic
    %plot([0 T(end)/24],[CLIM_m CLIM_m],'-k','DisplayName','CLIM')
    plot(T./24/365,C_m,'r-','LineWidth',1,'DisplayName','batch reaction');
    plot(mage/24/365,C_Qmstat,dotstyle,'MarkerSize',3,'Color',col,...
        'DisplayName','mixing');
    ylim([0 1.05*CLIM_m])
    xlim([0 1.05*agemax/365])
    %set(gca,'FontSize',17);
    xlabel('mean age [y]')
    ylabel('Concentration [\muM]')
    legend(gca,'Location','SE')
    
    s2=subplot(3,1,2);
    set(gca,'TickDir','out','box','on','NextPlot','add')
    title('Isotope ratio (\delta)')
    % plot([0 T(end)/24],[dLIM_af dLIM_af],'-k','DisplayName','\deltaLIM')
    plot(T./24/365,d_f,'r-','LineWidth',1,'DisplayName','batch reaction');
    plot(mage/24/365,d_Q,dotstyle,'MarkerSize',3,'Color',col,'DisplayName','mixing');
    %set(gca,'FontSize',17);
    xlabel('mean age [y]')
    ylabel(['\delta [',char(8240),']'])
    %ylim([-Inf 1.05*dLIM_af])
    xlim([0 1.05*agemax/365])
    %legend(gca,'Location','SE')
    
    s3=subplot(3,1,3);
    set(gca,'TickDir','out','box','on','NextPlot','add')
    title('\delta VS C')
    plot(C_m/CLIM_m,d_f,'r-','LineWidth',1,'DisplayName','batch reaction');
    plot(C_Qmstat/CLIM_m,d_Q,dotstyle,'MarkerSize',3,'Color',col,'DisplayName','mixing');
    xlim([0 1])
    %ylim([-Inf 1*dLIM_af])
    %set(gca,'FontSize',17);
    xlabel('C/CLIM [-]')
    ylabel(['\delta [',char(8240),']'])
    %legend(gca,'Location','NorthWest')
end

% PLOT with C_m and delta against Da number
if plot_Da == 1
    figure
    
    subplot(1,2,1)
    hold on
    set(gca,'TickDir','out','box','on','YLim',[0 1])
    title('Damkohler VS C')
    plot([1 1],[0 1],...
        'Color',[.2 .2 .2],'HandleVisibility','off') %to mark Da=1
    plot(Da,C_Qmstat/CLIM_m,dotstyle,'Color',col,'MarkerSize',3,'DisplayName','');
    xlim([0 10])
    xlabel('Da [-]')
    ylabel('C/CLIM [-]')
    
    subplot(1,2,2)
    hold on
    set(gca,'TickDir','out','box','on','YLim',[d_Q(1), d_Q(end)])
    title('Damkohler VS \delta')
    plot([1 1],[d_Q(1), d_Q(end)],...
        'Color',[.2 .2 .2],'HandleVisibility','off') %to mark Da=1
    plot(Da,d_Q,dotstyle,'Color',col,'MarkerSize',3,'DisplayName','');
    xlim([0 10])
    xlabel('Da [-]')
    ylabel(['\delta [',char(8240),']'])
    
end





