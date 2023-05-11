% this is no standalone script. It needs to be launched from the
% plot_results_isoweathering.m script

figure
sp1 = subplot(1,2,1);
set(gca,'TickDir','out','box','on','NextPlot','add')
axis square
title('\delta-C relationship')
xlim([0 1])
%set(gca,'FontSize',17);
xlabel('C/CLIM [-]')
ylabel(['\delta [',char(8240),']'])

if show_batch == 1
    
    % re-compute the batch reaction
    par.Clim = Pars([data.SASQl+data.SASETl+2, data.SASQl+data.SASETl+3]);
    par.lambda = Pars([data.SASQl+data.SASETl+4, data.SASQl+data.SASETl+5]);
    par.alph_f = [1, Pars(data.SASQl+data.SASETl+6)]; %when fractionation occurs for the second reaction
    par.rlim = Pars([data.SASQl+data.SASETl+7, data.SASQl+data.SASETl+8]);
    par.C0 = 0; %rainfall major isotope concentration
    par.r0=par.rlim(1); %rainfall isotope ratio
    par.dt=24*1; %time step [hours]
    N=200/par.lambda(2)/par.dt; %number of timesteps
    T=(0:par.dt:(N-1)*par.dt)'; %age vector in dt timesteps
    
    % define parameters of the gamma distributions to be used
    a = [0.25, 0.5, 1, 5]; %'shape' parameter
    Ng1 = length(a); %number of shape parameter values to be used
    Ng2 = 100; %number of scale parameter values to be used
    Ng22 = 40;
    gm = logspace(log10(24*1),log10(24*5000),Ng2);
    
    % compute and plot the batch reaction
    tmp1 = c_m(T,par);
    par.C_m = tmp1;
    tmp2 = data.frac_reaction(T,par);
    plot(tmp1/CLIM_m,(tmp2./tmp1-par.r0)/par.r0*1000,'r-','LineWidth',1,...
        'DisplayName','batch reaction'); %batch reaction
    
    % compute and plot the concentraton for static gamma distributions
    C_Qmstat = zeros(Ng2,Ng1,1);
    C_Qnstat = zeros(Ng2,Ng1,1);
    Da = zeros(Ng2,Ng1,1);
    mage = zeros(Ng2,Ng1,1);
    for i = 1:Ng1
        for j = 1:Ng2
            cumTTD=cdf('gamma',T,a(i),gm(j)/a(i)); %matlab built-in function 'cdf'
            TTD = diff([0;cumTTD]);
            C_Qmstat(j,i) = tmp1'*TTD;
            C_Qnstat(j,i) = tmp2'*TTD;
        end
    end
    R_Q = C_Qnstat./C_Qmstat;
    d_Q = (R_Q-par.r0)./par.r0*1000;
    plot(C_Qmstat/CLIM_m,d_Q,'--','MarkerSize',3,'Color',[.3 .3 .3],...
        'DisplayName','mixing','HandleVisibility','off');
    
end

% add all the new data
plot(C_Qm/CLIM_m,d_f,'o','MarkerSize',2,'Color',[.8 .8 .8],'HandleVisibility','off')

% identify and plot the events
if exist('event_sel','var')
    for i = event_sel
        per=data.dates>=datenum(events{i,1}) & data.dates<= datenum(events{i,2});
        plot(C_Qm(per)/CLIM_m,d_f(per),'-o','MarkerSize',3,...
            'Color',cols4ev(i,:),'DisplayName',sprintf('Event %d',i))
    end
end


%figure
sp2 = subplot(1,2,2);
set(gca,'TickDir','out','box','on','NextPlot','add','YLim',get(sp1,'YLim'))
axis square
title('\delta-f_{diss} relationship')
xlim([0 1])
xlabel('f_{diss} [-]')
ylabel(['\delta [',char(8240),']'])

if show_batch == 1
    % define a Rayleigh model here
    alphaR1 = 0.9983;
    alphaR2 = 0.9986;
    f_ray = linspace(0,1);
    d_ray1 = 0 - (1-f_ray).*(alphaR1-1).*1000;
    d_ray2 = (((par.r0.*f_ray.^(alphaR2-1))./par.r0)-1).*1000;
    d_ray3 = 0 - (1-f_ray).*(0.998-1).*1000;
    d_ray4 = (((par.r0.*f_ray.^(0.998-1))./par.r0)-1).*1000;
end

% plot the data
f_diss = C_Qm./C_Qm_NP;
plot(f_ray,d_ray1,'k-.','LineWidth',1.5,'DisplayName','batch reaction 1'); %batch reaction
plot(f_ray,d_ray2,'k-','LineWidth',1.5,'DisplayName','batch reaction 2'); %batch reaction
plot(f_ray,d_ray3,'r-.','LineWidth',1.5); %batch reaction
plot(f_ray,d_ray4,'r-','LineWidth',1.5); %batch reaction
plot(f_diss,d_f,'o','MarkerSize',2,'Color',[.8 .8 .8],'HandleVisibility','off')

ylim([-0.2 1.2])

% identify and plot the events
if exist('event_sel','var')
    for i = event_sel
        per=data.dates>=datenum(events{i,1}) & data.dates<= datenum(events{i,2});
        f_diss = C_Qm./C_Qm_NP;
        plot(f_diss(per),d_f(per),'-o','MarkerSize',3,...
            'Color',cols4ev(i,:),'DisplayName',sprintf('Event %d',i))
    end
end
%legend(gca,'Location','NW')
