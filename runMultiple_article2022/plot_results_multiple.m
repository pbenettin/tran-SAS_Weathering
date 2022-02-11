% this is in a separate script just to shorten the script with the results
% analysis

% the idea here is to explore the variability of the C-Q relationships for
% various parameter combinations
% re-define the events
events={'22-Jan-2020','03-Feb-2020';...
    '22-Jun-2020','01-Jul-2020';...
    '03-Dec-2020','19-Dec-2020';...
    '16-Jan-2021','30-Jan-2021';...
    '31-Oct-2021','21-Nov-2021';...
    };

% prepare for the loops
Ne=length(event_selection); %number of selected events
Nq=sum(q); %number of selected simulations
hF1 = zeros(size(events,1),1); %preallocate the handles of the figures
hF2 = zeros(size(events,1),1); %preallocate the handles of the figures
show_spinup = 0;

% set some colors
% cols = repmat([...
%     0    0.4470    0.7410;
%     0.8500    0.3250    0.0980;
%     0.9290    0.6940    0.1250;
%     0.4940    0.1840    0.5560;
%     0.4660    0.6740    0.1880;
%     0.3010    0.7450    0.9330;
%     0.6350    0.0780    0.1840;...
%     ],2,1);
% 
cols = [    0.8500    0.3250    0.0980;
    0    0.4470    0.7410;
    0.9290    0.6940    0.1250];

% choose the young water fraction to show
ywt_index = 3; %this should correspond to a ywt ~200 d

% FIGURE WITH INDIVIDUAL EVENTS and VARIOUS SIMULATIONS

% external loop on the events
load(sprintf('%s/%s/output_%05d',mainfolder,foldertoopen,1)) % load a first simulation just to identify the events
for j = event_selection
    
    % PREPARE the figures and axes
    
    % figure with the events timeseries
    hF1(j) = figure;
    set(hF1(j),'Units','Normalized','Position',[.1,.1,.3,.6],'Visible',false)
    
    % axes for discharge
    s1=subplot(2,2,[1,2],'Parent',hF1(j));
    set(s1,'TickDir','out','NextPlot','add',...
        'XTick',[ ],'YTick',[ ])
    title(sprintf('Event %d: %s to %s',j,events{j,1},events{j,2}))
    
    % axes for rainfall
    s11=axes('Parent',hF1(j),'Position',get(s1,'Position'),'NextPlot','add',...
        'YAxisLocation','left','TickDir','out','XTick',[ ],...
        'YTick',[ ],'XTickLabel','none','Color','none',...
        'YDir','reverse');
    
    % axes with major ion concentration timeseries
    s2=axes('Parent',hF1(j),'Position',get(s1,'Position'),'NextPlot','add',...
        'YAxisLocation','left','TickDir','out','XTick',[ ],...
        'XTickLabel',{},'Color','none');
    ylabel('C [umol/l]')
    
    % axes with isotope ratio timeseries
    s22=axes('Parent',hF1(j),'Position',get(s1,'Position'),'NextPlot','add',...
        'YAxisLocation','right','TickDir','out','XTick',[ ],...
        'XTickLabel',{},'Color','none');
    ylabel(['\delta [',char(8240),']'])
    
    % axes with Cm-Q relationship
    s3=subplot(2,2,3,'Parent',hF1(j));
    set(s3,'TickDir','out','NextPlot','add','XScale','log','box','on')
    xlabel('Q [mm/h]')
    ylabel('C/CLIM [-]')
    title('C-Q major isotope');
    axis square
    
    % axes with delta-Q relationship
    s5=subplot(2,2,4,'Parent',hF1(j));
    set(s5,'TickDir','out','NextPlot','add','XScale','log','box','on')
    xlabel('Q [mm/h]')
    ylabel('\delta [-]')
    title('\delta-Q');
    axis square
    
    % figure with timeseries of water age fraction
    hF2(j) = figure;
    set(hF2(j),'Visible',false)
    
    % axes for streamflow
    s10 = axes('Parent',hF2(j));
    set(s10,'TickDir','out','NextPlot','add',...
        'XTick',[ ],'YTick',[ ])
    title(sprintf('Event %d: %s to %s',j,events{j,1},events{j,2}))
    
    % axes for rainfall
    s110=axes('Parent',hF2(j),'Position',get(s10,'Position'),'NextPlot','add',...
        'YAxisLocation','left','TickDir','out','XTick',[ ],...
        'YTick',[ ],'XTickLabel','none','Color','none',...
        'YDir','reverse');
    
    % axes with age
    s20=axes('Parent',hF2(j),'Position',get(s10,'Position'),'NextPlot','add',...
        'YAxisLocation','left','TickDir','out','XTick',[ ],...
        'XTickLabel',{},'Color','none','YLim',[0 1]);
    ylabel('Fyw [-]')
    
    % LOOP THROUGH THE SIMULATIONS
    
    % start a loop over the queried simulations
    sim_list = T.ID(q);
    for jj = 1:length(sim_list)
        
        % load a simulation
        load(sprintf('%s/%s/output_%05d',mainfolder,foldertoopen,sim_list(jj)))
        if show_spinup==0
            C_Qm(1:data.ini_shift) = NaN; %discard results during spinup
            d_f(1:data.ini_shift) = NaN; %discard results during spinup
            Fbp(1:data.ini_shift) = NaN; %discard results during spinup
            Fyw(1:data.ini_shift,:) = NaN; %discard results during spinup
        end
        
        % identify the events
        per=data.dates>=datenum(events{j,1}) &...
            data.dates<= datenum(events{j,2});
        strt_ind=find(per==1,1,'first'); end_ind=find(per==1,1,'last');
        x = datetime(datestr(data.dates(per))); % x-axis dates, common to all timeseries
        
        % FILL the FIGURE WITH C and Q TIMESERIES
        lwd = 1; %linewidth
        
        % plot once those things that are always the same
        if jj == 1
            % plot flow
            p1=area(x,data.Q(per),...
                'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8],...
                'DisplayName','discharge','Parent',s1);
            plot(datetime(datestr(data.dates(strt_ind))),data.Q(strt_ind),'s','MarkerFaceColor','r',...
                'MarkerEdgeColor','k','Parent',s1); %marker for the beginning
            plot(datetime(datestr(data.dates(end_ind))),data.Q(end_ind),'d','MarkerFaceColor','g',...
                'MarkerEdgeColor','k','Parent',s1); %marker for the end

            % add rainfall on top
            p11=bar(x,data.J(per),...
                'FaceColor','w','EdgeColor',[.5 .5 .5],...
                'DisplayName','rainfall','Parent',s11);
        end
        
        % plots with major and rare isotope concentration timeseries
        plot(x,C_Qm(per),'-','LineWidth',lwd,'Color',cols(jj,:),...
            'MarkerSize',3,'DisplayName','major isotope','Parent',s2);
        plot(x,d_f(per),'o:','LineWidth',lwd,'Color',cols(jj,:),...
            'MarkerSize',3,'DisplayName','isotope ratio','Parent',s22)
        
        % plot with Cm-Q relationship
        plot(data.Q(per),C_Qm(per)/CLIM_m,'-','LineWidth',lwd,'Color',cols(jj,:),'MarkerSize',3,...
            'DisplayName','Major isotope','Parent',s3)
        plot(data.Q(strt_ind),C_Qm(strt_ind)/CLIM_m,'s',...
            'MarkerFaceColor','r','MarkerEdgeColor','k',...
            'HandleVisibility','off','Parent',s3); %marker for the beginning
        plot(data.Q(end_ind),C_Qm(end_ind)/CLIM_m,'d','MarkerFaceColor','g',...
            'MarkerEdgeColor','k','HandleVisibility','off','Parent',s3); %marker for the end
        
        % plot with delta-Q relationship
        plot(data.Q(per),d_f(per),'o:','LineWidth',lwd,'Color',cols(jj,:),'MarkerSize',3,...
            'DisplayName','fractionating rare isotope','Parent',s5)
        plot(data.Q(strt_ind),d_f(strt_ind),'s','MarkerFaceColor','r',...
            'MarkerEdgeColor','k','HandleVisibility','off','Parent',s5); %marker for the beginning
        plot(data.Q(end_ind),d_f(end_ind),'d','MarkerFaceColor','g',...
            'MarkerEdgeColor','k','HandleVisibility','off','Parent',s5); %marker for the end
        
        % some final things in each figure
        linkaxes([s1,s2,s22,s11],'x');
        set(s1,'XLim',[min(x), max(x)])
        set(s2,'XTickLabels',{})
        set(s11,'XTickLabels',{})
        set(s22,'XTickLabels',{})
        set(s1,'Ylim',[0 1.1*max(data.Q)])
        set(s11,'YLim',[0,5*max(data.J(per))])
        set(s2,'YLim',[0 CLIM_m])
        set(s22,'YLim',[-0.12 dMAX])
        set(s3,'XLim',[0.9*min(data.Q) 1.1*max(data.Q)],'YLim',[0 1]);
        set(s5,'XLim',[0.9*min(data.Q) 1.1*max(data.Q)],'YLim',[-0.12 dMAX]);
        
        % FILL the FIGURE WITH AGE METRICS and Q TIMESERIES
        
        % plot once those things that are always the same
        if jj == 1
            % plot flow
            p1=area(x,data.Q(per),...
                'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8],...
                'DisplayName','discharge','Parent',s10);
            
            % add rainfall on top
            bar(x,data.J(per),...
                'FaceColor',[.7 .7 .7],'EdgeColor',[.6 .6 .6],...
                'DisplayName','rainfall','Parent',s110);
        end
        
        % plot with age stats
        plot(x,Fyw(per,ywt_index),'-','Color',cols(jj,:),...
            'MarkerSize',3,'DisplayName','major isotope','Parent',s20);
        
        % some final things in each figure
        linkaxes([s10,s20,s110],'x');
        set(s20,'XTickLabels',{})
        set(s110,'XTickLabels',{})
        set(s10,'Ylim',[0 1.1*max(data.Q)])
        set(s110,'YLim',[0,10*max(data.J(per))])
        
    end
end

% make figures visible
for handle_i = hF1'
    if handle_i>0
        set(handle_i,'Visible',true)
    end
end
for handle_i = hF2'
    if handle_i>0
        set(handle_i,'Visible',true)
    end
end



% FIGURE WITH ALL DATA for VARIOUS SIMULATIONS

% some settings
show_batch = 1; % show delta-Q results on top of a batch-reactor simulation?

% prepare the figures

% figure 1: C-delta relationships
f1 = figure;
ax1 = subplot(1,2,1,'Parent',f1,'TickDir','out','box','on','NextPlot','add');
title('\delta-C relationship')
xlim([0 1])
ylim([-0.2 1])
axis square
xlabel('C/CLIM [-]')
ylabel(['\delta [',char(8240),']'])
ax2 = subplot(1,2,2,'Parent',f1,'TickDir','out','box','on','NextPlot','add');
title('\delta-f_{diss} relationship')
xlim([0 1.1])
%ylim([-0.2 1.2])
ylim([-0.2 1])
axis square
xlabel('f_{diss} [-]')
ylabel(['\delta [',char(8240),']'])

% figure 2: age-C and age-delta relationships
f2 = figure;
s1 = subplot(1,3,2,'Parent',f2,'TickDir','out','NextPlot','add','Box','on');
xlabel('Fyw [-]')
ylabel('C/CLIM [-]')
title('Age-C relationship')
xlim([0 1])
ylim([0 1])
axis square
s2 = subplot(1,3,3,'Parent',f2,'TickDir','out','NextPlot','add','Box','on');
xlabel('Fyw [-]')
ylabel(['\delta [',char(8240),']'])
title('Age-\delta relationship')
xlim([0 1])
set(s2,'YLim',[-0.12 dMAX]);
axis square

% figure 3: age-Q relationships
%f3 = figure;
% ax3 = axes('Parent',f3,'TickDir','out','NextPlot','add','XScale','log','box','on');
ax3 = subplot(1,3,1,'Parent',f2,'TickDir','out','NextPlot','add','XScale','log','Box','on');
xlabel('Q [mm/h]')
ylabel('Fyw [-]')
title('AgeFraction-Q relationships')
%xlim([0.02 5])
ylim([0 1])
axis square
% ylim([-0.8 0])

% figure 4: SAVI data
if exist(fullfile('..','SAVIdata'),'dir') == 7
    run(fullfile('..','SAVIdata','plot_SAVIdata.m')) %(this will create axes sp1 and sp2)
else
    figure
    set(gcf,'color','white');

    sp1 = subplot(1,2,2);
    set(sp1,'NextPlot','add','TickDir','out',...
        'XLim',[0 1.2],'box','on')
    axis square
    xlabel('f_{diss} Si')
    ylabel(['\delta ^{30}Si [',char(8240),']'])
    
    sp2 = subplot(1,2,1);
    set(sp2,'NextPlot','add','TickDir','out','XScale','log',...
        'XLim',[0.002 15],'box','on')
    axis square
    xlabel('runoff [mm/h]')
    ylabel(['\delta ^{30}Si [',char(8240),']'])
end


% LOOP THROUGH THE SIMULATIONS and populate the figures

% start a loop over the queried simulations
sim_list = T.ID(q);
for jj = 1:length(sim_list)
    
    % load a simulation
    load(sprintf('%s/%s/output_%05d',mainfolder,foldertoopen,sim_list(jj)))
    if show_spinup==0
        C_Qm(1:data.ini_shift) = NaN; %discard results during spinup
        d_f(1:data.ini_shift) = NaN; %discard results during spinup
        Fbp(1:data.ini_shift) = NaN; %discard results during spinup
        Fyw(1:data.ini_shift,:) = NaN; %discard results during spinup
    end
    
    % figure 1: C-delta relationships
    if show_batch == 1 && jj == 1 %just need to plot this once
        
        % re-compute the batch reaction
        par.Clim = Pars([data.SASQl+data.SASETl+2, data.SASQl+data.SASETl+3]);
        par.lambda = Pars([data.SASQl+data.SASETl+4, data.SASQl+data.SASETl+5]);
        par.alph_f = [1, Pars(data.SASQl+data.SASETl+6)]; %when fractionation occurs for the second reaction
        par.rlim = Pars([data.SASQl+data.SASETl+7, data.SASQl+data.SASETl+8]);
        par.C0 = 0; %rainfall major isotope concentration
        par.r0=par.rlim(1); %rainfall isotope ratio
        par.dt=24*1; %time step [hours]
        N=200/par.lambda(2)/par.dt; %number of timesteps
        TT=(0:par.dt:(N-1)*par.dt)'; %age vector in dt timesteps
        
        % define parameters of the gamma distributions to be used
        a = [0.25, 0.5, 1, 5]; %'shape' parameter
        Ng1 = length(a); %number of shape parameter values to be used
        Ng2 = 100; %number of scale parameter values to be used
        Ng22 = 40;
        gm = logspace(log10(24*1),log10(24*5000),Ng2);
        
        % compute and plot the batch reaction
        tmp1 = c_m(TT,par);
        par.C_m = tmp1;
        tmp2 = data.frac_reaction(TT,par);
        plot(tmp1/CLIM_m,(tmp2./tmp1-par.r0)/par.r0*1000,'r-','LineWidth',1,...
            'DisplayName','batch reaction','Parent',ax1); %batch reaction
        
        % compute and plot the concentraton for static gamma distributions
        C_Qmstat = zeros(Ng2,Ng1,1);
        C_Qnstat = zeros(Ng2,Ng1,1);
        Da = zeros(Ng2,Ng1,1);
        mage = zeros(Ng2,Ng1,1);
        for i = 1:Ng1
            for j = 1:Ng2
                cumTTD=cdf('gamma',TT,a(i),gm(j)/a(i)); %matlab built-in function 'cdf'
                TTD = diff([0;cumTTD]);
                C_Qmstat(j,i) = tmp1'*TTD;
                C_Qnstat(j,i) = tmp2'*TTD;
            end
        end
        R_Q = C_Qnstat./C_Qmstat;
        d_Q = (R_Q-par.r0)./par.r0*1000;
        plot(C_Qmstat/CLIM_m,d_Q,'--','MarkerSize',3,'Color',[.3 .3 .3],...
            'DisplayName','mixing','HandleVisibility','off','Parent',ax1);
    end
    
    % plot all the new data
    plot(C_Qm/CLIM_m,d_f,'o','MarkerSize',2,'Color',cols(jj,:),...
        'HandleVisibility','off','Parent',ax1)
    
    if show_batch == 1 && jj == 1 %just need to plot this once
        
        % define a Rayleigh model here
        alphaR1 = 0.9983;
        alphaR2 = 0.9986;
        f_ray = linspace(0,1);
        d_ray1 = 0 - (1-f_ray).*(alphaR1-1).*1000;
        d_ray2 = (((par.r0.*f_ray.^(alphaR2-1))./par.r0)-1).*1000;
        d_ray3 = 0 - (1-f_ray).*(0.998-1).*1000;
        d_ray4 = (((par.r0.*f_ray.^(0.998-1))./par.r0)-1).*1000;
        
        % add it to the plot
        plot(f_ray,d_ray1,'k-.','LineWidth',1.5,'DisplayName','batch reaction 1','Parent',ax2); %batch reaction
        plot(f_ray,d_ray2,'k-','LineWidth',1.5,'DisplayName','batch reaction 2','Parent',ax2); %batch reaction
        plot(f_ray,d_ray3,'r-.','LineWidth',1.5,'Parent',ax2); %batch reaction
        plot(f_ray,d_ray4,'r-','LineWidth',1.5,'Parent',ax2); %batch reaction        
    end
    
    % plot the data
    f_diss = C_Qm./C_Qm_NP;
    plot(f_diss,d_f,'o','MarkerSize',2,'Color',cols(jj,:),'HandleVisibility','off','Parent',ax2)
        
    
    % figure 2: age-C-delta relationships
    lwd = .75;
    mksz = 1;
    plot(Fyw(:,ywt_index),C_Qm/CLIM_m,'o','LineWidth',lwd,'MarkerSize',mksz,'Color',cols(jj,:),...
        'HandleVisibility','off','Parent',s1)
    plot(Fyw(:,ywt_index),d_f,'o','LineWidth',lwd,'MarkerSize',mksz,'Color',cols(jj,:),...
        'HandleVisibility','off','Parent',s2)
    
    % figure 3: age-Q relationship
    plot(data.Q,Fyw(:,ywt_index),'o','LineWidth',lwd,'MarkerSize',mksz,'Color',cols(jj,:),...
        'HandleVisibility','off','Parent',ax3)
    
    % figure 4: SAVI data
    mrkr = 'o';
    mksiz = 4;
    per2 = true(size(C_Qm));
    f_diss = C_Qm./C_Qm_NP;
    pm1 = plot(f_diss(per2),d_f(per2),'o','MarkerSize',4,...
        'MarkerFaceColor',cols(jj,:),...
        'Marker',mrkr,...
        'MarkerSize',mksiz,...
        'MarkerEdgeColor','none',...
        'LineWidth',.25,...
        'DisplayName','model','Parent',sp1);
    %'Color',[.1 .1 .1],'MarkerFaceColor',[.4 .4 .4],'MarkerEdgeColor','k',...
    uistack(pm1,'bottom')
    pm2 = plot(data.Q(per2),d_f(per2),'o','MarkerSize',4,...
        'MarkerFaceColor',cols(jj,:),...
        'Marker',mrkr,...
        'MarkerSize',mksiz,...
        'MarkerEdgeColor','none',...
        'LineWidth',.25,...
        'DisplayName','model','Parent',sp2);
    uistack(pm2,'bottom')
    
    
end




%}