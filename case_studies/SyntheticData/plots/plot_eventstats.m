% this is no standalone script. It needs to be launched from the
% plot_results_isoweathering.m script

% FIGURE WITH age-Q, C-Q and delta-Q relationships
figure

for jj = 1:length(ii_sel)
    
    ii = ii_sel(jj);
    
    % all age-Q relationships
    subplot(length(ii_sel),3,3*(jj-1)+1)
    set(gca,'TickDir','out','NextPlot','add','XScale','log','box','on')
    plot(data.Q,-Fyw(:,ii),'o','MarkerSize',2,'Color',[.8 .8 .8],'HandleVisibility','off')
    xlabel('Q [mm/h]')
    ylabel(sprintf('-Fyw [-] (ywt = %.0f d)',ywt(ii)))
    title('AgeFraction-Q relationships')
    xlim([0 5])
    %ylim([-0.9 0])
    
    % identify the events and plot them
    for i = event_sel
        per=data.dates>=datenum(events{i,1}) & data.dates<= datenum(events{i,2});
        plot(data.Q(per),-Fyw(per,ii),'-o','MarkerSize',3,...
            'Color',cols4ev(i,:),'DisplayName',sprintf('Event %d',i))
    end
    if jj == 1
        legend(gca,'Location','NW')
    end
    
    % all C-Q relationships
    subplot(length(ii_sel),3,3*(jj-1)+2)
    set(gca,'TickDir','out','NextPlot','add','XScale','log','box','on')
    plot(data.Q,C_Qm/CLIM_m,'o','MarkerSize',2,'Color',[.8 .8 .8],'HandleVisibility','off')
    xlabel('Q [mm/h]')
    ylabel('C/CLIM [-]')
    title('C-Q relationships')
    xlim([0 5])
    %ylim([-0.9 0])
    
    % identify the events and plot them
    for i = event_sel
        per=data.dates>=datenum(events{i,1}) & data.dates<= datenum(events{i,2});
        plot(data.Q(per),C_Qm(per)/CLIM_m,'-o','MarkerSize',3,...
            'Color',cols4ev(i,:),'DisplayName',sprintf('Event %d',i))
    end

    
    % all delta-Q relationships
    subplot(length(ii_sel),3,3*(jj-1)+3)
    set(gca,'TickDir','out','NextPlot','add','XScale','log','box','on')
    plot(data.Q,d_f,'o','MarkerSize',2,'Color',[.8 .8 .8],'HandleVisibility','off')
    xlabel('Q [mm/h]')
    ylabel(['\delta [',char(8240),']'])
    title('\delta-Q relationships')
    xlim([0 5])
    %ylim([-0.9 0])
    
    % identify the events and plot them
    for i = event_sel
        per=data.dates>=datenum(events{i,1}) & data.dates<= datenum(events{i,2});
        plot(data.Q(per),d_f(per),'-o','MarkerSize',3,...
            'Color',cols4ev(i,:),'DisplayName',sprintf('Event %d',i))
    end

    
end



% FIGURE WITH age-C and age-delta values
figure

% show them
for jj = 1:length(ii_sel)
    
    ii = ii_sel(jj);
    
    subplot(length(ii_sel),2,2*(jj-1)+1,'NextPlot','add','Box','on')

    plot(Fyw(:,ii),C_Qm/CLIM_m,'o','MarkerSize',2,'Color',[.8 .8 .8],'HandleVisibility','off')
    xlabel(sprintf('Fyw [-] (ywt = %.0f d)',ywt(ii)))
    ylabel('C/CLIM [-]')
    title('Age-C relationship')
    xlim([0 0.9])
    %ylim([0 1])
    
    % identify the events and plot them
    for i = event_sel
        per=data.dates>=datenum(events{i,1}) & data.dates<= datenum(events{i,2});
        plot(Fyw(per,ii),C_Qm(per)/CLIM_m,'-o','MarkerSize',3,...
            'Color',cols4ev(i,:),'DisplayName',sprintf('Event %d',i))
    end
    if jj ==1
        legend('show')
    end
    
    subplot(length(ii_sel),2,2*(jj-1)+2,'NextPlot','add','Box','on')
    plot(Fyw(:,ii),d_f,'o','MarkerSize',2,'Color',[.8 .8 .8],'HandleVisibility','off')
    xlabel(sprintf('Fyw [-] (ywt = %.0f d)',ywt(ii)))
    ylabel(['\delta [',char(8240),']'])
    title('Age-\delta relationship')
    xlim([0 0.9])
    ylim([0.3 0.6])
    
    % identify the events and plot them
    for i = event_sel
        per=data.dates>=datenum(events{i,1}) & data.dates<= datenum(events{i,2});
        plot(Fyw(per,ii),d_f(per),'-o','MarkerSize',3,...
            'Color',cols4ev(i,:),'DisplayName',sprintf('Event %d',i))
    end
end
