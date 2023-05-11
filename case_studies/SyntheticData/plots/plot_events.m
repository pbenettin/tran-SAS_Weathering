% this is no standalone script. It needs to be launched from the
% plot_results_isoweathering.m script

% some settings
ylim_d1 = [-0.12, dMAX];
ylim_d = [0.3, 0.6];
ylim_d1 = ylim_d1;
ylim_Cm = [0,1];
mks = 3; %markersize
lwd = 1; %linewidth

% make a figure for each event
for i = event_sel
    
    % identify the event period
    per=data.dates>=datenum(events{i,1}) & data.dates<= datenum(events{i,2});
    strt_ind=find(per==1,1,'first'); end_ind=find(per==1,1,'last');
    
    % if the timeseries figure still open, show the period
    if exist('fCoutax1','var')
        if isvalid(fCoutax1)
            aa = area([min(data.dates(per)),max(data.dates(per))],...
                [max(data.Q),max(data.Q)],...
                'FaceColor','none','EdgeColor','k','LineStyle',':',...
                'Parent',fCoutax1,'HandleVisibility','off');
            %uistack(aa,'bottom')
        end
    end
    
    %figure('Units','Centimeters','Position',[2,2,18,19])
    figure('Units','Normalized','Position',[.1,.1,.3,.6])
    
    % show C and Q timeseries
    s1=subplot(2,2,[1,2]);
    xtck = data.dates(strt_ind):3:data.dates(end_ind);
    set(s1,'TickDir','out','NextPlot','add',...
        'XTick',xtck,'YTick',[ ])
    %p1=plot(data.dates(per),data.Q(per),'-','Color',[.6 .6 .6],'DisplayName','discharge','Parent',s1);
    p1=area(data.dates(per),data.Q(per),...
        'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8],...
        'DisplayName','discharge','Parent',s1);
    plot(data.dates(strt_ind),data.Q(strt_ind),'s',...
        'MarkerFaceColor','r','MarkerEdgeColor','k'); %marker for the beginning
    plot(data.dates(end_ind),data.Q(end_ind),'d',...
        'MarkerFaceColor','g','MarkerEdgeColor','k'); %marker for the end
%     datetick('x','dd-mmm','keeplimits')
    ylim([0 1.1*max(data.Q)])
    %ylabel('Q mm/h')
    title(sprintf('Event %d: %s to %s',i,events{i,1},events{i,2}))
    
    % mark the time when TTDs were computed
    if show_TTDs == 1
        for j = 1:length(data.index_datesel)
            plot([data.dates(data.index_datesel(j)),data.dates(data.index_datesel(j))],...
                get(s1,'YLim'),'--','Color',[.6 .6 .6],'HandleVisibility','off');
            text(data.dates(data.index_datesel(j)),0.12,sprintf('%d',j),'Parent',s1)
        end
    end
    
    % add other data in the same figure on different axes
    if use_addaxis == 1
        addaxis(data.dates(per),C_Qm(per),'o-','Color',[0, 0.4470, 0.7410], 'MarkerSize',3,'DisplayName','concentration')
        addaxis(data.dates(per),d_f(per),'o-','Color',[0.8500, 0.3250, 0.0980],'MarkerSize',3,'DisplayName','\delta')
    else
        % add rainfall on top
        s11=axes('Position',get(s1,'Position'),'NextPlot','add',...
            'YAxisLocation','left','TickDir','out','XTick',[ ],...
            'YTick',[ ],'XTickLabel','none','Color','none',...
            'YDir','reverse','YLim',[0,5*max(data.J)]);
        p11=bar(data.dates(per),data.J(per),...
            'FaceColor','w','EdgeColor',[.5 .5 .5],...
            'DisplayName','rainfall','Parent',s11);
        
        % major ion
        s2=axes('Position',get(s1,'Position'),'NextPlot','add',...
            'YAxisLocation','left','TickDir','out','XTick',[ ],...
            'Color','none','YColor',cols(1,:));
        p2=plot(data.dates(per),C_Qm(per)/CLIM_m,'-',...
            'MarkerSize',mks,'LineWidth',lwd,'Color',cols(1,:),...
            'DisplayName','major isotope','Parent',s2);
        %ylabel('C [umol/l]')
        ylabel('C/CLIM [-]')
        set(s2,'YLim',ylim_Cm);
        
        % isotope ratio
%         s3=axes('Position',get(s1,'Position'),'NextPlot','add',...
%             'YAxisLocation','right','TickDir','out','XTick',[ ],...
%             'Color','none','YColor',cols(2,:));
%         p3=plot(data.dates(per),d_f(per),'o-',...
%             'MarkerSize',mks,'LineWidth',lwd,'Color',cols(2,:),...
%             'DisplayName','isotope ratio','Parent',s3);
        s3=axes('Position',get(s1,'Position'),'NextPlot','add',...
            'YAxisLocation','right','TickDir','out','XTick',[ ],...
            'Color','none','YColor',cols(1,:));
        p3=plot(data.dates(per),d_f(per),':o',...
            'MarkerSize',mks,'LineWidth',lwd,'Color',cols(1,:),...
            'DisplayName','isotope ratio','Parent',s3);

        ylabel(['\delta [',char(8240),']'])
        set(s3,'YLim',ylim_d); 
        
        xlim([data.dates(strt_ind) data.dates(end_ind)])
        linkaxes([s1,s11,s2,s3],'x')
        %datetick(s1,'x','dd-mmm','keeplimits')
        set(s1,'XTickLabel',datestr(xtck,'dd mmm'))
        legend(s1,[p11,p1,p2,p3],'Location','SE')
    end
    
    
    
    % show C_Qm vs Q
    s3=subplot(2,2,3);
    set(s3,'TickDir','out','NextPlot','add','XScale','log','box','on')
    plot(data.Q,C_Qm/CLIM_m,'o','MarkerSize',1,'Color',[.8 .8 .8],'HandleVisibility','off')
    plot(data.Q(per),C_Qm(per)/CLIM_m,'-','Color',cols(1,:),'LineWidth',lwd,'MarkerSize',3,'DisplayName','Major isotope')
    plot(data.Q(strt_ind),C_Qm(strt_ind)/CLIM_m,'s','MarkerFaceColor','r','MarkerEdgeColor','k','HandleVisibility','off'); %marker for the beginning
    plot(data.Q(end_ind),C_Qm(end_ind)/CLIM_m,'d','MarkerFaceColor','g','MarkerEdgeColor','k','HandleVisibility','off'); %marker for the end
    xlabel('Q [mm/h]')
    ylabel('C/CLIM [-]')
    %legend('show'); legend('boxoff')
    title('C-Q major isotope');
    axis square
    set(gca,'YLim',ylim_Cm);
    
    
    % show d_f vs Q
    s5=subplot(2,2,4);
    set(s5,'TickDir','out','NextPlot','add','XScale','log','box','on')
    plot(data.Q,d_f,'o','MarkerSize',1,'Color',[.8 .8 .8],'HandleVisibility','off')
    plot(data.Q(per),d_f(per),':o','LineWidth',lwd,'Color',cols(1,:),'MarkerSize',3,'DisplayName','fractionating rare isotope')
    plot(data.Q(strt_ind),d_f(strt_ind),'s','MarkerFaceColor','r','MarkerEdgeColor','k','HandleVisibility','off'); %marker for the beginning
    plot(data.Q(end_ind),d_f(end_ind),'d','MarkerFaceColor','g','MarkerEdgeColor','k','HandleVisibility','off'); %marker for the end
    xlabel('Q [mm/h]')
    ylabel(['\delta [',char(8240),']'])
    %legend('show'); legend('boxoff')
    title('\delta-Q');
    axis square
    set(gca,'YLim',ylim_d);
    
    
end
