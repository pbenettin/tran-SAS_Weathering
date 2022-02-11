% settings
use_addaxis = 0; %option to use the addaxis function
cols = lines(2);
mks = 3; %markersize
lwd = 1; %linewidth
xlm = [t_1,t_2]; %x-axis limits
%xlm = datenum({'01-Oct-2020','01-Oct-2021'}); %just one year

% figure with both discharge and C
figure

% plot 1: discharge
fCoutax1=axes('TickDir','out','XLim',xlm,'YLim',[0 1.1*max(data.Q)],...
    'NextPlot','add','YTick',[ ],'box','on');
p1=area(data.dates,data.Q,...
    'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8],...
    'DisplayName','discharge','Parent',fCoutax1);
%title('\bf full timeseries')
plot(1,1,'HandleVisibility','off')
%datetick(fCoutax1,'x','mmm-yy','keepticks')
datetick(fCoutax1,'x','mmm-yy','keeplimits')

% add rainfall on top
fCoutax2=axes('Position',get(fCoutax1,'Position'),'NextPlot','add',...
    'YAxisLocation','left','TickDir','out','XLim',xlm,'XTick',[ ],...
    'YTick',[ ],'XTickLabel','none','Color','none',...
    'YDir','reverse','YLim',[0,5*max(data.J)]);
p11=bar(data.dates,data.J,...
    'FaceColor','w','EdgeColor',[.5 .5 .5],...
    'DisplayName','rainfall','Parent',fCoutax2);


% add major isotope concentration and isotope ratio (delta value)
if use_addaxis == 1
    addaxis(data.dates,C_Qm,'-','Color',[0, 0.4470, 0.7410], 'MarkerSize',3,'DisplayName','major isotope')
    addaxis(data.dates,d_f,'-','Color',[0.8500    0.3250    0.0980],'DisplayName','\delta');
else
    ax2=axes('Position',get(fCoutax1,'Position'),'XLim',xlm,...
        'NextPlot','add','YAxisLocation','left','TickDir','out',...
        'Color','none','YColor',cols(1,:),'XTick',[ ]);
    p2=plot(data.dates,C_Qm,'-',...
        'MarkerSize',mks,'LineWidth',lwd,'Color',cols(1,:),...
        'DisplayName','major isotope','Parent',ax2);
    ylabel('C [umol/l]')
    set(ax2,'YLim',[0 CLIM_m]);
    
%     ax3=axes('Position',get(fCoutax1,'Position'),'XLim',xlm,'NextPlot','add',...
%         'YAxisLocation','right','TickDir','out','XTick',[ ],...
%         'Color','none','YColor',cols(2,:));
%     p3=plot(data.dates,d_f,'-',...
%         'MarkerSize',mks,'LineWidth',lwd,'Color',cols(2,:),...
%         'DisplayName','isotope ratio','Parent',ax3);
    ax3=axes('Position',get(fCoutax1,'Position'),'XLim',xlm,'NextPlot','add',...
        'YAxisLocation','right','TickDir','out','XTick',[ ],...
        'Color','none','YColor',cols(1,:));
    p3=plot(data.dates,d_f,':o',...
        'MarkerSize',mks,'LineWidth',lwd,'Color',cols(1,:),...
        'DisplayName','isotope ratio','Parent',ax3);
    ylabel(['\delta [',char(8240),']'])
    set(ax3,'YLim',[-0.12 dMAX]);
    %set(ax3,'YLim',[0.3 0.6]); 
    
    
    linkaxes([fCoutax1,fCoutax2,ax2,ax3],'x','keeplimits')
    legend(fCoutax1,[p11,p1,p2,p3],'Location','SE')
    
end