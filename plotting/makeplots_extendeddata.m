function [] = makeplots_extendeddata(plotnum);

clearvars -except plotnum;
close all;

% font preferences
fs = 7;
fn = 'helvetica-narrow';
fi = 'tex';
% plot colours
c1 = [0,0.447,0.741];
c2 = [0.85,0.325,0.098];
c3 = [0.929,0.694,0.125];
c4 = [0.494,0.184,0.556];
c5 = [0.466,0.674,0.188];
c6 = [0.301,0.745,0.933];
c7 = [0.635,0.078,0.184];
cols = [0.000,0.447,0.741;
        0.850,0.325,0.098;
        0.929,0.694,0.125;
        0.494,0.184,0.556;
        0.466,0.674,0.188;
        0.301,0.745,0.933;
        0.635,0.078,0.184];
% colour for plotting Q
qcol = c6;
% colour for plotting TF
tfcol = c7;
% colour for plotting both
bothcol = [0,0,0];

% melt rate parameterisation
ks = 0.142;
a = 0.31;
b = 1.19;

%% extended fig 1
if plotnum == 1,
    
    lspace = 0.01;
    hspace = 0.1;
    rspace = 0.02;
    bspace1 = 0.01;
    tspace1 = 0.01;
    bspace2 = 0.14;
    tspace2 = 0.05;
    ph1 = 1-bspace1-tspace1;
    ph2 = 1-bspace2-tspace2;
    pw1 = 0.2;
    pw2 = 1-lspace-rspace-hspace-pw1;
    ss = 10;
    figure();
    
    a1 = axes('position',[lspace,bspace1,pw1,ph1]);
    load ../../data/runoff/bedmachine450.mat
    terrain(terrain==0)=NaN;
    pcolor(x(1:ss:end),y(1:ss:end),terrain(1:ss:end,1:ss:end)); shading flat;
    xlim([min(x),max(x)]); ylim([min(y),max(y)]);
    colormap(a1,[0.6,0.3,0;1,1,0.8]);
    set(gca,'color','none','visible','off','box','off');
    
    a1b = axes('position',get(a1,'position'));    
    load ../data/drainagebasins.mat
    drainagebasins(drainagebasins~=1)=NaN;
    pcolor(x(1:ss:end),y(1:ss:end),drainagebasins(1:ss:end,1:ss:end)); shading flat;
    xlim([min(x),max(x)]); ylim([min(y),max(y)]);
    colormap(a1b,c1);
    set(gca,'box','on','xtick',[],'ytick',[],'color','none');
    xl = double(get(gca,'xlim')); yl = double(get(gca,'ylim'));
    tx = 0.05; ty = 0.95;
    text(xl(1)+tx*diff(xl),yl(1)+ty*diff(yl),'a','fontsize',fs+4,'fontname',fn,'interpreter',fi,'fontweight','bold');
    
    a2 = axes('position',[lspace+pw1+hspace,bspace2,pw2,ph2]); hold on;    
    load ../data/twglaciers.mat
    jid = 1;
    t = twglaciers(jid).runoff.RACMO.t;
    Q = twglaciers(jid).runoff.RACMO.Q;
    Qb = 0*Q+twglaciers(jid).basalmelt;
    plot(t,Q,'linewidth',0.5,'color',c1);
    plot(t,Qb,'linewidth',0.5,'color',c2);
    plot(t,Q+Qb,'linewidth',0.5,'color',c3);
    legend('surface runoff','basal melting','total','fontsize',fs,'fontname',fn,'interpreter',fi,'location','north','box','off','orientation','horizontal');
    xlim([1979 2018]); ylim([0 2500]);
    xlabel('year','fontsize',fs,'fontname',fn,'interpreter',fi);
    ylabel('monthly subglacial discharge Q (m^3/s)','fontsize',fs,'fontname',fn,'interpreter',fi);
    set(gca,'fontsize',fs,'fontname',fn,'ticklabelinterpreter',fi,'box','on');
    xl = double(get(gca,'xlim')); yl = double(get(gca,'ylim'));
    tx = 0.025; ty = 0.92;
    text(xl(1)+tx*diff(xl),yl(1)+ty*diff(yl),'b','fontsize',fs+4,'fontname',fn,'interpreter',fi,'fontweight','bold');
    
    fw = 17;
    fh = pw1*fw*(max(y)-min(y))/(ph1*(max(x)-min(x)));
    set(gcf,'renderer','painters');
    saveplot_pdf(fw,fh,300,'extfig1.pdf');
    close all;
    
end

%% extended fig 2
if plotnum == 2,
    
    lspace = 0.07;
    hspace = 0.1;
    rspace = 0.01;
    bspace = 0.04;
    vspace = 0.1;
    tspace = 0.04;
    pw = (1-lspace-rspace-hspace)/2;
    ph = (1-bspace-tspace-vspace)/2;
    
    glacid = 3;
    depths = [-500,-350];
    ms = 3;
    xlims = [1.5,5]*10^5;
    ylims = [-2.8,-2.55]*10^6;
    bm = load('../../data/runoff/bedmachine450.mat');
    bm.terrain(bm.terrain~=2)=NaN;
    load ../data/twglaciers.mat
    
    % get ORAS5 grid points and depths
    lon = ncread('~/Documents/ORAS5/votemper_ORAS5_1m_201812_grid_T_02.nc','nav_lon');
    lat = ncread('~/Documents/ORAS5/votemper_ORAS5_1m_201812_grid_T_02.nc','nav_lat');
    [x,y] = latlon2utm(lat,lon);
    oras5.x = x(:)';
    oras5.y = y(:)';
    T = ncread('~/Documents/ORAS5/votemper_ORAS5_1m_201812_grid_T_02.nc','votemper');
    Tdepth = reshape(T,size(T,1)*size(T,2),size(T,3));
    z = -ncread('~/Documents/ORAS5/votemper_ORAS5_1m_201812_grid_T_02.nc','deptht');
    for i=1:size(Tdepth,1),
        if ~isempty(find(~isnan(Tdepth(i,:)))),
            oras5.depth(i) = z(max(find(~isnan(Tdepth(i,:)))));
        else
            oras5.depth(i) = NaN;
        end
    end
    inds0 = find(~isnan(squeeze(Tdepth(:,1))));
    inds1 = find(oras5.depth<=depths(1));
    inds2 = find(oras5.depth<=depths(2));
    
    % contouring
    x_inds = find(bm.x<=xlims(2)&bm.x>=xlims(1));
    y_inds = find(bm.y<=ylims(2)&bm.y>=ylims(1));
    x0 = bm.x(x_inds);
    y0 = bm.y(y_inds);
    bed0 = bm.bed(y_inds,x_inds);
    bed1 = bed0; bed1(bed0<=depths(1))=1; bed1(bed0>depths(1))=0;
    bed2 = bed0; bed2(bed0<=depths(2))=1; bed2(bed0>depths(2))=0;
    mask0 = bm.terrain(y_inds,x_inds);
    cmap_bed = cptcmap('GMT_relief');
    
    figure();
    a1 = axes('position',[lspace,bspace+vspace+ph,pw,ph]); hold on;
    pcolor(x0,y0,bed0); shading flat;
    contour(x0,y0,bed1,[0,1],'r','linewidth',1);
    xlim(xlims); ylim(ylims); colormap(a1,cmap_bed); caxis([-1000 1000]);
    h = colorbar('southoutside');
    set(h,'position',[lspace,bspace+vspace+ph-0.12,pw/2,0.03]);
    set(h,'fontsize',fs,'fontname',fn,'ticklabelinterpreter',fi);
    xlabel(h,'bathymetry (m)','fontsize',fs,'fontname',fn,'interpreter',fi);
    set(gca,'box','off','color','none','visible','off');
    a1b = axes('position',get(a1,'position')); hold on;
    pcolor(x0,y0,mask0); shading flat;
    colormap(a1b,[1,1,0.8]);
    plot(twglaciers(glacid).x,twglaciers(glacid).y,'pk','markerfacecolor',c6,'markersize',ms+3);
    plot(oras5.x(inds0),oras5.y(inds0),'ko','markerfacecolor',c3,'markersize',ms);
    plot(oras5.x(inds1),oras5.y(inds1),'ko','markerfacecolor',c6,'markersize',ms);
    xlim(xlims); ylim(ylims);
    set(gca,'box','on','fontsize',fs,'fontname',fn,'ticklabelinterpreter',fi,'color','none');
    xlabel('x (km)','fontsize',fs,'fontname',fn,'interpreter',fi);
    ylabel('y (km)','fontsize',fs,'fontname',fn,'interpreter',fi);
    set(gca,'xtick',[xlims(1):1e5:xlims(2)],'xticklabel',{'0','100','200','300'});
    set(gca,'ytick',[ylims(1):1e5:ylims(2)],'yticklabel',{'0','100','200'});
    tx = -0.1; ty = 0.95;
    text(xlims(1)+tx*diff(xlims),ylims(1)+ty*diff(ylims),'a','fontsize',fs+4,'fontname',fn,'interpreter',fi,'fontweight','bold');
    
    % custom legend
    lx = 1.6e5;
    dx = 0.1e5;
    ly = -2.56e6;
    dy = 0.015e6;
    plot(lx,ly,'pk','markerfacecolor',c6,'markersize',ms+3);
    plot([lx-dx/2,lx+dx/2],[ly-dy,ly-dy],'r','linewidth',1);
    plot(lx,ly-2*dy,'ko','markerfacecolor',c3,'markersize',ms);
    plot(lx,ly-3*dy,'ko','markerfacecolor',c6,'markersize',ms);
    text(lx+dx,ly,'Helheim Glacier','fontsize',fs,'fontname',fn,'interpreter',fi,'verticalalignment','middle');
    text(lx+dx,ly-dy,'500 m contour','fontsize',fs,'fontname',fn,'interpreter',fi,'verticalalignment','middle');
    text(lx+dx,ly-2*dy,'ORAS5 point <500m','fontsize',fs,'fontname',fn,'interpreter',fi,'verticalalignment','middle');
    text(lx+dx,ly-3*dy,'ORAS5 point >500m','fontsize',fs,'fontname',fn,'interpreter',fi,'verticalalignment','middle');
    
    a2 = axes('position',[lspace+hspace+pw,bspace+vspace+ph,pw,ph]); hold on;
    pcolor(x0,y0,bed0); shading flat;
    contour(x0,y0,bed2,[0,1],'r','linewidth',1);
    xlim(xlims); ylim(ylims); colormap(a2,cmap_bed); caxis([-1000 1000]);
    set(gca,'box','off','color','none','visible','off');
    a2b = axes('position',get(a2,'position')); hold on;
    pcolor(x0,y0,mask0); shading flat;
    colormap(a2b,[1,1,0.8]);
    p(1) = plot(twglaciers(glacid).x,twglaciers(glacid).y,'pk','markerfacecolor',c6,'markersize',ms+3);
    p(2) = plot(oras5.x(inds0),oras5.y(inds0),'ko','markerfacecolor',c3,'markersize',ms);
    p(3) = plot(oras5.x(inds2),oras5.y(inds2),'ko','markerfacecolor',c6,'markersize',ms);
    plot(twglaciers(glacid).ocean.ORAS5.x,twglaciers(glacid).ocean.ORAS5.y,'kd','markerfacecolor',c6,'markersize',ms+3);
    plot(twglaciers(glacid).ocean.CHORE.x,twglaciers(glacid).ocean.CHORE.y,'kd','markerfacecolor',c5,'markersize',ms+3);
    plot(twglaciers(glacid).ocean.EN4.x,twglaciers(glacid).ocean.EN4.y,'kd','markerfacecolor',c7,'markersize',ms+3);
    plot(twglaciers(glacid).ocean.ASTE.x,twglaciers(glacid).ocean.ASTE.y,'kd','markerfacecolor',c3,'markersize',ms+3);
    xlim(xlims); ylim(ylims);
    set(gca,'box','on','fontsize',fs,'fontname',fn,'ticklabelinterpreter',fi,'color','none');
    xlabel('x (km)','fontsize',fs,'fontname',fn,'interpreter',fi);
    ylabel('y (km)','fontsize',fs,'fontname',fn,'interpreter',fi);
    set(gca,'xtick',[xlims(1):1e5:xlims(2)],'xticklabel',{'0','100','200','300'});
    set(gca,'ytick',[ylims(1):1e5:ylims(2)],'yticklabel',{'0','100','200'});
    tx = -0.1; ty = 0.95;
    text(xlims(1)+tx*diff(xlims),ylims(1)+ty*diff(ylims),'b','fontsize',fs+4,'fontname',fn,'interpreter',fi,'fontweight','bold');
    
    % custom legend
    plot(lx,ly,'pk','markerfacecolor',c6,'markersize',ms+3);
    plot([lx-dx/2,lx+dx/2],[ly-dy,ly-dy],'r','linewidth',1);
    plot(lx,ly-2*dy,'ko','markerfacecolor',c3,'markersize',ms);
    plot(lx,ly-3*dy,'ko','markerfacecolor',c6,'markersize',ms);
    plot(lx,ly-4*dy,'kd','markerfacecolor',c6,'markersize',ms+3);
    plot(lx,ly-5*dy,'kd','markerfacecolor',c5,'markersize',ms+3);
    plot(lx,ly-6*dy,'kd','markerfacecolor',c7,'markersize',ms+3);
    plot(lx,ly-7*dy,'kd','markerfacecolor',c3,'markersize',ms+3);
    text(lx+dx,ly,'Helheim Glacier','fontsize',fs,'fontname',fn,'interpreter',fi,'verticalalignment','middle');
    text(lx+dx,ly-dy,'350 m contour','fontsize',fs,'fontname',fn,'interpreter',fi,'verticalalignment','middle');
    text(lx+dx,ly-2*dy,'ORAS5 point <350m','fontsize',fs,'fontname',fn,'interpreter',fi,'verticalalignment','middle');
    text(lx+dx,ly-3*dy,'ORAS5 point >350m','fontsize',fs,'fontname',fn,'interpreter',fi,'verticalalignment','middle');
    text(lx+dx,ly-4*dy,'final ORAS5 point','fontsize',fs,'fontname',fn,'interpreter',fi,'verticalalignment','middle');
    text(lx+dx,ly-5*dy,'final CHORE point','fontsize',fs,'fontname',fn,'interpreter',fi,'verticalalignment','middle');
    text(lx+dx,ly-6*dy,'final EN4 point','fontsize',fs,'fontname',fn,'interpreter',fi,'verticalalignment','middle');
    text(lx+dx,ly-7*dy,'final ASTE point','fontsize',fs,'fontname',fn,'interpreter',fi,'verticalalignment','middle');
    
    a3 = axes('position',[lspace+pw/2+0.1,bspace+0.03,pw+hspace,ph-0.03]); hold on;
    plot(twglaciers(glacid).ocean.ORAS5.t,twglaciers(glacid).ocean.ORAS5.TF_GL,'linewidth',0.5,'color',c6);
    plot(twglaciers(glacid).ocean.CHORE.t,twglaciers(glacid).ocean.CHORE.TF_GL,'linewidth',0.5,'color',c5);
    plot(twglaciers(glacid).ocean.EN4.t,twglaciers(glacid).ocean.EN4.TF_GL,'linewidth',0.5,'color',c7);
    plot(twglaciers(glacid).ocean.ASTE.t,twglaciers(glacid).ocean.ASTE.TF_GL,'linewidth',0.5,'color',c3);
    plot(twglaciers(glacid).ocean.average.t,twglaciers(glacid).ocean.average.TF,'linewidth',1,'color','k');
    legend('ORAS5','CHORE\_RL','EN4','ASTE','mean','fontsize',fs,'fontname',fn,'interpreter',fi,'box','off','location','northwest');
    set(gca,'box','on','fontsize',fs,'fontname',fn,'ticklabelinterpreter',fi);
    xlabel('year','fontsize',fs,'fontname',fn,'interpreter',fi);
    ylabel(['monthly thermal forcing TF (',char(176),'C)'],'fontsize',fs,'fontname',fn,'interpreter',fi);
    xlim([1979,2018]); ylim([5,9]);
    tx = -0.1; ty = 1; xlims2 = get(gca,'xlim'); ylims2 = get(gca,'ylim');
    text(xlims2(1)+tx*diff(xlims2),ylims2(1)+ty*diff(ylims2),'c','fontsize',fs+4,'fontname',fn,'interpreter',fi,'fontweight','bold');
    
    fw = 17;
    fh = fw*pw*diff(ylims)/(ph*diff(xlims));
    saveplot_pdf(fw,fh,300,'extfig2.pdf');
    close all;
       
end

%% extended fig 4
if plotnum == 4,
    
    load ~/'OneDrive - University of Edinburgh'/greenland/write-up/draft3/modeloutput_methods_revised.mat
    cols = [0,0,0;c6;c7];
    tlims = [1979,2019];
    
    lspace = 0.055;
    hspace = 0.09;
    rspace = 0.01;
    bspace = 0.09;
    vspace = 0.09;
    tspace = 0.03;
    pw = (1-lspace-rspace-3*hspace)/4;
    ph = (1-bspace-tspace-vspace)/2;
    
    figure();
    
    % submarine melt rate
    axes('position',[lspace+0*(hspace+pw),bspace+1*(vspace+ph),pw,ph]); hold on;
    for jj=1:3,
        plot(eg_t,86400*eg_melt(jj,:),'color',cols(jj,:));
    end
    xlim(tlims);
    xlabel('year','fontname',fn,'fontsize',fs,'interpreter',fi);
    ylabel('submarine melt (m/d)','fontname',fn,'fontsize',fs,'interpreter',fi);
    set(gca,'box','on','fontname',fn,'fontsize',fs,'ticklabelinterpreter',fi);
    grid on;
    set(gca,'xtick',[1980:10:2010]); xtickangle(0);
    
    ylims = get(gca,'ylim');
    tdx = 0;
    tdy = -0.6;
    tddy = 0.08;
    text(tlims(1)+diff(tlims)*tdx,ylims(1)+(tdy-0*tddy)*diff(ylims),'atmosphere (Q) varies',...
        'fontsize',fs,'fontname',fn,'fontweight','bold','interpreter',fi,'color',c6,...
        'horizontalalignment','left');
    text(tlims(1)+diff(tlims)*tdx,ylims(1)+(tdy-1*tddy)*diff(ylims),'ocean (TF) varies',...
        'fontsize',fs,'fontname',fn,'fontweight','bold','interpreter',fi,'color',c7,...
        'horizontalalignment','left');
    text(tlims(1)+diff(tlims)*tdx,ylims(1)+(tdy-2*tddy)*diff(ylims),'both Q & TF vary',...
        'fontsize',fs,'fontname',fn,'fontweight','bold','interpreter',fi,'color','k',...
        'horizontalalignment','left');
    
    % ice flux
    axes('position',[lspace+1*(hspace+pw),bspace+1*(vspace+ph),pw,ph]); hold on;
    for jj=1:3,
        plot(eg_t,eg_iceflux(jj,:),'color',cols(jj,:),'linewidth',0.5);
        plot(eg_t,smooth(eg_iceflux(jj,:),200,'lowess'),'color',cols(jj,:),'linewidth',1);
    end
    xlim(tlims);
    ylabel('ice discharge (Gt/yr)','fontname',fn,'fontsize',fs,'interpreter',fi);
    set(gca,'box','on','fontname',fn,'fontsize',fs,'ticklabelinterpreter',fi);
    grid on;
    set(gca,'xtick',[1980:10:2010]); xtickangle(0);
    
    % normalised ice flux
    axes('position',[lspace+1*(hspace+pw),bspace+0*(vspace+ph),pw,ph]); hold on;
    for jj=1:3,
        D = smooth(eg_iceflux(jj,:),200,'lowess');
        plot(eg_t,D/D(1),'color',cols(jj,:),'linewidth',1);
    end
    xlim(tlims);
    xlabel('year','fontname',fn,'fontsize',fs,'interpreter',fi);
    ylabel('normalised ice discharge','fontname',fn,'fontsize',fs,'interpreter',fi);
    set(gca,'box','on','fontname',fn,'fontsize',fs,'ticklabelinterpreter',fi);
    grid on;
    set(gca,'xtick',[1980:10:2010]); xtickangle(0);
    
    % terminus position
    axes('position',[lspace+2*(hspace+pw),bspace+1*(vspace+ph),pw,ph]); hold on;
    for jj=1:3,
        plot(eg_t,eg_L(jj,:)/10^3,'color',cols(jj,:),'linewidth',0.5);
        plot(eg_t,smooth(eg_L(jj,:)/10^3,200,'lowess'),'color',cols(jj,:),'linewidth',1);
    end
    xlim(tlims);
    ylabel('terminus position (km)','fontname',fn,'fontsize',fs,'interpreter',fi);
    set(gca,'box','on','fontname',fn,'fontsize',fs,'ticklabelinterpreter',fi);
    grid on;
    set(gca,'xtick',[1980:10:2010]); xtickangle(0);
    
    % normalised terminus position
    axes('position',[lspace+2*(hspace+pw),bspace+0*(vspace+ph),pw,ph]); hold on;
    for jj=1:3,
        L = smooth(eg_L(jj,:)/10^3,200,'lowess');
        L = L-L(1);
        if jj==1, Lmax = L(end); end
        plot(eg_t,-L/Lmax,'color',cols(jj,:),'linewidth',1);
    end
    xlim(tlims); ylim([-1.2 0.2]);
    xlabel('year','fontname',fn,'fontsize',fs,'interpreter',fi);
    ylabel('normalised retreat','fontname',fn,'fontsize',fs,'interpreter',fi);
    set(gca,'box','on','fontname',fn,'fontsize',fs,'ticklabelinterpreter',fi);
    grid on;
    set(gca,'xtick',[1980:10:2010],'ytick',[-1:0.25:0]); xtickangle(0);
    
    % mass loss
    axes('position',[lspace+3*(hspace+pw),bspace+1*(vspace+ph),pw,ph]); hold on;
    for jj=1:3,
        plot(eg_t,cumtrapz(eg_t,eg_dVdt(jj,:)),'color',cols(jj,:),'linewidth',0.5);
        plot(eg_t,smooth(cumtrapz(eg_t,eg_dVdt(jj,:)),50),'color',cols(jj,:),'linewidth',1);
    end
    xlim(tlims);
    ylabel('cum. mass change (Gt)','fontname',fn,'fontsize',fs,'interpreter',fi);
    set(gca,'box','on','fontname',fn,'fontsize',fs,'ticklabelinterpreter',fi);
    grid on;
    set(gca,'xtick',[1980:10:2010]); xtickangle(0);
    
    % normalised mass loss
    axes('position',[lspace+3*(hspace+pw),bspace+0*(vspace+ph),pw,ph]); hold on;
    for jj=1:3,
        SL = smooth(cumtrapz(eg_t,eg_dVdt(jj,:)),50);
        if jj==1, SLmax = SL(end); end
        plot(eg_t,SL/SLmax,'color',cols(jj,:),'linewidth',1);
    end
    xlim(tlims); ylim([-0.1 1.1]);
    xlabel('year','fontname',fn,'fontsize',fs,'interpreter',fi);
    ylabel({'normalised dynamic';'sea level contribution'},'fontname',fn,'fontsize',fs,'interpreter',fi);
    set(gca,'box','on','fontname',fn,'fontsize',fs,'ticklabelinterpreter',fi);
    grid on;
    set(gca,'xtick',[1980:10:2010],'ytick',[0:0.25:1]); xtickangle(0);
    
    % subplot labels
    dtx = -0.005;
    dty = 0.305;
    letters1 = {'a','b','c','d'};
    letters2 = {'e','f','g'};
    for i=1:4,
        annotation('textbox','position',[lspace+(i-1)*(pw+hspace)+dtx,bspace+dty+ph+vspace,0.1,0.1],...
            'string',letters1{i},'fontsize',fs,'fontname',fn,'interpreter',fi,'fontweight','bold','edgecolor','none');
        if i~=4,
             annotation('textbox','position',[lspace+i*(pw+hspace)+dtx,bspace+dty,0.1,0.1],...
            'string',letters2{i},'fontsize',fs,'fontname',fn,'interpreter',fi,'fontweight','bold','edgecolor','none');
        end            
    end
    
    saveplot_pdf(17,10,300,'extfig4.pdf');
    close all;    
    
end

%% extended fig 5
if plotnum == 5,

load ../data/twglaciers.mat

Qmean = [twglaciers.Qmean]+[twglaciers.basalmelt];
TFmean = [twglaciers.TFmean];
zgl = [twglaciers.gldepth];
meanmelt = ks*Qmean.^a.*TFmean.^b;

% sector indices
inds(1).inds = find([twglaciers.ismip6_sector]<=3); % south
inds(2).inds = find([twglaciers.ismip6_sector]==4); % central-west
inds(3).inds = find([twglaciers.ismip6_sector]==6); % north-west
inds(4).inds = find([twglaciers.ismip6_sector]==5); % north-east
inds(5).inds = find([twglaciers.ismip6_sector]==7); % north-north

lspace = 0.07;
hspace = 0.05;
rspace = 0.02;
bspace = 0.12;
tspace = 0.1;
pw = (1-lspace-rspace-2*hspace)/3;
ph = 1-bspace-tspace;
ms = 8;

% Q vs gl depth
a1 = axes('position',[lspace+0*(pw+hspace),bspace,pw,ph]); hold on;
for l=1:5,
    sectorinds = inds(l).inds;
    xdata = Qmean(sectorinds); ydata = zgl(sectorinds);
    p = polyfit(ydata,log10(xdata),1);
    [r,pval] = corrcoef(ydata',xdata');
    pval = pval(1,2);    
    scatter(xdata,ydata,ms,'marker','s','markerfacecolor',cols(l,:),'markeredgecolor','none','markerfacealpha',0.8);
    if pval<0.05,
        plot(10.^(p(1)*[min(ydata),max(ydata)]+p(2)),[min(ydata),max(ydata)],'--','color',cols(l,:),'linewidth',0.5);
    end
end
xlabel('subglacial discharge (m^3/s)','fontsize',fs,'fontname',fn,'interpreter',fi);
ylabel('grounding line depth (m)','fontsize',fs,'fontname',fn,'interpreter',fi);
set(gca,'box','on','fontsize',fs,'fontname',fn,'ticklabelinterpreter',fi,'xscale','log');
xlim([2 250]); ylim([-900 0]);

% TF vs gl depth
a2 = axes('position',[lspace+1*(pw+hspace),bspace,pw,ph]); hold on;
for l=1:5,
    sectorinds = inds(l).inds;
    xdata = TFmean(sectorinds); ydata = zgl(sectorinds);
    p = polyfit(ydata,xdata,1);
    [r,pval] = corrcoef(ydata',xdata');
    pval = pval(1,2);
    scatter(xdata,ydata,ms,'marker','s','markerfacecolor',cols(l,:),'markeredgecolor','none','markerfacealpha',0.8);
    if pval<0.05,
        plot(p(1)*[min(ydata),max(ydata)]+p(2),[min(ydata),max(ydata)],'--','color',cols(l,:),'linewidth',0.5);
    end
end
xlabel(['thermal forcing (',char(176),'C)'],'fontsize',fs,'fontname',fn,'interpreter',fi);
set(gca,'box','on','fontsize',fs,'fontname',fn,'ticklabelinterpreter',fi,'yticklabel',[]);
xlim([0 8]); ylim([-900 0]);

% melt vs gl depth
a3 = axes('position',[lspace+2*(pw+hspace),bspace,pw,ph]); hold on;
for l=1:5,
    sectorinds = inds(l).inds;
    xdata = meanmelt(sectorinds); ydata = zgl(sectorinds);
    p = polyfit(ydata,xdata,1);
    [r,pval] = corrcoef(ydata',xdata');
    pval = pval(1,2);
    scatter(xdata,ydata,ms,'marker','s','markerfacecolor',cols(l,:),'markeredgecolor','none','markerfacealpha',0.8);
    if pval<0.05,
        plot(p(1)*[min(ydata),max(ydata)]+p(2),[min(ydata),max(ydata)],'--','color',cols(l,:),'linewidth',0.5);
    end
end
xlabel('submarine melt rate (m/d)','fontsize',fs,'fontname',fn,'interpreter',fi);
set(gca,'box','on','fontsize',fs,'fontname',fn,'ticklabelinterpreter',fi,'yticklabel',[]);
xlim([0 5]); ylim([-900 0]);

% custom legend
regions_name = {'South','Central-West','North-West','North-East','North'};
tx = [0.42,0.53,0.67,0.8,0.9];
ty = 0.9;
for l=1:5,
    t = annotation('textbox');
    t.Position = [tx(l),ty,0.1,0.1];
    t.String = regions_name{l};
    t.FontName = fn;
    t.FontSize = fs+2;
    t.FontWeight = 'bold';
    t.Interpreter = fi;
    t.Color = cols(l,:);
    t.EdgeColor = 'none';
    t.HorizontalAlignment = 'center';
    t.VerticalAlignment = 'middle';
end


% subplot labels
lbs = {'a','b','c'};
tx = 0.27;
ty = 0.81;
for l=1:3,
    t = annotation('textbox');
    t.Position = [tx+(l-1)*(hspace+pw),ty,0.1,0.1];
    t.String = lbs{l};
    t.FontName = fn;
    t.FontSize = fs+2;
    t.FontWeight = 'bold';
    t.Interpreter = fi;
    t.Color = 'k';
    t.EdgeColor = 'none';
    t.HorizontalAlignment = 'center';
    t.VerticalAlignment = 'middle';
end

set(gcf,'renderer','painters');
saveplot_pdf(17,8,300,'extfig5.pdf');
close all;

end

%% extended fig 6
if plotnum == 6,
    
    % load data
    load ../data/twglaciers.mat
    
    % time periods
    % normalisation period
    t0 = [1979,1993];
    t = twglaciers(1).submelt.t;
    inds0 = find(t>=t0(1) & t<=t0(end));
    % smoothing period (months)
    smt = 10*12;
    % smoothing method
    smm = 'lowess';
    
    % put required quantities into arrays
    prods = {'average'}; % can use this to put products into shading too
    i0 = 1;
    for i=1:length(twglaciers),
        for j=1:length(prods),
            
            % record sector
            secs(i0) = twglaciers(i).ismip6_sector;
        
            % Q
            Qa = interp1(twglaciers(i).runoff.RACMO.t,twglaciers(i).runoff.RACMO.Q,t)+twglaciers(i).basalmelt;
            Q = smooth(Qa,smt,smm);
            Q0 = mean(Q(inds0));
            Qnorm(i0,:) = 100*(Q/Q0-1);

            % TF
            TFa = interp1(twglaciers(i).ocean.(prods{j}).t,twglaciers(i).ocean.(prods{j}).TF,t);
            TFa(TFa<0)=0;
            TF = smooth(TFa,smt,smm);
            TF0 = mean(TF(inds0));
            TFnorm(i0,:) = 100*(TF/TF0-1);

            % melt rate
            ma = ks*Qa.^a.*TFa.^b;
            m = smooth(ma,smt,smm);
            m0 = mean(m(inds0));
            mnorm(i0,:) = 100*(m/m0-1);

            % melt rate Q varies
            % make TF0 a repeating version of the mean cycle during t0
            TFt = TFa(inds0);
            TF0 = NaN*TF';
            for i=1:12, TF0([i:12:end]) = mean(TFt([i:12:end])); end
            mQ = ks*Qa.^a.*TF0.^b;
            mQ = smooth(mQ,smt,smm);
            mQ0 = mean(mQ(inds0));
            mQnorm(i0,:) = 100*(mQ/mQ0-1);

            % melt rate TF varies
            % make Q0 a repeating version of the mean cycle during t0
            Qt = Qa(inds0);
            Q0 = NaN*Q';
            for i=1:12, Q0([i:12:end]) = mean(Qt([i:12:end])); end
            mTF = ks*Q0.^a.*TFa.^b;
            mTF = smooth(mTF,smt,smm);
            mTF0 = mean(mTF(inds0));
            mTFnorm(i0,:) = 100*(mTF/mTF0-1);
        
            i0=i0+1;
        end
    end
    
    % get sector indices    
    inds(1).inds = find([twglaciers.ismip6_sector]==5); inds(1).name = 'North-east (NE)';
    inds(2).inds = find([twglaciers.ismip6_sector]==7); inds(2).name = 'North (NO)';
    
    % plots
    lspace = 0.15;
    hspace = 0.05;
    rspace = 0.015;
    bspace = 0.085;
    vspace = 0.1;
    tspace = 0.04;
    pw = (1-lspace-rspace-hspace)/2;
    ph = (1-bspace-tspace-vspace)/2;
    
    ylims1 = [-20,150];
    ylims2 = [-25,130];
    tlims = [1979,2018];
    lw1 = 0.5;
    lw2 = 1.5;
    fa = 0.3;
    tdx1 = 0.05;
    tdy1 = 0.92;
    tdx2 = 0.05;
    tdy2 = 0.92;
    tddy = 0.1;
    lx1 = -0.35/5;
    lx2 = -0.3/5;
        
    figure();
    
    % forcings
    for i=1:2,
        
        a(i) = axes('position',[lspace+(i-1)*(hspace+pw),bspace,pw,ph]); hold on;
        plot(t,0*t,'k--','linewidth',lw1);
        patch([t,fliplr(t)],[prctile(Qnorm(inds(i).inds,:),25),fliplr(prctile(Qnorm(inds(i).inds,:),75))],...
            qcol,'facealpha',fa,'edgecolor','none');
        plot(t,prctile(Qnorm(inds(i).inds,:),50),'color',qcol,'linewidth',lw2);
        patch([t,fliplr(t)],[prctile(TFnorm(inds(i).inds,:),25),fliplr(prctile(TFnorm(inds(i).inds,:),75))],...
            tfcol,'facealpha',fa,'edgecolor','none');
        plot(t,prctile(TFnorm(inds(i).inds,:),50),'color',tfcol,'linewidth',lw2);
        xlim(tlims); ylim(ylims1);
        set(a(i),'box','on','ytick',[-40:40:120]); grid on;
        set(a(i),'fontsize',fs,'fontname',fn,'ticklabelinterpreter',fi);
        xlabel('year','fontsize',fs,'fontname',fn,'interpreter',fi);
        if i==1,
            text(tlims(1)+diff(tlims)*tdx1,ylims1(1)+(tdy1-0*tddy)*diff(ylims1),'subglacial discharge Q',...
                'fontsize',fs,'fontname',fn,'fontweight','bold','interpreter',fi,'color',qcol,'horizontalalignment','left');
            text(tlims(1)+diff(tlims)*tdx1,ylims1(1)+(tdy1-1*tddy)*diff(ylims1),'thermal forcing TF',...
                'fontsize',fs,'fontname',fn,'fontweight','bold','interpreter',fi,'color',tfcol,'horizontalalignment','left');
            ylabel({'% change in forcing';['versus ',num2str(t0(1)),'-',num2str(t0(end))]},'fontsize',fs,'fontname',fn,'interpreter',fi);
            set(gca,'ytick',[-40:40:120]);
        else
            set(gca,'ytick',[-40:40:120],'yticklabel',[]);
        end
        
    end
    
    % melt rates
    for i=1:2,
        
        a(i) = axes('position',[lspace+(i-1)*(hspace+pw),bspace+vspace+ph,pw,ph]); hold on;
        plot(t,0*t,'k--','linewidth',lw1);
        % Q
        patch([t,fliplr(t)],[prctile(mQnorm(inds(i).inds,:),25),fliplr(prctile(mQnorm(inds(i).inds,:),75))],...
            qcol,'facealpha',fa,'edgecolor','none');
        plot(t,prctile(mQnorm(inds(i).inds,:),50),'color',qcol,'linewidth',lw2);
        % TF
        patch([t,fliplr(t)],[prctile(mTFnorm(inds(i).inds,:),25),fliplr(prctile(mTFnorm(inds(i).inds,:),75))],...
            tfcol,'facealpha',fa,'edgecolor','none');
        plot(t,prctile(mTFnorm(inds(i).inds,:),50),'color',tfcol,'linewidth',lw2);
        % both
        patch([t,fliplr(t)],[prctile(mnorm(inds(i).inds,:),25),fliplr(prctile(mnorm(inds(i).inds,:),75))],...
            bothcol,'facealpha',fa,'edgecolor','none');
        plot(t,prctile(mnorm(inds(i).inds,:),50),'color',bothcol,'linewidth',lw2);
        
        xlim(tlims); ylim(ylims2);
        set(a(i),'box','on'); grid on;
        set(a(i),'fontsize',fs,'fontname',fn,'ticklabelinterpreter',fi);
        title(inds(i).name,'fontsize',fs,'fontname',fn,'interpreter',fi);
        if i==1,
            text(tlims(1)+diff(tlims)*tdx2,ylims2(1)+(tdy2-0*tddy)*diff(ylims2),'atmosphere (Q) varies',...
                'fontsize',fs,'fontname',fn,'fontweight','bold','interpreter',fi,'color',qcol,...
                'horizontalalignment','left');
            text(tlims(1)+diff(tlims)*tdx2,ylims2(1)+(tdy2-1*tddy)*diff(ylims2),'ocean (TF) varies',...
                'fontsize',fs,'fontname',fn,'fontweight','bold','interpreter',fi,'color',tfcol,...
                'horizontalalignment','left');
            text(tlims(1)+diff(tlims)*tdx2,ylims2(1)+(tdy2-2*tddy)*diff(ylims2),'both Q & TF vary',...
                'fontsize',fs,'fontname',fn,'fontweight','bold','interpreter',fi,'color',bothcol,...
                'horizontalalignment','left');
            ylabel({'% change in submarine melt rate';['versus ',num2str(t0(1)),'-',num2str(t0(end))]},'fontsize',fs,'fontname',fn,'interpreter',fi);
            set(gca,'ytick',[-40:40:120]);
        else
            set(gca,'ytick',[-40:40:120],'yticklabel',[]);
        end
        
    end
    
    % subplot labels
    dtx = -0.05;
    dty = 0.32;
    letters1 = {'a','b'};
    letters2 = {'c','d'};
    for i=1:2,
        annotation('textbox','position',[lspace+(i-1)*(pw+hspace)+dtx,bspace+1*(ph+vspace)+dty,0.1,0.1],...
            'string',letters1{i},'fontsize',fs,'fontname',fn,'interpreter',fi,'fontweight','bold','edgecolor','none');
        annotation('textbox','position',[lspace+(i-1)*(pw+hspace)+dtx,bspace+0*(ph+vspace)+dty,0.1,0.1],...
            'string',letters2{i},'fontsize',fs,'fontname',fn,'interpreter',fi,'fontweight','bold','edgecolor','none');
    end
    
    saveplot_pdf(9,10,300,'extfig6.pdf');
    close all;
    
end

%% extended fig 7
if plotnum == 7,
    
% run model
% [t,modeloutput] = glacier_response_everyglacier(ks,a,b);
load ../model/modeloutput.mat

lspace = 0.14;
hspace = 0.05;
rspace = 0.01;
bspace = 0.16;
tspace = 0.075;
pw = (1-lspace-rspace-hspace)/2;
ph = 1-tspace-bspace;
ylims = [-0.1,1.1];
tlims = [1979,2019];
tdx = 0.025;
tdy = 0.95;
tddy = 0.08;
fa = 0.3;
lw = 1.5;

casecols = [bothcol;qcol;tfcol];
inds = [4,5];
modeloutput(4).name = 'North-east (NE)';
modeloutput(5).name = 'North (NO)';

figure();
for i=1:2,
    a(i) = axes('position',[lspace+(i-1)*(pw+hspace),bspace,pw,ph]); hold on;
    for jj=3:-1:1,
        patch([t,fliplr(t)],[modeloutput(inds(i)).case(jj).SL25,...
            fliplr(modeloutput(inds(i)).case(jj).SL75)],...
            casecols(jj,:),'facealpha',fa,'edgecolor','none');
        plot(t,modeloutput(inds(i)).case(jj).SL50,'color',casecols(jj,:),'linewidth',lw);
    end
    set(a(i),'fontsize',fs,'fontname',fn,'ticklabelinterpreter',fi,'box','on','xtick',[1980:10:2020]); grid on;
    ylim(ylims); xlim(tlims);
    title(modeloutput(inds(i)).name,'fontsize',fs,'fontname',fn,'interpreter',fi);
    xlabel('year','fontsize',fs,'fontname',fn,'interpreter',fi);
    if i==1,
        ylabel({'normalised dynamic';'sea level contribution'},'fontsize',fs,'fontname',fn,'interpreter',fi);
        set(a(i),'ytick',[0:0.2:1]);
        text(tlims(1)+diff(tlims)*tdx,ylims(1)+(tdy-0*tddy)*diff(ylims),'atmosphere (Q) varies',...
            'fontsize',fs,'fontname',fn,'fontweight','bold','interpreter',fi,'color',qcol,...
            'horizontalalignment','left');
        text(tlims(1)+diff(tlims)*tdx,ylims(1)+(tdy-1*tddy)*diff(ylims),'ocean (TF) varies',...
            'fontsize',fs,'fontname',fn,'fontweight','bold','interpreter',fi,'color',tfcol,...
            'horizontalalignment','left');
        text(tlims(1)+diff(tlims)*tdx,ylims(1)+(tdy-2*tddy)*diff(ylims),'both Q & TF vary',...
            'fontsize',fs,'fontname',fn,'fontweight','bold','interpreter',fi,'color',bothcol,...
            'horizontalalignment','left');
    else
        set(a(i),'ytick',[0:0.2:1],'yticklabel',[]);
    end

end

% subplot labels
dtx = -0.05;
dty = 0.73;
letters = {'a','b'};
for i=1:2,
    annotation('textbox','position',[lspace+(i-1)*(pw+hspace)+dtx,bspace+dty,0.1,0.1],...
        'string',letters{i},'fontsize',fs,'fontname',fn,'interpreter',fi,'fontweight','bold','edgecolor','none');
end

saveplot_pdf(9,5,300,'extfig7.pdf');
close all;
    
end

%% extended fig 8
if plotnum == 8,
    
% run model
load ../model/modeloutput.mat

lspace = 0.09;
hspace = 0.04;
vspace = 0.09;
rspace = 0.01;
bspace = 0.09;
tspace = 0.04;
pw = (1-lspace-rspace-2*hspace)/3;
ph = (1-tspace-bspace-vspace)/2;
ylims = [-1.4,0.4];
tlims = [1979,2019];
tdx = 0.025;
tdy = 0.22;
tddy = 0.08;
fa = 0.3;
lw = 1.5;

% add abbreviations to name
modeloutput(1).name = 'South (SO)';
modeloutput(2).name = 'Central-west (CW)';
modeloutput(3).name = 'North-west (NW)';
modeloutput(4).name = 'North-east (NE)';
modeloutput(5).name = 'North (NO)';

casecols = [bothcol;qcol;tfcol];
inds = [1,2,3,4,5,6];
a(1) = axes('position',[lspace+0*(pw+hspace),bspace+1*(ph+vspace),pw,ph]);
a(2) = axes('position',[lspace+1*(pw+hspace),bspace+1*(ph+vspace),pw,ph]);
a(3) = axes('position',[lspace+2*(pw+hspace),bspace+1*(ph+vspace),pw,ph]);
a(4) = axes('position',[lspace+0*(pw+hspace),bspace+0*(ph+vspace),pw,ph]);
a(5) = axes('position',[lspace+1*(pw+hspace),bspace+0*(ph+vspace),pw,ph]);
a(6) = axes('position',[lspace+2*(pw+hspace),bspace+0*(ph+vspace),pw,ph]);

figure();
for i=1:6,
    axes(a(i)); hold on;
    for jj=3:-1:1,
        patch([t,fliplr(t)],[modeloutput(inds(i)).case(jj).L25,...
            fliplr(modeloutput(inds(i)).case(jj).L75)],...
            casecols(jj,:),'facealpha',fa,'edgecolor','none');
        plot(t,modeloutput(inds(i)).case(jj).L50,'color',casecols(jj,:),'linewidth',lw);
    end
    set(a(i),'fontsize',fs,'fontname',fn,'ticklabelinterpreter',fi,'box','on','xtick',[1980:10:2020]); grid on;
    ylim(ylims); xlim(tlims);
    title(modeloutput(inds(i)).name,'fontsize',fs,'fontname',fn,'interpreter',fi);
    if i==1 | i==4,
        ylabel('normalised retreat','fontsize',fs,'fontname',fn,'interpreter',fi);
        set(a(i),'ytick',[-1.2:0.2:0.2]);
    else
        set(a(i),'ytick',[-1.2:0.2:0.2],'yticklabel',[]);
    end
    if i==4,
        text(tlims(1)+diff(tlims)*tdx,ylims(1)+(tdy-2*tddy)*diff(ylims),'observed',...
            'fontsize',fs,'fontname',fn,'fontweight','bold','interpreter',fi,'color',c3,...
            'horizontalalignment','left');
        text(tlims(1)+diff(tlims)*tdx,ylims(1)+(tdy-1*tddy)*diff(ylims),'atmosphere (Q) varies',...
            'fontsize',fs,'fontname',fn,'fontweight','bold','interpreter',fi,'color',qcol,...
            'horizontalalignment','left');
        text(tlims(1)+diff(tlims)*tdx,ylims(1)+(tdy-0*tddy)*diff(ylims),'ocean (TF) varies',...
            'fontsize',fs,'fontname',fn,'fontweight','bold','interpreter',fi,'color',tfcol,...
            'horizontalalignment','left');
        text(tlims(1)+diff(tlims)*tdx,ylims(1)+(tdy+1*tddy)*diff(ylims),'both Q & TF vary',...
            'fontsize',fs,'fontname',fn,'fontweight','bold','interpreter',fi,'color',bothcol,...
            'horizontalalignment','left');
    end
    if i>3,
        xlabel('year','fontsize',fs,'fontname',fn,'interpreter',fi);
    end

end

% add king 2020 observed retreat
load ../../data/king_termpos/glaciers.mat
load ../data/ice_ocean_sectors.mat
tL_king = [1984:2019];
sec(1).inds = find(inpolygon([glaciers.x],[glaciers.y],regions(1).ocean.x,regions(1).ocean.y) | ...
                   inpolygon([glaciers.x],[glaciers.y],regions(2).ocean.x,regions(2).ocean.y) | ...
                   inpolygon([glaciers.x],[glaciers.y],regions(3).ocean.x,regions(3).ocean.y));
sec(2).inds = find(inpolygon([glaciers.x],[glaciers.y],regions(4).ocean.x,regions(4).ocean.y));
sec(3).inds = find(inpolygon([glaciers.x],[glaciers.y],regions(6).ocean.x,regions(6).ocean.y));
sec(4).inds = find(inpolygon([glaciers.x],[glaciers.y],regions(5).ocean.x,regions(5).ocean.y));
sec(5).inds = find(inpolygon([glaciers.x],[glaciers.y],regions(7).ocean.x,regions(7).ocean.y));
sec(6).inds = [1:length(glaciers)];

for j=1:6,
    i0 = 1;
    L_king =  [];
    for i=1:length(sec(j).inds),
        if ~isempty(glaciers(sec(j).inds(i)).termpos),
            L_king(i0,:) = interp1(glaciers(sec(j).inds(i)).termpos.t,glaciers(sec(j).inds(i)).termpos.L(:,1),tL_king);
            L_king(i0,:) = smooth(L_king(i0,:),4,'lowess');
            L_king(i0,:) = L_king(i0,:)-L_king(i0,1);
            L_king(i0,:) = -L_king(i0,:)/L_king(i0,end);
            i0 = i0+1;
        end
    end
    axes(a(j));
    plot(tL_king,trimmean(L_king,20),'-.','color',c3,'linewidth',1);
end


% subplot labels
dtx = -0.035;
dty = 0.32;
letters1 = {'a','b','c'};
letters2 = {'d','e','f'};
for i=1:3,
    annotation('textbox','position',[lspace+(i-1)*(pw+hspace)+dtx,bspace+vspace+ph+dty,0.1,0.1],...
        'string',letters1{i},'fontsize',fs,'fontname',fn,'interpreter',fi,'fontweight','bold','edgecolor','none');
    annotation('textbox','position',[lspace+(i-1)*(pw+hspace)+dtx,bspace+dty,0.1,0.1],...
        'string',letters2{i},'fontsize',fs,'fontname',fn,'interpreter',fi,'fontweight','bold','edgecolor','none');
end

set(gcf,'renderer','painters');
saveplot_pdf(13,10,300,'extfig8.pdf');
close all;
    
end

%% extended fig 9
if plotnum == 9,
    
% run model
% [t,modeloutput] = glacier_response_everyglacier(ks,a,b);
load ../model/modeloutput.mat

lspace = 0.09;
hspace = 0.04;
vspace = 0.09;
rspace = 0.01;
bspace = 0.09;
tspace = 0.04;
pw = (1-lspace-rspace-2*hspace)/3;
ph = (1-tspace-bspace-vspace)/2;
ylims = [0.98,1.16];
tlims = [1979,2019];
tdx = 0.025;
tdy = 0.95;
tddy = 0.08;
fa = 0.3;
lw = 1.5;

casecols = [bothcol;qcol;tfcol];
inds = [1,2,3,4,5,6];
a(1) = axes('position',[lspace+0*(pw+hspace),bspace+1*(ph+vspace),pw,ph]);
a(2) = axes('position',[lspace+1*(pw+hspace),bspace+1*(ph+vspace),pw,ph]);
a(3) = axes('position',[lspace+2*(pw+hspace),bspace+1*(ph+vspace),pw,ph]);
a(4) = axes('position',[lspace+0*(pw+hspace),bspace+0*(ph+vspace),pw,ph]);
a(5) = axes('position',[lspace+1*(pw+hspace),bspace+0*(ph+vspace),pw,ph]);
a(6) = axes('position',[lspace+2*(pw+hspace),bspace+0*(ph+vspace),pw,ph]);

% add abbreviations to name
modeloutput(1).name = 'South (SO)';
modeloutput(2).name = 'Central-west (CW)';
modeloutput(3).name = 'North-west (NW)';
modeloutput(4).name = 'North-east (NE)';
modeloutput(5).name = 'North (NO)';

figure();
for i=1:6,
    axes(a(i)); hold on;
    for jj=3:-1:1,
        patch([t,fliplr(t)],[modeloutput(inds(i)).case(jj).D25,...
            fliplr(modeloutput(inds(i)).case(jj).D75)],...
            casecols(jj,:),'facealpha',fa,'edgecolor','none');
        plot(t,modeloutput(inds(i)).case(jj).D50,'color',casecols(jj,:),'linewidth',lw);
    end
    set(a(i),'fontsize',fs,'fontname',fn,'ticklabelinterpreter',fi,'box','on','xtick',[1980:10:2020]); grid on;
    xlim(tlims); ylim(ylims);
    title(modeloutput(inds(i)).name,'fontsize',fs,'fontname',fn,'interpreter',fi);
    if i==1 | i==4,
        ylabel('normalised ice discharge','fontsize',fs,'fontname',fn,'interpreter',fi);
        set(a(i),'ytick',[0.95:0.05:1.15]);
    else
        set(a(i),'ytick',[0.95:0.05:1.15],'yticklabel',[]);
    end
    if i==4,
        text(tlims(1)+diff(tlims)*tdx,ylims(1)+(tdy-0*tddy)*diff(ylims),'atmosphere (Q) varies',...
            'fontsize',fs,'fontname',fn,'fontweight','bold','interpreter',fi,'color',qcol,...
            'horizontalalignment','left');
        text(tlims(1)+diff(tlims)*tdx,ylims(1)+(tdy-1*tddy)*diff(ylims),'ocean (TF) varies',...
            'fontsize',fs,'fontname',fn,'fontweight','bold','interpreter',fi,'color',tfcol,...
            'horizontalalignment','left');
        text(tlims(1)+diff(tlims)*tdx,ylims(1)+(tdy-2*tddy)*diff(ylims),'both Q & TF vary',...
            'fontsize',fs,'fontname',fn,'fontweight','bold','interpreter',fi,'color',bothcol,...
            'horizontalalignment','left');
    end
    if i>3,
        xlabel('year','fontsize',fs,'fontname',fn,'interpreter',fi);
    end

end

% subplot labels
dtx = -0.035;
dty = 0.32;
letters1 = {'a','b','c'};
letters2 = {'d','e','f'};
for i=1:3,
    annotation('textbox','position',[lspace+(i-1)*(pw+hspace)+dtx,bspace+vspace+ph+dty,0.1,0.1],...
        'string',letters1{i},'fontsize',fs,'fontname',fn,'interpreter',fi,'fontweight','bold','edgecolor','none');
    annotation('textbox','position',[lspace+(i-1)*(pw+hspace)+dtx,bspace+dty,0.1,0.1],...
        'string',letters2{i},'fontsize',fs,'fontname',fn,'interpreter',fi,'fontweight','bold','edgecolor','none');
end

set(gcf,'renderer','painters');
saveplot_pdf(13,10,300,'extfig9.pdf');
close all;
    
end

%% extended fig 10
if plotnum == 10,

av = load('../model/modeloutput.mat');
ORAS5 = load('../model/modeloutput_ORAS5.mat');
CHORE = load('../model/modeloutput_CHORE.mat');
EN4 = load('../model/modeloutput_EN4.mat');

% renormalise CHORE to take account of earlier end of dataset
ind = find(av.t==CHORE.t(end));
for ii=1:length(CHORE.modeloutput),
    for jj=1:3,
        CHORE.modeloutput(ii).case(jj).L50 = CHORE.modeloutput(ii).case(jj).L50*abs(av.modeloutput(ii).case(1).L50(ind));
        CHORE.modeloutput(ii).case(jj).SL50 = CHORE.modeloutput(ii).case(jj).SL50*abs(av.modeloutput(ii).case(1).SL50(ind));
    end
end

lspace = 0.08;
hspace = 0.04;
rspace = 0.015;
bspace = 0.085;
vspace = 0.1;
tspace = 0.04;
pw = (1-lspace-rspace-3*hspace)/4;
ph = (1-bspace-tspace-vspace)/2;
tlims = [1979,2019];
tdx = 0.025;
tdy = 0.95;
tddy = 0.08;
fa = 0.3;
lw = 1.5;

tfcol1 = tfcol;
tfcol2 = 'r';
tfcol3 = tfcol;
tfcol4 = 'r';

casecols = [bothcol;qcol;tfcol];
inds = [1,2,3,6];
plottext = {'South (SO)','Central-west (CW)','North-west (NW)','All Greenland'};

figure();

% retreat
ylims = [-1.4,0.4];
for i=1:4,
    a(i) = axes('position',[lspace+(i-1)*(pw+hspace),bspace+vspace+ph,pw,ph]); hold on;

    % TF
    plot(ORAS5.t,ORAS5.modeloutput(inds(i)).case(3).L50,'--','color',tfcol,'linewidth',lw/2);
    plot(CHORE.t,CHORE.modeloutput(inds(i)).case(3).L50,':','color',tfcol,'linewidth',lw/2);
    plot(EN4.t,EN4.modeloutput(inds(i)).case(3).L50,'-.','color',tfcol,'linewidth',lw/2);
    plot(av.t,av.modeloutput(inds(i)).case(3).L50,'-','color',tfcol,'linewidth',lw);
    % Q
    plot(ORAS5.t,ORAS5.modeloutput(inds(i)).case(2).L50,'--','color',qcol,'linewidth',lw/2);
    plot(CHORE.t,CHORE.modeloutput(inds(i)).case(2).L50,':','color',qcol,'linewidth',lw/2);
    plot(EN4.t,EN4.modeloutput(inds(i)).case(2).L50,'-.','color',qcol,'linewidth',lw/2);
    plot(av.t,av.modeloutput(inds(i)).case(2).L50,'-','color',qcol,'linewidth',lw);

    set(a(i),'fontsize',fs,'fontname',fn,'ticklabelinterpreter',fi,'box','on','xtick',[1980:10:2020]); grid on;
    ylim(ylims); xlim(tlims);
    title(plottext{i},'fontsize',fs,'fontname',fn,'interpreter',fi);
    if i==1,
        ylabel('normalised retreat','fontsize',fs,'fontname',fn,'interpreter',fi);
        set(a(i),'ytick',[-1.2:0.2:0.2]);
    else
        set(a(i),'ytick',[-1.2:0.2:0.2],'yticklabel',[]);
    end
end

% sea level
ylims = [-0.1,1.1];
for i=1:4,
    a(i) = axes('position',[lspace+(i-1)*(pw+hspace),bspace,pw,ph]); hold on;
    % TF
    plot(ORAS5.t,ORAS5.modeloutput(inds(i)).case(3).SL50,'--','color',tfcol,'linewidth',lw/2);
    plot(CHORE.t,CHORE.modeloutput(inds(i)).case(3).SL50,':','color',tfcol,'linewidth',lw/2);
    plot(EN4.t,EN4.modeloutput(inds(i)).case(3).SL50,'-.','color',tfcol,'linewidth',lw/2);
    plot(av.t,av.modeloutput(inds(i)).case(3).SL50,'-','color',tfcol,'linewidth',lw);
    % Q
    plot(ORAS5.t,ORAS5.modeloutput(inds(i)).case(2).SL50,'--','color',qcol,'linewidth',lw/2);
    plot(CHORE.t,CHORE.modeloutput(inds(i)).case(2).SL50,':','color',qcol,'linewidth',lw/2);
    plot(EN4.t,EN4.modeloutput(inds(i)).case(2).SL50,'-.','color',qcol,'linewidth',lw/2);
    plot(av.t,av.modeloutput(inds(i)).case(2).SL50,'-','color',qcol,'linewidth',lw);
    
    set(a(i),'fontsize',fs,'fontname',fn,'ticklabelinterpreter',fi,'box','on','xtick',[1980:10:2020]); grid on;
    ylim(ylims); xlim(tlims);
    xlabel('year','fontsize',fs,'fontname',fn,'interpreter',fi);
    if i==1,
        ylabel({'normalised dynamic';'sea level contribution'},'fontsize',fs,'fontname',fn,'interpreter',fi);
        set(a(i),'ytick',[0:0.2:1]);
    else
        set(a(i),'ytick',[0:0.2:1],'yticklabel',[]);
    end

    % custom legend
    x1 = 1982;
    x2 = 1990;
    y1 = 1.02;
    dy = 0.1;
    tx1 = 1992;
    if i==1,        
        plot([x1,x2],(y1-0*dy)*[1,1],'k-','linewidth',lw);
        plot([x1,x2],(y1-1*dy)*[1,1],'k--','linewidth',lw/2);
        plot([x1,x2],(y1-2*dy)*[1,1],'k:','linewidth',lw/2);
        plot([x1,x2],(y1-3*dy)*[1,1],'k-.','linewidth',lw/2);
        text(tx1,y1-0*dy,'average','fontsize',fs,'verticalalignment','middle');
        text(tx1,y1-1*dy,'ORAS5','fontsize',fs,'verticalalignment','middle');
        text(tx1,y1-2*dy,'CHORE','fontsize',fs,'verticalalignment','middle');
        text(tx1,y1-3*dy,'EN4','fontsize',fs,'verticalalignment','middle');
    end

end

% subplot labels
dtx = -0.025;
dty = 0.32;
letters1 = {'a','b','c','d'};
letters2 = {'e','f','g','h'};
for i=1:4,
    annotation('textbox','position',[lspace+(i-1)*(pw+hspace)+dtx,bspace+1*(ph+vspace)+dty,0.1,0.1],...
        'string',letters1{i},'fontsize',fs,'fontname',fn,'interpreter',fi,'fontweight','bold','edgecolor','none');
    annotation('textbox','position',[lspace+(i-1)*(pw+hspace)+dtx,bspace+0*(ph+vspace)+dty,0.1,0.1],...
        'string',letters2{i},'fontsize',fs,'fontname',fn,'interpreter',fi,'fontweight','bold','edgecolor','none');
end

set(gcf,'renderer','painters');
saveplot_pdf(17,10,300,'extfig10.pdf');
close all;
    
end
