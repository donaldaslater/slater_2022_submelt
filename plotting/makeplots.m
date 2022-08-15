function [] = makeplots(plotnum);

clearvars -except plotnum;
close all;

fs = 7;
fn = 'helvetica-narrow';
fi = 'tex';

cols = [0.000,0.447,0.741;
        0.850,0.325,0.098;
        0.929,0.694,0.125;
        0.494,0.184,0.556;
        0.466,0.674,0.188;
        0.301,0.745,0.933;
        0.635,0.078,0.184];
    
c1 = [0.000,0.447,0.741];
c2 = [0.850,0.325,0.098];
c3 = [0.929,0.694,0.125];
c4 = [0.494,0.184,0.556];
c5 = [0.466,0.674,0.188];
c6 = [0.301,0.745,0.933];
c7 = [0.635,0.078,0.184];

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

%% fig 2a
if plotnum == 2,
    
    % plot parameters
    lspace = 0.01;
    rspace = 0.01;
    bspace = 0.01;
    tspace = 0.01;
    ph = 1-bspace-tspace;
    pw = 1-lspace-rspace;
    
    figure();
    ax1 = axes('position',[lspace,bspace,pw,ph]);
    xlims = [-7.6,9.9]*10^5;
    ylims = [-3.42,-0.53]*10^6;
    tlims = [1918,2010];
    ss = 10;
    fa = 0.65;
    
    greenland_bg_plot(ax1,ss,xlims,ylims,fa); 

    % get submarine melt rates
    load ../data/twglaciers.mat
    for i=1:length(twglaciers),
        melt(i) = nanmean(twglaciers(i).submelt.m);
    end
    x = [twglaciers.x];
    y = [twglaciers.y];
    
    % sort so that largest is on top
    [melt,ind] = sort(melt);
    x = x(ind);
    y = y(ind);

    % lines of lat-lon
    lat0 = [55:10:80];
    lon0 = [-100:20:20];
    lwl = 0.25;
    lcl = 0.7*[1,1,1];
    nn = 500;
    for i=1:length(lat0),
        [lgrid(i).x,lgrid(i).y] = latlon2utm(lat0(i)*ones(1,nn),linspace(-180,180,nn));
        plot(lgrid(i).x,lgrid(i).y,'color',lcl,'linewidth',lwl);
    end
    for i=1:length(lon0),
        [lgrid(i).x,lgrid(i).y] = latlon2utm(linspace(50,90,nn),lon0(i)*ones(1,nn));
        plot(lgrid(i).x,lgrid(i).y,'color',lcl,'linewidth',lwl);
    end
    
    % add a few drainage basins
    ax3 = axes('position',[lspace,bspace,pw,ph]); hold on;
    dd = load('../data/drainagebasins.mat');
    for i=1:3,
        db = dd.drainagebasins;
        db(find(dd.drainagebasins~=i)) = i+1;
        contour(dd.x(1:ss:end),dd.y(1:ss:end),db(1:ss:end,1:ss:end),[i,i+1e-6],'k--','linewidth',0.5);
    end 
    
    % plot melt rates
    mms = 140;    
    scatter(x,y,(mms/max(sqrt(melt)))*sqrt(melt),melt,'filled','s','markerfacealpha',0.8,'markeredgecolor','k','markeredgealpha',0.8);
    
    % add shelf edge
    load ../data/bathy1000.mat
    plot(x1000,y1000,'k:','linewidth',0.5);  
        
    % submarine melting legend
    xc = 2e5;
    yc = -1.5e6;
    m_leg = [0.1,0.5,1,2,3,4,5];
    x_leg = xc*ones(1,length(m_leg));
    y_leg = linspace(yc-0.25e6,yc+0.25e6,length(m_leg));
    scatter(x_leg,y_leg,(mms/max(sqrt(melt)))*sqrt(m_leg),m_leg,'filled','s','markerfacealpha',0.8,'markeredgecolor','k','markeredgealpha',0.8);    
    dx = 40000;
    for i=1:length(m_leg),
        text(x_leg(i)+dx,y_leg(i),num2str(m_leg(i)),'fontsize',fs,'verticalalignment','middle','fontname',fn,'interpreter',fi);
    end
    text(x_leg(1)+4*dx,y_leg(4),['submarine melt rate (m',char(8201),'d^{-1})'],'fontsize',fs,'rotation',90,'horizontalalignment','center','fontname',fn,'interpreter',fi);    
    set(gca,'xtick',[],'ytick',[],'color','none','visible','off','box','on');
    xlim(xlims); ylim(ylims);
    colormap(ax3,parula(100)); caxis([0 5]);
    
    % ocean thermal forcing colorbar
    xc_ax = lspace+pw*(xc-xlims(1))/diff(xlims);
    yc_ax = bspace+ph*(yc-ylims(1))/diff(ylims);
    ax4 = axes('position',[lspace,bspace,pw,ph]);
    h = colorbar('position',[xc_ax-0.15,yc_ax-0.1,0.02,0.2],'axislocation','in','fontsize',fs,'ticklabelinterpreter',fi,'fontname',fn);
    ylabel(h,['ocean thermal forcing (',char(176),'C)'],'fontsize',fs,'interpreter',fi,'fontname',fn);
    set(h,'fontsize',fs);
    colormap(ax4,flipud(cbrewer('div','RdYlBu',100))); caxis([1 7]);
    set(ax4,'visible','off','color','none');

    fw = 8.75;
    fh = fw*(diff(ylims)/diff(xlims))*(pw/ph);
    saveplot_pdf(fw,fh,300,'fig2a.pdf');
    close all;
    
end

%% fig 2b
if plotnum == 2.1,

    lspace = 0.09;
    bspace = 0.14;
    tspace = 0.02;
    rspace = 0.16;
    cbw = 0.02;
    cbh = 0.6;
    cbspace = 0.04;
    pw = 1-lspace-rspace;
    ph = 1-bspace-tspace;
    figure();
    ax1 = axes('position',[lspace,bspace,pw,ph]); hold on;
    regions_name = {'South','Central-West','North-West','North-East','North'};
    ms1 = 3;
    ms2 = 7;
    ms3 = 8;
    fs1 = 6;
    fs2 = 10;
    cmap = parula(10);

    load ../data/twglaciers.mat
    Qmean = double([twglaciers.Qmean]+[twglaciers.basalmelt]);
    TFmean = double([twglaciers.TFmean]);

    % make plot background
    Q = logspace(0,3,1000);
    TF = linspace(0,8,1000);
    for i=1:length(Q),
        for j=1:length(TF),
            melt(i,j) = ks*Q(i)^a*TF(j)^b;
        end
    end
    p = pcolor(Q,TF,melt'); shading flat;
    set(p,'facealpha',0.5);
    set(gca,'box','on','fontsize',fs,'ticklabelinterpreter',fi,'fontname',fn);
    xlabel(['subglacial discharge (m^3',char(8201),'s^{-1})'],'fontsize',fs,'interpreter',fi,'fontname',fn);
    ylabel(['ocean thermal forcing (',char(176),'C)'],'fontsize',fs,'interpreter',fi,'fontname',fn);

    h = colorbar('position',[lspace+pw+cbspace,bspace+(ph-cbh)/2,cbw,cbh]);
    set(h,'fontsize',fs,'ticklabelinterpreter',fi,'fontname',fn);
    xlabel(h,['submarine melt rate (m',char(8201),'d^{-1})'],'fontsize',fs,'interpreter',fi,'fontname',fn);
    caxis([0 8]); colormap(cmap);

    % sector indices
    inds(1).inds = find([twglaciers.ismip6_sector]<=3); % south
    inds(2).inds = find([twglaciers.ismip6_sector]==4); % central-west
    inds(3).inds = find([twglaciers.ismip6_sector]==6); % north-west
    inds(4).inds = find([twglaciers.ismip6_sector]==5); % north-east
    inds(5).inds = find([twglaciers.ismip6_sector]==7); % north-north

    plot(Qmean,TFmean,'s','markerfacecolor','k','markeredgecolor','k','markersize',ms1);
    for l=1:5,
        plot(Qmean(inds(l).inds),TFmean(inds(l).inds),'s','markerfacecolor',cols(l,:),'markeredgecolor',cols(l,:),'markersize',ms1);
    end

    for l=1:5,
        sectorinds = inds(l).inds;
        % get stats
        xdata = Qmean(sectorinds); ydata = TFmean(sectorinds);
        x50 = prctile(xdata,50); y50 = prctile(ydata,50);
        x25 = prctile(xdata,25); y25 = prctile(ydata,25);
        x75 = prctile(xdata,75); y75 = prctile(ydata,75);
        errorbar(x50,y50,y50-y25,y75-y50,x50-x25,x75-x50,'s','markerfacecolor',cols(l,:),'markeredgecolor','k','markersize',ms2,'color','k','capsize',4)
        TFregionmean(l) = y50;
    end

    set(gca,'xscale','log'); xlim([2 250]); ylim([0 7.5]);

    fw = 8;
    fh = 7;
    saveplot_pdf(fw,fh,300,'fig2b.pdf');
    close all;

end

%% fig 2c
if plotnum == 2.2,

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

lspace = 0.13;
hspace = 0.05;
rspace = 0.02;
bspace = 0.12;
tspace = 0.02;
pw = 1-lspace-rspace;
ph = 1-bspace-tspace;
ms = 12;

% melt vs gl depth
a3 = axes('position',[lspace,bspace,pw,ph]); hold on;
for l=1:5,
    sectorinds = inds(l).inds;
    xdata = meanmelt(sectorinds); ydata = zgl(sectorinds);
    p = polyfit(ydata,xdata,1)
    [r,pval] = corrcoef(ydata',xdata');
    pval = pval(1,2);
    scatter(xdata,ydata,ms,'marker','s','markerfacecolor',cols(l,:),'markeredgecolor','none','markerfacealpha',0.8);
    if pval<0.05,
        plot(p(1)*[min(ydata),max(ydata)]+p(2),[min(ydata),max(ydata)],'--','color',cols(l,:),'linewidth',1);
    end
end

xlabel(['submarine melt rate (m/d)'],'fontsize',fs,'fontname',fn,'interpreter',fi);
ylabel('grounding line depth (m)','fontsize',fs,'fontname',fn,'interpreter',fi);
set(gca,'box','on','fontsize',fs,'fontname',fn,'ticklabelinterpreter',fi,...
    'ytick',[-800:200:0],'yticklabel',{'800','600','400','200','0'});
xlim([0 5]); ylim([-910 0]);

saveplot_pdf(8,6.8,300,'fig2c.pdf');
close all;

end



%% fig 3
if plotnum == 3,
    
    % load data
    load ../data/twglaciers.mat
    % time periods
    % normalisation period
    t0 = [1979,1993];
    t = twglaciers(1).submelt.t;
    yrs = unique(floor(t));
    inds0 = find(t>=t0(1) & t<=t0(end));
    inds1 = find(yrs>=t0(1) & yrs<=t0(end));
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

            % Q summer only
            for k=1:length(yrs),
                ids = find(t-yrs(k)>0.4 & t-yrs(k)<0.7);
                Qs(k) = mean(Qa(ids));
            end
            Qs = smooth(Qs,smt/12,smm);
            Qs0 = mean(Qs(inds1));
            Qsnorm(i0,:) = 100*(Qs/Qs0-1);

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

            % summer melt rate
            for k=1:length(yrs),
                ids = find(t-yrs(k)>0.4 & t-yrs(k)<0.7);
                ms(k) = mean(ma(ids));
            end
            ms = smooth(ms,smt/12,smm);
            ms0 = mean(ms(inds1));
            msnorm(i0,:) = 100*(ms/ms0-1);
        
            i0=i0+1;
        end
    end
    
    % get sector indices
    inds(1).inds = find(secs<=3); inds(1).name = 'South (SO)';
    inds(2).inds = find(secs==4); inds(2).name = 'Central-west (CW)';
    inds(3).inds = find(secs==6); inds(3).name = 'North-west (NW)';
    inds(4).inds = [1:length(secs)]; inds(4).name = 'All Greenland';
    
    % plots
    lspace = 0.08;
    hspace = 0.04;
    rspace = 0.015;
    bspace = 0.085;
    vspace = 0.1;
    tspace = 0.04;
    pw = (1-lspace-rspace-3*hspace)/4;
    ph = (1-bspace-tspace-vspace)/2;
    
    ylims1 = [-20,100];
    ylims2 = [-12,55];
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
    for i=1:4,
        
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
    for i=1:4,
        
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
            set(gca,'ytick',[-25:25:75]);
        else
            set(gca,'ytick',[-25:25:75],'yticklabel',[]);
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
    
    saveplot_pdf(17,10,300,'fig3.pdf');
    close all;
    
end

%% fig 4
if plotnum == 4,
    
% run model
% tic
% [t,modeloutput] = glacier_response_everyglacier(ks,a,b);
% toc
load ../model/modeloutput.mat

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

casecols = [bothcol;qcol;tfcol];
inds = [1,2,3,6];
plottext = {'South (SO)','Central-west (CW)','North-west (NW)','All Greenland'};

figure();

% retreat
ylims = [-1.4,0.4];
for i=1:4,
    a(i) = axes('position',[lspace+(i-1)*(pw+hspace),bspace+vspace+ph,pw,ph]); hold on;
    for jj=3:-1:1,
        patch([t,fliplr(t)],[modeloutput(inds(i)).case(jj).L25,...
            fliplr(modeloutput(inds(i)).case(jj).L75)],...
            casecols(jj,:),'facealpha',fa,'edgecolor','none');
        plot(t,modeloutput(inds(i)).case(jj).L50,'color',casecols(jj,:),'linewidth',lw);
    end
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
    for jj=3:-1:1,
        patch([t,fliplr(t)],[modeloutput(inds(i)).case(jj).SL25,...
            fliplr(modeloutput(inds(i)).case(jj).SL75)],...
            casecols(jj,:),'facealpha',fa,'edgecolor','none');
        plot(t,modeloutput(inds(i)).case(jj).SL50,'color',casecols(jj,:),'linewidth',lw);
    end
    set(a(i),'fontsize',fs,'fontname',fn,'ticklabelinterpreter',fi,'box','on','xtick',[1980:10:2020]); grid on;
    ylim(ylims); xlim(tlims);
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
saveplot_pdf(17,10,300,'fig4.pdf');
close all;
    
end