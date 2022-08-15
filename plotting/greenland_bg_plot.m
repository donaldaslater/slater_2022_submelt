% script to create greenland background plot
function [] = greenland_bg_plot(ax1,ss,xlims,ylims,fa);

% load bed topo and mask
bmfile = '~/Documents/BedMachineGreenland-2021-04-20.nc';
x = ncread(bmfile,'x');
y = ncread(bmfile,'y');
bed = ncread(bmfile,'bed')';
terrain = ncread(bmfile,'mask')';
% subsample
x = double(x(1:ss:end));
y = double(y(1:ss:end));
bed = double(bed(1:ss:end,1:ss:end));
bed0 = bed;
bed0(find(bed<0)) = NaN;
terrain = double(terrain(1:ss:end,1:ss:end));
terrain(terrain~=2) = NaN;

% if xlims and ylims exceed bedmachine range, use Rtopo2 to fill out
if xlims(1)<min(x) | xlims(2)>max(x) | ylims(1)<min(y) | ylims(2)>max(y),

    dx = abs(diff(x(1:2))); dy = abs(diff(y(1:2)));
    x_extend = [fliplr([min(x)-dx:-dx:xlims(1)-dx]),x',max(x)+dx:dx:xlims(2)+dx]';
    y_extend = [fliplr([max(y)+dy:dy:ylims(2)+dy]),y',min(y)-dy:-dy:ylims(1)-dy]';
    bed_extend = NaN(length(y_extend),length(x_extend));
    bed_extend(find(y_extend==y(1)):find(y_extend==y(end)),find(x_extend==x(1)):find(x_extend==x(end))) = bed;
    terrain_extend = NaN(length(y_extend),length(x_extend));
    terrain_extend(find(y_extend==y(1)):find(y_extend==y(end)),find(x_extend==x(1)):find(x_extend==x(end))) = terrain;

    rtopo = load('../../data/bathymetry/bathy.mat');
    [rtopo.LAT,rtopo.LON]=meshgrid(rtopo.lat,rtopo.lon);
    [rtopo.X,rtopo.Y]=latlon2utm(rtopo.LAT(:),rtopo.LON(:));
    rtopo.topo = rtopo.topo(:);
    inds = find(rtopo.X<=xlims(2)+20000 & rtopo.X>=xlims(1)-20000 & rtopo.Y<=ylims(2)+20000 & rtopo.Y>=ylims(1)-20000);
    rtopo.X = rtopo.X(inds);
    rtopo.Y = rtopo.Y(inds);
    rtopo.topo = rtopo.topo(inds);
    f = scatteredInterpolant(double(rtopo.X),double(rtopo.Y),double(rtopo.topo),'linear','none');

    for i=1:size(bed_extend,1),
        i/size(bed_extend,1)
        for j=1:size(bed_extend,2),
            if isnan(bed_extend(i,j)),
                bed_extend(i,j) = f(x_extend(j),y_extend(i));
            end
        end
    end

    x = x_extend;
    y = y_extend;
    bed = bed_extend;
    bed0 = bed;
    bed0(find(bed<0)) = NaN;
    terrain = terrain_extend;

end

% colormaps
cmap_bed = cptcmap('GMT_relief');
cmap_ice = [1,1,0.8];
cmap_ocean = flipud(cbrewer('div','RdYlBu',100));

% get ocean bottom thermal forcing from ORAS5
ocean = load('../data/ocean_ORAS5_plot.mat');
Tmean = ocean.potentialT;
Smean = ocean.practicalS;
zmean = repmat(ocean.z',1,length(ocean.x))';
l1 = -5.73e-2;
l2 = 8.32e-2;
l3 = 7.53e-4;
TFmean = Tmean-(l1*Smean+l2+l3*zmean);

% interpolate ocean thermal forcing onto bathymetry
X_TF = NaN(size(TFmean));
Y_TF = NaN(size(TFmean));
Z_TF = NaN(size(TFmean));
for ii=1:size(TFmean,1),
    X_TF(ii,:) = ocean.x(ii);
    Y_TF(ii,:) = ocean.y(ii);
end
for jj=1:size(TFmean,2),
    Z_TF(:,jj) = ocean.z(jj);
end
inds_to_remove = find(isnan(TFmean));
TFmean(inds_to_remove) = [];
X_TF(inds_to_remove) = [];
Y_TF(inds_to_remove) = [];
Z_TF(inds_to_remove) = [];
warning off;
F = scatteredInterpolant(double(X_TF(:)),double(Y_TF(:)),double(Z_TF(:)),double(TFmean(:)),'linear','nearest');
warning on;
Zq = bed;
Zq(Zq<-1000) = -1000; % TF depth maxes out at 1000 m depth
for ii=1:size(Zq,1),
    for jj=1:size(Zq,2),
        Xq(ii,jj) = x(jj);
        Yq(ii,jj) = y(ii);
    end
end
TF0 = F(Xq,Yq,Zq);
TF0(find(Zq>0)) = NaN;    

% plot ocean TF     
p1 = pcolor(ax1,x,y,TF0);
set(p1,'facealpha',fa);
shading flat;
colormap(ax1,cmap_ocean);
caxis([1 7]);
xlim(xlims); ylim(ylims);
set(ax1,'xtick',[],'ytick',[]);

% plot topo
ax2 = axes('position',get(ax1,'position')); hold on;
p2 = pcolor(ax2,x,y,bed0);
set(p2,'facealpha',fa);
shading flat;
colormap(ax2,cmap_bed);
caxis([-3000 3000]);
xlim(xlims); ylim(ylims);
set(ax2,'visible','off','color','none','box','off','xtick',[],'ytick',[]);

% plot ice
ax3 = axes('position',get(ax1,'position')); hold on;
p3 = pcolor(ax3,x,y,terrain);
shading flat;
colormap(ax3,cmap_ice);
xlim(xlims); ylim(ylims);
set(ax3,'visible','off','color','none','box','off','xtick',[],'ytick',[]);

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

end

