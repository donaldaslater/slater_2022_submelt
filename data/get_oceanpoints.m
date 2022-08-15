% script to get best point from which to sample ocean for each glacier
clear; close all;

% plotting parameters for illustration
fs = 8;
ms = 3;
lspace = 0.1;
bspace = 0.13;
tspace = 0.1;
rspace = 0.15;
cbspace = 0.05;
cbw = 0.02;
pw = 1-lspace-rspace;
ph = 1-bspace-tspace;
plotsavepath = '~/Google Drive/data/master2/plots/';
cmap_bed = cptcmap('GMT_relief');

% load glaciers structure so far
load glaciers.mat

% load bedmachine
bmfile = '~/Documents/BedMachineGreenland-2020-07-16.nc';
bm.x = double(ncread(bmfile,'x'));
bm.y = double(ncread(bmfile,'y'));
bm.bed = double(ncread(bmfile,'bed'));

% product to analyse
product = 'CHORE'; % ORAS5 or EN4 or ASTE or CHORE

if strcmp(product,'ORAS5'),
    
    % get ocean dataset bathymetry
    lon = ncread('~/Documents/ORAS5/votemper_ORAS5_1m_201812_grid_T_02.nc','nav_lon');
    lat = ncread('~/Documents/ORAS5/votemper_ORAS5_1m_201812_grid_T_02.nc','nav_lat');
    [x,y] = latlon2utm(lat,lon);
    x_oceanmodel = x(:)';
    y_oceanmodel = y(:)';
    T = ncread('~/Documents/ORAS5/votemper_ORAS5_1m_201812_grid_T_02.nc','votemper');
    Tdepth = reshape(T,size(T,1)*size(T,2),size(T,3));
    z = -ncread('~/Documents/ORAS5/votemper_ORAS5_1m_201812_grid_T_02.nc','deptht');
    
elseif strcmp(product,'EN4'),
    
    % get EN4 depths
    lon = ncread('~/Documents/EN4/EN.4.2.1.f.analysis.g10.202006.nc','lon');
    lat = ncread('~/Documents/EN4/EN.4.2.1.f.analysis.g10.202006.nc','lat');
    z = -ncread('~/Documents/EN4/EN.4.2.1.f.analysis.g10.202006.nc','depth');
    T = ncread('~/Documents/EN4/EN.4.2.1.f.analysis.g10.202006.nc','temperature');
    [lat,lon]=meshgrid(lat,lon);
    [x,y] = latlon2utm(lat,lon);
    Tdepth = reshape(T,size(T,1)*size(T,2),size(T,3));
    x_oceanmodel = x(:)';
    y_oceanmodel = y(:)';
    
elseif strcmp(product,'ASTE'),
    
    % get ASTE depths
    tiles = [5,11,12,14,15,27];
    for i=1:length(tiles),
        if tiles(i)<10, numstr = ['000',num2str(tiles(i))];
        else numstr = ['00',num2str(tiles(i))];
        end
        lon(i,:,:) = ncread(['~/Documents/ASTE/THETA.',numstr,'.nc'],'lon');
        lat(i,:,:) = ncread(['~/Documents/ASTE/THETA.',numstr,'.nc'],'lat');
        Ti = ncread(['~/Documents/ASTE/THETA.',numstr,'.nc'],'THETA');
        Ti = squeeze(Ti(:,:,:,1));
        land = ncread(['~/Documents/ASTE/THETA.',numstr,'.nc'],'land');
        Ti(find(land==0)) = NaN;
        T(i,:,:,:) = Ti;
        tileno(i,:,:) = 0*lon(i,:,:)+tiles(i);
    end
    lon = lon(:);
    lat = lat(:);
    tileno = tileno(:);
    z = -ncread(['~/Documents/ASTE/THETA.',numstr,'.nc'],'dep');
    Tdepth = reshape(T,size(T,1)*size(T,2)*size(T,3),size(T,4));
    [x,y] = latlon2utm(lat,lon);
    x_oceanmodel = x';
    y_oceanmodel = y';
    
elseif strcmp(product,'CHORE'),
    
    % get EN4 depths
    lon = ncread('~/Documents/CHORE/GRID.nc','lon');
    lat = ncread('~/Documents/CHORE/GRID.nc','lat');
    z = -ncread('~/Documents/CHORE/GRID.nc','dep');
    T = ncread('~/Documents/CHORE/GRID.nc','tmsk');
    [x,y] = latlon2utm(lat,lon);
    T(find(T==0))=NaN;
    Tdepth = reshape(T,size(T,1)*size(T,2),size(T,3));
    x_oceanmodel = x(:)';
    y_oceanmodel = y(:)';

end

% get product ocean bathmetry
for i=1:size(Tdepth,1),
    if ~isempty(find(~isnan(Tdepth(i,:)))),
        depth_oceanmodel(i) = z(max(find(~isnan(Tdepth(i,:)))));
    else
        depth_oceanmodel(i) = NaN;
    end
end

% loop over glaciers
for ii=1:length(glaciers),
    
    final_ind = [];
    x_glacier = glaciers(ii).x;
    y_glacier = glaciers(ii).y;
    gldepth = glaciers(ii).gldepth;
    
    if gldepth<0,
    
        % vertical search points
        z0 = [50*ceil(gldepth/50):50:-50,-25,0];

        % find distance to shelf break to set search radius
        load('bathy1000.mat');
        dist_to_shelf = sqrt((x1000-x_glacier).^2+(y1000-y_glacier).^2);
        searchradius = 50+min(dist_to_shelf)/1000;

            kk = 0;
            foundflag = 0;

            % restrict to searchradius to speed things up
            xlims = x_glacier + 1000*searchradius*[-1,1];
            ylims = y_glacier + 1000*searchradius*[-1,1];
            x_inds = find(bm.x<=xlims(2) & bm.x>=xlims(1));
            y_inds = find(bm.y<=ylims(2) & bm.y>=ylims(1));
            x = bm.x(x_inds);
            y = bm.y(y_inds);
            bed = bm.bed(x_inds,y_inds);

            % get indices of glacier on bed subset
            x_glacier_ind = max(find(x<=x_glacier));
            y_glacier_ind = min(find(y<=y_glacier));

            while foundflag == 0,

                % advance vertical search
                kk = kk+1;

                % connectivity maps
                % everything shallower than z0(kk) is assigned 0
                % everything deeper than or equal to z0(kk) is assigned 1
                bed0 = bed;
                bed0(bed<=z0(kk))=1;
                bed0(bed>z0(kk))=0;

                % get connectivity map of points assigned with 1
                a = bwconncomp(bed0);

                % array of connectivity groups
                conngroups = NaN(size(bed));
                for j=1:a.NumObjects,
                    conngroups(a.PixelIdxList{j}) = j;
                end

                % get connectivity group of glacier
                glacier_conn = interp2(x,y,conngroups',x_glacier,y_glacier,'nearest');

                % for plotting purposes, make modified conngroups
                connplot = conngroups;
                connplot(find(bed>0)) = 0;
                connplot(find(bed<=0)) = 1;
                connplot(find(conngroups == glacier_conn)) = 2;

                % find ocean model points which are sufficiently deep
                oceanmodel_deepenough = find(depth_oceanmodel<=z0(kk));

                % find ocean model points which are deep enough and connected
                % get connectivity group of chore points
                oceanmodel_conn = interp2(x,y,conngroups',x_oceanmodel,y_oceanmodel,'nearest');
                % get chore points in same connectivity group
                oceanmodel_glaciergroup = find(oceanmodel_conn==glacier_conn);

                % acceptable points are those which are deep enough and connected
                inds_acceptable = intersect(oceanmodel_deepenough,oceanmodel_glaciergroup);
                if ~isempty(inds_acceptable),
                    % calculate distance to all points in connectivity group
                    % following path through connectivity group
                    distmatrix = conngroups;
                    distmatrix(find(conngroups==glacier_conn))=1;
                    distmatrix(find(conngroups~=glacier_conn))=0;
                    dd = (x(2)-x(1))/1000; % resolution in km
                    D = dd*bwdistgeodesic(logical(distmatrix),y_glacier_ind,x_glacier_ind,'quasi-euclidean');
                    oceanmodel_dist = interp2(x,y,D',x_oceanmodel(inds_acceptable),y_oceanmodel(inds_acceptable),'nearest');   
                    % find closest eligible ocean model point to glacier
                    min_ind = find(oceanmodel_dist==min(oceanmodel_dist));
                    % if there is more than 1 point, take the first
                    final_ind = inds_acceptable(min_ind(1));
                else
                    final_ind = NaN;
                end

                % plot if there is an acceptable point
                if ~isnan(final_ind),

                    figure('visible','off');
                    a1 = axes('position',[lspace,bspace,pw,ph]); hold on;
                    imagesc(x,y,bed'); axis xy;
                    h = colorbar('position',[lspace+pw+cbspace,bspace,cbw,ph]);
                    colormap(cmap_bed); caxis([-500 500]);

                    plot(x_glacier,y_glacier,'rs','markerfacecolor','r','markersize',ms);
                    plot(x_oceanmodel(oceanmodel_deepenough),y_oceanmodel(oceanmodel_deepenough),'bo','markerfacecolor','b','markersize',ms);
                    plot(x_oceanmodel(intersect(oceanmodel_deepenough,oceanmodel_glaciergroup)),y_oceanmodel(intersect(oceanmodel_deepenough,oceanmodel_glaciergroup)),'go','markerfacecolor','g','markersize',ms);
                    plot(x_oceanmodel(final_ind),y_oceanmodel(final_ind),'ro','markersize',2*ms);
                    xlim(xlims); ylim(ylims);
                    set(gca,'fontsize',fs,'box','on');

                    % save number
                    if glaciers(ii).morlighem_number<10, savenum = ['00',num2str(glaciers(ii).morlighem_number)];
                    elseif glaciers(ii).morlighem_number<100, savenum = ['0',num2str(glaciers(ii).morlighem_number)];
                    else savenum = num2str(glaciers(ii).morlighem_number);
                    end

                    % save plot
                    fw = 15;
                    fh = fw*pw/ph;
                    saveplot(fw,fh,100,[plotsavepath,savenum,'_',product,'_',num2str(abs(z0(kk))),'.png']);
                    close all;

                    foundflag = 1;
                    searchradiusflag = 1;

                end

                if kk==length(z0), foundflag = 1; end

            end

        if isnan(final_ind),
            glaciers(ii).ocean.(product).x = NaN;
            glaciers(ii).ocean.(product).y = NaN;
            glaciers(ii).ocean.(product).ind = NaN;
            glaciers(ii).ocean.(product).effdepth = NaN;
        else
            glaciers(ii).ocean.(product).x = single(x_oceanmodel(final_ind));
            glaciers(ii).ocean.(product).y = single(y_oceanmodel(final_ind));
            glaciers(ii).ocean.(product).ind = single(final_ind);
            glaciers(ii).ocean.(product).effdepth = single(z0(kk));
        end
    
    else
        
        glaciers(ii).ocean.(product).x = NaN;
        glaciers(ii).ocean.(product).y = NaN;
        glaciers(ii).ocean.(product).ind = NaN;
        glaciers(ii).ocean.(product).effdepth = NaN;
    
    end

end

% save structure
save glaciers.mat glaciers

%% once have run all ocean products, make a .mat structure of indices

load glaciers.mat

for ii=1:length(glaciers),
    if ~isempty(glaciers(ii).morlighem_number),
        morlighem_n(ii) = glaciers(ii).morlighem_number;
    else
        morlighem_n(ii) = NaN;
    end
    ind_EN4(ii) = glaciers(ii).ocean.EN4.ind;
    ind_ORAS5(ii) = glaciers(ii).ocean.ORAS5.ind;
    ind_ASTE(ii) = glaciers(ii).ocean.ASTE.ind;
    ind_CHORE(ii) = glaciers(ii).ocean.CHORE.ind;
end

save oceanindices.mat morlighem_n ind_EN4 ind_ORAS5 ind_ASTE ind_CHORE