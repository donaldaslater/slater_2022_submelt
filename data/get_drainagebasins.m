% script to delineate hydrological drainage basins for
% greenland's tidewater glaciers
clear; close all;

%% load BedMachine
disp('Loading BedMachine...');

bmfile = '~/Documents/BedMachineGreenland-2020-07-16.nc';
x = double(ncread(bmfile,'x'));
y = double(ncread(bmfile,'y'));
b = double(ncread(bmfile,'bed'))';
h = double(ncread(bmfile,'thickness'))';
m = double(ncread(bmfile,'mask'))';
m(find(m==3))=0; % set floating ice to same mask value as ocean
bmres = 150; % resolution in meters

% trim to subset if wanted for testing
% x_inds = find(x>-0.4e6 & x<0e6);
% y_inds = find(y>-2.45e6 & y<-1.9e6);
% x = x(x_inds);
% y = y(y_inds);
% b = b(y_inds,x_inds);
% h = h(y_inds,x_inds);
% m = m(y_inds,x_inds);

%% flow routing
disp('Doing flow routing...');

% set path to routing routines
restoredefaultpath;
addpath('/Users/donaldslater/Google Drive/matlab_scripts/topotoolbox-2.2/');

% constants
rho_w = 1000; % meltwater density
rho_i = 917; % ice density
g = 9.81; % gravity
f = 1; % flotation fraction

% hydropotential phi
phi = rho_w*g*b + f*rho_i*g*h;
% mask out ocean, floating ice and non-greenland
% keep land in so that water can flow off land onto ice
phi(m==0 | m==3 | m==4) = NaN;

% set up for flow routing
[X,Y] = meshgrid(double(x),double(y)); % create coord arrays
PHI = GRIDobj(X,Y,double(phi)); % create grid object for topotoolbox
PHI = fillsinks(PHI); % fill holes in hydropotential
FR = FLOWobj(PHI); % make flow object for topotoolbox

% do flow routing
A = flowacc(FR); % upstream area
A = double(A.Z); % just keep area
[DB,outlet_x,outlet_y] = drainagebasins(FR); % basins and outlets
DB = double(DB.Z); % just keep basins

% set ocean values to NaN
DB(find(DB==0)) = NaN;

%% get ice-ocean boundary pixels and do manual fixes
disp('Getting ice-ocean boundary pixels...');

% get ice-ocean points
iceoceanmask = zeros(size(m));
for i=2:size(m,1)-1,
    for j=2:size(m,2)-1,
        if m(i,j) == 2,
            testmask = m(i-1:i+1,j-1:j+1);
            if any(ismember(testmask(:),0)),
                iceoceanmask(i,j) = 1;
            end
        end
    end
end

% manual fixes
% these are variously used to separate out calving fronts that are joined
% according to the bedmachine mask, and occasionally to join a calving
% front to a larger drainage basin if it missed by one pixel
% having to do this is annoying because it needs to be redone if the
% bedmachine mask changes, but i can't see a better way...
iceoceanmask(find(y==-3200375),find(x==-153125)) = 1;
iceoceanmask(find(y==-1642025),find(x==-324275)) = 0;
iceoceanmask(find(y==-1850975),find(x==-302225)) = 0;
iceoceanmask(find(y==-1827425),find(x==-304025)) = 0;
iceoceanmask(find(y==-1567925),find(x==-356225)) = 0;
iceoceanmask(find(y==-1584875),find(x==-361775)) = 0;
iceoceanmask(find(y==-1429775),find(x==-436475)) = 0;
iceoceanmask(find(y==-1406825),find(x==-486275)) = 0;
iceoceanmask(find(y==-1392425),find(x==-573875)) = 0;
iceoceanmask(find(y==-1388375),find(x==-565925)) = 0;
iceoceanmask(find(y==-1384175),find(x==-554075)) = 0;
iceoceanmask(find(y==-1414925),find(x==-475775)) = 0;
iceoceanmask(find(y==-909125),find(x==-2225)) = 0;
iceoceanmask(find(y==-835625),find(x==69025)) = 0;
iceoceanmask(find(y==-1080425),find(x==502225)) = 0;
iceoceanmask(find(y==-1340225),find(x==546325)) = 0;
iceoceanmask(find(y==-2295575),find(x==491125)) = 0;
iceoceanmask(find(y==-2297075),find(x==493525)) = 0;
iceoceanmask(find(y==-2404025),find(x==497425)) = 0;
iceoceanmask(find(y==-2685275),find(x==229075)) = 0;
iceoceanmask(find(y==-2804675),find(x==215875)) = 0;
iceoceanmask(find(y==-2840525),find(x==177325)) = 0;
iceoceanmask(find(y==-3097025),find(x==135925)) = 0;
iceoceanmask(find(y==-3182975),find(x==99925)) = 0;
iceoceanmask(find(y==-2911025),find(x==164575)) = 1;
iceoceanmask(find(y==-1602875),find(x==664375)) = 1;
iceoceanmask(find(y==-1971725),find(x==606025)) = 1;
iceoceanmask(find(y==-1971875),find(x==605875)) = 1;
iceoceanmask(find(y==-1453775),find(x==601975)) = 1;
iceoceanmask(find(y==-1453325),find(x==602875)) = 1;
iceoceanmask(find(y==-1453325),find(x==603025)) = 1;
iceoceanmask(find(y==-1305875),find(x==-509825)) = 1;
iceoceanmask(find(y==-1305875),find(x==-509675)) = 1;
iceoceanmask(find(y==-1305725),find(x==-509375)) = 1;
iceoceanmask(find(y==-1305725),find(x==-509225)) = 1;
iceoceanmask(find(y==-2599625),find(x==385525)) = 1;
iceoceanmask(find(y==-2720225),find(x==205525)) = 1;
iceoceanmask(find(y==-2720225),find(x==205675)) = 1;
iceoceanmask(find(y==-880175),find(x==277525)) = 0;

%% associate basins to glaciers
disp('Associating basins to glaciers...');

% load x,y locations for basin outlets
% i previously manually added them to the Morlighem spreadsheet
[aa,bb] = xlsread('Morlighem2017_majorglaciers.xlsx');
morlighem_lat = aa(:,3);
morlighem_lon = aa(:,4);
addpath('~/Google Drive/matlab_scripts/');
[morlighem_x,morlighem_y] = latlon2utm(morlighem_lat,morlighem_lon);
xvals = aa(:,13);
yvals = aa(:,14);

% check if every xval and yval is on an ice-ocean boundary pixel
% (updated bedmachine dataset may have changed things)
iceoceancheck = NaN*xvals;
for i=1:length(xvals),
    if ~isempty(find(x==xvals(i))) & ~isempty(find(y==yvals(i))),
        xvals_ind = find(x==xvals(i));
        yvals_ind = find(y==yvals(i));
        iceoceancheck(i) = iceoceanmask(yvals_ind,xvals_ind);
    end
end
if sum(iceoceancheck==0) == 0,
    disp('Ice-ocean check passed');
else
    disp('Warning: some xval/yval is not on the ice-ocean boundary');
end

% optional figures to help fill out xvals and yvals in spreadsheet
% figure(); hold on;
% imagesc(x,y,m); axis xy equal;
% plot(morlighem_x,morlighem_y,'ko');
% for i=1:length(morlighem_x),
%     text(morlighem_x(i),morlighem_y(i),num2str(i),'fontsize',18);
% end
% plot(xvals,yvals,'ro');
% figure();
% imagesc(m); axis equal;

% group ice-ocean pixels into calving fronts
calvingfronts = bwconncomp(iceoceanmask);

% linear indices of all outlets
for i=1:length(outlet_x),
    outlet_x_ind(i) = find(x==outlet_x(i));
    outlet_y_ind(i) = find(y==outlet_y(i));
end
outlets_linear = sub2ind(size(iceoceanmask),outlet_y_ind,outlet_x_ind);

% load glaciers structure if it exists
if isfile('glaciers.mat'),
    load glaciers.mat
end

% number drainage basins according to morlighem glacier number
% and start glaciers structure
drainagebasins = NaN*DB;
for i=1:length(xvals),
    i
    if ~isnan(xvals(i)),
        icf = find(x==xvals(i));
        jcf = find(y==yvals(i));
        if ~isempty(icf) & ~isempty(jcf),
            linindex = sub2ind(size(iceoceanmask),jcf,icf);
            % find linindex in calvingfronts
            clearvars cfgroup;
            for j=1:calvingfronts.NumObjects,
                if ismember(linindex,calvingfronts.PixelIdxList{j}),
                    cfgroup = j;
                end
            end
            % find all outlets that draining over the calving front
            cfoutlets = find(ismember(outlets_linear,calvingfronts.PixelIdxList{cfgroup}));
            % get indices of all points in drainage basin
            inds = [];
            for k=1:length(cfoutlets),
                inds = [inds;find(DB==cfoutlets(k))];
            end
            % number pixels in drainage basin
            drainagebasins(inds) = i;
            % start glaciers structure (make singles to save space)
            glaciers(i).morlighem_number = single(i);
            glaciers(i).name = bb{2+i,1};
            % grounding line depth
            cfinds = calvingfronts.PixelIdxList{cfgroup};
            [gldepth,id] = min(b(cfinds));
            glaciers(i).gldepth = single(gldepth);
            glaciers(i).x = single(X(cfinds(id)));
            glaciers(i).y = single(Y(cfinds(id)));
            [glaciers(i).lat,glaciers(i).lon] = polarstereo_inv(glaciers(i).x,glaciers(i).y);
        end
    end
end

% optional loop to check effect of new bedmachine mask
% just loops through glaciers to check if they look ok
% plotmask = m;
% for i=1:length(xvals),
%     if ~isnan(xvals(i)),
%         icf = find(x==xvals(i));
%         jcf = find(y==yvals(i));
%         if ~isempty(icf) & ~isempty(jcf),
%             linindex = sub2ind(size(iceoceanmask),jcf,icf);
%             % find linindex in calvingfronts
%             clearvars cfgroup;
%             for j=1:calvingfronts.NumObjects,
%                 if ismember(linindex,calvingfronts.PixelIdxList{j}),
%                     cfgroup = j;
%                 end
%             end
%             plotmask(calvingfronts.PixelIdxList{cfgroup}) = 5;
%             xcf = X(calvingfronts.PixelIdxList{cfgroup});
%             ycf = Y(calvingfronts.PixelIdxList{cfgroup});
%             pad = 5e3;
%             xinds = find(x>min(xcf)-pad & x<max(xcf)+pad);
%             yinds = find(y>min(ycf)-pad & y<max(ycf)+pad);
%             figure(); hold on;
%             imagesc(x(xinds),y(yinds),plotmask(yinds,xinds));
%             colormap([parula(4);1,0,0]); caxis([-0.5 5.5]);
%             plot(morlighem_x(i),morlighem_y(i),'ko');
%             text(morlighem_x(i),morlighem_y(i),num2str(i),'fontsize',18);
%             plot(xvals(i),yvals(i),'m*');
%             waitforbuttonpress;
%             close;
%         end
%     end
% end

%% save outputs
disp('Saving...');

x = single(x); y = single(y); drainagebasins = single(drainagebasins);
save drainagebasins.mat drainagebasins x y
save glaciers.mat glaciers