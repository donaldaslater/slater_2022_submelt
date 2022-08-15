% script to estimate basal melting for glaciers structure
% using dataset from Karlsson 2021
% https://www.nature.com/articles/s41467-021-23739-z
clear; close all;

% load dataset
karlsson = load('../basal_melt/basalmelt.mat');

% load drainage basins
load drainagebasins.mat

% interpolate basal melt rates onto BM3 grid
[Xbm3,Ybm3] = meshgrid(x,y);
Fk = scatteredInterpolant(karlsson.x(:),karlsson.y(:),karlsson.totmelt(:),'linear','none');
Mk = Fk(double(Xbm3),double(Ybm3));
% convert rates from m/yr to m/s
Mk = Mk/(365*86400);

% load glaciers structure
load glaciers.mat

% loop over glaciers summing basal melt for the basins
for i=1:length(glaciers),
    % get linear indices of basin
    inds = find(drainagebasins==i);    
    % sum basal melt over basin
    glaciers(i).basalmelt = 150*150*nansum(Mk(inds));  
end

% save structure
save glaciers.mat glaciers

