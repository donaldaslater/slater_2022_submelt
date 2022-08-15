% script to assign RACMO runoff to tidewater glaciers
clear; close all;

% load required datasets
load drainagebasins.mat
load glaciers.mat

% get list of RACMO files
fn = dir('~/Documents/RACMO');
fn0 = fn(1).folder;
i0 = 1;
for i=1:length(fn),
    fnn = fn(i).name;
    fnn = strsplit(fnn,'.');
    if strcmp(fnn(end),'nc'),
        f{i0} = [fn(i).folder,'/',fn(i).name];
        i0 = i0+1;
    end
end

% use first RACMO file to load RACMO grid
% note x is called lon and y is called lat in the netcdf files
xr = ncread(f{1},'lon');
yr = ncread(f{1},'lat');
[Xr,Yr] = meshgrid(xr,yr); Xr=Xr'; Yr=Yr';

% interpolate drainage basins onto RACMO grid
[X,Y] = meshgrid(x,y);
racmobasins = interp2(X,Y,drainagebasins,Xr,Yr,'nearest');

% create structure of linear indices for each basin
for i=1:length(glaciers),
    b(i).inds = find(racmobasins==i);
end

% initialise glaciers runoff field
for j=1:length(glaciers),
    glaciers(j).runoff.RACMO.Q = [];
    glaciers(j).runoff.RACMO.t = [];
end

% initialise monthly time vector
m = ([1:12]-0.5)/12;

% now loop over RACMO files and read into glacier structure
for i=1:length(f),
    
    Q = ncread(f{i},'runoffcorr');
    Q = reshape(Q,size(Q,1)*size(Q,2),size(Q,3));
    yr = strsplit(f{i},'.');
    yr = str2num(yr{2})
    
    for j=1:length(glaciers),
        if ~isempty(b(j).inds),
            glaciers(j).runoff.RACMO.Q = [glaciers(j).runoff.RACMO.Q,...
                1000*1000*nansum(Q(b(j).inds,:))./(1000*86400*eomday(yr,[1:size(Q,2)]))];
            glaciers(j).runoff.RACMO.t = [glaciers(j).runoff.RACMO.t,...
                single(yr+m(1:size(Q,2)))];
        end
    end
    
end

% save structure
save glaciers.mat glaciers