clear; close all;

load coastline.mat
coastlon = coastlon-360;

for i=1:29,
    
    if i<10, numstr = ['000',num2str(i)];
    else numstr = ['00',num2str(i)];
    end
    
    lon(i,:,:) = ncread(['~/Documents/ASTE/ETAN.',numstr,'.nc'],'lon');
    lat(i,:,:) = ncread(['~/Documents/ASTE/ETAN.',numstr,'.nc'],'lat');
    etai = ncread(['~/Documents/ASTE/ETAN.',numstr,'.nc'],'ETAN');
    eta(i,:,:) = squeeze(etai(:,:,1));

end

tile = NaN*lon;
for i=1:29,
    tile(i,:,:) = i;
end

figure(); hold on;
scatter(lon(:),lat(:),10,tile(:),'filled'); colorbar;
plot(coastlon,coastlat,'k');

% get tiles for greenland coastline
f = scatteredInterpolant(lon(:),lat(:),tile(:),'nearest','none');
tiles_needed = f(coastlon,coastlat);
tiles_needed = unique(tiles_needed);

% double check
figure(); hold on;
for i=1:length(tiles_needed),
    loni = squeeze(lon(tiles_needed(i),:,:));
    lati = squeeze(lat(tiles_needed(i),:,:));
    tilei = squeeze(tile(tiles_needed(i),:,:));
    scatter(loni(:),lati(:),10,tilei(:),'filled');
    colorbar;
end
plot(coastlon,coastlat,'k');
