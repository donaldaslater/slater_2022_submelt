% script to get ocean properties to glaciers
clear; close all;

oceanprods = {'ORAS5','EN4','ASTE','CHORE'};

% load glaciers structure to this point
load glaciers.mat

% loop over ocean products
for kk=1:length(oceanprods),

    product = oceanprods{kk}

    % load ocean output
    if strcmp(product,'ORAS5'),
        load ocean_ORAS5_glaciers.mat
    elseif strcmp(product,'EN4'),
        load ocean_EN4_glaciers.mat
    elseif strcmp(product,'ASTE'),
        load ocean_ASTE_glaciers.mat
    elseif strcmp(product,'CHORE'),
        load ocean_CHORE_glaciers.mat
    end

    % loop over glaciers and assign ocean properties
    for ii=1:length(glaciers),
    
        if glaciers(ii).ocean.(product).effdepth<-25,
            
            % get ocean properties for this glacier
            ind = find(morlighem_n==glaciers(ii).morlighem_number);        
            Tii = squeeze(potentialT(ind,:,:));
            Sii = squeeze(practicalS(ind,:,:));
            % create depth vector for ocean properties for this glacier
            z_inds = find(z>glaciers(ii).ocean.(product).effdepth);
            glaciers(ii).ocean.(product).z = single([z(z_inds)';glaciers(ii).ocean.(product).effdepth;glaciers(ii).gldepth]);
    
            for j=1:length(t),
                
                % get properties at effective depth
                interpinds = [z_inds';z_inds(end)+1];
                effdepth_potentialT = interp1(z(interpinds),squeeze(Tii(interpinds,j)),glaciers(ii).ocean.(product).effdepth,'linear',NaN);
                effdepth_practicalS = interp1(z(interpinds),squeeze(Sii(interpinds,j)),glaciers(ii).ocean.(product).effdepth,'linear',NaN);
                % create profile including extrapolation below eff depth
                glaciers(ii).ocean.(product).potentialT(:,j) = single([Tii(z_inds,j);effdepth_potentialT;effdepth_potentialT]);
                glaciers(ii).ocean.(product).practicalS(:,j) = single([Sii(z_inds,j);effdepth_practicalS;effdepth_practicalS]);
                % value at grounding line only
                glaciers(ii).ocean.(product).potentialT_GL(j) = single(glaciers(ii).ocean.(product).potentialT(end,j));
                glaciers(ii).ocean.(product).practicalS_GL(j) = single(glaciers(ii).ocean.(product).practicalS(end,j));
    
            end
            
            glaciers(ii).ocean.(product).t = single(t);
            
        else
            
            glaciers(ii).ocean.(product).z = NaN;
            glaciers(ii).ocean.(product).potentialT = NaN;
            glaciers(ii).ocean.(product).practicalS = NaN;
            glaciers(ii).ocean.(product).t = NaN;
            glaciers(ii).ocean.(product).potentialT_GL = NaN;
            glaciers(ii).ocean.(product).practicalS_GL = NaN;    
        
        end
    
    end

end

save glaciers.mat glaciers