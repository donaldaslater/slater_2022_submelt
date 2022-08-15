% script to finalise glaciers structure
clear; close all;
load glaciers.mat

% create mean ocean thermal forcing by taking average over products
product = {'ORAS5','EN4','ASTE','CHORE'};

% thermal forcing constants
l1 = -5.73e-2;
l2 = 8.32e-2;
l3 = 7.53e-4;

% melt rate parameterisation constants
ks = 0.142;
a = 0.31;
b = 1.19;

% common time period
t = single([1979+1/24:1/12:2019]);

% calculate thermal forcing
for i=1:length(glaciers),   
    for j=1:length(product),    
        if ~isnan(glaciers(i).ocean.(product{j}).z),
            glaciers(i).ocean.(product{j}).TF_GL = glaciers(i).ocean.(product{j}).potentialT_GL - ...
                (l1*glaciers(i).ocean.(product{j}).practicalS_GL + l2 + l3*glaciers(i).gldepth);
            % check not less than 0
            glaciers(i).ocean.(product{j}).TF_GL(find(glaciers(i).ocean.(product{j}).TF_GL<0)) = 0;
        else
            glaciers(i).ocean.(product{j}).TF_GL = NaN;
        end
    end
end

% thermal forcing over common time period
for i=1:length(glaciers),
    for j=1:length(product),
        if ~isnan(glaciers(i).ocean.ORAS5.z),
            glaciers(i).ocean.(product{j}).tc = t;
            glaciers(i).ocean.(product{j}).TFc = interp1(glaciers(i).ocean.(product{j}).t,...
                glaciers(i).ocean.(product{j}).TF_GL,t);
        end
    end
end

% extend CHORE and ASTE by mean of ORAS5 and EN4
for i=1:length(glaciers),   
    if sum(isnan(glaciers(i).ocean.ORAS5.z))==0,       
        % mean of ORAS5 and EN4
        TF1 = glaciers(i).ocean.ORAS5.TFc;
        TF2 = glaciers(i).ocean.EN4.TFc;
        TFmean = mean([TF1;TF2]);

        % difference with CHORE during overlapping period
        TFbias1 = mean(glaciers(i).ocean.CHORE.TF_GL-interp1(t,TFmean,glaciers(i).ocean.CHORE.t));

        % difference with ASTE during overlapping period
        TFbias2 = mean(glaciers(i).ocean.ASTE.TF_GL-interp1(t,TFmean,glaciers(i).ocean.ASTE.t));

        % extend CHORE by mean of ORAS5 and EN4
        TFextend = interp1(glaciers(i).ocean.CHORE.t,glaciers(i).ocean.CHORE.TF_GL,t);
        inds = find(isnan(TFextend));
        TFextend(inds) = TFmean(inds);
        glaciers(i).ocean.CHORE.tc = t;
        glaciers(i).ocean.CHORE.TFc = TFextend + TFbias1;

        % extend ASTE by mean of ORAS5 and EN4
        TFextend = interp1(glaciers(i).ocean.ASTE.t,glaciers(i).ocean.ASTE.TF_GL,t);
        inds = find(isnan(TFextend));
        TFextend(inds) = TFmean(inds);
        glaciers(i).ocean.ASTE.tc = t;
        glaciers(i).ocean.ASTE.TFc = TFextend + TFbias2;   
    end    
end

% mean thermal forcing over common time period
for i=1:length(glaciers),   
    for j=1:length(product),
        if ~isnan(glaciers(i).ocean.(product{j}).z),
            TF(j,:) = interp1(glaciers(i).ocean.(product{j}).t,glaciers(i).ocean.(product{j}).TF_GL,t);
        else
            TF(j,:) = NaN;
        end
    end
    glaciers(i).ocean.average.t = single(t);
    glaciers(i).ocean.average.TF = nanmean(TF);
end

% submarine melt rate based on ocean average
for i=1:length(glaciers),
    if ~isempty(glaciers(i).runoff.RACMO.Q) & ~isnan(glaciers(i).ocean.average.TF(1)),
        Q = interp1(glaciers(i).runoff.RACMO.t,glaciers(i).runoff.RACMO.Q,t)+glaciers(i).basalmelt;
        TF = interp1(glaciers(i).ocean.average.t,glaciers(i).ocean.average.TF,t);
        glaciers(i).submelt.t = single(t);
        glaciers(i).submelt.m = ks*Q.^a.*TF.^b;
    else
        glaciers(i).submelt.t = NaN;
        glaciers(i).submelt.m = NaN;
    end
end

% submarine melt rate based on individual ocean products
product = {'ORAS5','ASTE','EN4','CHORE'};
for j=1:length(product),
    for i=1:length(glaciers),
        if ~isempty(glaciers(i).runoff.RACMO.Q) & ~isnan(glaciers(i).ocean.average.TF(1)),
            Q = interp1(glaciers(i).runoff.RACMO.t,glaciers(i).runoff.RACMO.Q,t)+glaciers(i).basalmelt;
            TF = interp1(glaciers(i).ocean.(product{j}).t,glaciers(i).ocean.(product{j}).TF_GL,t);
            glaciers(i).submelt.t = single(t);
            glaciers(i).submelt.(['m_',product{j}]) = ks*Q.^a.*TF.^b;
        else
            glaciers(i).submelt.t = NaN;
            glaciers(i).submelt.(['m_',product{j}]) = NaN;
        end
    end
end

save glaciers.mat glaciers

%% filter for final set of glaciers

load glaciers.mat
disp(['Initially have ',num2str(length(glaciers)),' glaciers']);

% first remove glaciers that are barely below sea level
g2r = [];
for i=1:length(glaciers),
    if glaciers(i).gldepth>-50, g2r = [g2r,i]; end
end
disp(['Removing ',num2str(length(g2r)),' glaciers that do not have gldepth<-50']);
glaciers(g2r) = [];

% then remove glaciers without ocean properties
g2r = [];
for i=1:length(glaciers),
    ed = glaciers(i).ocean.ORAS5.effdepth;
    if isnan(ed) | ed>-50, g2r = [g2r,i]; end
end
disp(['Removing ',num2str(length(g2r)),' glaciers that do not have effdepth<-25']);
glaciers(g2r) = [];

% then remove glaciers with only tiny runoff
g2r = [];
for i=1:length(glaciers),
    inds = find(glaciers(i).runoff.RACMO.t>=1979 & glaciers(i).runoff.RACMO.t<=2019);
    if mean(glaciers(i).runoff.RACMO.Q(inds))<2.5, g2r = [g2r,i]; end
end
disp(['Removing ',num2str(length(g2r)),' glaciers with minimal runoff']);
glaciers(g2r) = [];

% final dataset
disp(['Finally have ',num2str(length(glaciers)),' glaciers']);

% add ismip6 sector number
load ice_ocean_sectors.mat
for i=1:length(regions),
    inds = find(inpolygon([glaciers.x],[glaciers.y],regions(i).ice.x,regions(i).ice.y));
    for j=1:length(inds),
        glaciers(inds(j)).ismip6_sector = i;
    end
end

% add mean Q and TF over 1979-2018
for i=1:length(glaciers),
    Qinds = find(glaciers(i).runoff.RACMO.t>=1979 & glaciers(i).runoff.RACMO.t<=2019);
    TFinds = find(glaciers(i).ocean.average.t>=1979 & glaciers(i).ocean.average.t<=2019);
    glaciers(i).Qmean = mean(glaciers(i).runoff.RACMO.Q(Qinds));
    glaciers(i).TFmean = mean(glaciers(i).ocean.average.TF(TFinds));
end

% rename and save
twglaciers = glaciers;
save twglaciers.mat twglaciers

    