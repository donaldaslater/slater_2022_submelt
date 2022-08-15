% function to investigate glacier response to melt rate
function [t,modeloutput] = glacier_response(ks,a,b);

%% set-up

% load dataset
load ~/'OneDrive - University of Edinburgh'/data/natgeos/data/twglaciers.mat

% melt rate time series by glacier
t0 = [1979,1993];
t = twglaciers(1).submelt.t;
inds0 = find(t>=t0(1) & t<=t0(end));

% sector indices
sec(1).inds = find([twglaciers.ismip6_sector]<=3); sec(1).name = 'South';
sec(2).inds = find([twglaciers.ismip6_sector]==4); sec(2).name = 'Central-west';
sec(3).inds = find([twglaciers.ismip6_sector]==6); sec(3).name = 'North-west';
sec(4).inds = find([twglaciers.ismip6_sector]==5); sec(4).name = 'North-east';
sec(5).inds = find([twglaciers.ismip6_sector]==7); sec(5).name = 'North';
sec(6).inds = [1:length(twglaciers)]; sec(6).name = 'All Greenland';

% get timeseries per glacier
for i=1:length(twglaciers),
    % Q
    Q = interp1(twglaciers(i).runoff.RACMO.t,twglaciers(i).runoff.RACMO.Q,t)+twglaciers(i).basalmelt;
    % TF
    TF = interp1(twglaciers(i).ocean.average.t,twglaciers(i).ocean.average.TF,t);
    % melt rate
    m(i,:) = ks*Q.^a.*TF.^b;
    m0(i) = mean(m(i,inds0));
    % melt rate Q varies
    TFt = TF(inds0);
    TF0 = NaN*TF;
    for k=1:12, TF0([k:12:end]) = mean(TFt([k:12:end])); end
    mQ(i,:) = ks*Q.^a.*TF0.^b;
    mQ0(i) = mean(mQ(i,inds0));
    % melt rate TF varies
    Qt = Q(inds0);
    Q0 = NaN*Q;
    for k=1:12, Q0([k:12:end]) = mean(Qt([k:12:end])); end
    mTF(i,:) = ks*Q0.^a.*TF.^b;
    mTF0(i) = mean(mTF(i,inds0));
end

%% model parameters

% constants
n = 3;
mp = 1/3;
alpha = 2*n+1;
gamma = n;
beta = 1+(mp+n+3)/(mp+1);
g = 9.81;
rho_w = 1027;
rho_i = 917;
lambda = rho_w/rho_i;
secs_in_year = 365*86400;

% time stepping
dt = 0.02*secs_in_year;
tm = [floor(t(1)):dt/secs_in_year:ceil(t(end))];
% tm = tm(1:end-1);

% choices
% calving multiplier
eta = [0,1];
% initial velocity to melt factor
vm_init = [1,5];
% initial frontal ice thickness (m)
h_init = [200,700];
% initial inland ice thickness (m)
H_init = [1000,3000];
% initial glacier length (m)
L_init = [50e3,200e3];
% bed slope (m/m)
dbdx = -[0.5e-3,2e-3];

% params to save for methods plot
i1save = 1;
i2save = 1;
i3save = 1;
i4save = 1;
i5save = 1;
i6save = 1;

% make linear parameter structure
ii0 = 1;
for i1=1:length(eta), for i2=1:length(vm_init), for i3=1:length(h_init),
for i4=1:length(H_init), for i5=1:length(L_init), for i6=1:length(dbdx),  
    par(ii0).eta = eta(i1);
    par(ii0).vm_init = vm_init(i2);
    par(ii0).h_init = h_init(i3);
    par(ii0).H_init = H_init(i4);
    par(ii0).L_init = L_init(i5);
    par(ii0).dbdx = dbdx(i6);
    ii0 = ii0+1;
    if i1==i1save & i2==i2save & i3==i3save & i4==i4save & ...
            i5==i5save & i6==i6save,
        iisave = ii0;
    end
end; end; end; end; end; end;

%% run model

% loop over glaciers
for k=1:length(twglaciers),
    
    disp(['Simulating glacier ',num2str(k),'/',num2str(length(twglaciers)),': ',twglaciers(k).name]);

    % loop over parameter combinations
    for ii=1:length(par),

        % loop over melt rate scenarios
        for jj=1:3,

            % melt rates in m/s
            % both vary
            if jj==1,
                melt = interp1(t,m(k,:),tm,'nearest','extrap')/86400;
                melt0 = m0(k)/86400;
            % Q varies
            elseif jj==2,
                melt = interp1(t,mQ(k,:),tm,'nearest','extrap')/86400;
                melt0 = mQ0(k)/86400;
            % TF varies
            elseif jj==3,
                melt = interp1(t,mTF(k,:),tm,'nearest','extrap')/86400;
                melt0 = mTF0(k)/86400;
            end
            
            % needed quantities
            Q_init = par(ii).vm_init*melt0*(1+par(ii).eta)*par(ii).h_init;
            omega = (Q_init - melt0*(1+par(ii).eta)*par(ii).h_init)/par(ii).h_init^beta;
            P_init = Q_init/par(ii).L_init;
            nu = Q_init*par(ii).L_init^gamma/par(ii).H_init^alpha;
            x = [0:1:1000]*10^3;
            b = -(1/lambda)*par(ii).h_init-abs(par(ii).dbdx)*(x-par(ii).L_init);

            % initial conditions
            xg = par(ii).L_init;
            h = par(ii).H_init;
            
            % do time-stepping
            for ti = 1:length(tm), 

                % bed topo and ice thickness at terminus
                bnow = interp1(x,b,xg,'linear',NaN);    
                hg = -lambda*bnow;

                % ice fluxes
                Q = nu * h^alpha / xg^gamma;   
                Qg = omega*(hg^beta) + melt(ti)*(1+par(ii).eta)*hg; 

                % tendencies
                dh_dt = P_init - (Qg/xg) - (h/(xg*hg))*(Q-Qg);
                dxg_dt = (Q-Qg)/hg;

                % step forwards in time
                h = h + dh_dt*dt;
                xg = xg + dxg_dt*dt;  

                % store outputs
                sol(k).case(jj).par(ii).L(ti) = xg;
                sol(k).case(jj).par(ii).H(ti) = h;
                sol(k).case(jj).par(ii).iceflux(ti) = rho_i*secs_in_year*5000*Q/10^12;
                sol(k).case(jj).par(ii).dVdt(ti) = rho_i*5000*(P_init*xg - Q)*secs_in_year/10^12;

            end
            
            % save for methods plot
            if ii==iisave & k==1,
                
                eg_t = tm;
                eg_melt(jj,:) = melt;
                eg_L(jj,:) = sol(k).case(jj).par(ii).L;
                eg_iceflux(jj,:) = sol(k).case(jj).par(ii).iceflux;
                eg_dVdt(jj,:) = sol(k).case(jj).par(ii).dVdt;
                
            end                

        end
        
    end

end

%% post-processing

disp('Post-processing...');

% loop over glaciers
for k=1:length(sol),    
    % loop over cases
    for jj=1:3,
        
        % loop over parameters
        for ii=1:length(par),
            
            % dynamic sea level contribution
            SL = cumtrapz(tm,sol(k).case(jj).par(ii).dVdt);
            % smooth out seasonality and reset to 0
            sol(k).case(jj).par(ii).SL = smooth(SL,round(1/(dt/secs_in_year)));
            sol(k).case(jj).par(ii).SL = sol(k).case(jj).par(ii).SL - sol(k).case(jj).par(ii).SL(1);
            % normalise to both varying case (jj=1)
            SLnorm(k,jj,ii,:) = sol(k).case(jj).par(ii).SL/sol(k).case(1).par(ii).SL(end);
            
            % smooth terminus position and set start to 0
            L = smooth(sol(k).case(jj).par(ii).L,4*round(1/(dt/secs_in_year)),'lowess');
            sol(k).case(jj).par(ii).L = L-L(1);
            % normalise to both varying case
            Lnorm(k,jj,ii,:) = -sol(k).case(jj).par(ii).L/sol(k).case(1).par(ii).L(end);
            
            % smooth ice flux
            sol(k).case(jj).par(ii).iceflux = smooth(sol(k).case(jj).par(ii).iceflux,4*round(1/(dt/secs_in_year)),'lowess');
            % normalise to start
            Dnorm(k,jj,ii,:) = sol(k).case(jj).par(ii).iceflux/sol(k).case(jj).par(ii).iceflux(1);
            
        end
        
    end

end

% get sector-wide statistics
for k=1:length(sec),    
    % loop over cases
    for jj=1:3,
        
        % restrict to sector in question
        SLsec = squeeze(SLnorm(sec(k).inds,jj,:,:));
        SLsec = reshape(SLsec,size(SLsec,1)*size(SLsec,2),size(SLsec,3));
        Lsec = squeeze(Lnorm(sec(k).inds,jj,:,:));
        Lsec = reshape(Lsec,size(Lsec,1)*size(Lsec,2),size(Lsec,3));
        Dsec = squeeze(Dnorm(sec(k).inds,jj,:,:));
        Dsec = reshape(Dsec,size(Dsec,1)*size(Dsec,2),size(Dsec,3));
        
        % get stats for plotting
        modeloutput(k).case(jj).meanSL = mean(SLsec);
        modeloutput(k).case(jj).maxSL = max(SLsec);
        modeloutput(k).case(jj).minSL = min(SLsec);
        modeloutput(k).case(jj).stdSL = std(SLsec);
        modeloutput(k).case(jj).SL50 = prctile(SLsec,50);
        modeloutput(k).case(jj).SL25 = prctile(SLsec,25);
        modeloutput(k).case(jj).SL75 = prctile(SLsec,75);
        
        modeloutput(k).case(jj).meanL = mean(Lsec);
        modeloutput(k).case(jj).maxL = max(Lsec);
        modeloutput(k).case(jj).minL = min(Lsec);
        modeloutput(k).case(jj).stdL = std(Lsec);
        modeloutput(k).case(jj).L50 = prctile(Lsec,50);
        modeloutput(k).case(jj).L25 = prctile(Lsec,25);
        modeloutput(k).case(jj).L75 = prctile(Lsec,75);
        
        modeloutput(k).case(jj).meanD = mean(Dsec);
        modeloutput(k).case(jj).maxD = max(Dsec);
        modeloutput(k).case(jj).minD = min(Dsec);
        modeloutput(k).case(jj).stdD = std(Dsec); 
        modeloutput(k).case(jj).D50 = prctile(Dsec,50);
        modeloutput(k).case(jj).D25 = prctile(Dsec,25);
        modeloutput(k).case(jj).D75 = prctile(Dsec,75);

    end
end

% save time axis and sector names
for k=1:length(sec),
    modeloutput(k).name = sec(k).name;
end
t = tm;

% save so that don't have to re-run every time
save modeloutput.mat t modeloutput
save modeloutput_methods.mat eg_t eg_melt eg_L eg_iceflux eg_dVdt sol

