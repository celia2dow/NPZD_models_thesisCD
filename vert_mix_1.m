% Adding temperature dependence to model from:
%   Zooplankton Mortality and the Dynamical Behaviour of Plankton 
%   Population Models
%   Edwards. A.M., & Brindley. J. 1999
% or
%   Oscillatory behaviour in a three-component plankton population model
%   Edwards. A.M., & Brindley. J. 1996
% using
%   1. Heaviside function for turning mixing on and off given the
%   temperature
% The model consists of three or four coupled ordinary differential 
% equations, describing changes in the concentrations of nutrient (N), 
% phytoplankton (P), zooplankton (Z), and detritus (D) in a physically 
% homogeneous oceanic mixed layer. 

% Citation for "derivative" 
% Tamas Kis (2023). Numerical Differentiation of Data (derivative) 
% (https://github.com/tamaskis/derivative-MATLAB/releases/tag/v3.0.3), 
% GitHub. Retrieved March 6, 2023. 

close all 
% IMPORT Sea-Surface Temperatures (SSTs)
load("2_years_deeper.mat","hawaii_soest_7e38_7a7b_afxhffc2")

% SCALES
lat_centres = hawaii_soest_7e38_7a7b_afxhffc2.latitude;
lat_h = lat_centres(2)-lat_centres(1);
num_lats = 1; %num_lats = length(lat_centres); % 
long_centres = hawaii_soest_7e38_7a7b_afxhffc2.longitude;
long_h = long_centres(2)-long_centres(1);
num_longs = 1; %num_longs = length(long_centres); % 
earth_SA = 510072000*1000^2;

% TOLERANCES
tol = 0.5; % tolerance in temperature gradient indicating the thermocline
params.Tc = 20; %21.9; % critical temperature above which stratification occurs
params.smooth = 1; % to switch off, =0
params.Pc = 0.025; % critical concentration below which Z graze on D

temps = hawaii_soest_7e38_7a7b_afxhffc2.water_temp(:,:,1:num_lats,1:num_longs);
%temps = temps(1:180,:,:,:); %%%%%%%%%%%%%%%%%%%%%% DELETE
dims = size(temps);
num_days = dims(1); 
num_depths = dims(2);
params.days = 0:num_days-1;
depths = zeros(num_days,num_lats,num_longs); % preallocate space
depths_LEVs = zeros(num_days,num_lats,num_longs); % preallocate space
temps_av = zeros(num_days,num_lats,num_longs); % preallocate space
for i = 1:num_lats
    for j = 1:num_longs
        for depth = 1:num_depths
            % For any temperatures that are missing (NaN), replace them 
            % with the average temperature of adjacent temperatures
            temps_ijdepth = squeeze(temps(:,depth,i,j));
            indices = find(isnan(temps_ijdepth));
            for index = indices
                if index == 1
                    temps_ijdepth(index) = temps_ijdepth(index+1);
                elseif index == num_days
                    temps_ijdepth(index) = temps_ijdepth(index-1);
                else
                    temps_ijdepth(index) = (temps_ijdepth(index-1)+temps_ijdepth(index+1))./2;
                end
            end
            temps(:,depth,i,j)=temps_ijdepth;
        end
        % Figure out depth of mixed layer
        for day = 1:num_days
            temp_smooth = smooth(temps(day,:,i,j));
            temp_grad = abs(temp_smooth(2:num_depths) - temp_smooth(1:num_depths-1));
            depth = min(find(temp_grad>tol));
            % The row
            depths(day,i,j)=depth; 
            % The actual depth in meters
            depths_LEVs(day,i,j)=hawaii_soest_7e38_7a7b_afxhffc2.LEV(depth); 
            % Find average temperature within the mixed layer (degrees
            % Celsius)
            temps_av(day,i,j)=squeeze(mean(temps(day,1:depth,i,j),2));
        end
        if params.smooth
            temps_av(:,i,j)= smooth(temps_av(:,i,j));
            depths_LEVs(:,i,j)=smooth(depths_LEVs(:,i,j));
        end
    end
end

% PAPER MODELLED
params.paper = 1996; % 1996 for quadratic (d)

% MODEL CHOICES: on == 1, off == 0
params.temp_vert_mix = 0;
params.temp_func = @(T) (tanh((params.Tc-T)/2)+1)/2; 
%params.temp_func = @(T) heaviside(params.Tc-T); 
params.detritus_layer = 0;
params.change_mix_depth = 1;
params.depth_dep_growth = 0;
% Palatability of detritus to zooplankton (in [0,1], if p2 = 0 then no
% consumption of detritus by zooplankton)
p2= 0; p1 = 1-p2; params.omega = p2/p1; % p2=0.3333

% PARAMETERS table 1
% Related to the mazimum growth rate of P
params.a = 0.2; % 1/(m day)
% Light attenuation by water
params.b = 0.2; % 1/m
% P self-shading coefficient
params.c = 0.4; % m^2/(g C)
% Higher predation on Z (for 1996 paper)
params.d = 1; %1; % m^3/(g C day)
% D concentration below mixed layer
params.D_0 = 0.2; % g C/(m^3)
% Half-saturation constant for N uptake
params.e = 0.03; %0.05; % g C/(m^3)
% Cross-thermocline exchange rate
params.k = 0.05; % 1/day
% Maximum phytoplankton growth rate under optimal light conditions (0.5-2)
params.Pmax = 2; % 1/(m day)
% Higher predation on Z (for 1999 paper)
params.q = 0.075; %0.11; % 1/day
% P respiration rate
params.r = 0.15; % 1/day
% P sinking loss rate
params.s = 0.04; % 1/day
% N concentration below mixed layer
params.N_0 = 0.6; % g C/(m^3)
% Z growth efficiency
params.alpha = 0.25;
% Z excretion fraction
params.beta = 0.33;
% Regeneration of Z predation excretion
params.gamma = 0.5;
% Maximum Z grazing rate
params.lambda = 0.6; % 1/day
% Z grazing half-saturation coefficient 
params.mu = 0.035; % g C/(m^3)
% D remineralisation rate
params.phi = 0.1; % 1/day
% D sinking loss rate
params.psi = 0.08; % 1/day

% MOVIE
visual_switch = 0; % 1=on, 0=off
params.vid_speed = 4; % frames per sec

% TIME AND LOCATION
year1 = 2013; % First year of imported data
month1 = 1; % 1: Jan, 2: Feb, ... , 12: Dec
day1 = 1; 
hemisphere = "N"; % N: northern, S: southern

months_label = ["Jan","Feb","Mar","Apr","May","Jun",...
    "Jul","Aug","Sep","Oct","Nov","Dec"];
season_label = ["Summer","Autumn","Winter","Spring"];
if hemisphere == "N"
    season_label = circshift(season_label,-2);
end
season_nums = ceil(mod(2:13,12)/3);
season_nums(season_nums==0) = 4;
params.month_season = cell(2,12);
for indx = 1:12
    params.month_season{1, indx} = months_label(indx);
    params.month_season{2, indx} = season_label(season_nums(indx));
end
%%
params.months_per_day = zeros(1,num_days);
months_days = [31,0,31,30,31,30,31,31,30,31,30,31];
months_nums = circshift(1:12,1-month1);
num_years = ceil(num_days/366);
old_day = 1;
for yr = 1:num_years
    if mod(year1 + yr,4)==0 && mod(year1 + yr,400)==0
        months_days(2)=29;
    elseif mod(year1 + yr,4)==0 && mod(year1 + yr,100)~=0
        months_days(2)=29;
    else
        months_days(2)=28;
    end
    for indx = 1:12
        mnth = months_nums(indx);
        if indx == 1
            new_day = old_day + months_days(mnth)-day1;
        else 
            new_day = min(old_day + months_days(mnth)-1, num_days);
        end
        params.months_per_day(old_day:new_day) = mnth;
        old_day = new_day+1;
        if new_day == num_days
            break
        end
    end
end
params.months_per_day(params.months_per_day==0)=month1;
%%

% FOLDER PATH
date_time = num2str(fix(clock));
folder_name = ['r' num2str(params.r) ...
    '_paper' num2str(params.paper) ...
    '_' date_time];
folder_name = folder_name(find(~isspace(folder_name)));
folder_path = [pwd '/' date '/' folder_name];
% The folder paths used agrees with a macOS system and may need to be
% flipped to '\' for a Windows system
if ~exist(folder_path, 'dir')
    mkdir(folder_path)
end

% PREALLOCATE SPACE
if params.detritus_layer % NPZD model
    record.soluts = zeros(num_lats,num_longs,4,num_days*100); 
else % NPZ model
    record.soluts = zeros(num_lats,num_longs,3,num_days*100); 
end
record.times = zeros(num_lats,num_longs,num_days*100);

lats = zeros(num_lats,2);
longs = zeros(num_longs,2);
params.volumes = zeros(num_lats,num_longs,num_days);
params.Vderiv = zeros(num_lats,num_longs,num_days); %params.Vderiv = @(t) sin(t);
for i = 1:num_lats
    for j = 1:num_longs
        y0 = [0.3,0.125,0.05];%[0.3,0.1,0.2]; % N0 P0 Z0, [0.3,0.125,0.05] 
        if params.detritus_layer
            y0 = [0.3,0.125,0.05,0.3];%[0.3,0.1,0.2,0.3]; % [0.3,0.125,0.05,0.3]
        end
        params.temps_ij = temps_av(:,i,j);
        params.LEV_ij = depths_LEVs(:,i,j);
        lats(i,:) = [lat_centres(i)-lat_h; lat_centres(i)+lat_h];
        longs(j,:) = [long_centres(j)-long_h; long_centres(j)+long_h];
        % Lattice site area is surface area of earth * fraction of sphere
        % covered
        area_grid_site = earth_SA * ...
            areaquad(lats(i,1),longs(j,1),lats(i,2),longs(j,2));
        params.volumes(i,j,:)=params.LEV_ij.*area_grid_site;
        params.vols_ij = squeeze(params.volumes(i,j,:));
        params.Vderiv(i,j,:) = derivative(params.days, params.vols_ij);
        params.Vderiv_ij = params.Vderiv(i,j,:);
        if params.change_mix_depth
            % Scale initial condition and nutrient rich N0 to be in mass
            % not concentration
            y0 = y0 .* params.vols_ij(1);
        end
        params.y0=y0;
        % Solve for ODE concentrations
        options = odeset('NonNegative',1);
        [t,y] = ode45(@(t,x) odefunc(t,x,params),[0 num_days-1],y0, options); % or ode15s
        record.soluts(i,j,1,1:length(t))=y(:,1); % Nutrient (mixed layer)
        record.soluts(i,j,2,1:length(t))=y(:,2); % Phytoplankton
        record.soluts(i,j,3,1:length(t))=y(:,3); % Zooplankton 
        if params.detritus_layer
            record.soluts(i,j,4,1:length(t))=y(:,4); % Detritus
        end
        record.times(i,j,1:length(t))=t; % time
    end
end
%%
% Go through each day in time and save the spatial distributuion of
% temperature, and phytoplankton
% Interpolate produced data to the days for which we have temperature data
if visual_switch
    % Prepare movie.
    ocean_grid = figure(1);
    axis tight manual % ensures getframe() returns a consistent size
    ocean_movie(num_days) = struct('cdata',[],'colormap',[]);
    ocean_grid.Visible = 'off';

    % Prepare gif 
    gifpath = [folder_path '/simulation.gif'];

    % Interpolate phytoplankton data
    phyto = zeros(num_lats,num_longs,num_days);
    for i = 1:num_lats
        for j = 1:num_longs
            non0times = record.times(i,j,record.times(i,j,:)>0);
            time_lim = length(squeeze(non0times));
            phyto(i,j,:) = interp1(squeeze(non0times),...
                squeeze(record.soluts(i,j,2,1:time_lim)),params.days);
            if params.change_mix_depth
                % convert from mass to concentration
                phyto(i,j,:) = squeeze(phyto(i,j,:)) ./ ...
                    squeeze(params.volumes(i,j,:));
            end
        end
    end
    % Find limits
    lims.min_T = min(temps_av(temps_av>0),[],'all');
    lims.max_T = max(temps_av,[],'all');
    lims.min_P = min(phyto,[],'all');
    lims.max_P = max(phyto,[],'all');
    
    for day = params.days
        day = day+1;
        clf
        % Save each figure as a frame in the movie
        draw_frame(params,phyto,temps_av,lims,day)
        frameN = getframe(ocean_grid);
        ocean_movie(day) = frameN;
        % Save the gif frame.
        imageN = frame2im(frameN);
        [imind,cm]=rgb2ind(imageN,256);
        if day == 1
            imwrite(imind,cm,gifpath,'gif','Loopcount',1)    
        else
            imwrite(imind,cm,gifpath,'gif','WriteMode','append')  
        end
    end

    % Play movie
    ocean_grid.Visible = 'on';
    % movie(ocean_grid, ocean_movie, 1, params.vid_speed);
end

months = record.times./(365.25/12);
fig2 = figure(2);
season1 = ceil(month1/3);
if hemisphere == "N"
    strt = 3;
else
    strt = 1;
end
year_add=0;
months_final = months(1,1,months(1,1,:)>0);
months_final = months_final(end);
tcks = 0:3:months_final;
num_ticks = length(tcks);
tck_seas = cell(1,num_ticks);
for tck = season1:num_ticks+season1
    indx = mod(tck,4);
    indx(indx == 0) = 4;
    if tck ==season1
        tck_seas{tck} = num2str(year1) + " " + season_label{indx};
    elseif indx == 1
        year_add=year_add+1;
        tck_seas{tck} = num2str(year1 + year_add) + " " + season_label{indx};
    else
        tck_seas{tck}=season_label{indx};
    end
end
for i = 1:num_lats
    for j = 1:num_longs
        % NPZ(D) plot
        condit = months(i,j,:)>0;
        non0times = record.times(i,j,record.times(i,j,:)>0);
        time_lim = length(squeeze(non0times));
        vols_ij = interp1(params.days, squeeze(params.volumes(i,j,:)),...
            squeeze(non0times));
        if params.change_mix_depth
            % Convert back to concentration from mass
            if params.detritus_layer
                lst = 4;
            else
                lst = 3;
            end
            for ind = 1:lst
                record.soluts(i,j,ind,2:time_lim+1) = ...
                    squeeze(record.soluts(i,j,ind,2:time_lim+1)) ./ vols_ij;
            end
        end
        months_ij = squeeze(months(i,j,condit));
        plot(months_ij,squeeze(record.soluts(i,j,1,condit)), 'b');
        hold on
        xlabel("Time (months)")
        ylabel("Concentration (g C/{m^3})")
        plot(months_ij,squeeze(record.soluts(i,j,2,condit)), 'k');
        plot(months_ij,squeeze(record.soluts(i,j,3,condit)), 'r');
        if params.detritus_layer
            plot(months_ij,squeeze(record.soluts(i,j,4,condit)), 'g');
        end
        legend('Mixed layer nutrient (N)', 'Phytoplankton (P)', ...
                'Zooplankton (Z)')
        xticks(tcks);
        xticklabels(tck_seas);
        xtickangle(90);
        %xticklabels({'x = 0','x = 5','x = 10'});
        if params.detritus_layer
            legend('Mixed layer nutrient (N)', 'Phytoplankton (P)', ...
                'Zooplankton (Z)', 'Detritus (D)')
        end
    end
end
hold off
savefig(fig2, [folder_path '/State_variables_over_time'], 'compact')
saveas(fig2, [folder_path '/State_variables_over_time'], 'png')

% Edited ODE system
function dxdt = odefunc(t,x,params)
    % Interpolate temperature, depth, volume
    Vderiv = interp1(params.days, squeeze(params.Vderiv_ij), t); %Vderiv = params.Vderiv(t);
    T = interp1(params.days, params.temps_ij, t);
    LEV = interp1(params.days, params.LEV_ij, t);
    V = interp1(params.days, squeeze(params.vols_ij), t);
    % Rename params
    b=params.b; c=params.c; d=params.d; e=params.e; q=params.q;
    N_0=params.N_0; D_0=params.D_0; alpha=params.alpha;
    beta=params.beta; gamma=params.gamma; lambda=params.lambda;
    mu=params.mu; r=params.r; omega = params.omega; phi = params.phi;
    

    % Rename state variables (state variable values remain unchanged)
    N = max(x(1),0);
    P = max(x(2),0);
    Z = max(x(3),0);
    if params.detritus_layer
        D = max(x(4),0);
    end

    % Switch for temperature-dependent vertical mixing/sinking
    if params.temp_vert_mix
        k = params.k*params.temp_func(T);
        s = params.s*params.temp_func(T);
        psi = params.psi*params.temp_func(T);
    else
        k=params.k;
        s=params.s;
        psi = params.psi;
    end
    
    % If letting the maximum growth rate averaged over the mixed layer
    if params.depth_dep_growth
        a=2.58*params.Pmax/LEV; 
    else
        a=params.a; 
    end

    % Change with mixed layer depth
    if params.change_mix_depth
        % Then N, P, Z, D have been given as masses not concentrations
        N = N/V;%params.vols_ij(1);
        P = P/V;%params.vols_ij(1);
        Z = Z/V;%params.vols_ij(1);
        if params.detritus_layer
            D = D/V;%params.vols_ij(1);
        end
        % k, s and psi need to be updated with the new depth
        k = k*12.5/LEV;%params.LEV_ij(1)/LEV; %12.5/LEV;
        s = s*12.5/LEV;%params.LEV_ij(1)/LEV;
        psi = psi*12.5/LEV;%params.LEV_ij(1)/LEV;
    end

    % ODE SYSTEM equations 1-3
    if params.paper == 1999
        h_Z = q*Z;
    elseif params.paper == 1996
        h_Z = d*Z^2;
    end

    % Switch for detritus layer inclusion
    if params.detritus_layer % NPZD model
        dNdt = -N*P*a/((e+N)*(b+c*P)) + beta*lambda*(P^2+omega*D^2)*Z/ ...
            (mu^2+P^2+omega*D^2) + gamma*h_Z + phi*D + k*(N_0-N);
        dPdt = N*P*a/((e+N)*(b+c*P)) - r*P - lambda*P^2*Z/ ...
            (mu^2+P^2+omega*D^2) - (s+k)*P;
        dZdt = alpha*lambda*(P^2+omega*D^2)*Z/(mu^2+P^2+omega*D^2) - h_Z;
        dDdt = r*P + lambda*((1-alpha-beta)*P^2-(alpha+beta)*omega*D^2)*Z/...
            (mu^2+P^2+omega*D^2) - (phi+psi+k)*D;
    else % NPZ model
        dNdt = -N*P*a/((e+N)*(b+c*P)) + r*P + beta*lambda*P^2*Z/...
            (mu^2+P^2) + gamma*h_Z + k*(N_0-N);
        dPdt = N*P*a/((e+N)*(b+c*P)) - r*P - lambda*P^2*Z/(mu^2+P^2) ...
             - (s+k)*P;
        dZdt = alpha*lambda*P^2*Z/(mu^2+P^2) - h_Z;
    end

    if params.change_mix_depth
        dNdt = V * dNdt + Vderiv * ...
            (N_0*heaviside(Vderiv) + N*heaviside(-Vderiv));
        dPdt = V * dPdt + Vderiv*P*heaviside(-Vderiv);
        dZdt = V * dZdt;
        if params.detritus_layer
            dDdt = V * dDdt + Vderiv * ...
                (D_0*heaviside(Vderiv) + D*heaviside(-Vderiv));
        end
    end
    
    dxdt = zeros(3,1);
    dxdt(1) = dNdt;
    dxdt(2) = dPdt;
    dxdt(3) = dZdt;
    if params.detritus_layer
        dxdt(4) = dDdt;
    end
    
end