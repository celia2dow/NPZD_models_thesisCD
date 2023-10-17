% INITIAL CONDITIONS MATCH THOSE OF THE EDWARDS & BRINDLEY PAPERS

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

close all 
clear
% IMPORT Sea-Surface Temperatures (SSTs)
load("2_years_deeper.mat","hawaii_soest_7e38_7a7b_afxhffc2")

% SCALES
lat_centres = hawaii_soest_7e38_7a7b_afxhffc2.latitude;
lat_h = lat_centres(2)-lat_centres(1);
lats = length(lat_centres); % lats = 1; %
long_centres = hawaii_soest_7e38_7a7b_afxhffc2.longitude;
long_h = long_centres(2)-long_centres(1);
longs = length(long_centres); % longs = 1; %

% TOLERANCES
tol = 0.5; % tolerance in temperature gradient indicating the thermocline
params.Tc = 21.9; % critical temperature above which stratification occurs
params.Pc = 0.025; % critical concentration below which Z graze on D

temps = hawaii_soest_7e38_7a7b_afxhffc2.water_temp(:,:,1:lats,1:longs);
temps = temps(1:180,:,:,:); %%%%%%%%%%%%%%%%%%%%%% DELETE
dims = size(temps);
num_days = dims(1);
num_depths = dims(2);
params.days = 1:num_days;
depths = zeros(num_days,lats,longs); % preallocate space
depths_LEVs = zeros(num_days,lats,longs); % preallocate space
temps_av = zeros(num_days,lats,longs); % preallocate space
for i = 1:lats
    for j = 1:longs
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
    end
end


% PAPER MODELLED
params.paper = 1996; % 1996 for quadratic dZ^2, or 1999 for linear qZ

% MODEL CHOICES: on == 1, off == 0
params.temp_vert_mix = 1;
params.detritus_layer = 0;
params.change_mix_depth = 0;
% Palatability of detritus to zooplankton (in [0,1], if p2 = 0 then no
% consumption of detritus by zooplankton)
p2= 1/3; p1 = 1-p2; params.omega = p2/p1;

% PARAMETERS table 1
% Related to the mazimum growth rate of P
params.a = 0.2; % 1/(m day)
% Light attenuation by water
params.b = 0.2; % 1/m
% P self-shading coefficient
params.c = 0.4; % m^2/(g C)
% Higher predation on Z (for 1996 paper)
params.d = 1; %1.5; % m^3/(g C day)
% D concentration below mixed layer
params.D_0 = 0.2; % g C/(m^3)
% Half-saturation constant for N uptake
params.e = 0.03; %0.03; % g C/(m^3)
% Cross-thermocline exchange rate
params.k = 0.05; % 1/day
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
visual_switch = 1; % 1=on, 0=off
params.vid_speed = 4; % frames per sec

% TIME AND LOCATION
year1 = 2013; % First year of imported data
month1 = 1; % 1: Jan, 2: Feb, ... , 12: Dec
day1 = 15; 
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


% FOLDER PATH
date_time = num2str(fix(clock));
folder_name = ['r' num2str(params.r) ...
    '_paper' num2str(params.paper) ...
    '_' date_time];
folder_name = folder_name(find(~isspace(folder_name)));
folder_path = [pwd '/Comparison_to_old_papers_' date '/' folder_name];
% The folder paths used agrees with a macOS system and may need to be
% flipped to '\' for a Windows system
if ~exist(folder_path, 'dir')
    mkdir(folder_path)
end

% PREALLOCATE SPACE
if params.detritus_layer % NPZD model
    record.soluts = zeros(lats,longs,4,num_days*100); 
else % NPZ model
    record.soluts = zeros(lats,longs,3,num_days*100); 
end
record.times = zeros(lats,longs,num_days*100);
% Initial condition
y0 = [0.4,0.1,0.05]; % N0 P0 Z0, [0.3,0.025,0.05] [0.4 0.1 0.05] 
if params.detritus_layer
    y0 = [0.4,0.1,0.05,0.08]; % N0 P0 Z0 D0
end

for i = 1:lats
    for j = 1:longs
        params.temps_ij = temps_av(:,i,j);
        params.LEV_ij = depths_LEVs(:,i,j);
        lats_i = [lat_centres(i)-lat_h,lat_centres(i)+lat_h];
        long_j = [long_centres(j)-long_h,long_centres(j)+long_h];
        params.area_grid_site = areaquad(lats_i(1),long_j(1),lats_i(2),long_j(2));
        params.tOLD = 0;
        params.VNEW = params.area_grid_site * params.LEV_ij(1);
        if params.change_mix_depth
            % Scale initial condition and nutrient rich N0 to be in mass
            % not concentration
            y0 = y0 * params.VNEW;
        end
        % Solve for ODE concentrations
        [t,y] = ode15s(@odefunc,[0 num_days],y0, [], params); % or ode45
        record.soluts(i,j,1,1:length(t))=y(:,1); % Nutrient (mixed layer)
        record.soluts(i,j,2,1:length(t))=y(:,2); % Phytoplankton
        record.soluts(i,j,3,1:length(t))=y(:,3); % Zooplankton 
        if params.detritus_layer
            record.soluts(i,j,4,1:length(t))=y(:,4); % Detritus
        end
        record.times(i,j,1:length(t))=t; % time
    end
end

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
    phyto = zeros(lats,longs,num_days);
    for i = 1:lats
        for j = 1:longs
            non0times = record.times(i,j,record.times(i,j,:)>0);
            time_lim = length(squeeze(non0times));
            phyto(i,j,:) = interp1(squeeze(non0times),...
                squeeze(record.soluts(i,j,2,1:time_lim)),params.days);
        end
    end
    % Find limits
    lims.min_T = min(temps_av(temps_av>0),[],'all');
    lims.max_T = max(temps_av,[],'all');
    lims.min_P = min(phyto,[],'all');
    lims.max_P = max(phyto,[],'all');
    
    for day = params.days
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
    % ocean_grid.Visible = 'on';
    % movie(ocean_grid, ocean_movie, 1, params.vid_speed);
end

months = record.times./(365.25/12);
fig2 = figure(2);
for i = 1:lats
    for j = 1:longs
        % NPZ(D) plot
        condit = months(i,j,:)>0;
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
    % Time-step size
    delta_t = t - params.tOLD;
    params.tOLD = t;

    % Interpolate temperature, depth, volume
    T = interp1(params.days, params.temps_ij, t);
    LEV = interp1(params.days, params.LEV_ij, t);
    V_old = params.VNEW;
    params.VNEW = params.area_grid_site * LEV;
    V_new = params.VNEW;
    params.VMID = (V_old+V_new)/2;
    V_mid = params.VMID;
    if delta_t == 0
        fprintf('\n delta_t=0')
        fprintf('\n delta_V=%d',V_new-V_old)
    end

    % Rename params
    a=params.a; b=params.b; c=params.c; d=params.d; D_0=params.D_0; 
    e=params.e; q=params.q; s=params.s; N_0=params.N_0; alpha=params.alpha;
    beta=params.beta; gamma=params.gamma; lambda=params.lambda;
    mu=params.mu; r=params.r; omega = params.omega; phi = params.phi;
    psi = params.psi;

    % Rename state variables (state variable values remain unchanged)
    N = x(1);
    P = x(2);
    Z = x(3);
    if params.detritus_layer
        D = x(4);
    end

    % Switch for temperature-dependent vertical mixing
    if params.temp_vert_mix
        k=params.k*heaviside(params.Tc-T);
    else
        k=params.k;
    end

    if params.change_mix_depth
        % Then N, P, Z, D have been given as masses not concentrations
        N = N/V_old;
        P = P/V_old;
        Z = Z/V_old;
        if params.detritus_layer
            D = D/V_old;
        end
        % k, s and psi need to be updated with the new depth
        k = k*params.LEV_ij(1)/LEV;
        s = s*params.LEV_ij(1)/LEV;
        psi = psi*params.LEV_ij(1)/LEV;
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
        dNdt = -N*a*P/((e+N)*(b+c*P)) + r*P + beta*lambda*P^2*Z/ ...
            (mu^2+P^2) + gamma*h_Z+ k*(N_0-N);
        dPdt = N*a*P/((e+N)*(b+c*P)) - r*P - lambda*P^2*Z/(mu^2+P^2) ...
            - (s+k)*P;
        dZdt = alpha*lambda*P^2*Z/(mu^2+P^2) - h_Z;
    end

    if params.change_mix_depth
        dNdt = V_mid * dNdt + ((V_new-V_old)/delta_t) * ...
            (N_0*heaviside(V_new-V_old) + N*heaviside(V_old-V_new));
        dPdt = V_mid * dPdt + ((V_new-V_old)/delta_t) * ...
            P*heaviside(V_old-V_new);
        dZdt = V_mid * dZdt;
        if params.detritus_layer
            dDdt = V_mid * dDdt + ((V_new-V_old)/delta_t) * ...
                (D_0*heaviside(V_new-V_old) + D*heaviside(V_old-V_new));
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