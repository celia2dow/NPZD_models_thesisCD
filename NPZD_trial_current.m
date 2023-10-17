% Adding temperature dependence to model from:
%   Zooplankton Mortality and the Dynamical Behaviour of Plankton 
%   Population Models
%   Edwards. A.M., & Brindley. J. 1999
% or:
%   Oscillatory behaviour in a three-component plankton population model
%   Edwards. A.M., & Brindley. J. 1996
% with MLD dynamics similar to the seasonal forcing seen in Chapter 8 of:
%   A Rational Dynamical-Systems Approach to Plankton Population Modelling
%   Edwards. A.M. 1997
% using
%   1. Heaviside/sigmoidal function for turning mixing on and off given the
%   temperature
%   2. MLD-dependent rates for maximum P growth rate, sinking and mixing 
% The model consists of three or four coupled ordinary differential 
% equations, describing changes in the concentrations of nutrient (N), 
% phytoplankton (P), zooplankton (Z), and detritus (D) in a physically 
% homogeneous oceanic mixed layer. 

% Citation for "derivative" 
% Tamas Kis (2023). Numerical Differentiation of Data (derivative) 
% (https://github.com/tamaskis/derivative-MATLAB/releases/tag/v3.0.3), 
% GitHub. Retrieved March 6, 2023. 

% Citation for "gif"
% Chad A. (June 2017). Greene of the University of Texas 
% Institute for Geophysics (UTIG). 
% (https://au.mathworks.com/matlabcentral/fileexchange/63239-gif)
% MathWorks. Retrieved June 26, 2023.

close all 
% IMPORT Sea-Surface Temperatures (SSTs)
%load("2_years_deeper.mat","hawaii_soest_7e38_7a7b_afxhffc2")
load('hawaii_soest_7e38_7a7b_afe2_4b48_981e_16ea.mat')

% SCALES
lat_centres = hawaii_soest_7e38_7a7b_afxhffc2.latitude;
lat_h = lat_centres(2)-lat_centres(1);
params.num_lats = length(lat_centres); % params.num_lats = 1; %
long_centres = hawaii_soest_7e38_7a7b_afxhffc2.longitude;
long_h = long_centres(2)-long_centres(1);
params.num_longs = length(long_centres); % params.num_longs = 1; %
earth_SA = 510072000*1000^2;

% TOLERANCES
tol = 0.2; % tolerance in temperature change indicating the thermocline (0.2)
params.Tc = 21.9; %21.9; %21 critical temperature above which stratification occurs
% params.Pc = 0.025; % critical concentration below which Z graze on D

temps = hawaii_soest_7e38_7a7b_afxhffc2.water_temp(:,:,1:params.num_lats,1:params.num_longs);
temps = temps(1:180,:,:,:); %%%%%%%%%%%%%%%%%%%%%% DELETE
dims = size(temps);
num_days = dims(1); 
num_depths = dims(2);
v_north = ones(num_days,num_depths,params.num_lats,params.num_longs); % m/s
v_east = ones(num_days,num_depths,params.num_lats,params.num_longs); % m/s
params.num_sites = params.num_lats*params.num_longs;
params.days = 0:num_days-1;

% PAPER MODELLED
params.paper = 1996; % 1996 for quadratic (d), 1999 for linear (q)

% MODEL CHOICES: on == 1, off == 0
smoothing = 0; % smoothing of temperature data: 0 if off, span of smoothing if on
xth = 1; % take every xth data point when constructing a cubic spline
params.temp_vert_mix = 1; % temperature dependent vertical mixing
params.temp_func = @(T) (tanh((params.Tc-T)/2)+1)/2; % function for above
%params.temp_func = @(T) heaviside(params.Tc-T); 
params.detritus_layer = 0; % inclusion of detritus (D) compartment in model
params.change_mix_depth = 1; % use of Edwards seasonally forced MLD model
% Palatability of detritus to zooplankton (in [0,1], if p2 = 0 then no
% consumption of detritus by zooplankton)
p2= 0; p1 = 1-p2; params.omega = p2/p1; % p2=0.3333
% Which method is used to calculate the thermocline depth: (0) some input 
% function of MLD depths, MLD() (1) the change in temperature gradient by 
% tol, or (2) the drop in near surface (10m depth) temperature by tol (3)
% cubic spline the temperatures at all missing depths and then use method 2
thermocline_calc = 3;
% Including current in the model
params.current = 1;

% PARAMETERS (table 1)
% Maximum phytoplankton growth rate under optimal light conditions (0.5-2)
params.Pmax = 2; % 1/(m day) 0.2*12.5/2.58
% Related to the mazimum growth rate of P
params.a = 0.2; % 1/(m day)
%params.a = 2.58*params.Pmax/12.5; % 1/(m day)
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
params.vid_speed = 24; % frames per sec

% TIME AND LOCATION
year1 = 2013; % First year of imported data
month1 = 1; % 1: Jan, 2: Feb, ... , 12: Dec
params.day1 = 1; 
params.hemisphere = "N"; % N: northern, S: southern

months_label = ["Jan","Feb","Mar","Apr","May","Jun",...
    "Jul","Aug","Sep","Oct","Nov","Dec"];
season_label = ["Summer","Autumn","Winter","Spring"];
if params.hemisphere == "N"
    season_label = circshift(season_label,-2);
end
season_nums = ceil(mod(2:13,12)/3);
season_nums(season_nums==0) = 4;
params.month_season = cell(2,12);
for indx = 1:12
    params.month_season{1, indx} = months_label(indx);
    params.month_season{2, indx} = season_label(season_nums(indx));
end

% FOR EACH GRID SITE, CALCULATE THE THERMOCLINE, THE MLD, THE LATTICE SITE
% AREA, AND FIND THE CUBIC SPLINES FOR TEMPERATURE AND THE MLD AS FUNCTIONS
% OF TIME.
depth_indices = 1:num_depths; % array of indices for each depth 
MLD_indices = zeros(num_days,params.num_lats,params.num_longs); 
depths_LEVs = zeros(num_days,params.num_lats,params.num_longs); 
temps_av = zeros(num_days,params.num_lats,params.num_longs);
v_north_av = zeros(num_days,params.num_lats,params.num_longs);
v_east_av = zeros(num_days,params.num_lats,params.num_longs);
area_grid_site = zeros(params.num_lats,params.num_longs); 
params.volumes = zeros(params.num_lats,params.num_longs,num_days);
lats = zeros(params.num_lats,2);
longs = zeros(params.num_longs,2);

for i = 1:params.num_lats
    for j = 1:params.num_longs% Find MLD for each lattice site for each day
        for day = 1:num_days
            % Calculate the area of the lattice site
            lats(i,:) = [lat_centres(i)-lat_h; lat_centres(i)+lat_h];
            longs(j,:) = [long_centres(j)-long_h; long_centres(j)+long_h];
            % Lattice site area is surface area of earth * fraction of sphere
            % covered
            area_grid_site(i,j) = earth_SA * ...
                areaquad(lats(i,1),longs(j,1),lats(i,2),longs(j,2));

            % Note any temperatures that are missing (NaN) to be later
            % ignored
            temps_ijday = squeeze(temps(day,:,i,j));
            NANindices = find(isnan(temps_ijday)); % the indices of depths that are missing
            GOODindices = setdiff(depth_indices,NANindices); % the indices of depths that aren't missing
            num_GOODdepths = length(GOODindices);

            % Find MLD
            if thermocline_calc == 0 % some input function of MLD depths, MLD()
                [depth_meters, ~] = MLD2(day-1, params);
                difs = hawaii_soest_7e38_7a7b_afxhffc2.LEV-depth_meters;
                depth_index = max(find(difs<=0));
                depths_LEVs(day,i,j) = depth_meters;
            elseif thermocline_calc == 1 % the change in temperature gradient by tol
                temp_smooth = smooth(temps(day,GOODindices,i,j));
                temp_grad = abs(temp_smooth(2:num_GOODdepths) - temp_smooth(1:num_GOODdepths-1));
                depth_index = min(find(temp_grad>tol));
            elseif thermocline_calc == 2 % the drop in near surface (10m depth or shallower) temperature by tol
                NSD = 6; % near-surface-depth index for 10m
                while ~ismember(NSD, GOODindices) 
                    NSD = NSD-1;
                    if NSD == 0
                        % warning("Warning: No data exists for depth 0-10m for i=%d, j=%d, day=%d \n",i,j,day)
                        depth_index = [];
                        break
                    end
                end
                if NSD>0
                    %if NSD<6
                    %    fprintf("Note: i=%d, j=%d, day=%d, near surface depth taken at %d m \n",...
                    %        i,j,day,hawaii_soest_7e38_7a7b_afxhffc2.LEV(NSD));
                    %end
                    near_surf_temp = temps(day,NSD,i,j);
                    thresh_temp = near_surf_temp - tol; % Threshold temperature
                    depth_index = min(find(temps(day,:,i,j)<thresh_temp));
                end
            elseif thermocline_calc == 3 % Cubic spline the temps over depths
                NSD = 6; % near-surface-depth index for 10m
                while ~ismember(NSD, GOODindices) 
                    NSD = NSD-1;
                    if NSD == 0
                        % warning("Warning: No data exists for depth 0-10m for i=%d, j=%d, day=%d \n",i,j,day)
                        depth_index = [];
                        break
                    end
                end
                if NSD>0
                    %if NSD<6
                    %    fprintf("Note: i=%d, j=%d, day=%d, near surface depth taken at %d m \n",...
                    %        i,j,day,hawaii_soest_7e38_7a7b_afxhffc2.LEV(NSD));
                    %end
                    spline_depths = 10:500;
                    cs_temp_v_depth = spline(...
                        hawaii_soest_7e38_7a7b_afxhffc2.LEV, ...
                        [0 temps(day,:,i,j) 0], spline_depths);
                    near_surf_temp = temps(day,NSD,i,j);
                    thresh_temp = near_surf_temp - tol; % Threshold temperature
                    depth_meters = spline_depths(min(find(cs_temp_v_depth<thresh_temp)));
                    difs = hawaii_soest_7e38_7a7b_afxhffc2.LEV-depth_meters;
                    depth_index = max(find(difs<=0));
                    depths_LEVs(day,i,j) = depth_meters;
                end
            end
            
            % Record the MLD and the average temperature within this layer
            if depth_index
                % The index
                MLD_indices(day,i,j)=depth_index;
                 % The actual depth in meters
                if depths_LEVs(day,i,j) == 0
                    depths_LEVs(day,i,j) = hawaii_soest_7e38_7a7b_afxhffc2.LEV(depth_index);
                end
                % Find average temperature within the mixed layer (degrees
                 % Celsius)
                GOODindices_inMLD = GOODindices(GOODindices<=depth_index);
                if thermocline_calc == 0
                    temps_av(day,i,j) = temp_fixed(day,params);
                    temps(day,:,i,j) = temp_fixed(day,params);
                else
                    temps_av(day,i,j) = squeeze(mean(temps(day,GOODindices_inMLD,i,j),2));
                    v_north_av(day,i,j) = squeeze(mean(v_north(day,GOODindices_inMLD,i,j),2));
                    v_east_av(day,i,j) = squeeze(mean(v_north(day,GOODindices_inMLD,i,j),2));
                end
            else
                % If there are no data points for lattice site ij on this
                % day, ignore this entire case
                % fprintf("Note: no data for i=%d, j=%d, day=%d \n",i,j,day);
                MLD_indices(day,i,j)=NaN;
                depths_LEVs(day,i,j)=NaN;
                temps_av(day,i,j)=NaN;
            end
        end

        % Smooth temperature and MLD records if specified
        days_indices = 1:num_days;
        GOODindices_days = days_indices(~isnan(temps_av(:,i,j)));
        if smoothing
            temps_av(GOODindices_days,i,j) = ...
                smooth(temps_av(GOODindices_days,i,j),smoothing);
            depths_LEVs(GOODindices_days,i,j) = ...
                smooth(depths_LEVs(GOODindices_days,i,j),smoothing);
        end

        % Record volume using (smoothed) MLD
        params.volumes(i,j,GOODindices_days)= ...
            depths_LEVs(GOODindices_days,i,j) .* area_grid_site(i,j);

        % Find a cubic spline to approximate the temperature and MLD
        SPLINEindices = 1:xth:num_days;
        if params.num_lats == 1
            temps_av = temps_av';
            depths_LEVs = depths_LEVs';
            NANindices_days = days_indices(isnan(temps_av(:,i,j)));
            SPLINEindices = setdiff(SPLINEindices,NANindices_days);
            cs_temp = spline(SPLINEindices-1, [0 temps_av(SPLINEindices) 0]);
            cs_MLD = spline(SPLINEindices-1, [0 depths_LEVs(SPLINEindices) 0]);
        else
            NANindices_days = days_indices(isnan(temps_av(:,i,j)));
            SPLINEindices = setdiff(SPLINEindices,NANindices_days);
            cs_temp = spline(SPLINEindices-1, [0 squeeze(temps_av(SPLINEindices,i,j))' 0]);
            cs_MLD = spline(SPLINEindices-1, [0 squeeze(depths_LEVs(SPLINEindices,i,j))' 0]);
        end
        name_site = sprintf('i%i_j%i',i,j);
        params.temp_CS.(name_site) = cs_temp;
        params.MLD_CS.(name_site) = cs_MLD;
    end
end

% PLOT (SMOOTHED) TEMPERATURE AND MLD DATA AGAINST THE CUBIC SPLINE
figure(123)
plot_days = 0:0.1:num_days;
plot(plot_days, ppval(params.temp_CS.i1_j1, plot_days),'k-', ...
    GOODindices_days, temps_av(GOODindices_days), 'r*')
legend("Cubic spline", "Data")
title("Temperature")
xlabel("Time (days)")
ylabel("Temperature (^{\circ} C)")

figure(124)
plot(plot_days, ppval(params.MLD_CS.i1_j1, plot_days),'k-',...
    GOODindices_days, depths_LEVs(GOODindices_days), 'r*')
legend("Cubic spline", "Data")
title("Mixed layer depth")
xlabel("Time (days)")
ylabel("MLD (m)")

% CREATE ARRAYS FOR THE MONTHS (1 to 12) ASSOCIATED WITH EACH DAY IN THE
% SIMULATION
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
            new_day = old_day + months_days(mnth)-params.day1;
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
    record.soluts = zeros(params.num_lats,params.num_longs,4,num_days*100); 
else % NPZ model
    record.soluts = zeros(params.num_lats,params.num_longs,3,num_days*100); 
end

% PREPARE INITIAL CONDITION
% N 0.3
params.Narray = zeros(params.num_lats,params.num_longs);
N0array = 0.4 .* ones(params.num_lats,params.num_longs);
N0vect = reshape(N0array,[1,params.num_sites]);
% P 0.1
params.Parray = zeros(params.num_lats,params.num_longs);
P0array = 0.1 .* ones(params.num_lats,params.num_longs);
P0vect = reshape(P0array,[1,params.num_sites]);
% Z 0.2
params.Zarray = zeros(params.num_lats,params.num_longs);
Z0array = 0.05 .* ones(params.num_lats,params.num_longs);
Z0vect = reshape(Z0array,[1,params.num_sites]);
if params.detritus_layer
    % D
    params.Darray = zeros(params.num_lats,params.num_longs);
    D0array = 0.3 .* ones(params.num_lats,params.num_longs);
    D0vect = reshape(D0array,[1,params.num_sites]);
    y0 = [N0vect, P0vect, Z0vect, D0vect];
else
    y0 = [N0vect, P0vect, Z0vect];
end

% Solve for ODE concentrations
%options = odeset('NonNegative',1);
[t,y] = ode45(@(t,x) odefunc(t,x,params),[0 num_days-1],y0); % or ode15s, with options

record.times=t; 
Narray = y(:, 1 : params.num_sites); 
Parray = y(:, params.num_sites+1 : 2*params.num_sites);
Zarray = y(:, 2*params.num_sites+1 : 3*params.num_sites);
if params.detritus_layer
    Darray = y(:, 3*params.num_sites+1 : end);
end

for time = 1:length(t)
    % N
    record.soluts(:,:,1,time) = ...
        reshape(Narray(time,:), [params.num_lats params.num_longs]); 
    % P
    record.soluts(:,:,2,time) = ...
        reshape(Parray(time,:), [params.num_lats params.num_longs]); 
    % Z
    record.soluts(:,:,3,time) = ...
        reshape(Zarray(time,:), [params.num_lats params.num_longs]); 
    if params.detritus_layer
        % D
        record.soluts(:,:,4,time) = ...
            reshape(Darray(time,:), [params.num_lats params.num_longs]); 
    end
end

%%
% Go through each day in time and save the spatial distributuion of
% temperature, and phytoplankton
% Interpolate produced data to the days for which we have temperature data
if visual_switch
    % Prepare gif
    ocean_grid = figure(1);
    axis tight manual % ensures getframe() returns a consistent size
    %ocean_grid.Visible = 'off';
    gifpath = [folder_path '/simulation.gif'];

    % Interpolate phytoplankton data
    phyto = zeros(params.num_lats,params.num_longs,num_days);
    for i = 1:params.num_lats
        for j = 1:params.num_longs
            non0times = record.times(record.times(:)>0);
            time_lim = length(squeeze(non0times));
            phyto(params.num_lats+1-i,j,:) = interp1(squeeze(non0times),...
                squeeze(record.soluts(i,j,2,1:time_lim)),params.days);
        end
    end
    % Find limits
    lims.min_T = min(temps_av(temps_av>0),[],'all');
    lims.max_T = max(temps_av,[],'all');
    lims.min_P = min(phyto,[],'all');
    lims.max_P = max(phyto,[],'all');
    
    % Preallocate
    plts = cell(1,2);
    axs = cell(1,2);
    for day = params.days
        day = day+1;
        % Save each figure as a frame in the gif
        [axs,plts] = draw_frame(params,phyto,temps_av,lims,day,axs,plts);
        % Save the gif frame.
        frame = getframe(ocean_grid);
        [A,map] = rgb2ind(frame2im(frame),256);
        if day == 1
            imwrite(A,map,gifpath,'gif',"LoopCount",Inf,"DelayTime",1/params.vid_speed);
        else
            imwrite(A,map,gifpath,'gif',"WriteMode","append","DelayTime",1/params.vid_speed);
        end
    end
end
%%
months = record.times./(365.25/12);
fig2 = figure(2);
season1 = ceil(month1/3);
if params.hemisphere == "N"
    strt = 3;
else
    strt = 1;
end
year_add=0;

months_final = months(months>0);
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

% NPZ(D) plot
for i = 1:params.num_lats
    for j = 1:params.num_longs
        condit = months>0;
        non0times = record.times(record.times(:)>0);
        time_lim = length(squeeze(non0times));
        vols_ij = interp1(params.days, squeeze(params.volumes(i,j,:)),...
            squeeze(non0times));
        months_ij = months(condit);
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
        ylim([0,0.6])
        %%%% COMMENT THIS OUT
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
    % We requires that solutions are ordered: long1 lat1, long1 lat2, ..., 
    % long1 latEnd, long2 lat1, long2 lat2, ..., long2 latEnd, ...
    lattice_size = [params.num_lats params.num_longs];

    % Convert state variables from vector-form to matrix-form
    % N
    Nvec = x(1 : params.num_sites);
    N = reshape(Nvec, lattice_size);
    % P
    Pvec = x(params.num_sites+1 : 2*params.num_sites);
    P = reshape(Pvec, lattice_size);
    % Z
    Zvec = x(2*params.num_sites+1 : 3*params.num_sites);
    Z = reshape(Zvec, lattice_size);
    if params.detritus_layer
        % D
        Dvec = x(3*params.num_sites+1 : end);
        D = reshape(Dvec, lattice_size);
    end
    
    % Rename params
    a=params.a; b=params.b; c=params.c; d=params.d; e=params.e; k=params.k; 
    q=params.q; s=params.s; N_0=params.N_0; D_0=params.D_0; 
    alpha=params.alpha; beta=params.beta; gamma=params.gamma; 
    lambda=params.lambda; mu=params.mu; r=params.r; omega = params.omega; 
    phi = params.phi; psi=params.psi;
    
    dNdt = zeros(params.num_lats, params.num_longs);
    dPdt = zeros(params.num_lats, params.num_longs);
    dZdt = zeros(params.num_lats, params.num_longs);
    dDdt = zeros(params.num_lats, params.num_longs);
    dxdt = zeros(length(x),1);
    for i = 1:params.num_lats
        for j = 1:params.num_longs
            name_site = sprintf('i%i_j%i',i,j);

            % Interpolate temperature (T), mixed-layer depth (LEV), and
            % gradient of depth (h)
            % [T, LEV, h] = interpFromCS(params,i,j,t); 
            T = CSfunction(params.temp_CS.(name_site),t);
            LEV = CSfunction(params.MLD_CS.(name_site),t);
            h = CSgradient(params.MLD_CS.(name_site),t);
            
            
            if params.temp_vert_mix && params.change_mix_depth == 1
                % a, k, s and psi need to be updated with the new depth
                % Switch for temperature-dependent vertical mixing/sinking
                a = params.a*12.5/LEV;
                k = params.k*params.temp_func(T)*12.5/LEV;
                s = params.s*params.temp_func(T)*12.5/LEV;
                psi = params.psi*params.temp_func(T)*12.5/LEV;
            elseif params.change_mix_depth == 1
                % Change with mixed layer depth
                % a, k, s and psi need to be updated with the new depth
                a = params.a*12.5/LEV;
                k = params.k*12.5/LEV;
                s = params.s*12.5/LEV;
                psi= params.psi*12.5/LEV;
            elseif params.temp_vert_mix
                % Switch for temperature-dependent vertical mixing/sinking
                k = params.k*params.temp_func(T);
                s = params.s*params.temp_func(T);
                psi = params.psi*params.temp_func(T);
            end

            % ODE SYSTEM equations 1-3
            if params.paper == 1999
                h_Z = q*Z(i,j);
            elseif params.paper == 1996
                h_Z = d*Z(i,j)^2;
            end

            % Switch for detritus layer inclusion
            if params.detritus_layer % NPZD model
                dNdt(i,j) = -N(i,j)*P(i,j)*a/((e+N(i,j))*(b+c*P(i,j))) + ...
                    beta*lambda*(P(i,j)^2+omega*D(i,j)^2)*Z(i,j)/...
                    (mu^2+P(i,j)^2+omega*D(i,j)^2) + gamma*h_Z + ...
                    phi*D(i,j) + k*(N_0-N(i,j));
                dPdt(i,j) = N(i,j)*P(i,j)*a/((e+N(i,j))*(b+c*P(i,j))) - ...
                    r*P(i,j) - lambda*P(i,j)^2*Z(i,j)/...
                    (mu^2+P(i,j)^2+omega*D(i,j)^2) - (s+k)*P(i,j);
                dZdt(i,j) = alpha*lambda*(P(i,j)^2+omega*D(i,j)^2)*Z(i,j)/...
                    (mu^2+P(i,j)^2+omega*D(i,j)^2) - h_Z;
                dDdt(i,j) = r*P(i,j) + ...
                    lambda*((1-alpha-beta)*P(i,j)^2-(alpha+beta)*omega*D(i,j)^2)*Z(i,j)/...
                    (mu^2+P(i,j)^2+omega*D(i,j)^2) - (phi+psi+k)*D(i,j);
            else % NPZ model
                dNdt(i,j) = -N(i,j)*P(i,j)*a/((e+N(i,j))*(b+c*P(i,j))) + ...
                    r*P(i,j) + beta*lambda*P(i,j)^2*Z(i,j)/(mu^2+P(i,j)^2) ...
                    + gamma*h_Z + k*(N_0-N(i,j));
                dPdt(i,j) = N(i,j)*P(i,j)*a/((e+N(i,j))*(b+c*P(i,j))) - ...
                    r*P(i,j) - lambda*P(i,j)^2*Z(i,j)/(mu^2+P(i,j)^2) - ...
                    (s+k)*P(i,j);
                dZdt(i,j) = alpha*lambda*P(i,j)^2*Z(i,j)/(mu^2+P(i,j)^2) - ...
                    h_Z;
            end
            
            % Change with mixed layer depth
            if params.change_mix_depth == 1
                dNdt(i,j) = dNdt(i,j) + max(h,0) * (N_0-N(i,j))/LEV;
                dPdt(i,j) = dPdt(i,j) - max(h,0) * P(i,j)/LEV;
                dZdt(i,j) = dZdt(i,j) - h*Z(i,j)/LEV;
                if params.detritus_layer
                    dDdt(i,j) = dDdt(i,j) - max(h,0) * D(i,j)/LEV;
                end
            end
        end
    end
    
    % Convert state variables from matrix-form to vector-form
    % N
    dNdt_vect = reshape(dNdt,[1,params.num_sites]);
    % P
    dPdt_vect = reshape(dPdt,[1,params.num_sites]);
    % Z
    dZdt_vect = reshape(dZdt,[1,params.num_sites]);
    if params.detritus_layer
        % D
        dDdt_vect = reshape(dDdt,[1,params.num_sites]);
        dxdt = [dNdt_vect, dPdt_vect, dZdt_vect, dDdt_vect]';
    else
        dxdt = [dNdt_vect, dPdt_vect, dZdt_vect]';
    end
end

function [depth, h] = MLD(t,params)
    if params.hemisphere == "S"
        t = params.day1+t+182.625;
    end
    t_day = mod(t,365.25);
    if 0<=t_day && t_day<80 % Winter
        depth = 93.575 + 0.705*t_day;
        h = 0.705;
    elseif 80<=t_day && t_day<130 % Spring
        depth = 150 - 2.75*(t_day-80);
        h = -2.75;
    elseif 130<=t_day && t_day<250 % Summer
        depth = 12.5; 
        h = 0;
    elseif 250<=t_day && t_day < 365.25 % Autumn
        depth = 12.5 + 0.705*(t_day-250);
        h = 0.705;
    end
end

function [depth, h] = MLD2(t,params)
    if params.hemisphere == "S"
        t = params.day1+t+182.625;
    end
    depth = 68.75*sin(2*pi*(t+10.385)/365)+81.25;
    h = 68.75*2*pi*cos(2*pi*(t+10.385)/365)/365;
end

function temp = temp_fixed(t,params)
    if params.hemisphere == "S"
        t = params.day1+t+182.625;
    end
    temp = -4*cos(2*pi*(t-50)/365)+21.5;
end

function [temp, MLD, gradMLD] = interpFromCS(params,i,j,day)
    % Interpolate the temperature, MLD and gradient of the MLD function on
    % a specific date at a specific lattice site (i,j)
    name_site = sprintf('i%i_j%i',i,j);
    CSfuncTemp = params.temp_CS.(name_site);
    CSfuncMLD = params.MLD_CS.(name_site);
    temp = ppval(CSfuncTemp,day);
    MLD = ppval(CSfuncMLD,day);
    MLD_neighbourhood = ppval(CSfuncMLD,[day-0.01, day, day+0.01]);
    gradMLDs = gradient(MLD_neighbourhood,0.01);
    gradMLD = gradMLDs(2);
end

function y = CSfunction(cubic_spline, x)
    % Takes in a cubic spline structure and a value x at which to evaluate
    % the cubic spline and returns the value of the cubic spline y
    difs = cubic_spline.breaks - x; 
    indx = min(sum(difs<=0),cubic_spline.pieces);
    x0 = cubic_spline.breaks(indx);
    c = cubic_spline.coefs;
    y = c(indx,1)*(x-x0)^3 + c(indx,2)*(x-x0)^2 + c(indx,3)*(x-x0) + c(indx,4);
end

function dydx = CSgradient(cubic_spline, x)
    % Takes in a cubic spline structure and a value x at which to evaluate
    % the cubic spline and returns the value of the cubic spline y
    difs = cubic_spline.breaks - x;
    indx = min(sum(difs<=0), cubic_spline.pieces);
    x0 = cubic_spline.breaks(indx);
    c = cubic_spline.coefs;
    dydx = 3*c(indx,1)*(x-x0)^2 + 2*c(indx,2)*(x-x0) + c(indx,3);
end