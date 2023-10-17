% Adding temperature dependence and advection to models from:
%   Zooplankton Mortality and the Dynamical Behaviour of Plankton 
%   Population Models
%   Edwards. A.M., & Brindley. J. 1999
% and:
%   Oscillatory behaviour in a three-component plankton population model
%   Edwards. A.M., & Brindley. J. 1996
% with MLD dynamics similar to the seasohnal forcing seen in Chapter 8 of:
%   A Rational Dynamical-Systems Approach to Plankton Population Modelling
%   Edwards. A.M. 1997
% using:
%   1. Heaviside/sigmoidal function for turning mixing on and off given the
%   temperature
%   2. MLD-dependent rates for maximum P growth rate, sinking and mixing
%   3. Finite difference approximation of current
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

% Citation for "@odeprog" and "@odeabort"
% Tim Franklin (2023). ODE Progress Bar and Interrupt 
% (https://www.mathworks.com/matlabcentral/fileexchange/9904-ode-progress-bar-and-interrupt), 
% MATLAB Central File Exchange. Retrieved September 10, 2023. 

% Citation for "multigradient"
% Laurens R Krol (2023). multigradient: custom gradient colormap 
% (https://github.com/lrkrol/multigradient), 
% GitHub. Retrieved September 10, 2023. 

%clear
%close all 

% IMPORT DATA
%load("2_years_deeper.mat","hawaii_soest_7e38_7a7b_afxhffc2")
load("canary_islands_at_collection_site_2011.mat",..."20_depths_of_current.mat",..."antarctic_peninsula_200m_010115_301215.mat",... "antarctic_peninsula_150m_010110_010111.mat", ...  "20_depths_of_current.mat",...
    "hawaii_soest_7e38_7a7b_afxhffc2") % "antarctic_peninsula_150m_010109_010111.mat",
data_file_temps = hawaii_soest_7e38_7a7b_afxhffc2;  % Sea-Surface Temperatures (SSTs)
%load("hawaii_soest_7e38_7a7b_afe2_bc82_1a36_5f92")
data_file_flow = hawaii_soest_7e38_7a7b_afxhffc2;   % Flow Field


% MODEL CHOICES: on == 1, off == 0
tol = 0.1;              % Tolerance in temperature change indicating the thermocline (0.2)
params.Tc = 21.9;       % Critical temperature above which stratification occurs % 21.9; %21 
params.Pc = 0.025;      % Critical concentration below which Z graze on D (unused)
%params.repeat_num = 4;
params.paper = 1996;                                    % Paper modelled: 1996 for quadratic mortality term (d), 1999 for linear mortality term (q)
%smoothing = 0;                                          % Smoothing of temperature data: 0 if off, span of smoothing if on (unused)
xth = 1;                                                % Take every xth data point when constructing a cubic spline (unused)
%params.temp_vert_mix = 0;                               % Temperature dependent vertical mixing
%params.temp_func = @(T) (tanh((params.Tc-T)/2)+1)/2;    % Function for above: alternative @(T) heaviside(params.Tc-T)
%params.detritus_layer = 0;                              % Inclusion of detritus (D) compartment in model
%params.change_mix_depth = 1;                            % Use of Edwards seasonally forced MLD model as base
%thermocline_calc = 0;               % Which method is used to calculate the thermocline depth: (0) some input 
                                    % function of MLD depths, MLD(), (1) the change in temperature gradient by 
                                    % tol, (2) the drop in near surface (10m depth) temperature by tol, or (3)
                                    % cubic spline the temperatures at all missing depths and then use method 2
params.canary = thermocline_calc;   % Whether or not Canary Islands data is being used
%params.current = 0;                 % Inclusion of current in model
%params.diffus = 0;                  % Inclusion of diffusion in model
%params.BC = 0;                      % Boundary conditions of spatial model: 0 periodic, 1 open

lat_centres = data_file_temps.latitude;     % Y-centres of lattice rows
lat_h = lat_centres(2)-lat_centres(1);      % Y-width of lattice rows
long_centres = data_file_temps.longitude;   % X-centres of lattice columns
long_h = long_centres(2)-long_centres(1);   % X-width of lattice columns
earth_SA = 510072000*1000^2;                % Surface area of Earth (unused)

if sum([params.temp_vert_mix, params.change_mix_depth, params.current, params.diffus]) == 0
    params.num_lats = 1;                            % Random site to model with no data dependence
    params.num_longs = 1;                       
    sites_to_visualise = [1,1];                     % List of sites [i1, j1; i2, j2; ...] 
    lat_centres = lat_centres(4);
    long_centres = long_centres(9);
%elseif thermocline_calc == 0
    %params.num_lats = 1;                            % Random site to model with no data dependence
    %params.num_longs = 1;                       
    %sites_to_visualise = [1,1];                     % List of sites [i1, j1; i2, j2; ...] 
    %lat_centres = lat_centres(4);
    %long_centres = long_centres(9);
elseif sum([params.current, params.diffus]) == 0
    params.num_lats = min(10,length(lat_centres));  % Random sites to model with no spatial flow 
    params.num_longs = min(10,length(long_centres));   
    sites_to_visualise = [4,9];                     % List of sites [i1, j1; i2, j2; ...] 
else
    params.num_lats = length(lat_centres);          % Number of latitudes in the data 
    params.num_longs = length(long_centres);        % Number of longitudes in the data 
    sites_to_visualise = [4,9];                     % List of sites [i1, j1; i2, j2; ...] 
end
temps = data_file_temps.water_temp(...
    1:end,:,1:params.num_lats,1:params.num_longs);      % Imported temperatures
v_north = data_file_flow.water_v(...
    1:end,:,1:params.num_lats,1:params.num_longs);      % Imported northerly currents 
v_east = data_file_flow.water_u(...
    1:end,:,1:params.num_lats,1:params.num_longs);      % Imported easterly currents 



%temps = zeros(1461,23,1,1); %%%%%%%%%%%%%%%%%%%%%% DELETE
%temps = temps(1:150,:,:,:); %%%%%%%%%%%%%%%%%%%%%% DELETE
dims = size(temps);                                     % Number of days, number of depths, number of rows, number of columns
num_days = dims(1);                                     % Number of distinct days of data
%num_days = 100; %%%%%%%%%%%%%%%%%%%%%%%% Delete
params.num_days = num_days;
num_depths = dims(2);                                   % Number of distinct depths of data
%v_north = zeros(1461,23,1,1); %%%%%%%%%%%%%%%%%%%%%% DELETE
%v_east = zeros(1461,23,1,1); %%%%%%%%%%%%%%%%%%%%%% DELETE
%v_north = ones(num_days,num_depths,params.num_lats,params.num_longs); % m/s %%%%%%%%%%%%%%%%%%%%%% DELETE
%v_east = ones(num_days,num_depths,params.num_lats,params.num_longs); % m/s %%%%%%%%%%%%%%%%%%%%%% DELETE
v_north = v_north .* 60*60*24;                          % Convert from m/s to m/day
v_east = v_east .* 60*60*24;                            % Convert from m/s to m/day
params.num_sites = params.num_lats*params.num_longs;    % Number of lattice sites
params.days = 0:num_days*params.repeat_num-1;                             % Array of day indices commencing at day 0

% PARAMETERS 
params.assumed_MLD = 12.5; % m                  % Assumed MLD in non-chanding model
params.assumed_Zeu = 100; % m                    % Assumed euphotic zone depth
params.Pmax = 2; % 1/(day) 0.2*12.5/2.58      % Maximum phytoplankton growth rate under optimal light conditions (0.5-2)
params.a = 0.184; % 1/(m day)                     % Related to the mazimum growth rate of P
%params.a = 2.58*params.Pmax/params.assumed_MLD; % 1/(m day)
params.b = 0.23; % 1/m                           % Light attenuation by water
%params.b = log(100)/params.assumed_Zeu; % 1/m
params.c = 0.4; % m^2/(g C)                     % P self-shading coefficient
%params.d = 1; %1; % m^3/(g C day)               % Higher predation on Z (for 1996 paper)
params.e = 0.05; %0.05; % g C/(m^3)             % Half-saturation constant for N uptake
params.k = 0.05; % 1/day                        % Cross-thermocline exchange rate
params.q = 0.075; %0.11; % 1/day                % Higher predation on Z (for 1999 paper)
params.r = 0.15; % 1/day                        % P respiration rate
params.s = 0.04; % 1/day                        % P sinking loss rate
params.N_0 = 0.8; % g C/(m^3)                   % N concentration below mixed layer
params.alpha = 0.25;                            % Z growth efficiency
params.beta = 0.33;                             % Z excretion fraction
params.gamma = 0.5;                             % Regeneration of Z predation excretion
params.lambda = 0.6; % 1/day                    % Maximum Z grazing rate
params.mu = 0.035; % g C/(m^3)                  % Z grazing half-saturation coefficient 
p2= 0; p1 = 1-p2; %p1 = 0.5; p2=0.25            % Palatability of detritus to zooplankton (in [0,1], if p2 = 0
%params.omega = 0.5; %p2/p1;                     % then no consumption of detritus by zooplankton) 
params.phi = 0.05; % 1/day                       % D remineralisation rate
params.psi = 0.085; % 1/day                      % D sinking loss rate
params.dif_N = 1 .* 60*60*24; % m^2/day         % ivity constant for nutrients
params.dif_P = 0.05.* 60*60*24; % m^2/day       % Diffusivity constant for phytoplankton
params.dif_Z = 0.5 .* 60*60*24; % m^2/day       % Diffusivity constant for zooplankton
params.dif_D = 0.05 .* 60*60*24; % m^2/day      % Diffusivity constant for detritus

% VISUALISATION
visual_switch = 1;                              % 1=on, 0=off
params.vid_speed = 24;                          % frames per sec
vars_to_visualise = {'N','P','Z'};                  % State variables to visualise: N, P, Z, D (min 1, max 4)
[num_sites,~] = size(sites_to_visualise);
%%
% COLOURS: uisetcolor
col.black = [0 0 0];
col.white = [1     1     1];
col.off_white = [0.9020    0.9020    0.9020];
col.ocean_blue = [0    0.3294    0.6196];
col.nutri_yellow_bright = [1 1 0];
col.nutri_yellow = [0.9294    0.6941    0.1255];
col.nutri_orange = [1.0000    0.4118    0.1608];
col.phyto_green = [0.1412    0.6392    0.4392];
col.phyto_green_mid = [0.0157    0.7804    0.4235]; 
col.phyto_green_bright = [0.0196    1.0000    0.5412]; %[0 1 0];
col.zoop_light = [0    0.8667    1.0000];
col.zoop_mid = [0.0275    0.4471    0.9294];
col.zoop_dark = [0.1490    0.1490    0.6392];
col.detri_purp = [0.6314    0.0549    0.5137];
col.detri_mid = [0.9098    0.1451    0.7569];
col.detri_light = [1.0000    0.4902    0.7961];
col.cold_blue = [0.0745    0.6235    1.0000];
col.hot_red = [1.0000    0.2392    0.2392];

% COLOUR MAPS
col.map.temp = parula; %[0 0 0; multigradient([col.cold_blue; col.hot_red], length = 256)];
col.map.nutri = [multigradient([col.white; col.nutri_yellow_bright; col.nutri_yellow; col.nutri_orange], length = 256)];
col.map.phyto = [multigradient([col.white; col.phyto_green_bright; col.phyto_green_mid; col.phyto_green], length = 256)];
col.map.zoop = [multigradient([col.white; col.zoop_light; col.zoop_mid; col.zoop_dark], length = 256)];
col.map.detri = [multigradient([col.white; col.detri_light; col.detri_mid; col.detri_purp], length = 256)];
col.map.MLD = bone; % [0 0 0; multigradient([col.white; col.black], length = 256)];
col.map.MLD = col.map.MLD(end:-1:1,:);
col.map.scale = copper;
%%
% TIME AND LOCATION
params.year1 = 2011;                                   % First calendar year of data %2013
month1 = 1;                                     % First calendar month of data: 1: Jan, 2: Feb, ... , 12: Dec
params.day1 = 1;                                % First day of the month of idata
params.hemisphere = "N";                        % Hemisphere of data: N: northern, S: southern

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

% IF VERTICAL STRATIFICATION IS AT PLAY AND DATA IS INFORMING IT, LOAD THE
% FILE WITH THE DECADE AVERAGE IN THE AREA
load("decade_TSI_values.mat");
for ii = 1:params.num_lats
    for jj =1:params.num_longs
        params.TSI_CS{ii,jj} = cubic_splines.TSI_CS{ii,jj};
        if params.num_lats == 1
            params.TSI_{ii,jj}=cubic_splines.TSI_CS{4,9};
        end
    end
end
params.TSI_mean = av_TSI_2006to2015;
params.TSI_std_dev = std_dev_TSI_2006to2015;
% Function to implement stratification: alternative @(T) heaviside(params.Tc-T)
%params.strat_func = @(x) (tanh((params.TSI_std_dev-(params.TSI_mean - x))/2)+1)/2; 
params.strat_func = @(x) ...
    1/(1+exp(-log(40000)*(x-params.TSI_mean)/params.TSI_std_dev ...
    + log(1/sqrt(40000))));
TSI_2011 = store_2011;


% FOR EACH GRID SITE, CALCULATE THE THERMOCLINE, THE MLD, THE LATTICE SITE
% AREA, AND FIND THE CUBIC SPLINES FOR TEMPERATURE, THE MLD AND CURRENTS 
% AS FUNCTIONS OF TIME.
thermocline_cubicSplines
%%
% PREALLOCATE SPACE FOR SOLUTIONS
if params.detritus_layer % NPZD model
    record.soluts = zeros(params.num_lats,params.num_longs,4,num_days*100*params.repeat_num); 
else % NPZ model
    record.soluts = zeros(params.num_lats,params.num_longs,3,num_days*100+params.repeat_num); 
end

% PREPARE INITIAL CONDITION
% N 0.3
N0array = 0.4 .* ones(params.num_lats,params.num_longs);
N0array(logical(params.landmass)) = 0;
N0vect = reshape(N0array,[1,params.num_sites]);
% P 0.1
P0array = 0.1 .* ones(params.num_lats,params.num_longs);
P0array(logical(params.landmass)) = 0;
P0vect = reshape(P0array,[1,params.num_sites]);
% Z 0.2
Z0array = 0.05 .* ones(params.num_lats,params.num_longs);
Z0array(logical(params.landmass)) = 0;
Z0vect = reshape(Z0array,[1,params.num_sites]);
if params.detritus_layer
    % D
    D0array = 0.08 .* ones(params.num_lats,params.num_longs);
    D0array(logical(params.landmass)) = 0;
    D0vect = reshape(D0array,[1,params.num_sites]);

    y0 = [N0vect, P0vect, Z0vect, D0vect];
else
    y0 = [N0vect, P0vect, Z0vect];
end

% SOLVE FOR ODE CONCENTRATIONS
%options = odeset('NonNegative',1);
maxstep = 1/...
    (max(abs(v_north_av),[],'all')/params.delta_y + ...
    max(abs(v_east_av),[],'all')/params.delta_x) - 0.01;
% Set maximise time-step size to avoid convergence errors
% Include progress bar and abort option
%options = odeset('MaxStep', maxstep,'OutputFcn', @odeprog,'Events', @odeabort);%,'Events', @myEvent);  
options = odeset('MaxStep', maxstep);%,'Events', @myEvent);      
[t,y] = ode45(@(t,x) odefunc(t,x,params),...
    [params.days(1) params.days(end)],y0,options); % or ode15s, with options

record.times=t; 
Narray = y(:, 1 : params.num_sites); 
Parray = y(:, params.num_sites+1 : 2*params.num_sites);
Zarray = y(:, 2*params.num_sites+1 : 3*params.num_sites);
if params.detritus_layer
    Darray = y(:, 3*params.num_sites+1 : end);
end

% CONVERT SOLUTIONS FROM VECTORS TO MATRICES
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
record.soluts = record.soluts(:,:,:,1:length(t));

% CREATE ARRAYS FOR THE MONTH AND SEASON ASSOCIATED WITH EACH DAY OF
% THE SIMULATION
month_season_arrays


% SAVE WORKSPACE
%save([folder_path '/Workspace.mat'])

%% VISUALISATION
%load('/Users/aliladearie/Documents/2023 sem 2/masters research/code/19-Sep-2023/r0.15_paper1996_202391917351/Workspace.mat')
%folder_path = '/Users/aliladearie/Documents/2023 sem 2/masters research/code/19-Sep-2023/r0.15_paper1996_202391917351';
if visual_switch
    % INTERPOLATE SOLUTIONS FOR EACH DAY AND SAVE MAPS IN A GIF
    %create_gif

    % PLOT PHASE PORTAIT
    %phase_portraits

    % PLOT TRAJECTORIES OF NPZ(D) CONCENTRATIONS FOR ALL LATTICE SITES
    %plot_all_NPZD_trajectories

    % PLOT TRAJECTORIES OF NPZ(D) CONCENTRATIONS FOR CHOSEN LATTICE SITES
    %plot_NPZD_trajectories
end