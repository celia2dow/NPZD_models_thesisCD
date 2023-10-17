% Adding temperature dependence to model from:
%   Zooplankton Mortality and the Dynamical Behaviour of Plankton 
%   Population Models
%   Edwards. A.M., & Brindley. J. 1999
% or
%   Oscillatory behaviour in a three-component plankton population model
%   Edwards. A.M., & Brindley. J. 1996
% using
%   Sequential variations of phytoplankton growth and mortality in an NPZ 
%   model: A remote-sensing-based assessment
%   Roy. S., Broomhead. D.S., Platt. R., Sathyendranath. S., & Ciavatta. S.
% The model consists of three coupled ordinary differential equations, 
% describing changes in the concentrations of nutrient (N), phytoplankton 
% (P) and zooplankton (Z) in a physically homogeneous oceanic mixed layer. 

close all 
% IMPORT Sea-Surface Temperatures (SSTs)
load("5_years.mat","hawaii_soest_7e38_7a7b_afxhffc2")

% SCALES
lats = length(hawaii_soest_7e38_7a7b_afxhffc2.latitude); %lats = 1; %
longs = length(hawaii_soest_7e38_7a7b_afxhffc2.longitude); %longs = 1; %
temps = hawaii_soest_7e38_7a7b_afxhffc2.water_temp(:,:,1:lats,1:longs);
temps = squeeze(mean(temps,2)); % average temp over all depths in degrees Celcius
num_days = length(temps);
params.days = 1:num_days;

% For any temperatures that are missing (NaN), replace them with the
% average temperature of adjacent temperatures
indices = find(isnan(temps));
for index = indices
    if index == 1
        temps(index) = temps(index+1);
    elseif index == num_days
        temps(index) = temps(index-1);
    else
        temps(index) = (temps(index-1)+temps(index+1))./2;
    end
end

% PAPER MODELLED
params.paper = 1996;

% PARAMETERS table 1
% Related to the mazimum growth rate of P
params.a = 0.2; % 1/(m day)
% Light attenuation by water
params.b = 0.2; % 1/m
% P self-shading coefficient
params.c = 0.4; % m^2/(g C)
% Higher predation on Z (for 1996 paper)
params.d = 1; %1; % m^3/(g C day)
% Half-saturation constant for N uptake
params.e = 0.05; %0.03; % g C/(m^3)
% Cross-thermocline exchange rate
params.k = 0.05; % 1/day
% Higher predation on Z (for 1999 paper)
params.q = 0.075; %0.11; % 1/day
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

% NEW/MODIFIED parameters
% P respiration rate
params.r = 0.07; % 1/day - 0.07, 0.15
% Recycling efficiency
params.mu_N = 0.2; % 0.2
% Specific mortality at 20 deg Cels *in this area*
mu_20 = 0.03; % 0.03 1/day
% Change in specific loss due to a temperature change of 10 deg Cels
Q_10 = 1.88;
% Specific growth rate *in this area*
beta_20 = 0.66; % 1/day

% MOVIE
visual_switch = 1; % 1=on, 0=off
params.vid_speed = 20; % frames per sec

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
record.beta_T = zeros(num_days,lats,longs);
record.mu_pT = zeros(num_days,lats,longs);
record.soluts = zeros(lats,longs,3,num_days*100);
record.times = zeros(lats,longs,num_days*100);
y0 = [0.3,0.025,0.05]; % N0 P0 Z0, [0.4 0.1 0.05] 

for i = 1:lats
    for j = 1:longs
        % PARAMETERS equation 8
        % Maximum specific growth rate of phytoplankton at any T (deg Cels)
        params.beta_T = beta_20 .* 1.066 .^ (temps(:,i,j) - 20); % 1/day
        record.beta_T(:,i,j)=params.beta_T;

        % PARAMETERS equation 9
        % Mortality rate of phytoplankton at any T (deg Cels)
        params.mu_pT = mu_20 .* Q_10 .^ ((temps(:,i,j)-20)./10); % 1/day
        record.mu_pT(:,i,j)=params.mu_pT;

        % Solve and plot ODE
        [t,y] = ode15s(@odefunc,[1 num_days],y0, [], params);
        record.soluts(i,j,1,1:length(t))=y(:,1); % Nutrient
        record.soluts(i,j,2,1:length(t))=y(:,2); % Phytoplankton
        record.soluts(i,j,3,1:length(t))=y(:,3); % Zooplankton
        record.times(i,j,1:length(t))=t; % time

    end
end

months = t./(365.25/12);

if lats == 1
    % NPZ plot
    figure(1)
    plot(months,y(:,1));
    xlabel("Time (months)")
    ylabel("Concentration (g C/{m^3})")
    hold on
    plot(months,y(:,2));
    plot(months,y(:,3));
    hold off
    legend('Nutrient (N)', 'Phytoplankton (P)', 'Zooplankton (Z)')
    
    % NPZ as well as temp, beta_T, mu_pT
    figure(2)
    Days = months;
    Nutrient = y(:,1);
    Phytoplankton = y(:,2);
    Zooplankton = y(:,3);
    Temperature = interp1(params.days,temps,t);
    beta_T = interp1(params.days, params.beta_T, t);
    mu_pT = interp1(params.days, params.mu_pT, t);
    data_plot = table(Days, Nutrient, Phytoplankton, Zooplankton, ...
        Temperature, beta_T, mu_pT);
    vars = ["Temperature","beta_T","mu_pT","Nutrient","Phytoplankton","Zooplankton"];
    stackedplot(data_plot,vars);

    % Plot mu_pT against beta_T
    figure(3)
    subplot(2,1,1)
    xlabel("days")
    yyaxis left
    plot(params.days, params.mu_pT)
    ylabel("per day")
    hold on
    yyaxis right
    plot(params.days, params.beta_T)
    ylabel("per day")
    hold off
    legend("\mu_{pT}","\beta_T")

    subplot(2,1,2)
    xlabel("months")
    yyaxis left
    plot(months, mu_pT)
    ylabel("per day")
    hold on
    yyaxis right
    plot(months, beta_T)
    plot(months, 1./(1+params.c.*y(:,2)./params.b))
    plot(months, y(:,1)./(params.e+y(:,1)))
    ylabel("per day")
    hold off
    legend("\mu_{pT}","\beta_T", "1/(1+cP/b)","N/(e+N)")
end

% Go through each day in time and save the spatial distributuion of
% temperature, and phytoplankton
% Interpolate produced data to the days for which we have temperature data

if visual_switch
    % Prepare movie.
    ocean_grid = figure;
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
    lims.min_T = min(temps(temps>0),[],'all');
    lims.max_T = max(temps,[],'all');
    lims.min_P = min(phyto,[],'all');
    lims.max_P = max(phyto,[],'all');


    
    for day = params.days
        clf
        % Save each figure as a frame in the movie
        draw_frame(params,phyto,temps,lims,day)
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
    movie(ocean_grid, ocean_movie, 1, params.vid_speed);
end

% Edited ODE system
function dxdt = odefunc(t,x,params)
    % Rename params
    a=params.a;
    b=params.b;
    c=params.c;
    d=params.d;
    e=params.e;
    k=params.k;
    q=params.q;
    s=params.s;
    N_0=params.N_0;
    alpha=params.alpha;
    beta=params.beta;
    gamma=params.gamma;
    lambda=params.lambda;
    mu=params.mu;
    r=params.r;
    mu_N=params.mu_N;
    
     % Interpolate betaT and mupT values
    beta_T = interp1(params.days, params.beta_T, t);
    mu_pT = interp1(params.days, params.mu_pT, t);
    % Rename state variables
    N = x(1);
    P = x(2);
    Z = x(3);

    % ODE SYSTEM equations 1-3
    if params.paper == 1999
        h_Z = q*Z;
    elseif params.paper == 1996
        h_Z = d*Z^2;
    end
    dNdt = -N*P*beta_T/((e+N)*(1+c*P/b)) + r*P + beta*lambda*P^2*Z/(mu^2+P^2) ...
        + gamma*h_Z + k*(N_0-N) + mu_N*mu_pT*P;
    dPdt = N*P*beta_T/((e+N)*(1+c*P/b)) - r*P - lambda*P^2*Z/(mu^2+P^2) - (s+k)*P ...
        - mu_pT*P;
    dZdt = alpha*lambda*P^2*Z/(mu^2+P^2) - h_Z;
    
    dxdt = zeros(3,1);
    dxdt(1) = dNdt;
    dxdt(2) = dPdt;
    dxdt(3) = dZdt;
end