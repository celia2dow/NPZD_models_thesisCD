% Attempt at recreating experiment C from:
%   Sequential variations of phytoplankton growth and mortality in an NPZ 
%   model: A remote-sensing-based assessment
%   Roy. S., Broomhead. D.S., Platt. T., Sathyendranath. S., % Ciavatta. S.
% Experiment C details: The reference form of the NPZ model is run over the 
% observation period keeping all parameters constant (as in Table 1) except 
% βT; this parameter is considered to be a function of temperature and is 
% determined by Eq. (8) (in this equation β20 is constant and the time 
% series of βT is generated using satellite-derived SSTs). This treatment 
% yields time series of nutrient, phytoplankton and zooplankton forced by a 
% growth parameter corrected for temperature effects.

close all 
% Import Sea-Surface Temperatures (SSTs)
load("5_years.mat","hawaii_soest_7e38_7a7b_afxhffc2")
temps = hawaii_soest_7e38_7a7b_afxhffc2.water_temp(:,:,1,1);
temps = mean(temps,2); % average temp over all depths in degrees Celcius

% Times
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

% PARAMETERS table 1
% Mortality rate of phytoplankton (per day) IF STATIC
params.mupT = 0.03; 
% Zooplankton searching rate ((mmol N m^(-3))^(-2) per day)
params.eta = 1.0;
% Zooplankton conversion efficiencty (-)
params.gamma = 0.75;
% Maximum rate of zooplankton grazing (per day)
params.alpha = 2.0;
% Nutrient input rate (per day)
params.s = 0.648;
% Half saturation of nutrient uptake (mmol N m^(-3))
params.kN = 0.5;
% Recylcing efficiency (-)
params.muN = 0.2;
% Zooplankton mortality ((mmol N m^(-3))^(-1) per day)
params.muZ = 0.2;
% Nutrient level below the mixed later (mmol N m^(-3))
params.N0 = 8.0;

% PARAMETERS equation 8
% Specific growth rate of phytoplankton (per day)
beta20 = 0.66; % at T=20
params.betaT = beta20 .* 1.066 .^ (temps - 20); % at any T

% PARAMETERS equation 9
% Mortality rate of phytoplankton (per day) IF DYNAMIC
mu20 = 0.03; % specific mortality at T=20
Q10 = 1.88; % change in specific loss rate due to a temperature change of 10 degrees Celcius
% params.mupT = mu20 .* Q10 .^ ((T-20)./10); % If dynamic

% Solve and plot ODE
% Initial N: from Fig 4a (mmol N m^(-3) per day)
% Initial P: from Fig 2a (mg m^(-3))
% Initial Z: from Fig 4d (mg m^(-3))
MolWeightN = 14.006747;
% y0 = [6.5,0.5/MolWeightN,1/MolWeightN];
 y0 = [6.5,0.5,1];
%y0 = [0.4/10,0.1/10,0.05/10];
[t,y] = ode15s(@odefunc,[1 num_days],y0, [], params);
months = t./(365.25/12);

subplot(3,1,1)
plot(months,y(:,1));
xlabel("Time (months)")
ylabel("Nutrient abundance (?)")

subplot(3,1,2)
plot(months,y(:,2));
xlabel("Time (months)")
ylabel("Phytoplankton abundance (?)")

subplot(3,1,3)
plot(months,y(:,3));
xlabel("Time (months)")
ylabel("Zooplankton abundance (?)")


function dxdt = odefunc(t,x,params)
    % Rename params
    mupT = params.mupT;
    eta = params.eta;
    gamma = params.gamma;
    alpha = params.alpha;
    s = params.s;
    kN = params.kN;
    muN = params.muN;
    muZ = params.muZ;
    N0 = params.N0;
    
    % Interpolate betaT and mupT values
    betaT = interp1(params.days, params.betaT, t);
    %mupT = interp1(params.days, params.mupT, t);

    % Rename state variables
    N = x(1);
    P = x(2);
    Z = x(3);

    % ODE SYSTEM equations 1-3
    dNdt = s*(N0-N) - betaT*N*P/(kN+N) + muN*((1-gamma)*alpha*eta*P^2*Z/...
        (alpha+eta*P^2) + mupT*P + muZ*Z^2);
    dPdt = betaT*N*P/(kN+N) - alpha*eta*P^2*Z/(alpha+eta*P^2) - mupT*P;
    dZdt = gamma*alpha*eta*P^2*Z/(alpha+eta*P^2) - muZ*Z^2;

    dxdt = zeros(3,1);
    dxdt(1) = dNdt;
    dxdt(2) = dPdt;
    dxdt(3) = dZdt;
end