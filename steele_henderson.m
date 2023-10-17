% Attempt at recreating model from:
%   The role of predation in plankton models
%   Steele. J.H., & Henderson. E.W.
% The model consists of three coupled ordinary differential equations, 
% describing changes in the concentrations of nutrient (N), phytoplankton 
% (P) and zooplankton (Z) in a physically homogeneous oceanic mixed layer. 

close all 

% Times
num_days = 200;
params.days = 1:num_days;

% Paper modelled
params.paper = 1996;

% PARAMETERS table 1
% Related to the mazimum growth rate of P
params.a = 0.2; % 1/(m day)
% Light attenuation by water
params.b = 0.2; % 1/m
% P self-shading coefficient
params.c = 0.4; % m^2/(g C)
% Higher predation on Z (for 1996 paper)
params.d = 1.5; %1; % m^3/(g C day)
% Half-saturation constant for N uptake
params.e = 0.03; % g C/(m^3)
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

% Solve and plot ODE
y0 = [0.4,0.1,0.05];
[t,y] = ode15s(@odefunc,[1 num_days],y0, [], params);

plot(t,y(:,1));
xlabel("Time (days)")
ylabel("Concentration (g C/{m^3})")
hold on
plot(t,y(:,2));
plot(t,y(:,3));
legend('Nutrient (N)', 'Phytoplankton (P)', 'Zooplankton (Z)')


function dxdt = odefunc(~,x,params)
    % Rename params
    a=params.a;
    b=params.b;
    c=params.c;
    d=params.d;
    e=params.e;
    k=params.k;
    q=params.q;
    r=params.r;
    s=params.s;
    N_0=params.N_0;
    alpha=params.alpha;
    beta=params.beta;
    gamma=params.gamma;
    lambda=params.lambda;
    mu=params.mu;
    
    % Rename state variables
    N = x(1);
    P = x(2);
    Z = x(3);

    % ODE SYSTEM equations 1-3
    dNdt = -N*a*P/((e+N)*(b+c*P)) + r*P + beta*lambda*P^2*Z/(mu^2+P^2) ...
        + gamma*q*Z + k*(N_0-N);
    dPdt = N*a*P/((e+N)*(b+c*P)) - r*P - lambda*P^2*Z/(mu^2+P^2) - (s+k)*P;
    if params.paper == 1999
        dZdt = alpha*lambda*P^2*Z/(mu^2+P^2) - q*Z;
    elseif params.paper == 1996
        dZdt = alpha*lambda*P^2*Z/(mu^2+P^2) - d*Z^2;
    end
    dxdt = zeros(3,1);
    dxdt(1) = dNdt;
    dxdt(2) = dPdt;
    dxdt(3) = dZdt;
end