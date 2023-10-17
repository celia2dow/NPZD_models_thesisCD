function [depth, h] = MLD(t,params)
% MLD from Edwards PhD paper
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