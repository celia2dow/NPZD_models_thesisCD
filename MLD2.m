function [depth, h] = MLD2(t,params)
% Approximately sinusoidal MLD
    if params.hemisphere == "S"
        t = params.day1+t+182.625;
    end
    depth = 68.75*sin(2*pi*(t+10.385)/365)+81.25;
    h = 68.75*2*pi*cos(2*pi*(t+10.385)/365)/365;
end