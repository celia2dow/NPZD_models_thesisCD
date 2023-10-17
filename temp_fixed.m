function temp = temp_fixed(t,params)
% Approximately sinusoidal temperature
    if params.hemisphere == "S"
        t = params.day1+t+182.625;
    end
    temp = -4*cos(2*pi*(t-50)/params.num_days)+21.5;
end