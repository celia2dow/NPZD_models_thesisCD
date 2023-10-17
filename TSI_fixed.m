function TSI = TSI_fixed(t,params)
% Approximately sinusoidal temperature
    if params.hemisphere == "S"
        t = params.day1+t+182.625;
    end
    TSI = 1.2*(-params.TSI_std_dev*cos(2*pi*(t+90)/params.num_days)...
        +params.TSI_mean);
end