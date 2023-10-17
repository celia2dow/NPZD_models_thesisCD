function [value, isterminal, direction] = myEvent(T, Y)
value      = (Y(1) >= 1000);
isterminal = 1;   % Stop the integration
direction  = 0;
end