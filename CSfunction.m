function y = CSfunction(cubic_spline, x)
% Takes in a cubic spline structure and a value x at which to evaluate
% the cubic spline and returns the value of the cubic spline y
    %indx = min(floor(x)+1,cubic_spline.pieces);
    difs = cubic_spline.breaks - x; 
    indx = min(sum(difs<=0),cubic_spline.pieces);
    x0 = cubic_spline.breaks(indx);
    c = cubic_spline.coefs;
    y = c(indx,1)*(x-x0)^3 + c(indx,2)*(x-x0)^2 + c(indx,3)*(x-x0) + c(indx,4);
end