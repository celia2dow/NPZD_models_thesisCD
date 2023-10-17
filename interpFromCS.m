function [temp, MLD, gradMLD] = interpFromCS(params,i,j,day)
% Interpolate the temperature, MLD and gradient of the MLD function on
% a specific date at a specific lattice site (i,j)
    name_site = sprintf('i%i_j%i',i,j);
    CSfuncTemp = params.temp_CS.(name_site);
    CSfuncMLD = params.MLD_CS.(name_site);
    temp = ppval(CSfuncTemp,day);
    MLD = ppval(CSfuncMLD,day);
    MLD_neighbourhood = ppval(CSfuncMLD,[day-0.01, day, day+0.01]);
    gradMLDs = gradient(MLD_neighbourhood,0.01);
    gradMLD = gradMLDs(2);
end