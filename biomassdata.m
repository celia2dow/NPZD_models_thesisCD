% Import data of zooplankton biomass and plot
data = readtable('krillbase_data.csv');
[rows, cols] = size(data);
inRegionOfInterest = zeros(rows,4);
inRegionOfInterest(:,1) = data.LATITUDE <-60;%<-62.5;   % Maximum latitude
inRegionOfInterest(:,2) = ones(rows,1);%data.LATITUDE >-63.5;   % Minimum latitude ones(rows,1);%
inRegionOfInterest(:,3) = data.LONGITUDE <-60;    % Maximum longitude
inRegionOfInterest(:,4) = ones(rows,1);%data.LONGITUDE >-63.5;  % Minimum longitude
isFoundInRegion_logical = min(inRegionOfInterest,[],2);
row_numbers = 1:rows;
isFoundInRegion_index = row_numbers(logical(isFoundInRegion_logical'));
dataInRegion = data(isFoundInRegion_index,:);
size(dataInRegion)

isFoundInYears_logical = ismember(dataInRegion.DATE.Year,[2013, 2014, 2015]);
isFoundinYears_index = row_numbers(logical(isFoundInYears_logical'));
dataInYears = dataInRegion(isFoundinYears_index,:);
size(dataInYears)
