% CALCULATE THE AVERAGE FOR A DECADE
tol = 1;              % Tolerance in temperature change indicating the thermocline base
TSI_per_daySite = zeros(366*10*50*50,1);
indx = 1;
store_2011 = zeros(366,50,50);
therm_base_2011 = zeros(366,50,50);
for hh = 6:15
    single_digit = hh<10;
    data_file_name = ['canary_islands_at_collection_site_' num2str(2000+hh)];
    load(data_file_name,"hawaii_soest_7e38_7a7b_afxhffc2") 
    data_file_temps = hawaii_soest_7e38_7a7b_afxhffc2;  % Sea-Surface Temperatures (SSTs)
    
    lat_centres = data_file_temps.latitude;     % Y-centres of lattice rows
    long_centres = data_file_temps.longitude;   % X-centres of lattice columns
    num_lats = length(lat_centres);          % Number of latitudes in the data 
    num_longs = length(long_centres);        % Number of longitudes in the data 
    temps = data_file_temps.water_temp(1:end,:,1:num_lats,1:num_longs); 
    dims = size(temps);                                 % Number of days, number of depths, number of rows, number of columns
    num_days = dims(1);                                 % Number of distinct days of data
    num_depths = dims(2);                               % Number of distinct depths of data
    depth_indices = 1:num_depths; 

    textprogressbar(['Year: ' num2str(2000+hh) '     ']);
    
    % Store TSIs and thermocline base if year is 2011
    if hh ==11
        store_2011 = store_2011(1:num_days,1:num_lats,1:num_longs);
        therm_base_2011 = therm_base_2011(1:num_days,1:num_lats,1:num_longs);
    end
    for ii = 1:num_lats
        textprogressbar(ii*100/num_lats)
        for jj = 1:num_longs 
            num_NaNs = 0;       % Count the number of NaN entries for site (i,j)
            for day = 1:num_days 

                % Find MLD andtemperature for lattice site (i,j) for each day
                % Note any entries that are missing (NaN) to be later ignored
                temps_ijday = squeeze(temps(day,:,ii,jj));
                if isnan(temps_ijday(1))
                    num_NaNs=num_NaNs+1;
                end
                
                GOODindicesT = find(~isnan(temps_ijday));       % The indices of depths with data
                num_GOODdepthsT = length(GOODindicesT);
                
                % Find base of thermocline using cubic splines of the 
                % temperatures with respect to depth
                if GOODindicesT
                    lowest_temp = temps_ijday(GOODindicesT(end));
                else
                    lowest_temp = [];
                end
                if lowest_temp
                    spline_depths = 1 : ...
                        data_file_temps.LEV(GOODindicesT(end)); % Possible depths 
                    cs_temp_v_depth = spline(data_file_temps.LEV(GOODindicesT), ...  
                        temps_ijday(GOODindicesT), spline_depths);
                    abs_dif_to_lowest_temp = abs(cs_temp_v_depth-lowest_temp);
                    depth_thermocline_base = spline_depths(max(...
                        find(abs_dif_to_lowest_temp>tol))); 

                    if isempty(depth_thermocline_base) && ~isempty(spline_depths)
                        depth_meters = spline_depths(end);  % Define the MLD in meters as deepest if threshold change isn't met
                    elseif isempty(depth_thermocline_base)                
                        depth_meters = [];                  % Define the MLD as empty if no data
                    else
                        depth_meters = depth_thermocline_base;
                    end
                end
                
                % Find the temperature stratification index (TSI) of this
                % site and day
                if depth_meters
                    indx_thermocline_base = find(spline_depths == depth_meters);
                    depth_av_temp = mean(cs_temp_v_depth(1:indx_thermocline_base));
                    TSI_ij_day = mean((cs_temp_v_depth(1:indx_thermocline_base) ...
                        -depth_av_temp) .* spline_depths(1:indx_thermocline_base));
                    TSI_per_daySite(indx) = TSI_ij_day;
                    indx = indx + 1;
                else
                    depth_meters = NaN;
                    TSI_ij_day = NaN;
                end

                % Store STI if year is 2011
                if hh ==11
                    store_2011(day,ii,jj) = TSI_ij_day;
                    therm_base_2011(day,ii,jj)= depth_meters;
                end
            end

            % Create cubic spline for TSIs at lattice site over time if
            % year is 2011
            if hh ==11
                days_indices = 1:num_days;
                GOODindices_days_TSI = days_indices(~isnan(store_2011(:,ii,jj)));
                cs_TSI = spline(GOODindices_days_TSI-1, ...
                    store_2011(GOODindices_days_TSI,ii,jj)); % Cubic spline of TSI w.r.t time
                cubic_splines.TSI_CS{ii,jj} = cs_TSI;
            end
        end
    end
    textprogressbar('done');
end
%%

TSI_per_daySite = TSI_per_daySite(1:indx-1);
av_TSI_2006to2015 = mean(TSI_per_daySite);
std_dev_TSI_2006to2015 = std(TSI_per_daySite);
save("decade_TSI_values.mat","therm_base_2011","TSI_per_daySite","store_2011","av_TSI_2006to2015","std_dev_TSI_2006to2015","cubic_splines")