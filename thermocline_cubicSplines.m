% Code calculating the thermocline and cubic splines for the MLD, 
% temperature, and eastern and northern currents. Called by 
% NPZD_model_forcing.m

% FOR EACH GRID SITE, CALCULATE THE THERMOCLINE, THE MLD, THE LATTICE SITE
% AREA, AND FIND THE CUBIC SPLINES FOR TEMPERATURE, THE MLD AND CURRENTS 
% AS FUNCTIONS OF TIME.
depth_indices = 1:num_depths;                                       % Array of indices for each depth 
MLD_indices = zeros(num_days,params.num_lats,params.num_longs);     % Initialise array for indices of MLD
depths_LEVs = zeros(num_days,params.num_lats,params.num_longs);     % Initialise array for MLDd in meters
temps_av = zeros(num_days,params.num_lats,params.num_longs);        % Initialise array for average temperature of mixed layer
area_grid_site = zeros(params.num_lats,params.num_longs);           % Initialise array for lattice site area
params.measures = zeros(2,params.num_lats,params.num_longs);        % Initialise array for eastern/northern distances of lattice sites in meters
v_north_av = zeros(num_days,params.num_lats,params.num_longs);      % Initialise array for average northerly current of mixed layer
v_east_av = zeros(num_days,params.num_lats,params.num_longs);       % Initialise array for average easterly current of mixed layer
params.landmass = zeros(params.num_lats,params.num_longs);          % Initialise array for sites that are landmasses
params.volumes = zeros(params.num_lats,params.num_longs,num_days);  % Initialise array for volume of mixed layer in cubed meters (unused)
lats = zeros(params.num_lats,2);                                    % Array of latitudes relevant to lattice site area being calculated
longs = zeros(params.num_longs,2);                                  % Array of longitudes relevant to lattice site area being calculated
wgs84 = wgs84Ellipsoid("m");                                        % World Geodetic System of 1984 (WGS84) 
params.temp_CS = cell(params.num_lats,params.num_longs);           % Initialise cell array of cubic splines for temperature
params.MLD_CS = cell(params.num_lats,params.num_longs);            % Initialise cell array of cubic splines for MLD
params.v_north_CS = cell(params.num_lats,params.num_longs);        % Initialise cell array of cubic splines for northerly current
params.v_east_CS= cell(params.num_lats,params.num_longs);          % Initialise cell array of cubic splines for easterly current

for ii = 1:params.num_lats
    for jj = 1:params.num_longs 
        % Calculate the area of the lattice site according to the World
        % Geodetic System 1984
        lats(ii,:) = [lat_centres(ii)-lat_h/2; lat_centres(ii)+lat_h/2];
        longs(jj,:) = [long_centres(jj)-long_h/2; long_centres(jj)+long_h/2];
        area_grid_site(ii,jj) = ...
            areaquad(lats(ii,1),longs(jj,1),lats(ii,2),longs(jj,2),wgs84);
        % Eastern distances
        params.measures(1,ii,jj) = distance(lat_centres(ii),longs(jj,1),...
            lat_centres(ii),longs(jj,2), wgs84);
        % Northern distances
        params.measures(2,ii,jj) = distance(lats(ii,1), long_centres(jj), ...
            lats(ii,2),long_centres(jj), wgs84);
        num_NaNs = 0;       % Count the number of NaN entries for site (i,j)
        for day = 1:num_days 
            % Find MLD, temperature, and northern + easterly currents for lattice site (i,j) for each day
            % Note any entries that are missing (NaN) to be later ignored
            temps_ijday = squeeze(temps(day,:,ii,jj));
            if isnan(temps_ijday(1))
                num_NaNs=num_NaNs+1;
            end
            if length(data_file_flow.LEV) ~=1
                v_north_ijday = squeeze(v_north(day,:,ii,jj));
                v_east_ijday = squeeze(v_east(day,:,ii,jj));
                GOODindicesCN = find(~isnan(v_north_ijday));    % The indices of depths with data
                num_GOODdepthsCN = length(GOODindicesCN);
                GOODindicesCE = find(~isnan(v_east_ijday));     % The indices of depths with data
                num_GOODdepthsCE = length(GOODindicesCE);
            else
                v_north_ijday = squeeze(v_north(day,1,ii,jj));
                v_east_ijday = squeeze(v_east(day,1,ii,jj));
            end
            GOODindicesT = find(~isnan(temps_ijday));           % The indices of depths with data
            num_GOODdepthsT = length(GOODindicesT);
            
            % Find MLD
            if thermocline_calc == 0                            % Some input function of MLD depths, MLD()
                [depth_meters_low, ~] = MLD(day-1, params);
                difs = data_file_temps.LEV-depth_meters_low;
                depth_index = max(find(difs<=0));               % Define the index of the MLD
                depths_LEVs(day,ii,jj) = depth_meters_low;            % Define the MLD in meters
                TSI_2011(day,ii,jj) = TSI_fixed(day,params);
            elseif thermocline_calc == 1                        % The change in temperature gradient by tol
                temp_smooth = smooth(temps(day,GOODindicesT,ii,jj));
                temp_grad = abs(temp_smooth(2:num_GOODdepthsT) - temp_smooth(1:num_GOODdepthsT-1));
                depth_index = min(find(temp_grad>tol));         % Define the index of the MLD
            elseif thermocline_calc == 2                        % The drop in near surface (10m depth or shallower) temperature by tol
                NSD = 6;                                        % Near-surface-depth index for 10m
                while ~ismember(NSD, GOODindicesT) 
                    NSD = NSD-1;
                    if NSD == 0
                        % warning("Warning: No data exists for depth 0-10m for i=%d, j=%d, day=%d \n",i,j,day)
                        depth_index = [];                       % Define the index of the MLD as empty if no data
                        break
                    end
                end
                if NSD>0
                    %if NSD<6
                    %    fprintf("Note: i=%d, j=%d, day=%d, near surface depth taken at %d m \n",...
                    %        i,j,day,data_file_temps.LEV(NSD));
                    %end
                    near_surf_temp = temps(day,NSD,ii,jj);
                    thresh_temp_low = near_surf_temp - tol;         % Threshold temperature
                    depth_index = find(temps(day,:,ii,jj)...
                        <thresh_temp_low, 1 );                      % Define the index of the MLD
                end
            elseif thermocline_calc == 3                        % Cubic spline the temperatures with respect to depth
                NSD = 6;                                        % Near-surface-depth index for 10m
                while ~ismember(NSD, GOODindicesT) 
                    NSD = NSD-1;
                    if NSD == 0
                        % warning("Warning: No data exists for depth 0-10m for i=%d, j=%d, day=%d \n",i,j,day)
                        depth_index = [];                       % Define the index of the MLD as empty if no data
                        break
                    end
                end
                if NSD>0
                    %if NSD<6
                    %    fprintf("Note: i=%d, j=%d, day=%d, near surface depth taken at %d m \n",...
                    %        i,j,day,data_file_temps.LEV(NSD));
                    %end
                    %spline_depths = min(10,data_file_temps.LEV(GOODindicesT(end))) : ...
                    %    data_file_temps.LEV(GOODindicesT(end)); % Possible depths below 10m %data_file_temps.LEV(end);
                    %cs_temp_v_depth = spline(data_file_temps.LEV(GOODindicesT), ...  
                    %    temps(day,GOODindicesT,i,j), spline_depths);
                    %near_surf_temp = temps(day,NSD,i,j);

                    spline_depths = 1 : ...
                        data_file_temps.LEV(GOODindicesT(end)); % Possible depths below 10m %data_file_temps.LEV(end);
                    cs_temp_v_depth = spline(data_file_temps.LEV(GOODindicesT), ...  
                        temps(day,GOODindicesT,ii,jj), spline_depths);
                    near_surf_temp = mean(cs_temp_v_depth(1:data_file_temps.LEV(NSD)),'all');
                    thresh_temp_low = near_surf_temp - tol;         % Threshold temperature below near surface temperature
                    thresh_temp_high = near_surf_temp + tol;        % Threshold temperature above near surface temperature
                    depth_meters_low = spline_depths(min(...
                        find(cs_temp_v_depth<thresh_temp_low)));    % Define the MLD in meters
                    depth_meters_high = spline_depths(min(...
                        find(cs_temp_v_depth>thresh_temp_high)));  
                    depth_meters = [depth_meters_low, depth_meters_high];
                    if isempty(depth_meters) && ~isempty(spline_depths)
                        depth_meters = spline_depths(end);      % Define the MLD in meters as deepest if threshold change isn't met
                        depth_index = GOODindicesT(end);
                    elseif isempty(depth_meters)                
                        depth_meters = [];                  % Define the MLD as empty if no data
                        depth_index = [];
                    else
                        depth_meters_max = max(depth_meters);
                        if sum(depth_meters_max == depth_meters_high) && day>170
                            %display(depth_meters_max)
                        end
                        difs = data_file_temps.LEV-depth_meters_max;
                        depth_index = max(find(difs<=0));       % Define the index of the MLD
                        depths_LEVs(day,ii,jj) = depth_meters_max;% Define the MLD in meters
                    end
                end
            end
            
            % Record the MLD, average temperature and current within this layer
            if depth_index
                % The index of the MLD
                MLD_indices(day,ii,jj)=depth_index;
                % The actual depth in meters
                if depths_LEVs(day,ii,jj) == 0
                    depths_LEVs(day,ii,jj) = data_file_temps.LEV(depth_index);
                end
                % Find average temperature within the mixed layer (degrees Celsius)
                GOODindices_inMLD_T = GOODindicesT(GOODindicesT<=depth_index);
                if length(data_file_flow.LEV) ~=1
                    GOODindices_inMLD_CN = GOODindicesCN(GOODindicesCN<=depth_index);
                    GOODindices_inMLD_CE = GOODindicesCE(GOODindicesCE<=depth_index);
                end
                if thermocline_calc == 0
                    temps_av(day,ii,jj) = temp_fixed(day,params);
                    temps(day,:,ii,jj) = temp_fixed(day,params);
                else
                    temps_av(day,ii,jj) = squeeze(mean(temps(day,GOODindices_inMLD_T,ii,jj),2));
                end
                % Find average current velocities within the mixed layer (degrees Celsius)
                if length(data_file_flow.LEV) ~=1
                    v_north_av(day,ii,jj) = squeeze(mean(v_north(day,GOODindices_inMLD_CN,ii,jj),2));
                    v_east_av(day,ii,jj) = squeeze(mean(v_east(day,GOODindices_inMLD_CE,ii,jj),2));
                else
                    v_north_av(day,ii,jj) = squeeze(mean(v_north(day,1,ii,jj),2));
                    v_east_av(day,ii,jj) = squeeze(mean(v_east(day,1,ii,jj),2));
                end
            else
                % If there are no data points for lattice site (i,j) on this
                % day, ignore this entire case
                % fprintf("Note: no data for i=%d, j=%d, day=%d \n",i,j,day);
                MLD_indices(day,ii,jj)=NaN;
                depths_LEVs(day,ii,jj)=NaN;
                temps_av(day,ii,jj)=NaN;
                v_north_av(day,ii,jj)=NaN;
                v_east_av(day,ii,jj)=NaN;
            end

            % 
        end

        % Smooth temperature, MLD and currents if specified
        days_indices = 1:num_days;
        GOODindices_days_T = days_indices(~isnan(temps_av(:,ii,jj)));
        GOODindices_days_CN = days_indices(~isnan(v_north_av(:,ii,jj)));
        GOODindices_days_CE = days_indices(~isnan(v_east_av(:,ii,jj)));
        if smoothing
            temps_av(GOODindices_days_T,ii,jj) = ...
                smooth(temps_av(GOODindices_days_T,ii,jj),smoothing);
            depths_LEVs(GOODindices_days_T,ii,jj) = ...
                smooth(depths_LEVs(GOODindices_days_T,ii,jj),smoothing);
            v_north_av(GOODindices_days_CN,ii,jj) = ...
                smooth(v_north_av(GOODindices_days_CN,ii,jj),smoothing);
            v_east_av(GOODindices_days_CE,ii,jj) = ...
                smooth(v_east_av(GOODindices_days_CE,ii,jj),smoothing);
        end

        % Record volume using (smoothed) MLD, assuming uniform lattice
        params.delta_x = mean(params.measures(1,:,:),"all");            % Average delta x
        params.delta_y = mean(params.measures(2,:,:),"all");            % Average delta y
        params.site_area = params.delta_x * params.delta_y;             % Average site area
        params.volumes(ii,jj,GOODindices_days_T)= ...
            depths_LEVs(GOODindices_days_T,ii,jj) .* params.site_area;    % Site volumes based on average site area
        %params.volumes(i,j,GOODindices_days_T)= ...
        %    depths_LEVs(GOODindices_days_T,i,j) .* area_grid_site(i,j);

        % Find a cubic spline to approximate the temperature, MLD, and
        % currents with respect to time
        SPLINEindices = 1:xth:num_days;
        if num_NaNs==num_days   % If there is no data for all days
            params.landmass(ii,jj)=1;    % Aassume there is a landmass here 
        else
            if params.num_lats == 1
                temps_av_rep_ij = temps_av';
                depths_LEVs_rep_ij = depths_LEVs';
                v_north_av_rep_ij = v_north_av';
                v_east_av_rep_ij = v_east_av';
            else
                temps_av_rep_ij = temps_av(:,ii,jj);
                depths_LEVs_rep_ij = depths_LEVs(:,ii,jj);
                v_north_av_rep_ij = v_north_av(:,ii,jj);
                v_east_av_rep_ij = v_east_av(:,ii,jj);
            end
            
            NANindices_days_T = days_indices(isnan(temps_av_rep_ij));               % Days with no temperature data
            SPLINEindices_T = setdiff(SPLINEindices,NANindices_days_T);             % Days with temperature data
            NANindices_days_CN = days_indices(isnan(v_north_av_rep_ij));            % Dito
            SPLINEindices_CN = setdiff(SPLINEindices,NANindices_days_CN);
            NANindices_days_CE = days_indices(isnan(v_east_av_rep_ij));
            SPLINEindices_CE = setdiff(SPLINEindices,NANindices_days_CE);
            cs_temp = spline(SPLINEindices_T-1, temps_av_rep_ij(SPLINEindices_T));         % Cubic spline of temperature w.r.t time
            cs_MLD = spline(SPLINEindices_T-1, depths_LEVs_rep_ij(SPLINEindices_T));       % Cubic spline of MLD w.r.t time
            cs_v_north = spline(SPLINEindices_CN-1, v_north_av_rep_ij(SPLINEindices_CN));  % Cubic spline of north current w.r.t time
            cs_v_east = spline(SPLINEindices_CE-1, v_east_av_rep_ij(SPLINEindices_CE));    % Cubic spline of east current w.r.t time
            if thermocline_calc == 0
                params.TSI_CS{ii,jj} = spline(GOODindices_days_T-1, ...
                    TSI_2011(GOODindices_days_T,ii,jj)); % Cubic spline of TSI w.r.t time
            end

            params.temp_CS{ii,jj} = cs_temp;
            params.MLD_CS{ii,jj} = cs_MLD;
            params.v_north_CS{ii,jj} = cs_v_north;
            params.v_east_CS{ii,jj} = cs_v_east;
        end

%        else
%             temps_av_rep_ij = repmat(temps_av(:,i,j),[repeat_num,1]);
%             depths_LEVs_rep_ij = repmat(depths_LEVs(:,i,j),[repeat_num,1]);
%             v_north_av_rep_ij = repmat(v_north_av(:,i,j),[repeat_num,1]);
%             v_east_av_rep_ij = repmat(v_east_av(:,i,j),[repeat_num,1]);

%             NANindices_days_T = days_indices(isnan(temps_av_rep_ij));
%             SPLINEindices_T = setdiff(SPLINEindices,NANindices_days_T);
%             NANindices_days_CN = days_indices(isnan(v_north_av_rep_ij));
%             SPLINEindices_CN = setdiff(SPLINEindices,NANindices_days_CN);
%             NANindices_days_CE = days_indices(isnan(v_east_av_rep_ij));
%             SPLINEindices_CE = setdiff(SPLINEindices,NANindices_days_CE);
%             cs_temp = spline(SPLINEindices_T-1, squeeze(temps_av(SPLINEindices_T,i,j))');
%             cs_MLD = spline(SPLINEindices_T-1, squeeze(depths_LEVs(SPLINEindices_T,i,j))');
%             cs_v_north = spline(SPLINEindices_CN-1, squeeze(v_north_av(SPLINEindices_CN,i,j))');
%             cs_v_east = spline(SPLINEindices_CE-1, squeeze(v_east_av(SPLINEindices_CE,i,j))');
%        end

%        if num_NaNs~=num_days   % If there is no landmass, save the cubic splines for site (i,j)
%             params.temp_CS{i,j} = cs_temp;
%             params.MLD_CS{i,j} = cs_MLD;
%             params.v_north_CS{i,j} = cs_v_north;
%             params.v_east_CS{i,j} = cs_v_east;
%        end
    end
end