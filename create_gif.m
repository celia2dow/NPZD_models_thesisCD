% Code creating the GIF of temperature, current, and desired variables.
% Called by NPZD_model_forcing.m

% Interpolate NPZ(D) data and rearrange temperature data for days after the
% burn in
condit_afterBI = record.times >= cut_off_day;
days_afterBI=record.times(condit_afterBI);
days_afterBI_shift = days_afterBI-days_afterBI(1);

plt.nutri = zeros(params.num_lats,params.num_longs,num_days);       % Initialise array for N
plt.phyto = zeros(params.num_lats,params.num_longs,num_days);       % Initialise array for P
plt.zoop = zeros(params.num_lats,params.num_longs,num_days);        % Initialise array for Z
if params.detritus_layer
    plt.detri = zeros(params.num_lats,params.num_longs,num_days);   % Initialise array for D
end
temp_rearr = zeros(num_days,params.num_lats,params.num_longs);      % Initialise array for temperatures
v_north_rearr = zeros(num_days,params.num_lats,params.num_longs);   % Initialise array for north currents
v_east_rearr = zeros(num_days,params.num_lats,params.num_longs);    % Initialise array for east currents
MLD_rearr = zeros(num_days,params.num_lats,params.num_longs);
TSI_rearr = zeros(num_days,params.num_lats,params.num_longs);
min_temp = min(temps_av,[],'all');
for i = 1:params.num_lats
    for j = 1:params.num_longs
        time_lim = length(squeeze(record.times));
        if ~params.landmass(i,j)        % If there is no landmass at (i,j)
            plt.nutri(params.num_lats+1-i,j,:) = interp1(days_afterBI_shift,...
                squeeze(record.soluts(i,j,1,condit_afterBI)),params.days);
            plt.phyto(params.num_lats+1-i,j,:) = interp1(days_afterBI_shift,...
                squeeze(record.soluts(i,j,2,condit_afterBI)),params.days);
            plt.zoop(params.num_lats+1-i,j,:) = interp1(days_afterBI_shift,...
                squeeze(record.soluts(i,j,3,condit_afterBI)),params.days);
            if params.detritus_layer
                plt.detri(params.num_lats+1-i,j,:) = interp1(days_afterBI_shift,...
                    squeeze(record.soluts(i,j,4,condit_afterBI)),params.days);
            end
            if ~(params.num_lats ==1)
                temp_rearr(:,params.num_lats+1-i,j) = spline(GOODindices_days_T, ...
                    temps_av(GOODindices_days_T,i,j), params.days(1:num_days));
                v_north_rearr(:,params.num_lats+1-i,j) = spline(GOODindices_days_CN, ...
                    v_north_av(GOODindices_days_CN,i,j), params.days(1:num_days));
                v_east_rearr(:,params.num_lats+1-i,j) = spline(GOODindices_days_CE, ...
                    v_east_av(GOODindices_days_CE,i,j), params.days(1:num_days));
                MLD_rearr(:,params.num_lats+1-i,j)= spline(GOODindices_days_T, ...
                    depths_LEVs(GOODindices_days_T,i,j), params.days(1:num_days));
                TSI_rearr(:,params.num_lats+1-i,j)= spline(GOODindices_days_T, ...
                    store_2011(GOODindices_days_T,i,j), params.days(1:num_days));
            else
                temp_rearr(:,params.num_lats+1-i,j) = spline(GOODindices_days_T, ...
                    temps_av(GOODindices_days_T), params.days(1:num_days));
                v_north_rearr(:,params.num_lats+1-i,j) = spline(GOODindices_days_CN, ...
                    v_north_av(GOODindices_days_CN), params.days(1:num_days));
                v_east_rearr(:,params.num_lats+1-i,j) = spline(GOODindices_days_CE, ...
                    v_east_av(GOODindices_days_CE), params.days(1:num_days));
                MLD_rearr(:,params.num_lats+1-i,j)= spline(GOODindices_days_T, ...
                    depths_LEVs(GOODindices_days_T), params.days(1:num_days));
                TSI_rearr(:,params.num_lats+1-i,j)= spline(GOODindices_days_T, ...
                    store_2011(GOODindices_days_T), params.days(1:num_days));
            end
        else                            % If there is a landmass at (i,j)
            plt.nutri(params.num_lats+1-i,j,:) = -0.01;
            plt.phyto(params.num_lats+1-i,j,:) = -0.01;
            plt.zoop(params.num_lats+1-i,j,:) = -0.01;
            if params.detritus_layer
                plt.detri(params.num_lats+1-i,j,:) = -0.01;
            end
            temp_rearr(:,params.num_lats+1-i,j) = min_temp - 0.01;
            v_north_rearr(:,params.num_lats+1-i,j) = 0;
            v_east_rearr(:,params.num_lats+1-i,j) = 0;
            MLD_rearr(:,params.num_lats+1-i,j)=0;
            TSI_rearr(:,params.num_lats+1-i,j)=0;
        end
    end
end
scales_rearr = params.strat_func(TSI_rearr);


% Find limits
lims.min_T = min(temps_av(temps_av~=0),[],'all');
lims.max_T = max(temps_av,[],'all');
lims.min_N = -0.01;
lims.max_N = max(plt.nutri,[],'all');
lims.min_P = -0.01;
lims.max_P = max(plt.phyto,[],'all');
lims.min_Z = -0.01;
lims.max_Z = max(plt.zoop,[],'all');
lims.max_MLD = max(MLD_rearr,[],'all');
lims.min_MLD = min(MLD_rearr,[],'all');
lim.max_scale = 1;
lim.min_scale = 0;
if params.detritus_layer
    lims.min_D = -0.01;
    lims.max_D = max(plt.detri,[],'all');
end

% Change temperature colour map
% [len_map,~] = size(col.map.temp);
% zero_index = ceil((0-lims.min_T) / (lims.max_T - lims.min_T) * len_map);
% if zero_index<=0
%     zero_index = 1;
%     [x_lm,y_lm] = ind2sub([params.num_lats,params.num_longs],find(params.landmass));
%     for i = x_lm
%         for j = y_lm
%             temp_rearr(:,i,j) = lims.min_T -1/100;
%         end
%     end
% end
% col.map.temp(zero_index,:) = [0 0 0];

% Preallocate
num_plots = length(vars_to_visualise)+1;
plts = cell(1,num_plots+1);
axs = cell(1,num_plots+1);

% Prepare gif
ocean_grid = figure(1);
ocean_grid.WindowState = 'maximized';
%axis tight manual % ensures getframe() returns a consistent size
%ocean_grid.Visible = 'off';
gifpath = [folder_path '/simulation.gif'];
for day = 1:min(params.days(end)+1,length(params.years_per_day))
    
    % Save each figure as a frame in the gif
    [axs,plts] = draw_frame(params,plt,temp_rearr,v_north_rearr,...
        v_east_rearr,lims,day,axs,plts,vars_to_visualise,...
        long_centres,lat_centres,col);
    % Save the gif frame.
    frame = getframe(ocean_grid);
    [A,map] = rgb2ind(frame2im(frame),256);
    if day == 1
        imwrite(A,map,gifpath,'gif',"LoopCount",Inf,"DelayTime",1/params.vid_speed);
    else
        imwrite(A,map,gifpath,'gif',"WriteMode","append","DelayTime",1/params.vid_speed);
    end
end