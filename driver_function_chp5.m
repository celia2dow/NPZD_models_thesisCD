% Code creating the GIF of temperature, current, and desired variables.
% Called by NPZD_model_forcing.m

clear
close all

params.current = 1;                 % Inclusion of current in model
params.diffus = 1;                  % Inclusion of diffusion in model
params.temp_vert_mix = 0;
params.change_mix_depth = 0;
params.detritus_layer = 1;
params.d = 1.5;
params.omega = 0; %0.5; 
params.repeat_num = 4;
smoothing = 7;          
thermocline_calc = 0;
params.BC = 1;                      % Boundary conditions of spatial model: 0 periodic, 1 open

lil_str = ['detritus' num2str(params.detritus_layer) 'd' num2str(params.d)];
%%
NPZD_model_forcing

%%
% FOLDER PATH
folder_name = 'Final_gifs';
folder_path_driver = [pwd '/' date '/' folder_name];
% The folder paths used agrees with a macOS system and may need to be
% flipped to '\' for a Windows system
if ~exist(folder_path_driver, 'dir')
    mkdir(folder_path_driver)
end

save([folder_path_driver '/Workspace_open_noStrat_noMLDchange_NPZ_d1.mat'])
% Interpolate NPZ(D) data and rearrange temperature, TSI and MLD data for 
% days after the burn in
%num_days = 100 ; %100; %
condit_afterBI = record.times >= cut_off_day; %0; %
days_afterBI=record.times(condit_afterBI);
days_afterBI_shift = days_afterBI-days_afterBI(1);

params.landmass_rearr = zeros(params.num_lats, params.num_longs);
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
                squeeze(record.soluts(i,j,1,condit_afterBI)),params.days(1:num_days));
            plt.phyto(params.num_lats+1-i,j,:) = interp1(days_afterBI_shift,...
                squeeze(record.soluts(i,j,2,condit_afterBI)),params.days(1:num_days));
            plt.zoop(params.num_lats+1-i,j,:) = interp1(days_afterBI_shift,...
                squeeze(record.soluts(i,j,3,condit_afterBI)),params.days(1:num_days));
            if params.detritus_layer
                plt.detri(params.num_lats+1-i,j,:) = interp1(days_afterBI_shift,...
                    squeeze(record.soluts(i,j,4,condit_afterBI)),params.days(1:num_days));
            end
            if ~(params.num_lats ==1)
                temp_rearr(:,params.num_lats+1-i,j) = spline(GOODindices_days_T-1, ...
                    temps_av(GOODindices_days_T,i,j), params.days(1:num_days));
                v_north_rearr(:,params.num_lats+1-i,j) = spline(GOODindices_days_CN-1, ...
                    v_north_av(GOODindices_days_CN,i,j), params.days(1:num_days));
                v_east_rearr(:,params.num_lats+1-i,j) = spline(GOODindices_days_CE-1, ...
                    v_east_av(GOODindices_days_CE,i,j), params.days(1:num_days));
                MLD_rearr(:,params.num_lats+1-i,j)= spline(GOODindices_days_T-1, ...
                    depths_LEVs(GOODindices_days_T,i,j), params.days(1:num_days));
                TSI_rearr(:,params.num_lats+1-i,j)= spline(GOODindices_days_T-1, ...
                    store_2011(GOODindices_days_T,i,j), params.days(1:num_days));
            else
                temp_rearr(:,params.num_lats+1-i,j) = spline(GOODindices_days_T-1, ...
                    temps_av(GOODindices_days_T), params.days(1:num_days));
                v_north_rearr(:,params.num_lats+1-i,j) = spline(GOODindices_days_CN-1, ...
                    v_north_av(GOODindices_days_CN), params.days(1:num_days));
                v_east_rearr(:,params.num_lats+1-i,j) = spline(GOODindices_days_CE-1, ...
                    v_east_av(GOODindices_days_CE), params.days(1:num_days));
                MLD_rearr(:,params.num_lats+1-i,j)= spline(GOODindices_days_T-1, ...
                    depths_LEVs(GOODindices_days_T), params.days(1:num_days));
                TSI_rearr(:,params.num_lats+1-i,j)= spline(GOODindices_days_T-1, ...
                    store_2011(GOODindices_days_T), params.days(1:num_days));
            end
        else                            % If there is a landmass at (i,j)
            params.landmass_rearr(params.num_lats+1-i,j) = 1;
        end
    end
end
scales_rearr = params.strat_func(TSI_rearr);
[x_landmass, y_landmass] = find(params.landmass_rearr);

% Find limits
lims.min_T = min(temps_av(temps_av~=0),[],'all') - 1;
lims.max_T = max(temps_av,[],'all');
lims.min_N = -0.01;
lims.max_N = 1; %max(plt.nutri,[],'all'); %
lims.min_P = -0.01;
lims.max_P = 1; %max(plt.phyto,[],'all'); %
lims.min_Z = -0.01;
lims.max_Z = max(plt.zoop,[],'all');
lims.max_MLD = max(MLD_rearr,[],'all');
lims.min_MLD = min(MLD_rearr,[],'all')-10;
lims.max_scale = 1;
lims.min_scale = -0.01;
if params.detritus_layer
    lims.min_D = -0.01;
    lims.max_D = 1; %max(plt.detri,[],'all');
end

for site = 1:length(x_landmass)
    plt.nutri(x_landmass(site),y_landmass(site),:) = -100;
    plt.phyto(x_landmass(site),y_landmass(site),:) = -100;
    plt.zoop(x_landmass(site),y_landmass(site),:) = -100;
    if params.detritus_layer
        plt.detri(x_landmass(site),y_landmass(site),:) = -100;
    end
    temp_rearr(:,x_landmass(site),y_landmass(site)) = -100;
    v_north_rearr(:,x_landmass(site),y_landmass(site)) = 0;
    v_east_rearr(:,x_landmass(site),y_landmass(site)) = 0;
    MLD_rearr(:,x_landmass(site),y_landmass(site))= -100;
    scales_rearr(:,x_landmass(site),y_landmass(site))=-100;
end

%%
%%% MAKE GIFS

close all
%create_enviro_gif

close all
create_NPZD_gif
%%

%%%% Plots of mean concentrations over time with burnin
x_vals1 = record.times/num_days; % In terms of years
condit_BI = x_vals1 <= cut_off_day;
condit_afterBI = record.times >= cut_off_day;
days_afterBI=record.times(condit_afterBI);
x_vals2 = days_afterBI-cut_off_day; % In terms of days
if params.detritus_layer
    num_plots = 4;
else
    num_plots = 3;
end

fig3000=figure(3000);
off_on = {'off','on'};
str_label2 = ['vertMix_' off_on{params.temp_vert_mix+1} ...
    '__changeMLD_' off_on{params.change_mix_depth+1}];
tiles = tiledlayout(num_plots,2,'TileSpacing','Compact','Padding','Compact');
axs = cell(1,num_plots*2);
cols = {col.nutri_orange, col.phyto_green, col.zoop_dark, col.detri_purp};
labels = {'\mu_N', '\mu_P', '\mu_Z', '\mu_D',...
    ['\mu_N ' char(177) ' ln(\sigma_N)'],...
    ['\mu_P ' char(177) ' ln(\sigma_P)'],...
    ['\mu_Z ' char(177) ' ln(\sigma_Z)'],...
    ['\mu_D ' char(177) ' ln(\sigma_D)']};

for state_var = 1:num_plots
    plot_burn_in_spatial
end

x0=10;
y0=10;
width=1050;
height=500;
set(gcf,'position',[x0,y0,width,height])
savefig(fig3000, [folder_path_driver '/Spatial_burn_in__' str_label2], 'compact')
saveas(fig3000, [folder_path_driver '/Spatial_burn_in__' str_label2], 'png')

% Plots of max and min concentrations over time in final year
% plot(x_vals1(condit_BI),c.max(condit_BI), 'LineStyle', ':', ...
%     'Color', col.nutri_orange, 'Parent', axs{state_var}); % Max value
% plot(x_vals1(condit_BI),c.max(condit_BI), 'LineStyle', ':', ...
%     'Color', col.nutri_orange, 'Parent', axs{state_var}); % Min value