% DRIVER FUNCTION FOR TEMP-RELATED PLOTS FOR CHAPTER 4
clear
close all

% FOLDER PATH
folder_name = 'ExtraThings_CHP4';
folder_path_driver = [pwd '/' date '/' folder_name];
% The folder paths used agrees with a macOS system and may need to be
% flipped to '\' for a Windows system
if ~exist(folder_path_driver, 'dir')
    mkdir(folder_path_driver)
end

params.temp_vert_mix = 0;
params.change_mix_depth = 0;
params.detritus_layer = 0;
params.d = 1;
params.omega = 0.5; 
params.repeat_num = 1;
smoothing = 1;          
thermocline_calc = 3;

NPZD_model_forcing

day = 261;
depths = data_file_temps.LEV;
temps_11_day = temps(day,:);
spline_11_day = spline(depths, temps_11_day);
temp_10m = temps_11_day(depths == 10);
all_depths = 0:depths(end);
temp_11_interp = ppval(spline_11_day, all_depths);
difs = abs(temp_10m - temp_11_interp);
poss_MLD = all_depths(difs>tol);
MLD = min(poss_MLD);
temp_MLD = temp_11_interp(all_depths == MLD);

% CALCULAIING THERMOCLINE GRAPH
fig18 = figure(18);
tiles = tiledlayout(3,2,'TileSpacing','Compact','Padding','Compact');
ax1 = nexttile(tiles,[3,1]);
plot(temp_11_interp, all_depths,'k-')
grid on
set(gca, 'YDir','reverse')
hold on
plot(temps_11_day, depths, 'Color', ...
    [0.4667, 0.6745, 0.1882], 'LineStyle', 'none', 'Marker', '.', 'LineWidth',2)
plot(temp_10m,10,'Color',[0.0275, 0.4471, 0.9294], ...
    'LineStyle', 'none', 'Marker', '+', 'LineWidth',2,  ...
    'MarkerSize',5,'HandleVisibility','off')
plot(temp_MLD,MLD,'Color',[0.8510, 0.3255, 0.0980], ...
    'LineStyle', 'none', 'Marker', '+',  'LineWidth',2,  ...
    'MarkerSize',5,'HandleVisibility','off')
x = [23.5, 24, 24, 23.5, 23.5];
y = [0, 0, 30, 30, 0];
plot(x, y, 'Color',[0.4667, 0.6745, 0.1882], 'LineWidth', 2);
hold off
xlim([min(temps_11_day)-0.5,max(temps_11_day)+0.5])
ylim([0,depths(end)])
set(gca,'XAxisLocation','top');
xlabel('Temperature (\circ C)')
ylabel('Depth (m)')
leg = legend("Cubic spline", "Data");
leg.Location = "southoutside";
leg.NumColumns = 2;

ax2 = nexttile(tiles,[2,1]);
plot(temp_11_interp, all_depths,'k-')
hold on
plot([23.5,temp_10m], [10,10], 'Color', ...
    [0.6510, 0.6510, 0.6510], 'LineStyle', '--', 'HandleVisibility','off') % The 10m depth line
text(23.55,8, '10 m', 'Color',[0.0275, 0.4471, 0.9294])
plot([temp_10m,temp_10m], [0,10], 'Color', ...
    [0.6510, 0.6510, 0.6510], 'LineStyle', '--', 'HandleVisibility','off') % The temp at 10m line
text(temp_10m-0.06,5,'T_{10 m}', 'Color',[0.0275, 0.4471, 0.9294])
plot([23.5,temp_MLD], [MLD,MLD], 'Color', ...
    [0.8510, 0.3255, 0.0980], 'LineStyle', '--', 'HandleVisibility','off') % The MLD depth line
text(23.55,MLD-2,'MLD', 'Color',[0.8510, 0.3255, 0.0980])
plot([temp_MLD,temp_MLD], [0,MLD], 'Color', ...
    [0.8510, 0.3255, 0.0980], 'LineStyle', '--', 'HandleVisibility','off') % The temp at 10m line
text(temp_MLD-0.13,5,'T_{10 m}+\Delta T','Color',[0.8510, 0.3255, 0.0980])
plot(temps_11_day, depths, 'Color', ...
    [0.4667, 0.6745, 0.1882], 'LineStyle', 'none', 'Marker', '.')
plot(temp_10m,10,'Color',[0.0275, 0.4471, 0.9294], ...
    'LineStyle', 'none', 'Marker', '+', 'HandleVisibility','off', ...
    'LineWidth',2,  'MarkerSize',10)
plot(temp_MLD,MLD,'Color',[0.8510, 0.3255, 0.0980], ...
    'LineStyle', 'none', 'Marker', '+', 'HandleVisibility','off', ...
    'LineWidth',2, 'MarkerSize',10)
ax2.YColor = [0.4667, 0.6745, 0.1882];
ax2.LineWidth = 2;
ax2.XColor = [0.4667, 0.6745, 0.1882];
hold off
xlim([23.5,24])
grid on
set(gca, 'YDir','reverse')
set(gca,'XAxisLocation','top','YAxisLocation','right');
xlabel('Temperature (\circ C)')
ylabel('Depth (m)')

%%%%%%%%%%% SAVE FIGURE %%%%%%%%%%%%%%%%%
savefig(fig18, [folder_path_driver '/temps_per_depth'], 'compact')
saveas(fig18, [folder_path_driver '/temps_per_depth'], 'png')




% VISUALISING MLD ON COLORMAP 
fig19=figure(19);
more_temps = zeros(length(all_depths),num_days);
for day = params.days+1
    temps_11_day = temps(day,:);
    spline_11_day = spline(depths, temps_11_day);
    more_temps(:,day) = ppval(spline_11_day,all_depths);
end
x0=50;
y0=50;
width=1250;
height=1250;
set(gcf,'position',[x0,y0,width,height])
ax = gca;

imagesc(more_temps);
colorbar;
a=colorbar;
a.Label.String = 'Temperature (\circ C)';
a.Label.FontSize = 10;
hold on


yticks(depths);
ax.YColor = [0.4667, 0.6745, 0.1882];
axMLD= axes('Visible','off','HandleVisibility','off'); 
plot(params.days, ppval(params.MLD_CS{1,1},params.days),'Color','k','LineWidth',2)
plot(params.days, therm_base_2011(:,1,1),'Color','w','LineWidth',2)
ax.FontSize = 6;
ylabel('Depth (m)','FontSize',10)
xlabel('Time (days)','FontSize',10)
text(225,45,'MLD','FontWeight','bold')
text(225,125,'Thermocline base','FontWeight','bold','Color','w')
hold off
coords = [lat_centres(1),long_centres(1)];
if lat_centres(1)<0
    NS = 'S';
else
    NS = 'N';
end
if long_centres(1)<0
    EW = 'W';
else
    EW = 'E';
end
addon = [num2str(abs(floor(100*coords(1)))) NS num2str(abs(floor(100*coords(2)))) EW];

%%%%%%%%%%% SAVE FIGURE %%%%%%%%%%%%%%%%%
savefig(fig19, [folder_path_driver '/temps_per_day' addon], 'compact')
saveas(fig19, [folder_path_driver '/temps_per_day' addon], 'png')



% VISUALISING NOISY MLD VS SMOOTH MLD
fig20=figure(20);
plot(params.days, ppval(params.MLD_CS{1,1},params.days),'Color','k','LineStyle','-.')
hold on
plot(params.days, smooth(depths_LEVs(GOODindices_days_T,1,1),7),'Color','k')
set(gca, 'YDir','reverse')
ylabel('Depth (m)')
xlabel('Time (days)')
leg = legend('Original calculated MLD','7-day moving average MLD');
leg.Location = 'southoutside';
leg.NumColumns = 2;
hold off

%%%%%%%%%%% SAVE FIGURE %%%%%%%%%%%%%%%%%
savefig(fig20, [folder_path_driver '/noisy_Vs_smooth_MLD' addon], 'compact')
saveas(fig20, [folder_path_driver '/noisy_Vs_smooth_MLD' addon], 'png')



% VISUALISING NOISY TEMP AND TSI AND THE SIGMOIDAL SCALING
tiles1 = tiledlayout(3,4,'TileSpacing','compact','Padding','Compact');
nexttile(tiles1,[1,3])
plot(params.days, ppval(params.temp_CS{1,1},params.days),'Color',[0.6350 0.0780 0.1840])
ylabel('Temperature (\circ C)')
xlim([params.days(1),params.days(end)])
set(gca,'XTick',[])

nexttile(tiles1,[1,3])
plot(params.days, ppval(params.TSI_CS{1,1},params.days),'Color','r')
hold on
plot([params.days(1),params.days(end)], (params.TSI_mean-params.TSI_std_dev).*ones(1,2),...
    'Color','r','LineStyle', '--')
text(220,params.TSI_mean-params.TSI_std_dev+60,'\mu_{TSI,10yr} - \sigma_{TSI,10yr}', 'Color', 'r')
hold off
ylabel("TSI (\circ C m)")
low_TSI = min(ppval(params.TSI_CS{1,1},params.days));
high_TSI = max(ppval(params.TSI_CS{1,1},params.days));
xlim([params.days(1),params.days(end)])
ylim([low_TSI,high_TSI])
xlim([0,num_days])
num_seasons = num_years*4;
xticks(tcks(1:num_seasons+1)*num_days/12);
xticklabels(tck_seas(end-num_seasons-1:end));
set(gca,'XAxisLocation','top');

nexttile(tiles1,8)
TSI_vals = low_TSI:0.01:high_TSI;
plot(params.strat_func(TSI_vals),TSI_vals,'Color',[0.8510, 0.3255, 0.0980])
xlim([-0.2,1.2])
ylim([low_TSI,high_TSI])
text(0.1,high_TSI-20,['no' newline 'stratification'], 'Color', [0.8510, 0.3255, 0.0980])
text(0.1,low_TSI+20,['strong' newline 'stratification'], 'Color', [0.8510, 0.3255, 0.0980])
xlabel("Scaling")
set(gca,'YTick',[])

nexttile(tiles1,[1,3])
plot(params.days, params.strat_func(ppval(params.TSI_CS{1,1},params.days)),'Color',[0.8510, 0.3255, 0.0980])
xlim([params.days(1),params.days(end)])
xlim([params.days(1),params.days(end)])
ylim([-0.2,1.2])
ylabel("Scaling")
xlabel('Time (days)')

set(gcf,'position',[x0,y0,600,600])


%%%%%%%%%%% SAVE FIGURE %%%%%%%%%%%%%%%%%
savefig(fig20, [folder_path_driver '/TSI_plots' addon], 'compact')
saveas(fig20, [folder_path_driver '/TSI_plots' addon], 'png')