% Plot time series trajectories with burn-in
fig3000=figure(3000);
tiles = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
axs1 = nexttile(tiles);


% PLOT WITH BURN-IN PERIOD

% burn in year plots
plot(x_vals1(condit_BI),squeeze(record.soluts(1,1,1,condit_BI)), ...
    'Color', col.nutri_orange);              % Nutrient
hold on
plot(x_vals1(condit_BI),squeeze(record.soluts(1,1,2,condit_BI)), ...
    'Color', col.phyto_green);               % Phytoplankton
plot(x_vals1(condit_BI),squeeze(record.soluts(1,1,3,condit_BI)), ...
    'Color', col.zoop_dark);                 % Zooplankton
% final year plots
plot(x_vals1(condit_afterBI),squeeze(record.soluts(1,1,1,condit_afterBI)), ...
    'Color', col.nutri_orange, 'LineWidth', 2);              
plot(x_vals1(condit_afterBI),squeeze(record.soluts(1,1,2,condit_afterBI)), ...
    'Color', col.phyto_green, 'LineWidth', 2);               
plot(x_vals1(condit_afterBI),squeeze(record.soluts(1,1,3,condit_afterBI)), ...
    'Color', col.zoop_dark, 'LineWidth', 2);   
if params.detritus_layer
    plot(x_vals1(condit_BI),squeeze(record.soluts(1,1,4,condit_BI)), ...
        'Color', col.detri_purp);           % Detritus
    plot(x_vals1(condit_afterBI),squeeze(record.soluts(1,1,4,condit_afterBI)), ...
        'Color', col.detri_purp, 'LineWidth', 2);  
end

ylabel("Concentration (g C/{m^3})")
xlabel("Time (years)")
xlim([0,num_years*params.repeat_num])

%num_seasons = num_years*4;
%xticks(tcks(1:num_seasons+1));
%xticklabels(tck_seas(end-num_seasons-1:end));
%xtickangle(90);


% PLOT WITHOUT BURN-IN PERIOD
axs2 = nexttile(tiles);

plot(x_vals2,squeeze(record.soluts(1,1,1,condit_afterBI)), ...
    'Color', col.nutri_orange, 'LineWidth', 2); % Nutrient
hold on
plot(x_vals2,squeeze(record.soluts(1,1,2,condit_afterBI)), ...
    'Color', col.phyto_green, 'LineWidth', 2);  % Phytoplankton
plot(x_vals2,squeeze(record.soluts(1,1,3,condit_afterBI)), ...
    'Color', col.zoop_dark, 'LineWidth', 2);    % Zooplankton
if params.detritus_layer
    plot(x_vals2,squeeze(record.soluts(1,1,4,condit_afterBI)), ...
        'Color', col.detri_purp, 'LineWidth', 2); % Detritus           
    leg = legend('Mixed layer nutrient (N)', 'Phytoplankton (P)', ...
        'Zooplankton (Z)', 'Detritus (D)');
    num_cols = 4;
else
    leg= legend('Mixed layer nutrient (N)', 'Phytoplankton (P)', ...
        'Zooplankton (Z)');
    num_cols = 4;
end

leg.NumColumns = num_cols;
leg.Location= 'southoutside';

ylabel("Concentration (g C/{m^3})")
xlabel("Time (days)")
xlim([0,num_days])

x0=100;
y0=100;
width=1050;
height=400;
set(gcf,'position',[x0,y0,width,height])

savefig(fig3000, [folder_path_driver '/Example_burn_in__' str_label2], 'compact')
saveas(fig3000, [folder_path_driver '/Example_burn_in__' str_label2], 'png')