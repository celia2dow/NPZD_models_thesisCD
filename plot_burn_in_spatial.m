% Plot time series trajectories with burn-in
c.min = squeeze(min(record.soluts(:,:,state_var,:),[],[1,2]));
c.max = squeeze(max(record.soluts(:,:,state_var,:),[],[1,2]));
c.mean = squeeze(mean(record.soluts(:,:,state_var,:),[1,2]));
c.std = squeeze(std(record.soluts(:,:,state_var,:),0,[1,2]));

% PLOT WITH BURN-IN PERIOD

% burn in year plots
axs{(state_var-1)*2+1} = nexttile(tiles,(state_var-1)*2+1);
plot(x_vals1(condit_BI),c.mean(condit_BI), 'LineStyle', '-', ...
    'Color', cols{state_var}); % Mean
hold on
% plot(x_vals1(condit_BI),c.mean(condit_BI)+log(c.std(condit_BI)), 'LineStyle', '--', ...
%     'Color', cols{state_var}); % Mean + stand. dev.
% plot(x_vals1(condit_BI),c.mean(condit_BI)-log(c.std(condit_BI)), 'LineStyle', '--', ...
%     'Color', cols{state_var}, 'HandleVisibility','off'); % Mean - stand. dev.
% final year plots
plot(x_vals1(condit_afterBI),c.mean(condit_afterBI), 'LineStyle', '-', ...
    'Color', cols{state_var}, 'LineWidth', 2); % Mean
% plot(x_vals1(condit_afterBI),c.mean(condit_afterBI)+log(c.std(condit_afterBI)), 'LineStyle', '--', ...
%     'Color', cols{state_var}, 'LineWidth', 2); % Mean + stand. dev.
% plot(x_vals1(condit_afterBI),c.mean(condit_afterBI)-log(c.std(condit_afterBI)), 'LineStyle', '--', ...
%     'Color', cols{state_var}, 'LineWidth', 2, 'HandleVisibility','off'); % Mean - stand. dev.
hold off
if state_var ~=num_plots
    set(gca,'XTick',[])
else
    xlabel("Time (years)");
end

ylabel("Concentration (g C/{m^3})")
xlim([0,num_years*params.repeat_num])

%num_seasons = num_years*4;
%xticks(tcks(1:num_seasons+1));
%xticklabels(tck_seas(end-num_seasons-1:end));
%xtickangle(90);


% PLOT WITHOUT BURN-IN PERIOD
axs{state_var*2} = nexttile(tiles,state_var*2);
plot(x_vals2,c.mean(condit_afterBI), 'LineStyle', '-', ...
    'Color', cols{state_var}, 'LineWidth', 2); % Mean
% hold on
% plot(x_vals2,c.mean(condit_afterBI)+log(c.std(condit_afterBI)), 'LineStyle', '--', ...
%     'Color', cols{state_var}, 'LineWidth', 2); % Mean + stand. dev.
% plot(x_vals2,c.mean(condit_afterBI)-log(c.std(condit_afterBI)), 'LineStyle', '--', ...
%     'Color', cols{state_var}, 'LineWidth', 2, 'HandleVisibility','off'); % Mean - stand. dev.
% hold off
if state_var ~=num_plots
    set(gca,'XTick',[])
else
    xlabel("Time (days)");
end
leg = legend(labels{state_var});%,labels{state_var+num_plots});
leg.Location= 'best';

xlim([0,num_days])

