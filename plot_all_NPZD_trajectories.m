% Code plotting all trajectories of NPZ(D) concentrations. Called by 
% NPZD_model_forcing.m

fig2 = figure(2);

for i = 1:params.num_lats
    for j = 1:params.num_longs
        condit = months>0;
        non0times = record.times(record.times(:)>0);
        time_lim = length(squeeze(non0times));
        %vols_ij = interp1(params.days(1:num_days), squeeze(params.volumes(i,j,:)),...
        %    squeeze(non0times));
        months_ij = months(condit);
        if ~params.temp_vert_mix && ~params.current
            months_ij = record.times(condit);
        end
        plot(months_ij,squeeze(record.soluts(i,j,1,condit)),...
            'Color', col.nutri_orange);              % Nutrient
        hold on
        plot(months_ij,squeeze(record.soluts(i,j,2,condit)), ...
            'Color', col.phyto_green);               % Phytoplankton
        plot(months_ij,squeeze(record.soluts(i,j,3,condit)), ...
            'Color', col.zoop_dark);                 % Zooplankton
        if params.detritus_layer
            plot(months_ij,squeeze(record.soluts(i,j,4,condit)), ...
                'Color', col.detri_purp);            % Detritus
        end
    end
end
if params.temp_vert_mix || params.current
    xticks(tcks);
    xticklabels(tck_seas);
    xtickangle(90);
    xlabel("Time (months)")
else
    xlabel("Time (days)")
end
ylabel("Concentration (g C/{m^3})")
ylim([0,0.6]) %%%% COMMENT THIS OUT
hold off
legend('Mixed layer nutrient (N)', 'Phytoplankton (P)', ...
    'Zooplankton (Z)')
if params.detritus_layer
    legend('Mixed layer nutrient (N)', 'Phytoplankton (P)', ...
        'Zooplankton (Z)', 'Detritus (D)')
end
savefig(fig2, [folder_path '/State_variables_over_time_ALLsites'], 'compact')
saveas(fig2, [folder_path '/State_variables_over_time_ALLsites'], 'png')