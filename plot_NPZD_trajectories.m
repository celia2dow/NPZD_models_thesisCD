% Code plotting chosen trajectories of NPZ(D) concentrations. Called by 
% NPZD_model_forcing.m

for site = 1:num_sites
    fig = figure(3 + site);
    i = sites_to_visualise(site,1);
    j = sites_to_visualise(site,2);
    name_site = sprintf('i%i_j%i',i,j);     % The plot for site (i,j)

    condit = months>0;
    non0times = record.times(record.times(:)>0);
    time_lim = length(squeeze(non0times));
    %vols_ij = interp1(params.days, squeeze(params.volumes(i,j,:)),...
    %    squeeze(non0times));
    months_ij = months(condit);
    if ~params.temp_vert_mix && ~params.current
        months_ij = record.times(condit);
    end
    plot(months_ij,squeeze(record.soluts(i,j,1,condit)), ...
        'Color', col.nutri_orange);              % Nutrient
    hold on
    ylabel("Concentration (g C/{m^3})")
    plot(months_ij,squeeze(record.soluts(i,j,2,condit)), ...
        'Color', col.phyto_green);               % Phytoplankton
    plot(months_ij,squeeze(record.soluts(i,j,3,condit)), ...
        'Color', col.zoop_dark);                 % Zooplankton
    if params.detritus_layer
        plot(months_ij,squeeze(record.soluts(i,j,4,condit)), ...
            'Color', col.detri_purp);            % Detritus
        legend('Mixed layer nutrient (N)', 'Phytoplankton (P)', ...
            'Zooplankton (Z)', 'Detritus (D)')
    else
        legend('Mixed layer nutrient (N)', 'Phytoplankton (P)', ...
            'Zooplankton (Z)')
    end
    hold off
    if ~params.temp_vert_mix && ~params.current
        xlabel("Time (days)")
        xlim([0,num_days*params.repeat_num])
    else
        xticks(tcks);
        xticklabels(tck_seas);
        xtickangle(90);
        xlabel("Time")
    end
    % ylim([0,0.6]) %%%% COMMENT THIS OUT
    if sum([params.temp_vert_mix,params.change_mix_depth])~=0
        sgtitle(sprintf("Solutions for %4.2f^{\\circ}N, %4.2f^{\\circ}E",...
            lat_centres(i),long_centres(j)))
    end
    savefig(fig, [folder_path '/State_variables_over_time_' name_site], 'compact')
    saveas(fig, [folder_path '/State_variables_over_time_' name_site], 'png')
end