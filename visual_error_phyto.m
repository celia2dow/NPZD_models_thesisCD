% Go through each day in time and save the spatial distributuion of
% temperature, and phytoplankton
% Interpolate produced data to the days for which we have temperature data
visual_switch=0;
if visual_switch
    % Prepare gif
    ocean_grid = figure(1);
    axis tight manual % ensures getframe() returns a consistent size
    %ocean_grid.Visible = 'off';
    gifpath = [folder_path '/simulation.gif'];
    phyto = zeros(params.num_lats,params.num_longs,num_days);
   
    % Find limits
    lims.min_T = min(temps_av(temps_av>0),[],'all');
    lims.max_T = max(temps_av,[],'all');
    lims.min_P = 0;
    lims.max_P = 1;
    
    % Preallocate
    plts = cell(1,3);
    axs = cell(1,3);
    for day = params.days
        day = day+1;
        % Save each figure as a frame in the gif
        [axs,plts] = draw_frame(params,phyto,temps_av,v_north_av,v_east_av,lims,day,axs,plts);
        % Save the gif frame.
        frame = getframe(ocean_grid);
        [A,map] = rgb2ind(frame2im(frame),256);
        if day == 1
            imwrite(A,map,gifpath,'gif',"LoopCount",Inf,"DelayTime",1/params.vid_speed);
        else
            imwrite(A,map,gifpath,'gif',"WriteMode","append","DelayTime",1/params.vid_speed);
        end
    end
end

% NPZ(D) plot
months = record.times./(365.25/12);
fig2 = figure(2);
season1 = ceil(month1/3);
if params.hemisphere == "N"
    strt = 3;
else
    strt = 1;
end
year_add=0;

months_final = months(months>0);
months_final = months_final(end);
tcks = 0:3:months_final;
num_ticks = length(tcks);
tck_seas = cell(1,num_ticks);
for tck = season1:num_ticks+season1
    indx = mod(tck,4);
    indx(indx == 0) = 4;
    if tck ==season1
        tck_seas{tck} = num2str(year1) + " " + season_label{indx};
    elseif indx == 1
        year_add=year_add+1;
        tck_seas{tck} = num2str(year1 + year_add) + " " + season_label{indx};
    else
        tck_seas{tck}=season_label{indx};
    end
end 

for i = 1:params.num_lats
    for j = 1:params.num_longs
        condit = months>0;
        non0times = record.times(record.times(:)>0);
        time_lim = length(squeeze(non0times));
        vols_ij = interp1(params.days, squeeze(params.volumes(i,j,:)),...
            squeeze(non0times));
        months_ij = months(condit);
        plot(months_ij,squeeze(record.soluts(i,j,1,condit)), 'b');
        hold on
        xlabel("Time (months)")
        ylabel("Concentration (g C/{m^3})")
        plot(months_ij,squeeze(record.soluts(i,j,2,condit)), 'k');
        plot(months_ij,squeeze(record.soluts(i,j,3,condit)), 'r');
        if params.detritus_layer
            plot(months_ij,squeeze(record.soluts(i,j,4,condit)), 'g');
        end
        legend('Mixed layer nutrient (N)', 'Phytoplankton (P)', ...
                'Zooplankton (Z)')
        xticks(tcks);
        xticklabels(tck_seas);
        xtickangle(90);
        %xticklabels({'x = 0','x = 5','x = 10'});
        ylim([0,0.6])
        %%%% COMMENT THIS OUT
        if params.detritus_layer
            legend('Mixed layer nutrient (N)', 'Phytoplankton (P)', ...
                'Zooplankton (Z)', 'Detritus (D)')
        end
    end
end
hold off
savefig(fig2, [folder_path '/State_variables_over_time'], 'compact')
saveas(fig2, [folder_path '/State_variables_over_time'], 'png')