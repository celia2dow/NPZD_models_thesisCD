function [axs,plts] = draw_frame2(params,phyto,temps,lims,day,axs,plts)
% DRAW_FRAME Save each figure as a frame in the movie illustrating the 
% density of phytoplankton in a region of the ocean

% Month (season) label
month_num = params.months_per_day(day);
month_label = params.month_season(1,month_num);
season_label = params.month_season(2,month_num);
str_label = month_label + " (" + season_label + ")";

% Draw a map of the temperature distribution
if day == 1
    axs{1}=subplot(1,2,1);
    hold on
    if params.num_lats == 1
        plts{1} = imagesc(squeeze(temps(day)),[lims.min_T lims.max_T]);
    else
        plts{1} = imagesc(squeeze(temps(day,:,:)),[lims.min_T lims.max_T]);
    end
    colormap(axs{1},cool);
    a=colorbar;
    a.Label.String = 'Temperature (\circ C)';
    title(str_label);
    hold off

    axs{2}=subplot(1,2,2);
    hold on
    plts{2} = imagesc(squeeze(phyto(:,:,day)), [lims.min_P lims.max_P]);  
    colormap(axs{2},winter);
    b=colorbar;
    b.Label.String = 'Phytoplankton concentration (g C m^{-3})';
    hold off
else
    if params.num_lats == 1
        set(plts{1}, 'CData', squeeze(temps(day)));
    else
        set(plts{1}, 'CData', squeeze(temps(day,:,:)));
    end
    set(axs{1}.Title, 'String', str_label);
    set(plts{2}, 'CData', squeeze(phyto(:,:,day)));
end

sgtitle(sprintf('Day %d',day));
drawnow;
