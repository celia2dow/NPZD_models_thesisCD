function [axs,plts] = draw_frame_NPZD(params,plt,temps,v_north,v_east,...
    lims,day,axs,plts,vars_to_visualise,long_centres,lat_centres,col)
% DRAW_FRAME Save each figure as a frame in the gif illustrating the
% density of phytoplankton in a region of the ocean

num_uniq_plots = length(plts); % Number of unique plots including temperature and flow field

% Month (season) label
month_num = params.months_per_day(day);
month_label = params.month_season(1,month_num);
season_label = params.month_season(2,month_num);
year = params.years_per_day(params.num_days*(params.repeat_num-1)+1);
%if year<params.year1
%    str_label = "BI year " + num2str(year) + ", " + month_label + " (" + season_label + ")" + "     Day " + num2str(day);
%else
%    str_label = num2str(year) + ", " + month_label + " (" + season_label + ")"+ "     Day " + num2str(day);
%end
str_label = num2str(year) + ", " + month_label + " (" + season_label + ")"+ "     Day " + num2str(day);
day_edited = mod(day,params.num_days);
if day_edited==0
    day_edited = params.num_days;
end

% Flow field prep
x = 1:params.num_longs;
y = 1:params.num_lats;
scaleFactor = 0.5; 

if day == 1
    % Construct tiled layout
    t1 = tiledlayout(1,4,'TileSpacing','Compact','Padding','Compact');

    % Zooplankton
    axs{1} = nexttile(t1);
    plts{1} = imagesc(squeeze(plt.zoop(:,:,day)), [lims.min_Z lims.max_Z]);  
    hold on
    colormap(axs{1},col.map.zoop);
    % Add invisible quiver plot to set the same aspect ratio
    axs{2} = axes('Visible','off','HandleVisibility','off'); 
    % Plot the quiver data
    quiver ( x , y , 0 .* squeeze(v_east(day,:,:)), ...
        0 .* squeeze(v_north(day,:,:)), scaleFactor, 'LineWidth', 2);
    hold off
    % Set properties for main axis
    axis equal
    a=colorbar;
    a.Location ='southoutside';
    a.Label.String = 'Zooplankton concentration (g C m^{-3})';
    set(a, 'ylim', [lims.min_Z+0.01 lims.max_Z])
    %title(str_label);
    xlim([0.5,params.num_longs+0.5])
    ylim([0.5,params.num_lats+0.5])
    y_ticks = 1:params.num_lats;
    yticks(y_ticks(5:5:end))
    yticklabels(flip(lat_centres(5:5:end)) + "\circ N")
    xticklabels(abs(long_centres) + "\circ W")

    % Nutrient
    axs{3} = nexttile(t1);
    plts{2} = imagesc(squeeze(plt.nutri(:,:,day)), [lims.min_N lims.max_N]);  
    hold on
    colormap(axs{3},col.map.nutri);
    b=colorbar;
    b.Location ='southoutside';
    b.Label.String = 'Nutrient concentration (g C m^{-3})';
    set(b, 'ylim', [lims.min_N+0.01 lims.max_N])
    % Add invisible quiver plot to set the same aspect ratio
    axs{2} = axes('Visible','off','HandleVisibility','off'); 
    % Plot the quiver data
    quiver ( x , y , 0 .* squeeze(v_east(day_edited,:,:)), ...
        0 .* squeeze(v_north(day_edited,:,:)), scaleFactor, 'LineWidth', 2);
    hold off
    % Set properties for main axis
    axis equal
    xlim([0.5,params.num_longs+0.5])
    ylim([0.5,params.num_lats+0.5])
    yticks(y_ticks(5:5:end))
    yticklabels(flip(lat_centres(5:5:end)) + "\circ N")
    xticklabels(abs(long_centres) + "\circ W")

    % Phytoplankton
    axs{4} = nexttile(t1);
    plts{3} = imagesc(squeeze(plt.phyto(:,:,day)), [lims.min_P lims.max_P]);  
    hold on
    colormap(axs{4},col.map.phyto);
    c=colorbar;
    c.Location ='southoutside';
    c.Label.String = 'Phytoplankton concentration (g C m^{-3})';
    set(c, 'ylim', [lims.min_P+0.01 lims.max_P])
    % Add invisible quiver plot to set the same aspect ratio
    axs{2} = axes('Visible','off','HandleVisibility','off'); 
    % Plot the quiver data
    quiver ( x , y , 0 .* squeeze(v_east(day_edited,:,:)), ...
        0 .* squeeze(v_north(day_edited,:,:)), scaleFactor, 'LineWidth', 2);
    hold off
    % Set properties for main axis
    axis equal
    xlim([0.5,params.num_longs+0.5])
    ylim([0.5,params.num_lats+0.5])
    yticks(y_ticks(5:5:end))
    yticklabels(flip(lat_centres(5:5:end)) + "\circ N")
    xticklabels(abs(long_centres) + "\circ W")

    if params.detritus_layer 
        % Detritus
        axs{5} = nexttile(t1);
        plts{4} = imagesc(squeeze(plt.detri(:,:,day)), [lims.min_D lims.max_D]);  
        hold on
        colormap(axs{5},col.map.detri);
        cc=colorbar;
        cc.Location ='southoutside';
        cc.Label.String = 'Detritus concentration (g C m^{-3})';
        set(cc, 'ylim', [lims.min_D+0.01 lims.max_D])
        % Add invisible quiver plot to set the same aspect ratio
        axs{2} = axes('Visible','off','HandleVisibility','off'); 
        % Plot the quiver data
        quiver ( x , y , 0 .* squeeze(v_east(day_edited,:,:)), ...
            0 .* squeeze(v_north(day_edited,:,:)), scaleFactor, 'LineWidth', 2);
        hold off
        % Set properties for main axis
        axis equal
        xlim([0.5,params.num_longs+0.5])
        ylim([0.5,params.num_lats+0.5])
        yticks(y_ticks(5:5:end))
        yticklabels(flip(lat_centres(5:5:end)) + "\circ N")
        xticklabels(abs(long_centres) + "\circ W")
    end
elseif day ~= 1
    %set(axs{1}.Title, 'String', str_label);                     % Update date title map
    set(plts{1}, 'CData', squeeze(plt.zoop(:,:,day)));          % Update zooplankton map
    set(plts{2}, 'CData', squeeze(plt.nutri(:,:,day)));          % Update zooplankton map
    set(plts{3}, 'CData', squeeze(plt.phyto(:,:,day)));          % Update phytoplankton map
    if params.detritus_layer
        set(plts{4}, 'CData', squeeze(plt.detri(:,:,day)));          % Update detritus map
    end
end
sgtitle(str_label)
%sgtitle(sprintf('Day %d',day));
drawnow;