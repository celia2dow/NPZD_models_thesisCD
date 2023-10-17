function [axs,plts] = draw_frame_enviro(params,MLD,scales_tsi,temps,v_north,v_east,...
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
    t1 = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
    
    % Temperature
    axs{1} = nexttile(t1);
    if params.temp_vert_mix
        if params.num_lats == 1
            plts{1} = imagesc(squeeze(temps(day_edited)),[lims.min_T lims.max_T]);
        else
            plts{1} = imagesc(squeeze(temps(day_edited,:,:)),[lims.min_T lims.max_T]);
        end
        colormap(axs{1},col.map.temp);
    end
    hold on

    % Flow Field
    % Add a 2nd invisible axis to host the 2nd colormap and colorbar
    % HandleVisibility is set to off to avoid accidentally accessing those axes.
    axs{2} = axes('Visible','off','HandleVisibility','off'); 
    % Plot the quiver data
    plts{2} = quiver ( x , y , squeeze(v_east(day_edited,:,:)), ...
        -squeeze(v_north(day_edited,:,:)), scaleFactor, 'LineWidth', 2);
    hold off
    % Set properties for main axis
    axis equal
    if params.temp_vert_mix
        a=colorbar;
        a.Location ='southoutside';
        a.Label.String = 'Temperature (\circ C)';
        set(a, 'ylim', [lims.min_T+0.1 lims.max_T])
    end
    %title(str_label);
    xlim([0.5,params.num_longs+0.5])
    ylim([0.5,params.num_lats+0.5])
    %yticks(1:params.num_lats)
    y_ticks = 1:params.num_lats;
    yticks(y_ticks(5:5:end))
    yticklabels(flip(lat_centres(5:5:end)) + "\circ N")
    %xticks(1:params.num_longs)
    xticklabels(abs(long_centres) + "\circ W")

    % MLD
    if params.change_mix_depth
        axs{3} = nexttile(t1);
        plts{3} = imagesc(squeeze(MLD(day_edited,:,:)),[lims.min_MLD lims.max_MLD]);
        hold on
        colormap(axs{3},col.map.MLD);
        % Add invisible quiver plot to set the same aspect ratio
        axs{2} = axes('Visible','off','HandleVisibility','off'); 
        % Plot the quiver data
        quiver ( x , y , 0 .* squeeze(v_east(day,:,:)), ...
            0 .* squeeze(v_north(day,:,:)), scaleFactor, 'LineWidth', 2);
        hold off
        % Set properties for main axis
        axis equal
        b=colorbar;
        b.Location ='southoutside';
        b.Label.String = 'MLD (m)';
        set(b, 'ylim', [lims.min_MLD + 10 lims.max_MLD])
        xlim([0.5,params.num_longs+0.5])
        ylim([0.5,params.num_lats+0.5])
        yticks(y_ticks(5:5:end))
        yticklabels(flip(lat_centres(5:5:end)) + "\circ N")
        xticklabels(abs(long_centres) + "\circ W")
    end

    % SCALING
    if params.temp_vert_mix
        axs{4} = nexttile(t1);
        plts{4} = imagesc(squeeze(scales_tsi(day_edited,:,:)),[-0.5 lims.max_scale]);
        hold on
        colormap(axs{4},col.map.scale);
        % Add invisible quiver plot to set the same aspect ratio
        axs{2} = axes('Visible','off','HandleVisibility','off'); 
        % Plot the quiver data
        quiver ( x , y , 0 .* squeeze(v_east(day,:,:)), ...
            0 .* squeeze(v_north(day,:,:)), scaleFactor, 'LineWidth', 2);
        hold off
        % Set properties for main axis
        axis equal
        c=colorbar;
        c.Location ='southoutside';
        c.Label.String = 'Scaling';
        c.Ticks = [0.1,0.9] ; 
        c.TickLabels = {'Strong stratification','No stratification'} ;   
        set(c, 'ylim', [lims.min_scale+0.01 lims.max_scale])
        xlim([0.5,params.num_longs+0.5])
        ylim([0.5,params.num_lats+0.5])
        yticks(y_ticks(5:5:end))
        yticklabels(flip(lat_centres(5:5:end)) + "\circ N")
        xticklabels(abs(long_centres) + "\circ W")
    end

elseif day ~= 1
    if params.num_lats == 1
        if params.temp_vert_mix
            set(plts{1}, 'CData', squeeze(temps(day_edited)));             % Update temperature map
        end
    elseif day_edited == 1
        if params.temp_vert_mix
            set(plts{1}, 'CData', squeeze(temps(day_edited,:,:)));         % Update temperature map
        end
        if squeeze(temps(day_edited,:,:)) ~= squeeze(temps(params.num_days,:,:))
            set(plts{2}, 'UData', squeeze(v_east(day_edited,:,:)), ...
                'VData', -squeeze(v_north(day_edited,:,:)))            % Update flow field map
        end
    else
        if params.temp_vert_mix
            set(plts{1}, 'CData', squeeze(temps(day_edited,:,:)));         % Update temperature map
        end
        set(plts{2}, 'UData', squeeze(v_east(day_edited,:,:)), ...
            'VData', -squeeze(v_north(day_edited,:,:)))            % Update flow field map
    end
    %set(axs{1}.Title, 'String', str_label);                     % Update date title map
    if params.change_mix_depth
        set(plts{3}, 'CData', squeeze(MLD(day,:,:)));               % Update MLD map
    end
    if params.temp_vert_mix
        set(plts{4}, 'CData', squeeze(scales_tsi(day,:,:)));        % Update TSI scaling map
    end
end

sgtitle(str_label)
%sgtitle(sprintf('Day %d',day));
drawnow;