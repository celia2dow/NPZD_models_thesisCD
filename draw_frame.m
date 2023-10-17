function [axs,plts] = draw_frame(params,plt,temps,v_north,v_east,...
    lims,day,axs,plts,vars_to_visualise,long_centres,lat_centres,col)
% DRAW_FRAME Save each figure as a frame in the gif illustrating the
% density of phytoplankton in a region of the ocean

num_uniq_plots = length(plts); % Number of unique plots including temperature and flow field

% Month (season) label
month_num = params.months_per_day(day);
month_label = params.month_season(1,month_num);
season_label = params.month_season(2,month_num);
year = params.years_per_day(day);
if year<params.year1
    str_label = "BI year " + num2str(year) + ", " + month_label + " (" + season_label + ")" + "     Day " + num2str(day);
else
    str_label = num2str(year) + ", " + month_label + " (" + season_label + ")"+ "     Day " + num2str(day);
end
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
    if num_uniq_plots <= 4
        t2 = tiledlayout(1,num_uniq_plots-1,'TileSpacing','Compact','Padding','Compact');
    elseif num_uniq_plots >= 5
        t1 = tiledlayout(2,1,'TileSpacing','compact');
        t2 = tiledlayout(t1,'flow',"TileSpacing",'Compact');
        t3 = tiledlayout(t1,'flow',"TileSpacing",'Compact');
        t3.Layout.Tile = 2;
    end

    % Temperature
    axs{1} = nexttile(t2);
    if params.num_lats == 1
        plts{1} = imagesc(squeeze(temps(day_edited)),[lims.min_T lims.max_T]);
    else
        plts{1} = imagesc(squeeze(temps(day_edited,:,:)),[lims.min_T lims.max_T]);
    end
    hold on
    colormap(axs{1},col.map.temp);

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
    a=colorbar;
    a.Label.String = 'Temperature (\circ C)';
    title(str_label);
    xlim([0.5,params.num_longs+0.5])
    ylim([0.5,params.num_lats+0.5])
    %yticks(1:params.num_lats)
    yticklabels(flip(lat_centres))
    %xticks(1:params.num_longs)
    xticklabels(long_centres)

    % Zooplankton
    axs{3} = nexttile(t2);
    plts{3} = imagesc(squeeze(plt.zoop(:,:,day)), [lims.min_Z lims.max_Z]);  
    hold on
    colormap(axs{3},col.map.zoop);
    % Add invisible quiver plot to set the same aspect ratio
    axs{2} = axes('Visible','off','HandleVisibility','off'); 
    % Plot the quiver data
    quiver ( x , y , 0 .* squeeze(v_east(day,:,:)), ...
        0 .* squeeze(v_north(day,:,:)), scaleFactor, 'LineWidth', 2);
    hold off
    % Set properties for main axis
    axis equal
    b=colorbar;
    b.Label.String = 'Zooplankton concentration (g C m^{-3})';
    xlim([0.5,params.num_longs+0.5])
    ylim([0.5,params.num_lats+0.5])
    yticks(1:params.num_lats)
    yticklabels(flip(lat_centres))
    xticklabels(long_centres)
    
    checkN = 1;
    checkP = 1;
    for ax_num = 4:num_uniq_plots
        if num_uniq_plots == 4
            axs{ax_num}=nexttile(t2);   % Plot in-line with temperature, current and zooplankton
        else
            axs{ax_num}=nexttile(t3);   % Plot on next line
        end

        if checkN && any(cellfun(@isequal, vars_to_visualise, ...
                repmat({'N'}, size(vars_to_visualise)))) 
            % Nutrient
            plts{ax_num} = imagesc(squeeze(plt.nutri(:,:,day)), [lims.min_N lims.max_N]);  
            hold on
            colormap(axs{ax_num},col.map.nutri);
            b=colorbar;
            b.Label.String = 'Nutrient concentration (g C m^{-3})';
            checkN = 0;                 % N has now been plotted
        elseif checkP && any(cellfun(@isequal, vars_to_visualise, ...
                repmat({'P'}, size(vars_to_visualise))))
            % Phytoplankton
            plts{ax_num} = imagesc(squeeze(plt.phyto(:,:,day)), [lims.min_P lims.max_P]);  
            hold on
            colormap(axs{ax_num},col.map.phyto);
            b=colorbar;
            b.Label.String = 'Phytoplankton concentration (g C m^{-3})';
            checkP = 0;                 % P has now been plotted
        elseif any(cellfun(@isequal, vars_to_visualise, ...
                repmat({'D'}, size(vars_to_visualise)))) 
            % Detritus
            plts{ax_num} = imagesc(squeeze(plt.detri(:,:,day)), [lims.min_D lims.max_D]);  
            hold on
            colormap(axs{ax_num},col.map.detri);
            b=colorbar;
            b.Label.String = 'Detritus concentration (g C m^{-3})';
        end
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
        %yticks(1:params.num_lats)
        yticklabels(flip(lat_centres))
        xticklabels(long_centres)
    end
elseif day ~= 1
    if params.num_lats == 1
        set(plts{1}, 'CData', squeeze(temps(day_edited)));             % Update temperature map
    elseif day_edited == 1
        set(plts{1}, 'CData', squeeze(temps(day_edited,:,:)));         % Update temperature map
        if squeeze(temps(day_edited,:,:)) ~= squeeze(temps(params.num_days,:,:))
            set(plts{2}, 'UData', squeeze(v_east(day_edited,:,:)), ...
                'VData', -squeeze(v_north(day_edited,:,:)))            % Update flow field map
        end
    else
        set(plts{1}, 'CData', squeeze(temps(day_edited,:,:)));         % Update temperature map
        if squeeze(temps(day_edited,:,:)) ~= squeeze(temps(day_edited-1,:,:))
            set(plts{2}, 'UData', squeeze(v_east(day_edited,:,:)), ...
                'VData', -squeeze(v_north(day_edited,:,:)))            % Update flow field map
        end
    end
    set(axs{1}.Title, 'String', str_label);                     % Update date title map
    set(plts{3}, 'CData', squeeze(plt.zoop(:,:,day)));          % Update zooplankton map
    checkN = 1;
    checkP = 1;
    for ax_num = 4:num_uniq_plots
        if checkN && any(cellfun(@isequal, vars_to_visualise, ...
                repmat({'N'}, size(vars_to_visualise)))) 
            set(plts{ax_num}, 'CData', ...
                squeeze(plt.nutri(:,:,day)));                   % Update nutrient map
            checkN = 0;                                         % N has now been updated
        elseif checkP && any(cellfun(@isequal, vars_to_visualise, ...
                repmat({'P'}, size(vars_to_visualise))))
            set(plts{ax_num}, 'CData', ...
                squeeze(plt.phyto(:,:,day)));                   % Update phytoplankton map
            checkP = 0;                                         % P has now been updated
        elseif any(cellfun(@isequal, vars_to_visualise, ...
                repmat({'D'}, size(vars_to_visualise))))
            set(plts{ax_num}, 'CData', ...
                squeeze(plt.detri(:,:,day)));                   % Update detritus map
        end
    end
end

%sgtitle(sprintf('Day %d',day));
drawnow;