% Run set for producing numerous figures for chapter 3
close all
clear

params.temp_vert_mix = 0; % SWITCH FOR VERTICAL STRATIFICATION
params.change_mix_depth = 0; % SWITCH FOR CHANGING MLD
params.current = 0;                 % Inclusion of current in model
params.diffus = 0;                  % Inclusion of diffusion in model
params.BC = 0;                      % Boundary conditions of spatial model: 0 periodic, 1 open
params.repeat_num = 1;
smoothing = 0;                      % Smoothing of temperature data: 0 if off, span of smoothing if on (unused)
thermocline_calc = 0;               % Which method is used to calculate the thermocline depth: (0) some input 
                                    % function of MLD depths, MLD(), (1) the change in temperature gradient by 
                                    % tol, (2) the drop in near surface (10m depth) temperature by tol, or (3)
                                    % cubic spline the temperatures at all missing depths and then use method 2

% FOLDER PATH
folder_name = 'E&B_CHP3';
folder_path_driver = [pwd '/' date '/' folder_name];
% The folder paths used agrees with a macOS system and may need to be
% flipped to '\' for a Windows system
if ~exist(folder_path_driver, 'dir')
    mkdir(folder_path_driver)
end
    
% SWITCH FOR DETRITUS LAYER
for j = 0:1:2
    params.detritus_layer = min(j,1);
    if j == 2
        params.omega = 0.5; % Grazing on detritus
    else
        params.omega = 0; % No grazing on detritus
    end

    % PLOT COMBINATION OF FIGURES FOR GIVEN MODEL
    fig299=figure(299);
    x0=100;
    y0=100;
    width=1050;
    height=400;
    set(gcf,'position',[x0,y0,width,height])
    tiles = tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');
    
    % SWITCH FOR HIGHER PREDATION RATE
    for d = [2.2, 2.3]  
        params.d = d; % m^3/(g C day)

        opts = {' model', 'D model', 'D model with grazing'};
        str_label1 = ['NPZ' opts{j+1}];
        str_label2 = "d = " + num2str(params.d) + " m^3/(g C day)";

        % GET NUMERICAL RESULTS
        NPZD_model_forcing
        x_vals = record.times; % In terms of days

        %%%%%%%%%%% TILE 1 %%%%%%%%%%%%%%%%%
        ax1 = nexttile([1,2]);

        fig299;
        max_y = max(record.soluts(1,1,:,:),[],'all');
        plot(x_vals,squeeze(record.soluts(1,1,1,:)), ...
            'Color', col.nutri_orange); % Nutrient
        hold on
        plot(x_vals,squeeze(record.soluts(1,1,2,:)), ...
            'Color', col.phyto_green);  % Phytoplankton
        plot(x_vals,squeeze(record.soluts(1,1,3,:)), ...
            'Color', col.zoop_dark);    % Zooplankton
        if params.detritus_layer
            plot(x_vals,squeeze(record.soluts(1,1,4,:)), ...
                'Color', col.detri_purp); % Detritus 
        end
        hold off
        
        ylabel("Concentration (g C/{m^3})")
        xlim([0,num_days])
        %text(150,2*max_y/3,str_label2) 

        %%%%%%%%%%% TILE 2 %%%%%%%%%%%%%%%%%
        ax2 = nexttile([1,1]);
        phase_portraits
        %num_seasons = num_years*4;
        %xticks(tcks(1:num_seasons+1));
        %xticklabels(tck_seas(end-num_seasons-1:end));
        %xtickangle(90);
    end
    ax1.XLabel.String = 'Time (days)';
    if params.detritus_layer
        leg = legend(ax1,{'Mixed layer nutrient (N)', 'Phytoplankton (P)', ...
            'Zooplankton (Z)', 'Detritus (D)'});
        leg.NumColumns = 4;
    else
        leg = legend(ax1,{'Mixed layer nutrient (N)', 'Phytoplankton (P)', ...
            'Zooplankton (Z)'});
        leg.NumColumns = 3;
    end
    leg.Location = 'southoutside';
    
    %%%%%%%%%%% SAVE FIGURE %%%%%%%%%%%%%%%%%
    savefig(fig299, [folder_path_driver '/timeSeriesPlots_and_PhasePortraits__' str_label1], 'compact')
    saveas(fig299, [folder_path_driver '/timeSeriesPlots_and_PhasePortraits__' str_label1], 'png')
end
