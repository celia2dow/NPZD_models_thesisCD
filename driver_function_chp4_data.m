% Run set for producing numerous figures for chapter 4 - Canary Islands
% data
clear
close all

% SWITCH FOR VERTICAL STRATIFICATION
for h = 0:1
    params.temp_vert_mix = h;

    params.omega = 0.5; % For zooplankton grazing preference
    params.repeat_num = 4;
    smoothing = 7;          % Smoothing of temperature data: 0 if off, span of smoothing if on 
    thermocline_calc = 3;   % Which method is used to calculate the thermocline depth: (0) some input 
                            % function of MLD depths, MLD(), (1) the change in temperature gradient by 
                            % tol, (2) the drop in near surface (10m depth) temperature by tol, or (3)
                            % cubic spline the temperatures at all missing depths and then use method 2

    % FOLDER PATH
    folder_name = 'Canary_CHP4';
    folder_path_driver = [pwd '/' date '/' folder_name];
    % The folder paths used agrees with a macOS system and may need to be
    % flipped to '\' for a Windows system
    if ~exist(folder_path_driver, 'dir')
        mkdir(folder_path_driver)
    end
    
    % SWITCH FOR CHANGING MLD
    for i = 0:1     
        params.change_mix_depth = i;
        if h+i==0 
            continue   % Ignore case with no seasonal forcing
        end
        
        % SWITCH FOR DETRITUS LAYER
        for j = 0:1 
            params.detritus_layer = j;
            
            % SWITCH FOR HIGHER PREDATION RATE
            for d = [1, 1.5]  
                params.d = d; % m^3/(g C day)

                opts = [" model", "D model with grazing"];
                str_label1 = "NPZ" + opts(params.detritus_layer+1) + ...
                    newline + "d = " + num2str(params.d) + " m^3/(g C day)";
                off_on = {'off','on'};
                str_label2 = ['vertMix_' off_on{params.temp_vert_mix+1} ...
                    '__changeMLD_' off_on{params.change_mix_depth+1}];

                % GET NUMERICAL RESULTS
                NPZD_model_forcing

                x_vals1 = record.times/num_days; % In terms of years
                condit_BI = x_vals1 <= cut_off_day;
                condit_afterBI = record.times >= cut_off_day;
                days_afterBI=record.times(condit_afterBI);
                x_vals2 = days_afterBI-cut_off_day; % In terms of days
                
               
                if j + d == 1
                    % PLOT EXAMPLE WITH BURN IN
                    plot_burn_in

                    % PLOT COMBINATION OF FIGURES FOR GIVEN SEASONAL FORCING
                    fig3001=figure(3001);
                    x0=100;
                    y0=100;
                    width=1050;
                    height=1050;
                    set(gcf,'position',[x0,y0,width,height])
                    tiles = tiledlayout(5,1,'TileSpacing','Compact','Padding','Compact');

                    %%%%%%%%%%% TILE 1 %%%%%%%%%%%%%%%%%
                    axs1 = nexttile(tiles);
                    
                    plot_days = 0:0.1:num_days;
                    if i && h
                        colororder({'k','r'})
                        yyaxis left
                        plot(plot_days, ppval(params.MLD_CS{1,1}, plot_days),...
                            'Color','k','Parent',axs1)
                        set(gca, 'YDir','reverse')
                        ylabel("MLD (m)")

                        yyaxis right
                        plot(plot_days, ppval(params.TSI_CS{1,1}, plot_days),...
                            'Color','r','LineStyle', '-', 'Parent',axs1)
                        hold on
                        plot(plot_days, (params.TSI_mean-params.TSI_std_dev).*ones(1,length(plot_days)),...
                            'Color','r','LineStyle', '--')
                        text(260,params.TSI_mean-params.TSI_std_dev+30,'\mu_{TSI,10yr} - \sigma_{TSI,10yr}', 'Color', 'r')
                        hold off
                        ylabel("TSI (\circ C m)")
                    elseif i
                        plot(plot_days, ppval(params.MLD_CS{1,1}, plot_days),...
                            'Color','k','Parent',axs1)
                        set(gca, 'YDir','reverse')
                        ylabel("MLD (m)")
                    elseif h
                        plot(plot_days, ppval(params.TSI_CS{1,1}, plot_days),...
                            'Color','r','Parent',axs1)
                        hold on
                        plot(plot_days, (params.TSI_mean-params.TSI_std_dev).*ones(1,length(plot_days)),...
                            'Color','r','LineStyle', '--')
                        text(260,params.TSI_mean-params.TSI_std_dev+30,'\mu_{TSI,10yr} - \sigma_{TSI,10yr}', 'Color', 'r')
                        hold off
                        axs1.YColor = 'r';
                        ylabel("TSI (\circ C m)")
                    end
                    xlim([0,num_days])
                    num_seasons = num_years*4;
                    %xticks(tcks(1:num_seasons+1));
                    xticks(tcks(1:num_seasons+1)*num_days/12);
                    xticklabels(tck_seas(end-num_seasons-1:end));
                    %xtickangle(90);
                    set(gca,'XAxisLocation','top');
                else
                    figure(3001)
                end
                
                %%%%%%%%%%% NEXT TILE %%%%%%%%%%%%%%%%%
                nexttile(tiles);
                colororder([col.zoop_dark;col.nutri_orange])
                max_y = max(record.soluts(1,1,2:end,condit_afterBI),[],'all');

                yyaxis right
                plot(x_vals2,squeeze(record.soluts(1,1,1,condit_afterBI)), ...
                    'Color', col.nutri_orange); % Nutrient
                ylabel("N (g C/{m^3})")
                if h == 1 && i == 0
                    ylim([0,0.4])
                elseif h == 1 && i == 1
                    ylim([0,0.65])
                elseif h == 0 && i == 1
                    ylim([0,0.65 ...
                        ])
                end
                %ylim([0 0.65])

                yyaxis left
                if params.detritus_layer
                    plot(x_vals2,squeeze(record.soluts(1,1,4,condit_afterBI)), ...
                        'Color', col.detri_purp, 'LineStyle', '-'); % Detritus 
                    ylabel("P, Z, D (g C/{m^3})")
                else
                    ylabel("P, Z (g C/{m^3})")
                end
                hold on
                plot(x_vals2,squeeze(record.soluts(1,1,2,condit_afterBI)), ...
                    'Color', col.phyto_green, 'LineStyle', '-');  % Phytoplankton
                plot(x_vals2,squeeze(record.soluts(1,1,3,condit_afterBI)), ...
                    'Color', col.zoop_dark, 'LineStyle', '-');    % Zooplankton
                if h == 1 && i == 0
                    ylim([0,0.15])
                elseif h == 1 && i == 1
                    ylim([0,0.2])
                elseif h == 0 && i == 1
                    ylim([0,0.08])
                end
                hold off
                if j+d ~=2.5
                    set(gca,'XTick',[])
                end
                %ylim([0, 0.2])

                xlim([0,num_days])
                %text(20,max_y/2,str_label1)
                
                
            end
        end
        xlabel("Time (days)")
        leg = legend( "Detritus (D)", "Phytoplankton (P)", ...
            "Zooplankton (Z)", "Mixed layer nutrient (N)");
        leg.NumColumns = 4;
        leg.Location= 'southoutside';
    
        %%%%%%%%%%% SAVE FIGURE %%%%%%%%%%%%%%%%%
        savefig(fig3001, [folder_path_driver '/timeSeriesPlots__' str_label2], 'compact')
        saveas(fig3001, [folder_path_driver '/timeSeriesPlots__' str_label2], 'png')
    
        close all
    end
    clear
end