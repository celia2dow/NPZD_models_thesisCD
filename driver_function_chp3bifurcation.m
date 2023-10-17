% Sampling of different $d$ values for 1D model without forcing

% Run set for producing numerous figures for chapter 3
close all
clear

params.temp_vert_mix = 0; % SWITCH FOR VERTICAL STRATIFICATION
params.change_mix_depth = 0; % SWITCH FOR CHANGING MLD
params.current = 0;                 % Inclusion of current in model
params.diffus = 0;                  % Inclusion of diffusion in model
params.BC = 0;                      % Boundary conditions of spatial model: 0 periodic, 1 open
params.repeat_num = 6;
smoothing = 0;                      % Smoothing of temperature data: 0 if off, span of smoothing if on (unused)
thermocline_calc = 0;               % Which method is used to calculate the thermocline depth: (0) some input 
                                    % function of MLD depths, MLD(), (1) the change in temperature gradient by 
                                    % tol, (2) the drop in near surface (10m depth) temperature by tol, or (3)
                                    % cubic spline the temperatures at all missing depths and then use method 2
labels = {'N^* (g C/{m^3})', 'P^*(g C/{m^3})', 'Z^* (g C/{m^3})', ...
    'D^* (g C/{m^3})'};
    

% FOLDER PATH
folder_name = 'E&B_CHP3';
folder_path_driver = [pwd '/' date '/' folder_name];
% The folder paths used agrees with a macOS system and may need to be
% flipped to '\' for a Windows system
if ~exist(folder_path_driver, 'dir')
    mkdir(folder_path_driver)
end
    
% SWITCH FOR DETRITUS LAYER
for j = 0 % 0:1:2
    close all
    params.detritus_layer = min(j,1);
    if j == 2
        params.omega = 0.5; % Grazing on detritus
    else
        params.omega = 0; % No grazing on detritus
    end
  
    % RANGE OF d VALUES
    d_vals = 0.5:0.01:3.5;
    nd = length(d_vals);
    Xstar.steady = zeros(2,nd);
    if params.detritus_layer
        Xstar.lims = zeros(4,3,nd);
        num_plots  =4;
    else
        Xstar.lims = zeros(3,3,nd);
        num_plots =3;
    end

    opts = {' model', 'D model', 'D model with grazing'};
    str_label1 = ['NPZ' opts{j+1}];
    % SWITCH FOR HIGHER PREDATION RATE
    textprogressbar(['Model: ' num2str(j+1) '     ']);
    for i = 1:nd
        d = d_vals(i);
        textprogressbar(i*100/nd)
        params.d = d; % m^3/(g C day)

        % GET NUMERICAL RESULTS
        NPZD_model_forcing
        cols = {col.nutri_orange, col.phyto_green, col.zoop_dark, col.detri_purp};

        % IS IT A STABLE STEADY STATE OR A LIMIT CYCLE?
        test_steadiness = record.soluts(1,1,1,end-100:end)-record.soluts(1,1,1,end);
        test_steadiness = sum(test_steadiness <= 0.0001);
        if test_steadiness == 101
            Xstar.steady(1,i) = 1;
            Xstar.lims(:,1,i) = squeeze(record.soluts(1,1,:,end));
            Xstar.lims(:,2,i) = Xstar.lims(:,1,i);
            Xstar.lims(:,3,i) = Xstar.lims(:,1,i);
        else
            % Find the period of the stable limit cycle
            [peaks,peak_locations] = findpeaks(squeeze(record.soluts(1,1,1,:)));
            periods = diff(record.times(peak_locations));
            include_peaks = abs(periods-periods(end))<0.2;
            check_for_damp = abs(peaks-peaks(end))<0.001;
            include_peak_locations = peak_locations(logical([0,include_peaks']));
            peaks_begin = include_peak_locations(1);
            if sum(check_for_damp)/length(check_for_damp)>=0.9
                % If there is no damping
                Xstar.steady(2,i) = mean(periods(include_peaks));
                % Find the min and max of the cycle
                min_peaks = min(squeeze(record.soluts(...
                    1,1,:,peaks_begin:end)),[],2);
                max_peaks = max(squeeze(record.soluts(...
                    1,1,:,peaks_begin:end)),[],2);
                mean_peaks = mean(squeeze(record.soluts(...
                    1,1,:,peaks_begin:end)),2);
                Xstar.lims(:,1,i) = min_peaks;
                Xstar.lims(:,2,i) = mean_peaks;
                Xstar.lims(:,3,i) = max_peaks;
            else
                Xstar.steady(1,i) = 1;
                Xstar.lims(:,1,i) = mean(squeeze(record.soluts(...
                    1,1,:,peaks_begin:end)),2);
                Xstar.lims(:,2,i) = Xstar.lims(:,1,i);
                Xstar.lims(:,3,i) = Xstar.lims(:,1,i);
            end
        end
    end
    textprogressbar('done');
 %%
    % PLOT COMBINATION OF FIGURES FOR GIVEN MODEL
    fig299=figure(299);
    x0=100;
    y0=100;
    width=1050;
    height=400;
    set(gcf,'position',[x0,y0,width,height])
    tiles = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
    
    steady_condit = logical(Xstar.steady(1,:));
    
    for k = 1:num_plotsx
        nexttile(k);
        not_steady_index = find(~steady_condit);
        ns_1 = not_steady_index(1)+1;
        ns_end = not_steady_index(end)+1;
        % Min
        plot(d_vals,squeeze(Xstar.lims(k,1,:)), 'Color', cols{k}, 'Marker', '.', 'LineWidth', 4);
        hold on
        plot(d_vals,squeeze(Xstar.lims(k,2,:)), 'Color', cols{k}, 'Marker', '.', 'LineWidth', 2);
        plot(d_vals,squeeze(Xstar.lims(k,3,:)), 'Color', cols{k}, 'Marker', '.', 'LineWidth', 4);
        hold off
%         plot([d_vals(1:ns_1),d_vals(ns_end:end)],...
%             [squeeze(Xstar.lims(k,1,1:ns_1)),squeeze(Xstar.lims(k,1,ns_end:end))],...
%             'Color', cols{k}, 'LineStyle', '-', 'LineWidth', 2);
%         if sum(~steady_condit)
%             hold on
%             plot(d_vals(ns_1:ns_end), squeeze(Xstar.lims(k,2,ns_1:ns_end)),...
%                 'Color',cols{k}, 'LineStyle', ':', 'LineWidth', 2)
%             plot(d_vals(ns_1:ns_end), squeeze(Xstar.lims(k,1,ns_1:ns_end)),...
%                 'Color',cols{k}, 'LineStyle', '-', 'LineWidth', 2)
%             plot(d_vals(ns_1:ns_end), squeeze(Xstar.lims(k,3,ns_1:ns_end)),...
%                 'Color',cols{k}, 'LineStyle', '-', 'LineWidth', 2)
%             hold off
%         end
        xlim([d_vals(1),d_vals(end)])
        ylabel(labels{k})
        xlabel('d (m^3/{g C day}')
    end
   
    %%%%%%%%%%% SAVE FIGURE %%%%%%%%%%%%%%%%%
    savefig(fig299, [folder_path_driver '/bifurcation_plots_' str_label1], 'compact')
    saveas(fig299, [folder_path_driver '/bifurcation_plots_' str_label1], 'png')

    % Plot periods
    fig300 = figure(300);
    plot(d_vals, Xstar.steady(2,:),'k.', 'LineWidth', 2)
    %%%%%%%%%%% SAVE FIGURE %%%%%%%%%%%%%%%%%
    savefig(fig300, [folder_path_driver '/Periods_' str_label1], 'compact')
    saveas(fig300, [folder_path_driver '/Periods_' str_label1], 'png')
end
