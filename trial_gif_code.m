%%
% Go through each day in time and save the spatial distributuion of
% temperature, and phytoplankton
% Interpolate produced data to the days for which we have temperature data
if visual_switch
    % Prepare gif
    ocean_grid = figure(1);
    axis tight manual % ensures getframe() returns a consistent size
    %ocean_grid.Visible = 'off';
    gifpath = [folder_path '/simulation.gif'];

    % Interpolate phytoplankton data
    phyto = zeros(params.num_lats,params.num_longs,num_days);
    for i = 1:params.num_lats
        for j = 1:params.num_longs
            non0times = record.times(record.times(:)>0);
            time_lim = length(squeeze(non0times));
            phyto(i,j,:) = interp1(squeeze(non0times),...
                squeeze(record.soluts(i,j,2,1:time_lim)),params.days);
        end
    end
    % Find limits
    lims.min_T = min(temps_av(temps_av>0),[],'all');
    lims.max_T = max(temps_av,[],'all');
    lims.min_P = min(phyto,[],'all');
    lims.max_P = max(phyto,[],'all');
    
    % Preallocate
    plts = cell(1,2);
    axs = cell(1,2);
    for day = params.days
        day = day+1;
        % Save each figure as a frame in the gif
        draw_frame2(params,phyto,temps_av,lims,day,axs,plts);
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


% 
% % Go through each day in time and save the spatial distributuion of
% % temperature, and phytoplankton
% % Interpolate produced data to the days for which we have temperature data
% if visual_switch
%     % Prepare gif
%     ocean_grid = figure(1);
%     axis tight manual % ensures getframe() returns a consistent size
%     %ocean_grid.Visible = 'off';
%     gifpath = [folder_path '/simulation.gif'];
% 
%     % Interpolate phytoplankton data
%     phyto = zeros(params.num_lats,params.num_longs,num_days);
%     for i = 1:params.num_lats
%         for j = 1:params.num_longs
%             non0times = record.times(record.times(:)>0);
%             time_lim = length(squeeze(non0times));
%             phyto(i,j,:) = interp1(squeeze(non0times),...
%                 squeeze(record.soluts(i,j,2,1:time_lim)),params.days);
%         end
%     end
%     % Find limits
%     lims.min_T = min(temps_av(temps_av>0),[],'all');
%     lims.max_T = max(temps_av,[],'all');
%     lims.min_P = min(phyto,[],'all');
%     lims.max_P = max(phyto,[],'all');
%     
%     % Preallocate
%     plts = cell(1,2);
%     axs = cell(1,2);
%     for day = params.days
%         day = day+1;
%         % Save each figure as a frame in the gif
%         draw_frame2(params,phyto,temps_av,lims,day,axs,plts);
%         % Save the gif frame.
%         frame = getframe(ocean_grid);
%         [A,map] = rgb2ind(frame2im(frame),256);
%         if day == 1
%             imwrite(A,map,gifpath,'gif',"LoopCount",Inf,"DelayTime",1/params.vid_speed);
%         else
%             imwrite(A,map,gifpath,'gif',"WriteMode","append","DelayTime",1/[arams.vid_speed]);
%         end
%     end
% end