% Create gifs for temperature, current, and TSI scaling
num_plots = 3;
plts = cell(1,num_plots+1);
axs = cell(1,num_plots+1);

if params.detritus_layer == 0
    % Prepare environmental gif
    ocean_grid_1 = figure(1);
    ocean_grid_1.WindowState = 'maximized';
    %axis tight manual % ensures getframe() returns a consistent size
    %ocean_grid.Visible = 'off';
    gifpath = [folder_path_driver '/simulation_enviro.gif'];
    for day = 1:num_days
        
        % Save each figure as a frame in the gif
        [axs,plts] = draw_frame_enviro(params,MLD_rearr,scales_rearr,temp_rearr,v_north_rearr,...
            v_east_rearr,lims,day,axs,plts,vars_to_visualise,...
            long_centres,lat_centres,col);
        % Save the gif frame.
        frame = getframe(ocean_grid_1);
        [A,map] = rgb2ind(frame2im(frame),256);
        if day == 1
            imwrite(A,map,gifpath,'gif',"LoopCount",Inf,"DelayTime",1/params.vid_speed);
        else
            imwrite(A,map,gifpath,'gif',"WriteMode","append","DelayTime",1/params.vid_speed);
        end
    end
end