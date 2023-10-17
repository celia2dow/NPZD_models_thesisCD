% Create gif for NPZD
num_plots = 4;
plts = cell(1,num_plots+2);
axs = cell(1,num_plots+2);

% Prepare environmental gif
ocean_grid_2 = figure(1);
ocean_grid_2.WindowState = 'maximized';
%axis tight manual % ensures getframe() returns a consistent size
%ocean_grid.Visible = 'off';
gifpath = [folder_path_driver '/simulation_final_' lil_str '.gif'];
for day = 1:num_days
    
    % Save each figure as a frame in the gif
    [axs,plts] = draw_frame_NPZD(params,plt,temp_rearr,v_north_rearr,...
        v_east_rearr,lims,day,axs,plts,vars_to_visualise,...
        long_centres,lat_centres,col);
    % Save the gif frame.
    frame = getframe(ocean_grid_2);
    [A,map] = rgb2ind(frame2im(frame),256);
    if day == 1
        imwrite(A,map,gifpath,'gif',"LoopCount",Inf,"DelayTime",1/params.vid_speed);
    else
        imwrite(A,map,gifpath,'gif',"WriteMode","append","DelayTime",1/params.vid_speed);
    end
end