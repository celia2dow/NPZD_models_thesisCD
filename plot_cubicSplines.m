% Code plotting cubic splines against data. Called by NPZD_model_forcing.m

for site = 1:num_sites
    if params.num_lats==1
        figure(124)
        i = 1;
        j = 1;
        
        plot_days = 0:0.1:num_days;
        %subplot(1,2,1)
        %plot(plot_days, ppval(params.temp_CS{i,j}, plot_days),'k-', ...
        %    GOODindices_days_T-1, temps_av(GOODindices_days_T), 'r*')
        %legend("Cubic spline", "Data")
        %title("Temperature of MLD")
        %xlabel("Time (days)")
        %ylabel("Temperature (^{\circ} C)")

        %subplot(1,2,2)
        name_site = sprintf('i%i_j%i',i,j);
        plot(plot_days, ppval(params.MLD_CS{i,i}, plot_days),'k-',...
            GOODindices_days_T-1, depths_LEVs(GOODindices_days_T), 'r*')
        title("Mixed layer depth")
        xlabel("Time (days)")
        ylabel("MLD (m)")

        savefig(figure(124), [folder_path '/Cubic_splines_' name_site], 'compact')
        saveas(figure(124), [folder_path '/Cubic_splines_' name_site], 'png')
        break
    end
    figure(123 + site) % The plot for site (i,j)
    i = sites_to_visualise(site,1);
    j = sites_to_visualise(site,2); 
    
    plot_days = 0:0.1:num_days;
    subplot(2,2,1)
    plot(plot_days, ppval(params.temp_CS{i,j}, plot_days),'k-', ...
        GOODindices_days_T-1, temps_av(GOODindices_days_T,i,j), 'r*')
    legend("Cubic spline", "Data")
    title("Temperature")
    xlabel("Time (days)")
    ylabel("Temperature (^{\circ} C)")
    
    subplot(2,2,2)
    plot(plot_days, ppval(params.MLD_CS{1,1}, plot_days),'k-',...
        GOODindices_days_T-1, depths_LEVs(GOODindices_days_T,1,1), 'r*')
    title("Mixed layer depth")
    xlabel("Time (days)")
    ylabel("MLD (m)")
    
    subplot(2,2,3)
    plot(plot_days, ppval(params.v_north_CS{i,j}, plot_days),'k-',...
        GOODindices_days_CN-1, v_north_av(GOODindices_days_CN,i,j), 'r*')
    title("Northerly current")
    xlabel("Time (days)")
    ylabel("Current (m/day)")
    
    subplot(2,2,4)
    plot(plot_days, ppval(params.v_east_CS{i,j}, plot_days),'k-',...
        GOODindices_days_CE-1, v_east_av(GOODindices_days_CE,i,j), 'r*')
    title("Easterly current")
    xlabel("Time (days)")
    ylabel("Current (m/day)")
    
    name_site = sprintf('i%i_j%i',i,j);
    sgtitle(sprintf("Data for %4.2f^{\\circ}N, %4.2f^{\\circ}E",...
        lat_centres(i),long_centres(j)))
    savefig(figure(123 + site), [folder_path '/Cubic_splines_' name_site], 'compact')
    saveas(figure(123 + site), [folder_path '/Cubic_splines_' name_site], 'png')
end

if thermocline_calc == 0
    figure(4000)
    plot(plot_days(1:3650), ppval(params.MLD_CS{i,i}, plot_days(1:3650)),'k-')
    title("Mixed layer depth")
    xlabel("Time (days)")
    ylabel("MLD (m)")

    savefig(figure(4000), [folder_path '/MLD'], 'compact')
    saveas(figure(4000), [folder_path '/MLD'], 'png')
end