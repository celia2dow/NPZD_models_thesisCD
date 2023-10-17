% Visualise temperature gradient and MLD
figure(56)
subplot(1,2,1)
flipT = temps(GOODindices_days_T,:,4,34)';
imagesc(flipT,[lims.min_T lims.max_T]);
colorbar
yticks(1:length(data_file_temps.LEV))
yticklabels(data_file_temps.LEV)
xlabel("Time (days)")

subplot(1,2,2)
plot(GOODindices_days_T-1, depths_LEVs(GOODindices_days_T,4,34), 'r*');
ax = gca;
ax.YDir = 'reverse';
title("Mixed layer depth")
xlabel("Time (days)")
ylabel("MLD (m)")