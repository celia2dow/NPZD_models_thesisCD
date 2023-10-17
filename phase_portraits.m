% plot phase portait
%pp = figure(55);

%axh = axes;

plot3(squeeze(record.soluts(:,:,2,:)), ...
    squeeze(record.soluts(:,:,1,:)), ...
    squeeze(record.soluts(:,:,3,:)), 'b')
hold on
dxdt = (record.soluts(:,:,2,3)-record.soluts(:,:,2,2))/(t(3)-t(2));
dydt = (record.soluts(:,:,1,3)-record.soluts(:,:,1,2))/(t(3)-t(2));
dzdt = (record.soluts(:,:,3,3)-record.soluts(:,:,3,2))/(t(3)-t(2));
q = quiver3(record.soluts(:,:,2,2), ...
    record.soluts(:,:,1,2), ...
    record.soluts(:,:,3,2), ...
    dxdt, dydt, dzdt, 'b');
q.MaxHeadSize = 2;
hold off


xlim([0,max(0.4, max(record.soluts(:,:,2,:),[],"all"))])
ylim([0,max(0.5, max(record.soluts(:,:,1,:),[],"all"))])
zlim([0,max(0.2, max(record.soluts(:,:,3,:),[],"all"))])
grid on
ylabel('Nutrient')
xlabel('Phytoplankton')
zlabel('Zooplankton')
set(gca, 'YDir','reverse')

% azimuth = -45;
% elevation = 15.264;
% % Isometric view, c. f. https://en.wikipedia.org/wiki/Isometric_projection
% view(axh,azimuth,elevation);
% camproj % returns 'orthographic'
% 
% unitx = [1;0;0];
% unity = [0;1;0];
% unitz = [0;0;1];
% projectedunitx = rotx(elevation) * rotz(-azimuth) * unitx;
% projectedunity = rotx(elevation) * rotz(-azimuth) * unity;
% xlabelangle = atan2d(projectedunitx(3),projectedunitx(1)) %#ok
% ylabelangle = -(180 - atan2d(projectedunity(3),projectedunity(1))) %#ok
% xlabelhandle = axh.XLabel;
% ylabelhandle = axh.YLabel;
% xlabelhandle.Rotation = xlabelangle;
% ylabelhandle.Rotation = ylabelangle;
% xlimits = xlim(axh);
% ylimits = ylim(axh);
% zlimits = zlim(axh);
% xmean = mean(xlimits);
% ymean = mean(ylimits);
% xbottom = xlimits(1);
% ybottom = ylimits(1);
% zbottom = zlimits(1);
% xlabelhandle.Position = [xmean ybottom zbottom];
% ylabelhandle.Position = [xbottom ymean zbottom];
% axis equal

%savefig(pp, [folder_path '/Phase_portrait'], 'compact')
%saveas(pp, [folder_path '/Phase_portrait'], 'png')

