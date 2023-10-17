close

rrr = [1 1 3; 2 2 2; 1 2 3];
axs = cell(1,2);
plts = cell(1,2);
gif1 = figure(11);
axs{1} = subplot(1,2,1);
hold on
plts{1}= imagesc(rrr,[0 3]);
colormap(axs{1},cool);
a=colorbar;
a.Label.String = 'rrr';
title('rrr');
hold off

axs{2} = subplot(1,2,2);
hold on
plts{2} = imagesc(rrr,[0,3]);
colormap(axs{2},winter);
b=colorbar;
b.Label.String = 'rrr og';
hold off
sgtitle(sprintf('Day 1'));
drawnow;
frame = getframe(gif1);
[A,map] = rgb2ind(frame2im(frame),256);
imwrite(A,map,'practice.gif','gif',"LoopCount",Inf,"DelayTime",1);

set(plts{1}, 'CData', [3 3 3; 3 3 3; 3 3 3]);
drawnow;
frame = getframe(gif1);
[A,map] = rgb2ind(frame2im(frame),256);
imwrite(A,map,'practice.gif','gif',"WriteMode","append","DelayTime",1);

set(plts{1},'CData', [0 0 0; 0 0 0; 0 0 0]);
drawnow;
frame = getframe(gif1);
[A,map] = rgb2ind(frame2im(frame),256);
imwrite(A,map,'practice.gif','gif',"WriteMode","append","DelayTime",1);


