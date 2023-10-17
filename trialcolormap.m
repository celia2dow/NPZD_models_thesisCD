% Create sample data in the range -90 to +90
dblImage = 180*rand(50,50)-90;
% Set some rows to black
dblImage(20:25, :) = 0;
% Define a custom colormap - use jet.
myColorMap = jet(256);
% Find out what row a value of 0 would be in.
minValue = min(dblImage(:))
maxValue = max(dblImage(:))
zeroIndex = ceil((0-minValue) / (maxValue - minValue) * 256)
% Make that row of the colormap be zero.
myColorMap(zeroIndex, :) = [0,0,0]
% Display it
imshow(dblImage, [], 'InitialMagnification', 1000);
% Apply our custom colormap.
colormap(myColorMap);
colorbar;