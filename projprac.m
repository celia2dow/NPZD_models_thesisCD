 % Overlay landmarks on 'boston.tif'.
    % Includes material (c) GeoEye, all rights reserved.
 
    % Specify latitude and longitude coordinates of landmarks in Boston.
    lat = [42.3604 42.3691 42.3469 42.3480 42.3612]; 
    lon = [-71.0580 -71.0710 -71.0623 -71.0968 -71.0941];
 
    % Obtain the projcrs object corresponding to 'boston.tif'.
    info = georasterinfo('boston.tif');
    proj = info.CoordinateReferenceSystem;
 
    % Project the landmarks.
    [x, y] = projfwd(proj, lat, lon);
 
    % Read the 'boston.tif' image.
    [RGB, R] = readgeoraster('boston.tif');
 
    % Display the image and projected coordinates.
    figure
    mapshow(RGB, R)
    mapshow(x,y,'DisplayType','point','Marker','o', ...
        'MarkerFaceColor','y','MarkerEdgeColor','none')
    xlabel('Easting in Survey Feet')
    ylabel('Northing in Survey Feet')