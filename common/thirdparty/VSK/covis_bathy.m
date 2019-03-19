function [xb,yb,zb] = covis_bathy(covis, grd)
%
% Load bathy data set and regrid it onto the grid defined with
% the input grid (grd)
%
% ----------
% This program is free software distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY. You can redistribute it and/or modify it.
% Any modifications of the original software must be distributed in such a 
% manner as to avoid any confusion with the original work.
% 
% Please acknowledge the use of this software in any publications arising  
% from research that uses it.
% 
% ----------
%  Version 1.0 - 10/2010, cjones@apl.washington.edu
%

% translate and rotate bathy to fit covis data
% these are set by manually trans and rotating the bathy
coff = 12.0 * (pi/180); % compass offset
xoff = -1;
yoff = 2;
zoff = 1;

% load bathy
load covis_bathy_1.mat

% position of Covis in UTM coords
if(isfield(covis.sonar.position, 'easting'))
    utm_x0 = covis.sonar.position.easting;
    utm_y0 = covis.sonar.position.northing;
else
    fprintf('No sonar positions available\n');
    return;
end

alt = covis.sonar.position.altitude;
if(isfield(covis.sonar.position, 'depth'))
   depth0 = covis.sonar.position.depth;
else
   depth0 = griddata(utm_x, utm_y, depth, utm_x0, utm_y0);
end
x = utm_x - utm_x0;
y = utm_y - utm_y0;
z = depth - depth0;


R = sqrt(x.^2 + y.^2);
theta = atan2(y, x) + coff;
x = R.*cos(theta) + xoff;
y = R.*sin(theta) + yoff;
z = z + zoff;

%rotate(hb,[0,0,1],theta_off,[0,0,0]);
%rotate(hc,[0,0,1],theta_off,[0,0,0]);

% resample bathy onto a uniform rectangular grid
% this is necessary because utm_x, and utm_y are not necessarily uniform
xmin = grd.bounds.xmin;
xmax = grd.bounds.xmax;
ymin = grd.bounds.ymin;
ymax = grd.bounds.ymax;
dx = grd.spacing.dx;
dy = grd.spacing.dy;

xb = xmin:dx:xmax;
yb = ymin:dy:ymax;
[xb,yb] = meshgrid(xb,yb);
zb = griddata(x,y,z,xb,yb);

end

