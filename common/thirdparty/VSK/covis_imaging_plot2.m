function [h, covis] = covis_imaging_plot2(covis, json_file,fig_num)
%
% Plot covis image grid as isosurface
%
% Plotting parameters are defined with the corresponding 
% JSON parameter files (json_file).  If no json file is supplied, 
% the default covis_image_plot.json file is used.
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
% ---------------------------
% This is a modified version of Dr. Christopher Jones's covis_imaging_plot  
% cjones@apl.washington.edu
% l.liu6819@gmail.com
%
% parsing the json file
% which contains all the user supplied parameters


if(isempty(json_file) | (json_file==0)) 
   json_file = fullfile('input','covis_image_plot.json')
end
json_str = fileread(json_file);
input = parse_json(json_str);

if(~isfield(covis,'grid'))
   return;
end

if(~isfield(input,'plot_bathy'))
    plot_bathy = 1;
end
plot_bathy = input.plot_bathy;

if(~isfield(input,'isosurface'))
   return;
end
isosurf = input.isosurface;

if(~isfield(input,'name')) 
   input.name = covis.grid.name;
end

% figure type, must be one of formats that the 'saveas' function accepts
if(~isfield(input,'format'))
   input.format = 'fig';
end
type = input.format;

h = figure(fig_num); clf;

% make local copies of the grids
grd = covis.grid; % intensity grid

% grid values
vg = grd.v;
xg = grd.x;
yg = grd.y;
zg = grd.z;

% load bathy
[xb,yb,zb] = covis_bathy(covis, grd);

% plot the bathymetry map
hold on
hb = surf(xb, yb, zb);
shading interp;
colormap('Summer')
[c,hc] = contour3(xb, yb, zb, 50,'k');
hold on
%hold off
%colorbar
%title('COVIS Bathymetry Map for Grotto');
%axis equal;

mask = zeros(size(vg));

% mask the data grid using the user inputs
if(isfield(input,'mask'))
    for n=1:length(input.mask)
        switch input.mask{n}.type
            case'cylinder'
                x0 = input.mask{n}.x;
                y0 = input.mask{n}.y;
                R = input.mask{n}.radius;
                r = sqrt((xg-x0).^2 + (yg-y0).^2);
                mask(find(r <= R)) = 1;
            case'cone'
                x0 = input.mask{n}.x;
                y0 = input.mask{n}.y;
                z0 = input.mask{n}.z;
                r = sqrt((xg-x0).^2 + (yg-y0).^2);
                R = input.mask{n}.radius;
                H = input.mask{n}.height;
                for m = 1:size(zg,3)
                    zslice = mask(:,:,m);
                    rc = (zg(1,1,m)-z0)*(R/H);
                    zslice(find(r(:,:,m) <= rc)) = 1;
                    mask(:,:,m) = zslice;
                end
            case'parabola'
                x0 = input.mask{n}.x;
                y0 = input.mask{n}.y;
                r0 = input.mask{n}.radius;
                r = sqrt((xg-x0).^2 + (yg-y0).^2);
                mask(find(r <= r0)) = 1;
        end
    end
end

% if mask is empty then unmask everything
if(isempty(find(mask)))
    mask = ones(size(vg));
end

% mask data below the bathy
mask_bathy = 1;
if(mask_bathy)
    dz = 1;
    for m = 1:size(zg,3)
        zslice = ones(size(zg,1),size(zg,2));
        k = find(zg(:,:,m) < (zb+dz));
        zslice(k) = 0;
        mask(:,:,m) = mask(:,:,m) .* zslice;
    end
end

vg(find(~mask)) = nan;

% plot isosurfaces
%figure(fig_num); 
%hold on;
for n=1:length(isosurf)
   v = vg;
   if(strcmp(isosurf{n}.units, 'db'))
      eps = nan;
      m = find(v==0);
      v(m) = eps;  % remove zeros
      v = 10*log10(v);
   end
   surf_value = isosurf{n}.value;
   surf_color = cell2mat(isosurf{n}.color);
   surf_alpha = isosurf{n}.alpha;
   p = patch(isosurface(xg, yg, zg, v, surf_value));
   isonormals(xg, yg, zg, v, p);
   %set(p,'FaceColor','red','EdgeColor','none');
   set(p,'EdgeColor','none','FaceColor',surf_color,'FaceAlpha',surf_alpha);
end
daspect([1 1 1])
hold off;
grid on

axis([-40 10 -40 10 0 40]);

xlabel('meters')
ylabel('meters')
zlabel('meters')

if(isfield(grd,'name'))
   str = grd.name;
   str(strfind(str,'_'))='-';
   title(str);
end

% set the view
if(isfield(input,'view'))
   view(input.view.azimuth, input.view.elevation);
else
   view([1,1,1]);
end


end

