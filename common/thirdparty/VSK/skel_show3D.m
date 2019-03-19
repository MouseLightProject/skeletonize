function skel_show3D(X,Y,Z,data,th,fig_number,AZ,EL,a)
%
% Compute distances between skeleton points found out by Distance Matrix 
% Transformation and determine the connection of each pair of skeleton
% points by using minimum spanning tree algorithm.
%
%
% input:
% X, Y, Z are coordinates of a 3D rectangular grid with the same size of
% the data for displaying
% data:        data for displaying
% th:          a value telling apart object and background 
% fig_number:  index of figure
% rows,cols,slices are 3 dimensions of voxel
% AZ and EL are both view angles of which AZ is the azimuth or horizontal
% rotation and EL is the vertical elevation (both in degrees).
% a:           alpha properties for objects in the current Axis
%
% ---------------------------
% written by Li Liu in 11/01/2012 
% l.liu6819@gmail.com
%

figure (fig_number)

p0 = patch(isosurface(X,Y,Z,data,th));
isonormals(X,Y,Z,data,p0);

if fig_number==1
    set(p0,'FaceColor','c','EdgeColor','none');
end

if fig_number==2
    set(p0,'FaceColor','c','EdgeColor','none');
end

if fig_number==3
    set(p0,'FaceColor','c','EdgeColor','none');
end

if fig_number==4
    set(p0,'FaceColor',[0,0.5,1],'EdgeColor','none');
end

if fig_number==5
    set(p0,'FaceColor','m','EdgeColor','none');
end

if fig_number>5
    set(p0,'FaceColor','c','EdgeColor','none');
end
 
view(AZ,EL); 
daspect;
daspect([1 1 1]);
axis tight
camlight; 
%lighting gouraud;
lighting phong
alpha(a);
% if a==1
    grid on
% else
%     box on
% end