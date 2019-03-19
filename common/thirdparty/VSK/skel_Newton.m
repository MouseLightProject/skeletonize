function [coordinate,inside] = skel_Newton(start_point, boundary_pos, boundary_num, xmin, ymin, zmin, xmax, ymax, zmax, RPF_x, RPF_y, RPF_z, RPForder, NewtonTH)
%
% Determine if a critical point lines in a candidate grid cell.
% A candidate grid cell is a cube where all three field vectors change sign
% at some of its vertex.
% 
%
% input:
% start_point:   starting point for Newton iterative search
% boundary_num:  number of boundary voxels
% xmin, ymin, zmin, xmax, ymax, zmax are the size of grid cell
% RPF_x, RPF_y, and RPF_z are x,y,and z component of Repulsive Potential
% boundary_pos:  coordinates of boundary points
% boundary_num:  number of boundary voxels
% i, j, and k are coordinates of inspecting point
% RPF_order:     order of RPF
% NewtonTH:      convergence criteria of Newton's method
%
% output:
% coordinate:    coordinates of the critical point
% inside:        a boolean value indicating if the inspecting grid cell
%                contain a critical point 
%
% ---------------------------
% written by Li Liu in 11/16/2012 
% l.liu6819@gmail.com
%

coordinate=[0,0,0];
inside=0;
P=start_point';
[rx, ry, rz]=skel_computeRPF(start_point, boundary_pos, boundary_num, RPForder); 
time=0;

% sum_old=rx*rx+ry*ry+rz*rz;
% sum_new=sum_old;
% rate=sum_new/sum_old;

while rx*rx+ry*ry+rz*rz>NewtonTH
% while rate>th
    
    time=time+1;
    F=[rx, ry, rz]';
    
    xx=floor(P(1)+0.5);
    yy=floor(P(2)+0.5);
    zz=floor(P(3)+0.5);
    
    if xx+1>size(RPF_x,1)
        xx=size(RPF_x,1)-1;
    end
    
    if yy+1>size(RPF_x,2)
        yy=size(RPF_x,2)-1;
    end
    
    if zz+1>size(RPF_x,3)
        zz=size(RPF_x,3)-1;
    end
    
    
    dRXdx = (RPF_x(xx+1, yy, zz)-rx) / (xx+1-P(1));
    dRXdy = (RPF_x(xx, yy+1, zz)-rx) / (yy+1-P(2));
    dRXdz = (RPF_x(xx, yy, zz+1)-rx) / (zz+1-P(3));
    
    dRYdx = (RPF_y(xx+1, yy, zz)-ry) / (xx+1-P(1));
    dRYdy = (RPF_y(xx, yy+1, zz)-ry) / (yy+1-P(2));
    dRYdz = (RPF_y(xx, yy, zz+1)-ry) / (zz+1-P(3));
    
    dRZdx = (RPF_z(xx+1, yy, zz)-rz) / (xx+1-P(1));
    dRZdy = (RPF_z(xx, yy+1, zz)-rz) / (yy+1-P(2));
    dRZdz = (RPF_z(xx, yy, zz+1)-rz) / (zz+1-P(3));
    
    J=[dRXdx dRXdy dRXdz
       dRYdx dRYdy dRYdz
       dRZdx dRZdy dRZdz];
   
    dP=-F./(J*P);
   
    P=P+dP;
    
    if xmin>P(1) || ymin>P(2) || zmin>P(3) || P(1)>xmax || P(2)>ymax || P(3)>zmax
        inside=0;
        break
    end
    
    [rx, ry, rz]=skel_computeRPF(P', boundary_pos, boundary_num, RPForder); 
%     
%     sum_new=rx*rx+ry*ry+rz*rz;
%     rate=sum_new/sum_old;   
%     
    if time>30
        break
    end
    
end

if rx*rx+ry*ry+rz*rz<=NewtonTH
% if rate<=th
    inside=1;
    coordinate=P';
end

