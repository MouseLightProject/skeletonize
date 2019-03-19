function [rx, ry, rz]=skel_computeRPF(object_pos, boundary_pos, boundary_num, order) 
%
% Compute Repulsive Potential Force (RPF) for each interior voxel by charge of
% whole boundary. RPF is defined as a force pushes away the object in the
% field away from the charge.
%
%
% input:
% object_pos:     coordinates of interior points  
% boundary_pos:   coordinates of boundary points   
% boundary_num:   number of boundary voxels
% order:          order of the Repulsive Potential Field
%
% output:
% rx, ry, and rz are three components of RPF in each interior voxels
%
% ---------------------------
% written by Li Liu in 11/12/2012 
% l.liu6819@gmail.com
%

temp = boundary_pos - repmat(object_pos, boundary_num, 1);
x=temp(:,1);
y=temp(:,2);
z=temp(:,3);
distance = sqrt(x.^2 + y.^2 + z.^2);
r_n=power(distance, order);

rx = sum(-x./r_n);
ry = sum(-y./r_n);
rz = sum(-z./r_n);



