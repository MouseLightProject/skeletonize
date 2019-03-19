function DisMatrix = skel_computeDisMatrix(data, boundary_pos, object_pos, boundary_num, object_num)
%
% Find minimum distance from each interior voxel to the boundary and save 
% into a 3D matrix of the same size of data. All voxels belong to boundary 
% or background have 0 value in the matrix.
%
%
% input:
% data:           original volume data
% boundary_pos:   coordinates of boundary points   
% object_pos:     coordinates of interior points    
% boundary_num:   number of boundary voxels
% object_num:     number of interior voxels
%
% output:
% DisMatrix is a 3D matrix of minimum distances between every interior voxel
% to boundary
%
% ---------------------------
% written by Li Liu in 10/05/2012 
% l.liu6819@gmail.com
%

[rows, cols, slices]=size(data);
Dist=ComputeDist(boundary_pos, object_pos, object_num, boundary_num, rows, cols, slices);
DisMatrix=reshape(Dist,[rows cols slices]);



