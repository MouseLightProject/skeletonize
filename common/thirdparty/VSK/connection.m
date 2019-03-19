function T=connection(voxel, ridge_num, ridge_pos, rows, cols, slices)
%
% Compute distances between skeleton points found out by Distance Matrix 
% Transformation and determine the connection of each pair of skeleton
% points by using minimum spanning tree algorithm.
%
%
% input:
% voxel:     a binary 3D matrix thresholded from original data
% ridge_num: number of skeleton points 
% ridge_pos: coordinates of each skeleton points
% rows,cols,slices are 3 dimensions of voxel
%
% output:
% T is a 2D matrix of which each dimension equals to the number of skeleton
% points. The items of T could have 3 different values: 0, 1, Inf.
%
% ---------------------------
% written by Li Liu in 11/16/2012 
% l.liu6819@gmail.com
%

Dist=Connect(voxel,ridge_pos,ridge_num,rows, cols, slices);                 % Compute distances between each pair of skeleton points

% If a connection line traverse background voxels, the distance between the two nodes the line 
% connected with equals to Inf.
Dist(Dist==100000)=inf;
Dist=reshape(Dist,[ridge_num ridge_num]);

% Apply minimum spanning tree algorithm to determine the connection of each pair of skeleton points
[T, cost] = minispanningtree(Dist);
