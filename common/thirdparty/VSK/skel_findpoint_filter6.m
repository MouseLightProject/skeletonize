function [num_ridge, ridge, ridge_pos] = skel_findpoint_filter6(voxel, DisMatrix, rows, cols, slices, th)

ridge=Filter6(voxel, DisMatrix, rows, cols, slices, th);
ridge=reshape(ridge,[rows, cols, slices]);

b=reshape(ridge,prod(size(ridge)),1);
[indx3,indy3,indz3]=ind2sub(size(ridge),find(b==255));
ridge_pos=[indx3,indy3,indz3];
num_ridge=size(ridge_pos,1);