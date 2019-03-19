function [num_ridge, ridge, ridge_pos] = skel_findpoints_filter18(voxel, DisMatrix, rows, cols, slices, th)

contain = input('Please set the maximum number of voxels inside the 3*3*3 filtering kernel.([]=default: 4): ');

if isempty(contain)
    contain = 4;
end


ridge=Filter18(voxel, DisMatrix, rows, cols, slices, th, contain);
ridge=reshape(ridge,[rows, cols, slices]);

b=reshape(ridge,prod(size(ridge)),1);
[indx3,indy3,indz3]=ind2sub(size(ridge),find(b==255));
ridge_pos=[indx3,indy3,indz3];
num_ridge=size(ridge_pos,1);

