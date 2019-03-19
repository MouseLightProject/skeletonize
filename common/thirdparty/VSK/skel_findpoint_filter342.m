function [num_ridge, ridge, ridge_pos] = skel_findpoint_filter342(voxel, DisMatrix, rows, cols, slices, th)

contain = input('Please set the maximum number of voxels inside the 7*7*7 filtering kernel.([]=default: 20): ');

if isempty(contain)
    contain = 20;
end


ridge=Filter342(voxel, DisMatrix, rows, cols, slices, th, contain);
ridge=reshape(ridge,[rows, cols, slices]);

b=reshape(ridge,prod(size(ridge)),1);
[indx3,indy3,indz3]=ind2sub(size(ridge),find(b==255));
ridge_pos=[indx3,indy3,indz3];
num_ridge=size(ridge_pos,1);
