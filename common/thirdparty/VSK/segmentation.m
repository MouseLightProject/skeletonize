function [object_pos, boundary_pos, object, boundary] = segmentation(voxel,rows, cols, slices)
%
% Classify the object in a binary volume data into interior voxels and
% boundary voxels.
%
%
% input:
% voxel: a binary 3D matrix thresholded from original data
% rows,cols,slices are 3 dimensions of voxel
%
% output:
% object_pos:   coordinates of all interior voxels
% boundary_pos: coordinates of all boundary voxels
% object:       a 3D matrix of the same size of voxel in which all voxels 
%               equal to 0 execept for those belong to interior part
% boundary:     a 3D matrix of the same size of voxel in which all voxels 
%               equal to 0 execept for those belong to boundary part
%
% ---------------------------
% written by Li Liu in 10/05/2012 
% l.liu6819@gmail.com
%

object=voxel;
boundary=Segment(voxel,rows, cols, slices);
boundary=reshape(boundary,[rows, cols, slices]);

b=reshape(boundary,prod(size(boundary)),1);
[indx1,indy1,indz1]=ind2sub(size(boundary),find(b==1));
object_pos=[indx1,indy1,indz1];

[indx2,indy2,indz2]=ind2sub(size(boundary),find(b==255));
boundary_pos=[indx2,indy2,indz2];

boundary(boundary==1)=0;
object(boundary==255)=0;



