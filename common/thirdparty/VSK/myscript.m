testvol = load('U:\CODE\MATLAB\common\thirdparty\Skeleton3D\testvol.mat')
data = double(testvol.testvol);
th = .5;
AZ=-120;
EL=15;
voxel=zeros(size(data));
voxel(data>th)=1;  

[rows,cols,slices] = size(data);
[X,Y,Z] = meshgrid(1:cols, 1:rows, 1:slices);
skel_show3D(X,Y,Z,voxel,0.5,4,AZ,EL, 1);
title('original data');
rev = 0
show = voxel;

filter=fspecial3('average', 3);
voxel2=imfilter(voxel, filter);
voxel=voxel2;
voxel((voxel>0))=1;
skel_show3D(X,Y,Z,voxel,0.5,4,AZ,EL, 1);    
%%

method = 2
switch method
    case 1,
        disp('You have chosen Distance Matrix Method.');
        skel_Distmethod(voxel, show, rows, cols, slices, X, Y, Z, AZ, EL, rev);
    case 2,
        disp('You have chosen Repulsive Potential Field Method.');
        skel_RPFmethod(voxel, show, rows, cols, slices, X, Y, Z, AZ, EL, rev);
    case 3,
        disp('You have chosen Thinning Method.');
        skel_thinningmethod(voxel, show, rows, cols, slices, X, Y, Z, AZ, EL, rev);
end