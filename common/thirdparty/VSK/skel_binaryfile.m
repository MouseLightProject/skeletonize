function skel_binaryfile(data)
%
% If the data file you selected in "main function" is binary (only has 0
% and 1 value), this function will be called.
%
% input:
% data is the volume data file selected by user in "skel_main" function.
%
% ---------------------------
% written by Li Liu in 01/03/2013 
% l.liu6819@gmail.com
%

%%
foldername='Skeleton_results';
if(~exist(foldername,'dir'))
    mkdir(foldername);
end

his = input('Do you want to show the histogram of the data?  ([]=no, other = yes) ','s');
if isempty(his)==0
    figure
    hist(data);   
end

disp(' ');
th = input('Please set a threshold to the data ([]=default: (minimum value + maximum value)/2 )');
large = input('Is the value of interior voxel larger than that of background voxel?  ([] = yes, other = no) ','s');

if isempty(th)
    th = ( min(min(min(data))) + max(max(max(data))) )/2 ;
end

voxel=zeros(size(data));
if isempty(large)
    voxel(data>th)=1;  
else
    voxel(data<th)=1;
end

[rows, cols, slices] = size(data);
[X,Y,Z] = meshgrid(1:cols, 1:rows, 1:slices);
AZ = input('Please input the horizontal rotation angle w.r.t origin: ([]:default=-120): ');
if isempty(AZ)
    AZ=-120;
end

EL = input('Please input the vertical rotation angle w.r.t origin: ([]:default=15): ');
if isempty(EL)
    EL=15;
end

skel_show3D(X,Y,Z,voxel,0.5,3,AZ,EL,1);
title('original data');

disp(' ');
reverse = input('Do you want to reverse the data view?  ([]=no, other = yes) ','s');

if isempty(reverse)
    rev=0;   
else
    set(gca, 'ZDir','reverse');
    disp('The data view has been reversed!');
    rev=1;
end
str='original data';
filename = fullfile(foldername, [str '.' 'fig']);  
saveas(gcf, filename);
filename = fullfile(foldername, [str '.' 'jpg']);  
print(gcf, '-djpeg', filename);
  
show=voxel;

%%
decision = input('Do you want to use 3D filter to filter your data?  ([]=yes, other = no) ','s');

if isempty(decision)
    disp('Please choose a type of the filter from the following:');
    disp(' ');
    disp('******************************************************');
    disp('average     averaging filter');
    disp('ellipsoid   ellipsoidal averaging filter');
    disp('gaussian    Gaussian lowpass filter');
    disp('laplacian   Laplacian operator');
    disp('log         Log of Gaussian filter');
    disp('******************************************************');
    disp(' ');
    filter_type = input('Please choose the type of the filter: ','s');
    filter_size = input('Please set a size to your filter (usually 3~9): ');
    filter=fspecial3(filter_type, filter_size);
    voxel2=imfilter(voxel, filter);
    voxel=voxel2;
    voxel((voxel>0))=1;
    skel_show3D(X,Y,Z,voxel,0.5,4,AZ,EL, 1);    
    title('Filtered data');
    if rev==1
        set(gca, 'ZDir','reverse');
    end
    str='Filtered data';
    filename = fullfile(foldername, [str '.' 'fig']);  
    saveas(gcf, filename);
    filename = fullfile(foldername, [str '.' 'jpg']);  
    print(gcf, '-djpeg', filename);  
end

%%
disp(' ');
disp('Please choose a method to compute the skeleton:');
disp(' ');
disp('******************************************************');
disp('1     Distance Matrix Method');
disp('2     Repulsive Potential Field Method');
disp('3     Thinning Method');
disp('******************************************************');
disp(' ');

method = input('Please choose a method (1~3): ');

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