


%
% Display volume data and its histogram. 
%
% ---------------------------
% written by Li Liu in 01/07/2013 
% l.liu6819@gmail.com
%


clear all
clc

disp(' ');
disp('******************************************************');
disp('****** Welcome to use Volume Skeleton toolbox! *******');
disp('******************************************************');

disp(' ');
disp('Please open your volume data.');
disp(' ');

[FileName,PathName] = uigetfile('*.mat', '*.vol', 'Select data file');
str=[PathName, FileName];

[TYPE, OPENCMD, LOADCMD, DESCR] = finfo(str);

if strcmp(TYPE, 'mat')==1
    load([PathName FileName]);
    data = covis.grid.v;
    disp(' ');
    disp('Data loading is completed!');
    disp(' ');
    [rows, cols, slices] = size(data);
    data(isnan(data))=0; 
%     figure (2)
%     hist(data);
%     title('histogram of the data');
    disp(' ');
    th = input('Please set a threshold to the data ([]=default: 5*1e-9): ');
    
    if isempty(th)
        th = 5*1e-9;
    end

    th1=2*th;

    AZ = input('Please input the horizontal rotation angle w.r.t origin: ([]:default=305): ');
    
    if isempty(AZ)
        AZ=305;
    end

    EL = input('Please input the vertical rotation angle w.r.t origin: ([]:default=15): ');
    
    if isempty(EL)
        EL=15;
    end

    [XX,YY,ZZ] = meshgrid(1:cols, 1:rows, 1:slices);
    skel_show3D(XX,YY,ZZ,data,th1,3,AZ,EL,1);
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
end


if strcmp(TYPE, 'vol')==1 || strcmp(TYPE, 'img')==1 || strcmp(TYPE, 'unknown')==1
    fid=fopen([PathName FileName]);
    disp(' ');
    disp('Data loading is completed!');
    disp(' ');
    disp(' ');
    disp('******************************************************');
    disp('1     uint8');
    disp('2     uint16');
    disp('3     uint32');
    disp('4     double');
    disp('******************************************************');
    disp(' ');

    datatype = input('Please choose a the datatype of your data (1~4): ');

    switch datatype
        case 1,
            src=fread(fid);
        case 2,
            src=fread(fid,'uint16');
        case 3,
            src=fread(fid,'uint32');
        case 4,
            src=fread(fid,'double');
    end

    disp(' ');
    rows = input('Please input the rows of the data: ');
    cols = input('Please input the columns of the data: ');
    slices = input('Please input the slices of the data: ');
    data=reshape(src, [rows cols slices]);
    disp(' ');
    disp('Data loading is completed!');
    figure (2)
    hist(data);
    title('histogram of the data');
    disp(' ');
    th = input('Please set a threshold to the plume data ([]=default: (minimum value + maximum value)/2: ');
    
    if isempty(th)
        th = ( min(min(min(data))) + max(max(max(data))) )/2 ;
    end

    voxel=zeros(size(data));
    voxel(data>th)=1;  
    
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

end