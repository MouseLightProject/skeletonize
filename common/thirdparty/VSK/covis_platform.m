%
%
%This is a special module of VSK toolbox used for group processing COVIS
%acoustic data. For each data, this program first computes skeleton of the 
%plumes, then fits a line for each plume skeleton by using Least Square Method 
%and providing the bending angles of the lines. After processing all data 
%selected, the program will make group analysis. 
%
%data for processing:
%The data for processing should be COVIS imaging, doppler, or diffuse acoustic
%data that has been transformed into MAT format. 
%(e.g. APLUWCOVISMBSONAR001_20101001T000511.281Z-IMAGING_masked.mat)
%
%The statiscal analysis need at least two data files be selected and the data  
%files all recorded at the same time (month/day/hour). If only one data be 
%selected, only skeleton and fitting line will be displayed.
%
%
% --- ------------------------
% written by Li Liu in 01/08/2013 
% l.liu6819@gmail.com
%

clear all
close all
clc

disp(' ');
disp('***************************************************************');
disp('****** Welcome to COVIS Plume Data Processing Platform! *******');
disp('***************************************************************');

%%
%%%%%%%%%%%%%%%%%%%% select data files and load parametres from JSON file %%%%%%%%%%%%%%%%%%%

[FileName,PathName] = uigetfile( {'*imaging*.mat',  'imaging file (*IMAGING*.mat)';...
                                  '*dopper*.mat', 'doppler file (*DOPPLER*.mat)';...
                                  '*diffuse*.mat', 'diffuse file (*DIFFUSE*.mat)'}, 'MultiSelect', 'on');

if iscell(FileName)==0                              
    while FileName==0
         [FileName,PathName] = uigetfile( {'*imaging*.mat',  'imaging file (*IMAGING*.mat)';...
                                           '*dopper*.mat', 'doppler file (*DOPPLER*.mat)';...
                                           '*diffuse*.mat', 'diffuse file (*DIFFUSE*.mat)'}, 'MultiSelect', 'on');
    end
end
                              
if iscell(FileName)==0                                    % determine the number of selected files
    filecount=1;
else
    filecount=size(FileName,2);
end

fprintf('\nYou chose %d COVIS data sets.\n\n',filecount);

tic

if filecount==1                                           % save the path of the data, if more than one data selected, save them all into a cell
    filepath = fullfile(PathName, FileName);
else
    for ii = 1:filecount                                      
        filepath{1,ii} = fullfile(PathName, FileName{ii});
    end
end


json_str = fileread('covis_skeleton.json');               % parameters are defined with the corresponding JSON parameter files.
input = parse_json(json_str);                             % parses a JSON string and returns parametres as cell arrays 

AZ=input.view.azimuth;                                    % azimuth or horizontal rotation degree
EL=input.view.elevation;                                  % vertical elevation degree
th=input.threshold;                                       % threhold of transforming original data to binary data
ker=input.morph_kernel;                                   % size of structure element for morphological operation (e.g. dilation, erosion)

mindist=input.skel_point.min_dist;                        % the minimum distance a skeleton point should be away from boundary
maxnum=input.skel_point.max_num;                          % the maximum number of local maximum point allowed within a 3*3*3 structure element                
sample=input.sample_num;                                  % the maximum number of samples between two skeleton points
line_lower=input.lsm_range.min;                           % the lower bounding of skeleton for fitting a line
line_upper=input.lsm_range.max;                           % the upper bounding of skeleton for fitting a line

folder_base=input.result_save.folder_base;                % the base name of the folder to be created 
saved_file1=input.result_save.file1;
saved_file2=input.result_save.file2;
saved_file3=input.result_save.file3;
saved_file4=input.result_save.file4;
saved_file5=input.result_save.file5;
saved_file6=input.result_save.file6;


Foldername=strcat(folder_base, date);
mkdir(Foldername);

if(~exist(Foldername,'dir'))
    mkdir(Foldername);
end


%%%%%%%%%%%%%%%%%%%%%%%%% begin to process each selected data file %%%%%%%%%%%%%%%%%%%%%%%%%

for data_num=1:filecount
    
    close all
    fprintf('Begin to process %dth data...\n',data_num);
    
    if filecount==1
       load(filepath);
    else
       load(filepath{1,data_num});
    end
    
    xmin=covis.grid.bounds.xmin;
    ymin=covis.grid.bounds.ymin;
    zmin=covis.grid.bounds.zmin;
    xmax=covis.grid.bounds.xmax;
    ymax=covis.grid.bounds.ymax;
    zmax=covis.grid.bounds.zmax;
    
    scale=(covis.grid.size(3)-1)/(zmax-zmin);                   % the inverse of the voxel edge length in meter
    
    Name = covis.grid.name;
    deliminater1 = strfind(Name,'_');
    deliminater2 = strfind(Name,'.');
    date_info =  Name( deliminater1(1)+1 : deliminater2-1);
    
    month{data_num}=date_info(5:6);                             % extract date information from file name
    day{data_num}=date_info(7:8);
    hour{data_num}=date_info(10:11);
    Month=str2num(month{data_num});
    Day=str2num(day{data_num});
    Hour=str2num(hour{data_num});
    
    datevalue(data_num)=covis_date2number(Month,Day,Hour);      % transform the input date (month/day/hour) into a numerical value   
    
    NameGroup{data_num}=Name;
    Name(strfind(Name,'_'))='-';
    foldername=fullfile(Foldername, Name);                      % create a folder to save results for current processing data file 
    mkdir(foldername);   
    
    data = covis.grid.v;
    [rows, cols, slices] = size(data);
    data(isnan(data))=0; 

    XX=covis.grid.x;
    YY=covis.grid.y;
    ZZ=covis.grid.z;
    
    %%
    %%%%%%%%%% data pre-processing %%%%%%%%%%  
    
    voxel=zeros(size(data));    
    voxel(data>th)=1;                                           % threshold original data and get a binary data matrix voxel in which 0 indicates a voxel 
                                                                % belongs to background and 1 for object                                                        
 
    voxel_dilated=imdilate(voxel, ones(ker,ker,ker));           % morphologically dilate the voxel with structure element of size ker*ker*ker
    voxel_eroded=imerode(voxel_dilated, ones(ker,ker,ker));     % morphologically erode the voxel with structure element of size ker*ker*ker
    v=voxel_eroded;
    voxel=v;
    
    [object_pos, boundary_pos, object, boundary] = segmentation(voxel, rows, cols, slices);  % classify the object in a binary volume data into 
                                                                                             % interior voxels and boundary voxels
    num_object=size(object_pos, 1);
    num_boundary=size(boundary_pos, 1);
    
    %%
    %%%%%%%%%% find skeleton point by using Distance Matrix Method %%%%%%%%%%
    
    Dist=ComputeDist(boundary_pos, object_pos, num_object, num_boundary, rows, cols, slices);   % Compute minimum distance from each interior voxel to the boundary
    DisMatrix=reshape(Dist,[rows cols slices]);
    
    ridge=Filter26(voxel, DisMatrix, rows, cols, slices, mindist, maxnum);                   % find out the skeleton points defined as local maxima of DisMatrix                
    ridge=reshape(ridge,[rows, cols, slices]);  
    b=reshape(ridge,prod(size(ridge)),1);
    [indx3,indy3,indz3]=ind2sub(size(ridge),find(b==255));
    ridge_pos=[indx3,indy3,indz3];
    num_ridge=size(ridge_pos,1);                                                             % count the number of skeleton points
    
    %%
    %%%%%%%%%% connect skeleton points to get skeleton %%%%%%%%%%
    
    T=connection(voxel, num_ridge, ridge_pos, rows, cols, slices);          % compute distances between skeleton points and determine the connection 
                                                                            % of each pair of skeleton points by using minimum spanning tree algorithm
    [m,n]=size(T);
    [h, covis] = covis_imaging_plot2(covis, 'covis_image_plot.json', 3);    % plot plumes and their bath
    hold on

    con=zeros(size(ridge_pos,1));
    
    for i=1:m                                                               % plot skeleton points and skeleton
        for j=i:n
            if T(i,j)==1               
                con(i)=1;
                con(j)=1;
                start_x=ridge_pos(i,2);
                start_y=ridge_pos(i,1);
                start_z=ridge_pos(i,3);
                end_x=ridge_pos(j,2);
                end_y=ridge_pos(j,1);
                end_z=ridge_pos(j,3);
                
                start_x=round( (start_x-1)/scale+xmin );
                start_y=round( (start_y-1)/scale+ymin );
                start_z=round( (start_z-1)/scale+zmin );
                end_x=round( (end_x-1)/scale+xmin );
                end_y=round( (end_y-1)/scale+ymin );
                end_z=round( (end_z-1)/scale+zmin );
            
                plot3([start_x,end_x],[start_y,end_y],[start_z,end_z], '-rs','LineWidth',1, 'MarkerEdgeColor','k', 'MarkerFaceColor','y',  'MarkerSize',3);
                hold on
            end
        end
    end
    hold off   
    
    str=strcat(saved_file1, Name);
    str(strfind(str,'_'))='-';
    title(str);
    filename = fullfile(foldername, [str '.' 'fig']);  
    saveas(gcf, filename);

    connect_num=sum(sum(T))/2;
    [sk_X,sk_Y,sk_Z]=GetSkel(voxel,T,ridge_pos,connect_num,rows,cols,slices,sample);     % get skeleton coordinates by interpolating points between each
                                                                                   % each pair connected skeleton points (sample is the maximum number of
                                                                                   % interpolated points for any connection line)
    skel_X=reshape(sk_X, [sample, connect_num]);
    skel_Y=reshape(sk_Y, [sample, connect_num]);
    skel_Z=reshape(sk_Z, [sample, connect_num]);
    skelet_X=skel_X';
    skelet_Y=skel_Y';
    skelet_Z=skel_Z';
    skeleton_X = (skelet_X-ones(size(skelet_X)))./scale + ymin*ones(size(skelet_X));    % transform the skeleton coordinates from voxel to meter
    skeleton_Y = (skelet_Y-ones(size(skelet_Y)))./scale + xmin*ones(size(skelet_Y));
    skeleton_Z = (skelet_Z-ones(size(skelet_Z)))./scale + zmin*ones(size(skelet_Z));
    
    str=strcat(saved_file2, Name);                                                 % save skeleton coordinates into MAT file
    str(strfind(str,'_'))='-'; 
    filename = fullfile(foldername, [str '.' 'mat']);   
    s=struct('Name',Name,'X',skeleton_X,'Y',skeleton_Y,'Z',skeleton_Z,'AZ',AZ,'EL',EL,'Reverse',0);
    save(filename, '-struct', 's');
   
    %%
    %%%%%%%%%% find a line to fit skeleton points by using Least Square Method %%%%%%%%%%
    
    count1=0;
    count2=0;
    rig1=[];
    rig2=[];   
    x1=[];
    x2=[];
    y1=[];
    y2=[];
    z1=[];
    z2=[];
                             
    for i=1:num_ridge                      
        if ridge_pos(i,2)>(cols+1)/2 && con(i)==1                           % tell the two plume skeletons apart
            count1=count1+1;
            rig1(count1,:)=ridge_pos(i,:);
        end
        if ridge_pos(i,2)<=(cols+1)/2 && con(i)==1 && ridge_pos(i,3)>=((line_lower-zmin)*scale+1) && ridge_pos(i,3)<=((line_upper-zmin)*scale+1)  % set range of big plume skeleton to fit 
            count2=count2+1;
            rig2(count2,:)=ridge_pos(i,:);
        end
    end
    
    
    if count1>=2                                                            % if there are at least two connected skeleton points for small plume, fit the line
        [a1,b1,c1,d1]=covis_fitline3D(rig1);
        x1 = a1*rig1(:,3) + b1;
        y1 = c1*rig1(:,3) + d1;
        z1 = rig1(:,3);
        x1 = (x1-ones(size(x1)))./scale + ymin*ones(size(x1));
        y1 = (y1-ones(size(y1)))./scale + xmin*ones(size(y1));
        z1 = (z1-ones(size(z1)))./scale + zmin*ones(size(z1));
        rig1(:,1) = (rig1(:,1)-ones(size(rig1(:,1))))./scale + ymin*ones(size(rig1(:,1)));
        rig1(:,2) = (rig1(:,2)-ones(size(rig1(:,2))))./scale + xmin*ones(size(rig1(:,2)));
        rig1(:,3) = (rig1(:,3)-ones(size(rig1(:,3))))./scale + zmin*ones(size(rig1(:,3)));
    end
    
    if count2>=2                                                            % if there are at least two connected skeleton points for big plume, fit the line
        [a2,b2,c2,d2]=covis_fitline3D(rig2);
        x2 = a2*rig2(:,3) + b2;
        y2 = c2*rig2(:,3) + d2;
        z2 = rig2(:,3);    
        x2 = (x2-ones(size(x2)))./scale + ymin*ones(size(x2));
        y2 = (y2-ones(size(y2)))./scale + xmin*ones(size(y2));
        z2 = (z2-ones(size(z2)))./scale + zmin*ones(size(z2));   
        rig2(:,1) = (rig2(:,1)-ones(size(rig2(:,1))))./scale + ymin*ones(size(rig2(:,1)));
        rig2(:,2) = (rig2(:,2)-ones(size(rig2(:,2))))./scale + xmin*ones(size(rig2(:,2)));
        rig2(:,3) = (rig2(:,3)-ones(size(rig2(:,3))))./scale + zmin*ones(size(rig2(:,3)));
    end

    [h, covis] = covis_imaging_plot2(covis, 'covis_image_plot.json', 5);    
    hold on

    if count1>=2                                                            % plot fitting line
        plot3(y1,x1,z1,'k','linewidth', 1.5);
        hold on
    end
    if count2>=2
        plot3(y2,x2,z2,'k','linewidth', 1.5);
        hold on
    end
    hold off
    
    str=strcat(saved_file3, Name);
    str(strfind(str,'_'))='-'; 
    title(str);
    filename = fullfile(foldername, [str '.' 'fig']);
    saveas(gcf, filename); 

    %%
    %%%%%%%%%% compute the bending angle for fitting line %%%%%%%%%%
    
    if count1>=2
        [TH,R,Z] = cart2pol(a1,c1,1);
        if TH>=0 && TH<pi/2
            angle_horizontal1(data_num) = (pi/2-TH)*180/pi;
        elseif TH>=pi/2 && TH<pi
            angle_horizontal1(data_num) = (2*pi+pi/2-TH)*180/pi;
        elseif TH>=-pi && TH<-pi/2
            angle_horizontal1(data_num) = (pi/2-TH)*180/pi;
        elseif TH>=-pi/2 && TH<0
            angle_horizontal1(data_num) = (pi/2-TH)*180/pi;
        end
        angle_vertical1(data_num)= acos(1/sqrt(a1*a1+c1*c1+1))*180/pi;
    else
        angle_horizontal1(data_num) = 0;
        angle_vertical1(data_num) = 0;
    end
         
    if count2>=2
        [TH,R,Z] = cart2pol(a2,c2,1);
        if TH>=0 && TH<pi/2
            angle_horizontal2(data_num) = (pi/2-TH)*180/pi;
        elseif TH>=pi/2 && TH<pi
            angle_horizontal2(data_num) = (2*pi+pi/2-TH)*180/pi;
        elseif TH>=-pi && TH<-pi/2
            angle_horizontal2(data_num) = (pi/2-TH)*180/pi;
        elseif TH>=-pi/2 && TH<0
            angle_horizontal2(data_num) = (pi/2-TH)*180/pi;
        end
        angle_vertical2(data_num) = acos(1/sqrt(a2*a2+c2*c2+1))*180/pi;
    else
        angle_horizontal2(data_num) = 0;
        angle_vertical2(data_num) = 0;
    end
    
end

time=toc;

disp(' ');
disp('Data Processing is completed!');
disp(' ');
fprintf('Elapsed time is %f seconds.\n',time);

%%
%%%%%%%%%% compute the mean and standard deviation of the bending angle %%%%%%%%%%

if filecount>1
    
    for i=1:filecount
        Mean_Horizontal1=mean(angle_horizontal1);
        Mean_Horizontal2=mean(angle_horizontal2);
        Mean_Vertical1=mean(angle_vertical1);
        Mean_Vertical2=mean(angle_vertical2);
        Std_Horizontal1=std(angle_horizontal1);
        Std_Horizontal2=std(angle_horizontal2);
        Std_Vertical1=std(angle_vertical1);
        Std_Vertical2=std(angle_vertical2);  
    end

%%
%%%%%%%%%% display and save statistical analysis results %%%%%%%%%%

    Start=find(min(datevalue));
    End=find(max(datevalue));
    figure('Name','Bending Angle of small Plume','NumberTitle','off')
    subplot(211);
    for i=1:filecount
        plot(datevalue(i), angle_horizontal1(i), 'r*');
        hold on
    end
    hold off
    m1=min(datevalue);
    m2=max(datevalue);
    axis([m1-str2num(hour{Start}) m2 0 400]);
    h=((m1-str2num(hour{Start})): 24 : m2);
    set(gca, 'XTick', h);

    for i=1:size(h,2)
        stringx{i}=covis_number2date(h(i));
    end
    
    set(gca,'XTickLabel',stringx);
    xlabel('Date');
    ylabel('Degree');
    title('Bending angle with respect to Positive Y-axis');

    subplot(212);
    for i=1:filecount
        plot(datevalue(i), angle_vertical1(i), 'b*');
        hold on
    end
    hold off
    axis([m1-str2num(hour{Start}) m2 0 90]);
    set(gca, 'XTick', h);
    set(gca,'XTickLabel',stringx);
    xlabel('Date');
    ylabel('Degree');
    title('Bending angle with respect to Z=0 plane');

    str=strcat('Bending_small plume_', datestr(now,30));
    filename = fullfile(Foldername, [str '.' 'jpg']);   
    print(gcf, '-djpeg', filename);
 
    figure('Name','Bending Angle of Big Plume','NumberTitle','off')
    subplot(211);
    for i=1:filecount
        plot(datevalue(i), angle_horizontal2(i), 'r*');
        hold on
    end
    hold off
    axis([m1-str2num(hour{Start}) m2 0 400]);
    set(gca, 'XTick', h);
    set(gca,'XTickLabel',stringx);
    xlabel('Date');
    ylabel('Degree');
    title('Bending angle with respect to Positive Y-axis');
    
    subplot(212);
    for i=1:filecount
        plot(datevalue(i), angle_vertical2(i), 'b*');
        hold on
    end
    
    hold off
    axis([m1-str2num(hour{Start}) m2 0 90]);
    set(gca, 'XTick', h);
    set(gca,'XTickLabel',stringx);
    xlabel('Date');
    ylabel('Degree');
    title('Bending angle with respect to Z=0 plane');

    str=strcat('Bending_big plume_', datestr(now,30));
    filename = fullfile(Foldername, [str '.' 'jpg']);   
    print(gcf, '-djpeg', filename);

    str=strcat('proecssing report_', datestr(now,30));                      % generate the processing report
    filename = fullfile(Foldername, [str '.' 'txt']); 
    fid = fopen(filename,'wt+');
    fprintf(fid,'\nDate: %s.\n\n',datestr(now,31));
    fprintf(fid,'Elapsed time: %f seconds.\n\n',time);
    fprintf(fid,'Data number: %d.\n\n\n\n',filecount);
    
    for i=1:ceil((size(NameGroup{1},2)+4)/8)+1
        fprintf(fid,'\t');
    end
    fprintf(fid,'%s\t','vertical_big');
    % fprintf(fid,'\t\t\t\t\t\t\t\t\t%s\t','vertical_big');
    fprintf(fid,'\t%s\t','horizontal_big');
    fprintf(fid,'\t%s\t','vertical_small');
    fprintf(fid,'\t%s\n\n','horizontal_small');

    for i = 1:filecount
        fprintf(fid,'%d: %s\t\t',i,NameGroup{i}); 
        fprintf(fid,'%f\t\t',angle_vertical2(i));
        fprintf(fid,'%f\t\t',angle_horizontal2(i));
        fprintf(fid,'%f\t\t',angle_vertical1(i));
        fprintf(fid,'%f\n',angle_horizontal1(i));
    end

    fprintf(fid,'\n\t\t\t\t%s\t\t\t\t\t','mean');
    fprintf(fid,'%f\t\t',Mean_Vertical2);
    fprintf(fid,'%f\t\t',Mean_Horizontal2);
    fprintf(fid,'%f\t\t',Mean_Vertical1);
    fprintf(fid,'%f\n',Mean_Horizontal1);
    fprintf(fid,'\t\t\t\t%s\t\t\t\t\t','std');
    fprintf(fid,'%f\t\t',Std_Vertical2);
    fprintf(fid,'%f\t\t',Std_Horizontal2);
    fprintf(fid,'%f\t\t',Std_Vertical1);
    fprintf(fid,'%f\n',Std_Horizontal1);

    fclose(fid);

end

