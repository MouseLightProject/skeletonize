%
% Runs the Volume Skeleton Toolbox.
%
% ---------------------------
% This toolbox is free software distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY. You can redistribute it and/or modify it.
% Any modifications of the original software must be distributed in such a 
% manner as to avoid any confusion with the original work.
% 
% Please acknowledge the use of this software in any publications arising  
% from research that uses it.

% ---------------------------
% Author: Li Liu, Vizlab, Rutgers University (http://coewww.rutgers.edu/www2/vizlab/home) 
% l.liu6819@gmail.com
% 


% clear all
close all
clc

cell_list = {};
fig_number = 1;

title_figure = 'Volume Skeleton Toolbox. (Vizlab, Rutgers University)';

cell_list{1,1} = {'Compute Volume Skeleton','skel_main;'};
cell_list{2,1} = {'COVIS Data Processing','covis_platform;'};
cell_list{3,1} = {'View Volume Data','skel_displaydata;'};
cell_list{4,1} = {'View Volume Skeleton','skel_displayskeleton;'};
cell_list{5,1} = {'Exit',['disp(''Bye. To run again, type vsk.''); close(' num2str(fig_number) ');']};

show_window(cell_list,fig_number,title_figure,400,30,10,'clean',15);
