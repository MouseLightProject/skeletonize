function [permit,critical] = skel_computecritical(RPF_x, RPF_y, RPF_z, boundary_pos, boundary_num, i, j, k, RPForder, NewtonTH)
%
% Determine if a point is critical point. Critical point is defined as a
% point where all components of the field simultaneously vanish.
%
%
% input:
% RPF_x, RPF_y, and RPF_z are x,y,and z component of Repulsive Potential
% boundary_pos:  coordinates of boundary points
% boundary_num:  number of boundary voxels
% i, j, and k are coordinates of inspecting point
% RPF_order:     order of RPF
% NewtonTH:      convergence criteria of Newton's method
%
% output:
% attracting:   coordinates of the attracting skeleton points
% repelling:    coordinates of the repelling skeleton points
% saddle1:      coordinates of the saddle skeleton points of which Jacobian 
%               matrix have two positive real eigen values
% saddle2:      coordinates of the saddle skeleton points of which Jacobian 
%               matrix have two negtive real eigen values
% attr_num:     number of attracting point
% saddle1_num:  number of the saddle1 point
% saddle2_num:  number of the saddle2 point
% pathline1:    eigen vector of attracting points       
% pathline2a:   eigen vector of saddle1 points
% pathline2b:   eigen vector of saddle2 points
%
% ---------------------------
% written by Li Liu in 11/15/2012 
% l.liu6819@gmail.com
%


            sign_x=0;
            sign_y=0;
            sign_z=0;
            critical=[0,0,0];
            permit=0;

            p1_x = RPF_x(i,j,k);
            p1_y = RPF_y(i,j,k);
            p1_z = RPF_z(i,j,k);
            
            p2_x = RPF_x(i+1,j,k);
            p2_y = RPF_y(i+1,j,k);
            p2_z = RPF_z(i+1,j,k);
            
            p3_x = RPF_x(i+1,j+1,k);
            p3_y = RPF_y(i+1,j+1,k);
            p3_z = RPF_z(i+1,j+1,k);
            
            p4_x = RPF_x(i,j+1,k);
            p4_y = RPF_y(i,j+1,k);
            p4_z = RPF_z(i,j+1,k);
            
            p5_x = RPF_x(i,j,k+1);
            p5_y = RPF_y(i,j,k+1);
            p5_z = RPF_z(i,j,k+1);
            
            p6_x = RPF_x(i+1,j,k+1);
            p6_y = RPF_y(i+1,j,k+1);
            p6_z = RPF_z(i+1,j,k+1);
            
            p7_x = RPF_x(i+1,j+1,k+1);
            p7_y = RPF_y(i+1,j+1,k+1);
            p7_z = RPF_z(i+1,j+1,k+1);
            
            p8_x = RPF_x(i,j+1,k+1);
            p8_y = RPF_y(i,j+1,k+1);
            p8_z = RPF_z(i,j+1,k+1);
            
            a=[p1_x, p2_x, p3_x, p4_x, p5_x, p6_x, p7_x, p8_x];
            b=[p1_y, p2_y, p3_y, p4_y, p5_y, p6_y, p7_y, p8_y];
            c=[p1_z, p2_z, p3_z, p4_z, p5_z, p6_z, p7_z, p8_z];
            
            if isempty(find(a>0))==0 && isempty(find(a<=0))==0
                sign_x=1;
            end
            
            if isempty(find(b>0))==0 && isempty(find(b<=0))==0
                sign_y=1;
            end
            
            if isempty(find(c>0))==0 && isempty(find(c<=0))==0
                sign_z=1;
            end
            
            
            
            if sign_x==1 && sign_y==1 && sign_z==1 
                
                p12_x=(p1_x+p2_x)/2;
                p23_x=(p2_x+p3_x)/2;
                p34_x=(p3_x+p4_x)/2;
                p41_x=(p1_x+p4_x)/2;
                p56_x=(p5_x+p6_x)/2;
                p67_x=(p6_x+p7_x)/2;
                p78_x=(p7_x+p8_x)/2;
                p85_x=(p5_x+p8_x)/2;
                p15_x=(p1_x+p5_x)/2;
                p26_x=(p2_x+p6_x)/2;
                p37_x=(p3_x+p7_x)/2;
                p48_x=(p4_x+p8_x)/2;
                p_down_x=(p1_x+p2_x+p3_x+p4_x)/4;
                p_up_x=(p5_x+p6_x+p7_x+p8_x)/4;
                p_left_x=(p1_x+p4_x+p5_x+p8_x)/4;
                p_right_x=(p2_x+p3_x+p6_x+p7_x)/4;
                p_fore_x=(p3_x+p4_x+p7_x+p8_x)/4;
                p_rear_x=(p1_x+p2_x+p5_x+p6_x)/4;
                p_center_x=(p1_x+p2_x+p3_x+p4_x+p5_x+p6_x+p7_x+p8_x)/8;
                
                sub1_x=[p1_x, p12_x, p15_x, p41_x, p_rear_x, p_left_x, p_down_x, p_center_x];
                sub2_x=[p2_x, p12_x, p26_x, p23_x, p_rear_x, p_right_x, p_down_x, p_center_x];
                sub3_x=[p3_x, p34_x, p37_x, p23_x, p_fore_x, p_right_x, p_down_x, p_center_x];
                sub4_x=[p4_x, p34_x, p48_x, p41_x, p_fore_x, p_left_x, p_down_x, p_center_x];
                sub5_x=[p5_x, p56_x, p15_x, p85_x, p_rear_x, p_left_x, p_up_x, p_center_x];
                sub6_x=[p6_x, p56_x, p26_x, p67_x, p_rear_x, p_right_x, p_up_x, p_center_x];
                sub7_x=[p7_x, p78_x, p37_x, p67_x, p_fore_x, p_right_x, p_up_x, p_center_x];
                sub8_x=[p8_x, p78_x, p48_x, p85_x, p_fore_x, p_left_x, p_up_x, p_center_x];
                
                
                
                p12_y=(p1_y+p2_y)/2;
                p23_y=(p2_y+p3_y)/2;
                p34_y=(p3_y+p4_y)/2;
                p41_y=(p1_y+p4_y)/2;
                p56_y=(p5_y+p6_y)/2;
                p67_y=(p6_y+p7_y)/2;
                p78_y=(p7_y+p8_y)/2;
                p85_y=(p5_y+p8_y)/2;
                p15_y=(p1_y+p5_y)/2;
                p26_y=(p2_y+p6_y)/2;
                p37_y=(p3_y+p7_y)/2;
                p48_y=(p4_y+p8_y)/2;             
                p_down_y=(p1_y+p2_y+p3_y+p4_y)/4;
                p_up_y=(p5_y+p6_y+p7_y+p8_y)/4;
                p_left_y=(p1_y+p4_y+p5_y+p8_y)/4;
                p_right_y=(p2_y+p3_y+p6_y+p7_y)/4;
                p_fore_y=(p3_y+p4_y+p7_y+p8_y)/4;
                p_rear_y=(p1_y+p2_y+p5_y+p6_y)/4;
                p_center_y=(p1_y+p2_y+p3_y+p4_y+p5_y+p6_y+p7_y+p8_y)/8;
                
                sub1_y=[p1_y, p12_y, p15_y, p41_y, p_rear_y, p_left_y, p_down_y, p_center_y];
                sub2_y=[p2_y, p12_y, p26_y, p23_y, p_rear_y, p_right_y, p_down_y, p_center_y];
                sub3_y=[p3_y, p34_y, p37_y, p23_y, p_fore_y, p_right_y, p_down_y, p_center_y];
                sub4_y=[p4_y, p34_y, p48_y, p41_y, p_fore_y, p_left_y, p_down_y, p_center_y];
                sub5_y=[p5_y, p56_y, p15_y, p85_y, p_rear_y, p_left_y, p_up_y, p_center_y];
                sub6_y=[p6_y, p56_y, p26_y, p67_y, p_rear_y, p_right_y, p_up_y, p_center_y];
                sub7_y=[p7_y, p78_y, p37_y, p67_y, p_fore_y, p_right_y, p_up_y, p_center_y];
                sub8_y=[p8_y, p78_y, p48_y, p85_y, p_fore_y, p_left_y, p_up_y, p_center_y];
                
                
                
                p12_z=(p1_z+p2_z)/2;
                p23_z=(p2_z+p3_z)/2;
                p34_z=(p3_z+p4_z)/2;
                p41_z=(p1_z+p4_z)/2;
                p56_z=(p5_z+p6_z)/2;
                p67_z=(p6_z+p7_z)/2;
                p78_z=(p7_z+p8_z)/2;
                p85_z=(p5_z+p8_z)/2;
                p15_z=(p1_z+p5_z)/2;
                p26_z=(p2_z+p6_z)/2;
                p37_z=(p3_z+p7_z)/2;
                p48_z=(p4_z+p8_z)/2;                
                p_down_z=(p1_z+p2_z+p3_z+p4_z)/4;
                p_up_z=(p5_z+p6_z+p7_z+p8_z)/4;
                p_left_z=(p1_z+p4_z+p5_z+p8_z)/4;
                p_right_z=(p2_z+p3_z+p6_z+p7_z)/4;
                p_fore_z=(p3_z+p4_z+p7_z+p8_z)/4;
                p_rear_z=(p1_z+p2_z+p5_z+p6_z)/4;
                p_center_z=(p1_z+p2_z+p3_z+p4_z+p5_z+p6_z+p7_z+p8_z)/8;
                
                sub1_z=[p1_z, p12_z, p15_z, p41_z, p_rear_z, p_left_z, p_down_z, p_center_z];
                sub2_z=[p2_z, p12_z, p26_z, p23_z, p_rear_z, p_right_z, p_down_z, p_center_z];
                sub3_z=[p3_z, p34_z, p37_z, p23_z, p_fore_z, p_right_z, p_down_z, p_center_z];
                sub4_z=[p4_z, p34_z, p48_z, p41_z, p_fore_z, p_left_z, p_down_z, p_center_z];
                sub5_z=[p5_z, p56_z, p15_z, p85_z, p_rear_z, p_left_z, p_up_z, p_center_z];
                sub6_z=[p6_z, p56_z, p26_z, p67_z, p_rear_z, p_right_z, p_up_z, p_center_z];
                sub7_z=[p7_z, p78_z, p37_z, p67_z, p_fore_z, p_right_z, p_up_z, p_center_z];
                sub8_z=[p8_z, p78_z, p48_z, p85_z, p_fore_z, p_left_z, p_up_z, p_center_z];
                
                
                sx1=0;
                sy1=0;
                sz1=0;
                sx2=0;
                sy2=0;
                sz2=0;
                sx3=0;
                sy3=0;
                sz3=0;
                sx4=0;
                sy4=0;
                sz4=0;
                sx5=0;
                sy5=0;
                sz5=0;
                sx6=0;
                sy6=0;
                sz6=0;
                sx7=0;
                sy7=0;
                sz7=0;
                sx8=0;
                sy8=0;
                sz8=0;
                subcell1=0;
                subcell2=0;
                subcell3=0;
                subcell4=0;
                subcell5=0;
                subcell6=0;
                subcell7=0;
                subcell8=0;
                
                
                if isempty(find(sub1_x>0))==0 && isempty(find(sub1_x<=0))==0
                   sx1=1;
                end
                if isempty(find(sub1_y>0))==0 && isempty(find(sub1_y<=0))==0
                   sy1=1;
                end               
                if isempty(find(sub1_z>0))==0 && isempty(find(sub1_z<=0))==0
                   sz1=1;
                end
                
                if sx1==1 && sy1==1 && sz1==1
                   start_point=[i+1/4, j+1/4, k+1/4];
                   xmin=i;
                   ymin=j; 
                   zmin=k;
                   xmax=i+1/2; 
                   ymax=j+1/2;
                   zmax=k+1/2;
                   [coordinate,inside] = skel_Newton(start_point, boundary_pos, boundary_num, xmin, ymin, zmin, xmax, ymax, zmax, RPF_x, RPF_y, RPF_z, RPForder, NewtonTH);
                   if inside==1
                       critical_pos=coordinate;
                       subcell1=1;
                   end
                end
                
                
                
                if isempty(find(sub2_x>0))==0 && isempty(find(sub2_x<=0))==0
                   sx2=1;
                end                
                if isempty(find(sub2_y>0))==0 && isempty(find(sub2_y<=0))==0
                   sy2=1;
                end                
                if isempty(find(sub2_z>0))==0 && isempty(find(sub2_z<=0))==0
                   sz2=1;
                end
                
                if sx2==1 && sy2==1 && sz2==1
                   start_point=[i+3/4, j+1/4, k+1/4];
                   xmin=i+1/2;
                   ymin=j; 
                   zmin=k;
                   xmax=i+1; 
                   ymax=j+1/2;
                   zmax=k+1/2;
                   [coordinate,inside] = skel_Newton(start_point, boundary_pos, boundary_num, xmin, ymin, zmin, xmax, ymax, zmax, RPF_x, RPF_y, RPF_z, RPForder, NewtonTH);
                   if inside==1
                       critical_pos=coordinate;
                       subcell2=1;
                   end
                end
                
                
                
                if isempty(find(sub3_x>0))==0 && isempty(find(sub3_x<=0))==0
                   sx3=1;
                end               
                if isempty(find(sub3_y>0))==0 && isempty(find(sub3_y<=0))==0
                   sy3=1;
                end                
                if isempty(find(sub3_z>0))==0 && isempty(find(sub3_z<=0))==0
                   sz3=1;
                end
                
                if sx3==1 && sy3==1 && sz3==1                  
                   start_point=[i+3/4, j+3/4, k+1/4];
                   xmin=i+1/2;
                   ymin=j+1/2; 
                   zmin=k;
                   xmax=i+1; 
                   ymax=j+1;
                   zmax=k+1/2;
                   [coordinate,inside] = skel_Newton(start_point, boundary_pos, boundary_num, xmin, ymin, zmin, xmax, ymax, zmax, RPF_x, RPF_y, RPF_z, RPForder, NewtonTH);
                   if inside==1
                       critical_pos=coordinate;
                       subcell3=1;
                   end
                end
                
                
                
                if isempty(find(sub4_x>0))==0 && isempty(find(sub4_x<=0))==0
                   sx4=1;
                end                
                if isempty(find(sub4_y>0))==0 && isempty(find(sub4_y<=0))==0
                   sy4=1;
                end                
                if isempty(find(sub4_z>0))==0 && isempty(find(sub4_z<=0))==0
                   sz4=1;
                end
                
                if sx4==1 && sy4==1 && sz4==1                   
                   start_point=[i+1/4, j+3/4, k+1/4];
                   xmin=i;
                   ymin=j+1/2; 
                   zmin=k;
                   xmax=i+1/2; 
                   ymax=j+1;
                   zmax=k+1/2;
                   [coordinate,inside] = skel_Newton(start_point, boundary_pos, boundary_num, xmin, ymin, zmin, xmax, ymax, zmax, RPF_x, RPF_y, RPF_z, RPForder, NewtonTH);
                   if inside==1
                       critical_pos=coordinate;
                       subcell4=1;
                   end
                end
                
                
                
                if isempty(find(sub5_x>0))==0 && isempty(find(sub5_x<=0))==0
                   sx5=1;
                end                
                if isempty(find(sub5_y>0))==0 && isempty(find(sub5_y<=0))==0
                   sy5=1;
                end                
                if isempty(find(sub5_z>0))==0 && isempty(find(sub5_z<=0))==0
                   sz5=1;
                end
                
                if sx5==1 && sy5==1 && sz5==1
                   start_point=[i+1/4, j+1/4, k+3/4];
                   xmin=i;
                   ymin=j; 
                   zmin=k+1/2;
                   xmax=i+1/2; 
                   ymax=j+1/2;
                   zmax=k+1;
                   [coordinate,inside] = skel_Newton(start_point, boundary_pos, boundary_num, xmin, ymin, zmin, xmax, ymax, zmax, RPF_x, RPF_y, RPF_z, RPForder, NewtonTH);
                   if inside==1
                       critical_pos=coordinate;
                       subcell5=1;
                   end
                end
                
                
                
                if isempty(find(sub6_x>0))==0 && isempty(find(sub6_x<=0))==0
                   sx6=1;
                end                
                if isempty(find(sub6_y>0))==0 && isempty(find(sub6_y<=0))==0
                   sy6=1;
                end                
                if isempty(find(sub6_z>0))==0 && isempty(find(sub6_z<=0))==0
                   sz6=1;
                end
                
                if sx6==1 && sy6==1 && sz6==1
                   start_point=[i+3/4, j+1/4, k+3/4];
                   xmin=i+1/2;
                   ymin=j; 
                   zmin=k+1/2;
                   xmax=i+1; 
                   ymax=j+1/2;
                   zmax=k+1;
                   [coordinate,inside] = skel_Newton(start_point, boundary_pos, boundary_num, xmin, ymin, zmin, xmax, ymax, zmax, RPF_x, RPF_y, RPF_z, RPForder, NewtonTH);
                   if inside==1
                       critical_pos=coordinate;
                       subcell6=1;
                   end
                end
                
                
                
                if isempty(find(sub7_x>0))==0 && isempty(find(sub7_x<=0))==0
                   sx7=1;
                end              
                if isempty(find(sub7_y>0))==0 && isempty(find(sub7_y<=0))==0
                   sy7=1;
                end                
                if isempty(find(sub7_z>0))==0 && isempty(find(sub7_z<=0))==0
                   sz7=1;
                end
                
                if sx7==1 && sy7==1 && sz7==1
                   start_point=[i+3/4, j+3/4, k+3/4];
                   xmin=i+1/2;
                   ymin=j+1/2; 
                   zmin=k+1/2;
                   xmax=i+1; 
                   ymax=j+1;
                   zmax=k+1;
                   [coordinate,inside] = skel_Newton(start_point, boundary_pos, boundary_num, xmin, ymin, zmin, xmax, ymax, zmax, RPF_x, RPF_y, RPF_z, RPForder, NewtonTH);
                   if inside==1
                       critical_pos=coordinate;
                       subcell7=1;
                   end
                end
                
                
                
                if isempty(find(sub8_x>0))==0 && isempty(find(sub8_x<=0))==0
                   sx8=1;
                end                
                if isempty(find(sub8_y>0))==0 && isempty(find(sub8_y<=0))==0
                   sy8=1;
                end                
                if isempty(find(sub8_z>0))==0 && isempty(find(sub8_z<=0))==0
                   sz8=1;
                end
                
                if sx8==1 && sy8==1 && sz8==1
                   start_point=[i+1/4, j+3/4, k+3/4];
                   xmin=i;
                   ymin=j+1/2; 
                   zmin=k+1/2;
                   xmax=i+1/2; 
                   ymax=j+1;
                   zmax=k+1;
                   [coordinate,inside] = skel_Newton(start_point, boundary_pos, boundary_num, xmin, ymin, zmin, xmax, ymax, zmax, RPF_x, RPF_y, RPF_z, RPForder, NewtonTH);
                   if inside==1
                       critical_pos=coordinate;
                       subcell8=1;
                   end
                end
                
                temp=[subcell1,subcell2,subcell3,subcell4,subcell5,subcell6,subcell7,subcell8];
                num=size(find(temp==1),2);
                
                if num==1
                    critical=critical_pos;
                    permit=1;
                end
                
                if num>1
                   start_point=[i+1/2, j+1/2, k+1/2];
                   xmin=i;
                   ymin=j; 
                   zmin=k;
                   xmax=i+1; 
                   ymax=j+1;
                   zmax=k+1;
                   [coordinate,inside] = skel_Newton(start_point, boundary_pos, boundary_num, xmin, ymin, zmin, xmax, ymax, zmax, RPF_x, RPF_y, RPF_z, RPForder, NewtonTH);
                   if inside==1
                       critical_pos=coordinate;
                   end
                   critical=critical_pos;
                   permit=1;
                end
                
            end
                