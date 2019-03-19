function [XXX1, YYY1, ZZZ1, XXX2a, YYY2a, ZZZ2a, XXX2b, YYY2b, ZZZ2b, full_skeleton] = skel_tracefield(voxel, object, boundary_pos, num_boundary, saddle1, saddle2, attr_num, saddle1_num, saddle2_num, pathline1, pathline2a, pathline2b, RPF_order)

step = input('Please set the step of particle tracing in Vector Fields: ([]:default=0.25) ');
if isempty(step)
    step=0.25;
end

maxtime = input('Please set the maximum time of particle tracing in Vector Fields: ([]:default=200) ');
if isempty(maxtime)
    maxtime=200;
end

XXX1=[];
YYY1=[];
ZZZ1=[];
XXX2a=[];
YYY2a=[];
ZZZ2a=[];
XXX2b=[];
YYY2b=[];
ZZZ2b=[];

full_skeleton=zeros(size(voxel));

for time=1:saddle1_num
    px=saddle1(time,1);
    py=saddle1(time,2);
    pz=saddle1(time,3);
    
    px=px-0.5*pathline1(time,1);
    py=py-0.5*pathline1(time,2);
    pz=pz-0.5*pathline1(time,3);
    
    number=0;
    
    sample_x=floor(px+0.5);
    sample_y=floor(py+0.5);
    sample_z=floor(pz+0.5);
    
    if sample_x-1<=0
        sample_x=2;
    end    
    if sample_y-1<=0
        sample_y=2;
    end
    if sample_z-1<=0
        sample_z=2;
    end
    if sample_x+1>size(voxel,1)
        sample_x=size(voxel,1)-1;
    end
    if sample_y+1>size(voxel,2)
        sample_y=size(voxel,2)-1;
    end
    if sample_z+1>size(voxel,3)
        sample_z=size(voxel,3)-1;
    end
    
    full_skeleton(sample_x, sample_y, sample_z)=1;
    box=object(sample_x-1:sample_x+1, sample_y-1:sample_y+1, sample_z-1:sample_z+1);
    
    if  sum(sum(sum(box)))<=9
        XXX1(1:maxtime,time)=px;
        YYY1(1:maxtime,time)=py;
        ZZZ1(1:maxtime,time)=pz;
    end
    
    while sum(sum(sum(box)))>9
        number = number + 1;
        point=[px, py, pz];
        [rx, ry, rz]=skel_computeRPF(point, boundary_pos, num_boundary, RPF_order); 
        AE=sqrt(rx^2+ry^2+rz^2);
        px=px+step*(rx/AE);
        py=py+step*(ry/AE);
        pz=pz+step*(rz/AE);
        
        if px<=size(voxel,1) && py<=size(voxel,2) && pz<=size(voxel,3) && px>=1 && py>=1 && pz>=1
            if object(floor(px+0.5),floor(py+0.5),floor(pz+0.5))==0
                XXX1(number:maxtime,time)=px;
                YYY1(number:maxtime,time)=py;
                ZZZ1(number:maxtime,time)=pz;
                break
            end
        end
        
        if px>size(voxel,1) || py>size(voxel,2) || pz>size(voxel,3) || px<=0 || py<=0 || pz<=0
            XXX1(number:maxtime,time)=px;
            YYY1(number:maxtime,time)=py;
            ZZZ1(number:maxtime,time)=pz;
            break
        end
        
        if number>maxtime
            break
        end        
        
        XXX1(number,time)=px;
        YYY1(number,time)=py;
        ZZZ1(number,time)=pz;
        
        sample_x=floor(px+0.5);
        sample_y=floor(py+0.5);
        sample_z=floor(pz+0.5);
        
         if sample_x-1<=0
             sample_x=2;
         end    
         if sample_y-1<=0
             sample_y=2;
         end
         if sample_z-1<=0
             sample_z=2;
         end
         if sample_x+1>size(voxel,1)
             sample_x=size(voxel,1)-1;
         end
         if sample_y+1>size(voxel,2)
             sample_y=size(voxel,2)-1;
         end
         if sample_z+1>size(voxel,3)
             sample_z=size(voxel,3)-1;
         end
         
        full_skeleton(sample_x, sample_y, sample_z)=1;
        box=object(sample_x-1:sample_x+1, sample_y-1:sample_y+1, sample_z-1:sample_z+1);
        
        if  sum(sum(sum(box)))<=9
            XXX1(number:maxtime,time)=px;
            YYY1(number:maxtime,time)=py;
            ZZZ1(number:maxtime,time)=pz;
        end
        
    end
end


for time=1:saddle1_num
    px=saddle1(time,1);
    py=saddle1(time,2);
    pz=saddle1(time,3);
    
    px=px+0.5*pathline1(time,1);
    py=py+0.5*pathline1(time,2);
    pz=pz+0.5*pathline1(time,3);
    
    number=0;
    
    sample_x=floor(px+0.5);
    sample_y=floor(py+0.5);
    sample_z=floor(pz+0.5);
    
    
         if sample_x-1<=0
             sample_x=2;
         end    
         if sample_y-1<=0
             sample_y=2;
         end
         if sample_z-1<=0
             sample_z=2;
         end
         if sample_x+1>size(voxel,1)
             sample_x=size(voxel,1)-1;
         end
         if sample_y+1>size(voxel,2)
             sample_y=size(voxel,2)-1;
         end
         if sample_z+1>size(voxel,3)
             sample_z=size(voxel,3)-1;
         end
    
    full_skeleton(sample_x, sample_y, sample_z)=1;
    box=object(sample_x-1:sample_x+1, sample_y-1:sample_y+1, sample_z-1:sample_z+1);
    
    if  sum(sum(sum(box)))<=9
        XXX1(1:maxtime,time+saddle1_num)=px;
        YYY1(1:maxtime,time+saddle1_num)=py;
        ZZZ1(1:maxtime,time+saddle1_num)=pz;
    end

    while sum(sum(sum(box)))>9
        number = number + 1;
        point=[px, py, pz];
        [rx, ry, rz]=skel_computeRPF(point, boundary_pos, num_boundary, RPF_order); 
        AE=sqrt(rx^2+ry^2+rz^2);
        px=px+step*(rx/AE);
        py=py+step*(ry/AE);
        pz=pz+step*(rz/AE);
        
        if px<=size(voxel,1) && py<=size(voxel,2) && pz<=size(voxel,3) && px>=1 && py>=1 && pz>=1
            if object(floor(px+0.5),floor(py+0.5),floor(pz+0.5))==0
                XXX1(number:maxtime,time+saddle1_num)=px;
                YYY1(number:maxtime,time+saddle1_num)=py;
                ZZZ1(number:maxtime,time+saddle1_num)=pz;
                break
            end
        end
        
        if px>size(voxel,1) || py>size(voxel,2) || pz>size(voxel,3) || px<=0 || py<=0 || pz<=0
            XXX1(number:maxtime,time+saddle1_num)=px;
            YYY1(number:maxtime,time+saddle1_num)=py;
            ZZZ1(number:maxtime,time+saddle1_num)=pz;
            break
        end
        
        if number>maxtime
            break
        end  
                
        XXX1(number,time+saddle1_num)=px;
        YYY1(number,time+saddle1_num)=py;
        ZZZ1(number,time+saddle1_num)=pz;
        
        sample_x=floor(px+0.5);
        sample_y=floor(py+0.5);
        sample_z=floor(pz+0.5);
        
        if sample_x-1<=0
             sample_x=2;
         end    
         if sample_y-1<=0
             sample_y=2;
         end
         if sample_z-1<=0
             sample_z=2;
         end
         if sample_x+1>size(voxel,1)
             sample_x=size(voxel,1)-1;
         end
         if sample_y+1>size(voxel,2)
             sample_y=size(voxel,2)-1;
         end
         if sample_z+1>size(voxel,3)
             sample_z=size(voxel,3)-1;
         end
        
        full_skeleton(sample_x, sample_y, sample_z)=1;
        box=object(sample_x-1:sample_x+1, sample_y-1:sample_y+1, sample_z-1:sample_z+1);
        
        if  sum(sum(sum(box)))<=9
            XXX1(number:maxtime,time+saddle1_num)=px;
            YYY1(number:maxtime,time+saddle1_num)=py;
            ZZZ1(number:maxtime,time+saddle1_num)=pz;
        end
        
    end
end




for time=1:saddle2_num
    px=saddle2(time,1);
    py=saddle2(time,2);
    pz=saddle2(time,3);
    
    px=px+0.5*real(pathline2a(time,1));
    py=py+0.5*real(pathline2a(time,2));
    pz=pz+0.5*real(pathline2a(time,3));
        
    number=0;
    
    sample_x=floor(px+0.5);
    sample_y=floor(py+0.5);
    sample_z=floor(pz+0.5);
    
    if sample_x-1<=0
        sample_x=2;
    end    
    if sample_y-1<=0
        sample_y=2;
    end
    if sample_z-1<=0
        sample_z=2;
    end
    if sample_x+1>size(voxel,1)
        sample_x=size(voxel,1)-1;
    end
    if sample_y+1>size(voxel,2)
        sample_y=size(voxel,2)-1;
    end
    if sample_z+1>size(voxel,3)
        sample_z=size(voxel,3)-1;
    end
    
    full_skeleton(sample_x, sample_y, sample_z)=1;
    
    if  object(sample_x, sample_y, sample_z)==0
        XXX2a(1:maxtime,time)=px;
        YYY2a(1:maxtime,time)=py;
        ZZZ2a(1:maxtime,time)=pz;
    end  
        
    
    while object(sample_x, sample_y, sample_z)==1
        number=number+1;
        point=[px, py, pz];
        [rx, ry, rz]=skel_computeRPF(point, boundary_pos, num_boundary, RPF_order); 
        AE=sqrt(rx^2+ry^2+rz^2);
        px=px+step*(rx/AE);
        py=py+step*(ry/AE);
        pz=pz+step*(rz/AE);
        
        if px<=size(voxel,1) && py<=size(voxel,2) && pz<=size(voxel,3) && px>=1 && py>=1 && pz>=1
            if object(floor(px+0.5),floor(py+0.5),floor(pz+0.5))==0
                XXX2a(number:maxtime,time)=px;
                YYY2a(number:maxtime,time)=py;
                ZZZ2a(number:maxtime,time)=pz;
                break
            end
        end
        
        if px>size(voxel,1) || py>size(voxel,2) || pz>size(voxel,3) || px<=0 || py<=0 || pz<=0
            XXX2a(number:maxtime,time)=px;
            YYY2a(number:maxtime,time)=py;
            ZZZ2a(number:maxtime,time)=pz;
            break
        end
        
        if number>maxtime
            break
        end
                
        XXX2a(number,time)=px;
        YYY2a(number,time)=py;
        ZZZ2a(number,time)=pz;
        
        sample_x=floor(px+0.5);
        sample_y=floor(py+0.5);
        sample_z=floor(pz+0.5);
        
        if sample_x-1<=0
            sample_x=2;
        end    
        if sample_y-1<=0
            sample_y=2;
        end
        if sample_z-1<=0
            sample_z=2;
        end
        if sample_x+1>size(voxel,1)
            sample_x=size(voxel,1)-1;
        end
        if sample_y+1>size(voxel,2)
            sample_y=size(voxel,2)-1;
        end
        if sample_z+1>size(voxel,3)
            sample_z=size(voxel,3)-1;
        end
        
        full_skeleton(sample_x, sample_y, sample_z)=1;
        
        if  object(sample_x, sample_y, sample_z)==0
            XXX2a(number:maxtime,time)=px;
            YYY2a(number:maxtime,time)=py;
            ZZZ2a(number:maxtime,time)=pz;
        end          
    end
    
end



for time=1:saddle2_num
    px=saddle2(time,1);
    py=saddle2(time,2);
    pz=saddle2(time,3);
    
    px=px-0.5*real(pathline2a(time,1));
    py=py-0.5*real(pathline2a(time,2));
    pz=pz-0.5*real(pathline2a(time,3));
    
    number=0;
    
    sample_x=floor(px+0.5);
    sample_y=floor(py+0.5);
    sample_z=floor(pz+0.5);
    
    if sample_x-1<=0
        sample_x=2;
    end    
    if sample_y-1<=0
        sample_y=2;
    end
    if sample_z-1<=0
        sample_z=2;
    end
    if sample_x+1>size(voxel,1)
        sample_x=size(voxel,1)-1;
    end
    if sample_y+1>size(voxel,2)
        sample_y=size(voxel,2)-1;
    end
    if sample_z+1>size(voxel,3)
        sample_z=size(voxel,3)-1;
    end
    
    full_skeleton(sample_x, sample_y, sample_z)=1;
    
    if  object(sample_x, sample_y, sample_z)==0
        XXX2a(1:maxtime,time+saddle2_num)=px;
        YYY2a(1:maxtime,time+saddle2_num)=py;
        ZZZ2a(1:maxtime,time+saddle2_num)=pz;
    end
    
    while object(sample_x, sample_y, sample_z)==1
        number=number+1;
        point=[px, py, pz];
        [rx, ry, rz]=skel_computeRPF(point, boundary_pos, num_boundary, RPF_order); 
        AE=sqrt(rx^2+ry^2+rz^2);
        px=px+step*(rx/AE);
        py=py+step*(ry/AE);
        pz=pz+step*(rz/AE);
        
        if px<=size(voxel,1) && py<=size(voxel,2) && pz<=size(voxel,3) && px>=1 && py>=1 && pz>=1
            if object(floor(px+0.5),floor(py+0.5),floor(pz+0.5))==0
                XXX2a(number:maxtime,time+saddle2_num)=px;
                YYY2a(number:maxtime,time+saddle2_num)=py;
                ZZZ2a(number:maxtime,time+saddle2_num)=pz;
                break
            end
        end
        
        if px>size(voxel,1) || py>size(voxel,2) || pz>size(voxel,3) || px<=0 || py<=0 || pz<=0
            XXX2a(number:maxtime,time+saddle2_num)=px;
            YYY2a(number:maxtime,time+saddle2_num)=py;
            ZZZ2a(number:maxtime,time+saddle2_num)=pz;
            break
        end
        
        if number>maxtime
            break
        end
           
        XXX2a(number,time+saddle2_num)=px;
        YYY2a(number,time+saddle2_num)=py;
        ZZZ2a(number,time+saddle2_num)=pz;
        
        sample_x=floor(px+0.5);
        sample_y=floor(py+0.5);
        sample_z=floor(pz+0.5);
        
        if sample_x-1<=0
             sample_x=2;
        end    
        if sample_y-1<=0
             sample_y=2;
        end
        if sample_z-1<=0
             sample_z=2;
        end
        if sample_x+1>size(voxel,1)
             sample_x=size(voxel,1)-1;
        end
        if sample_y+1>size(voxel,2)
             sample_y=size(voxel,2)-1;
        end
        if sample_z+1>size(voxel,3)
             sample_z=size(voxel,3)-1;
        end
        
        full_skeleton(sample_x, sample_y, sample_z)=1;
        
        if  object(sample_x, sample_y, sample_z)==0
            XXX2a(number:maxtime,time+saddle2_num)=px;
            YYY2a(number:maxtime,time+saddle2_num)=py;
            ZZZ2a(number:maxtime,time+saddle2_num)=pz;
        end         
    end
    
end


for time=1:saddle2_num
    px=saddle2(time,1);
    py=saddle2(time,2);
    pz=saddle2(time,3);
    
    px=px+0.5*real(pathline2b(time,1));
    py=py+0.5*real(pathline2b(time,2));
    pz=pz+0.5*real(pathline2b(time,3));   
    
    number=0;
    
    sample_x=floor(px+0.5);
    sample_y=floor(py+0.5);
    sample_z=floor(pz+0.5);
    
    if sample_x-1<=0
        sample_x=2;
    end    
    if sample_y-1<=0
        sample_y=2;
    end
    if sample_z-1<=0
        sample_z=2;
    end
    if sample_x+1>size(voxel,1)
        sample_x=size(voxel,1)-1;
    end
    if sample_y+1>size(voxel,2)
        sample_y=size(voxel,2)-1;
    end
    if sample_z+1>size(voxel,3)
        sample_z=size(voxel,3)-1;
    end
    
    full_skeleton(sample_x, sample_y, sample_z)=1;
    
    if  object(sample_x, sample_y, sample_z)==0
        XXX2b(1:maxtime,time)=px;
        YYY2b(1:maxtime,time)=py;
        ZZZ2b(1:maxtime,time)=pz;
    end
     
    while object(sample_x, sample_y, sample_z)==1
        number=number+1;
        point=[px, py, pz];
        [rx, ry, rz]=skel_computeRPF(point, boundary_pos, num_boundary, RPF_order); 
        AE=sqrt(rx^2+ry^2+rz^2);
        px=px+step*(rx/AE);
        py=py+step*(ry/AE);
        pz=pz+step*(rz/AE);
        
        if px<=size(voxel,1) && py<=size(voxel,2) && pz<=size(voxel,3) && px>=1 && py>=1 && pz>=1
            if object(floor(px+0.5),floor(py+0.5),floor(pz+0.5))==0
                XXX2b(number:maxtime,time)=px;
                YYY2b(number:maxtime,time)=py;
                ZZZ2b(number:maxtime,time)=pz;
                break
            end
        end
        
        if px>size(voxel,1) || py>size(voxel,2) || pz>size(voxel,3) || px<=0 || py<=0 || pz<=0
            XXX2b(number:maxtime,time)=px;
            YYY2b(number:maxtime,time)=py;
            ZZZ2b(number:maxtime,time)=pz;
            break
        end
        
        if number>maxtime
            break
        end 
        
        XXX2b(number,time)=px;
        YYY2b(number,time)=py;
        ZZZ2b(number,time)=pz;
        
        sample_x=floor(px+0.5);
        sample_y=floor(py+0.5);
        sample_z=floor(pz+0.5);
        
        if sample_x-1<=0
             sample_x=2;
        end    
        if sample_y-1<=0
             sample_y=2;
        end
        if sample_z-1<=0
             sample_z=2;
        end
        if sample_x+1>size(voxel,1)
             sample_x=size(voxel,1)-1;
        end
        if sample_y+1>size(voxel,2)
             sample_y=size(voxel,2)-1;
        end
        if sample_z+1>size(voxel,3)
             sample_z=size(voxel,3)-1;
        end
        
        full_skeleton(sample_x, sample_y, sample_z)=1;
        
        if object(sample_x, sample_y, sample_z)==0
            XXX2b(number:maxtime,time)=px;
            YYY2b(number:maxtime,time)=py;
            ZZZ2b(number:maxtime,time)=pz;
        end      
    end
    
end



for time=1:saddle2_num
    px=saddle2(time,1);
    py=saddle2(time,2);
    pz=saddle2(time,3);
    
    px=px-0.5*real(pathline2b(time,1));
    py=py-0.5*real(pathline2b(time,2));
    pz=pz-0.5*real(pathline2b(time,3));
       
    number=0;
    sample_x=floor(px+0.5);
    sample_y=floor(py+0.5);
    sample_z=floor(pz+0.5);
    
    if sample_x-1<=0
        sample_x=2;
    end    
    if sample_y-1<=0
        sample_y=2;
    end
    if sample_z-1<=0
        sample_z=2;
    end
    if sample_x+1>size(voxel,1)
        sample_x=size(voxel,1)-1;
    end
    if sample_y+1>size(voxel,2)
        sample_y=size(voxel,2)-1;
    end
    if sample_z+1>size(voxel,3)
        sample_z=size(voxel,3)-1;
    end
    
    full_skeleton(sample_x, sample_y, sample_z)=1;
    
    if  object(sample_x, sample_y, sample_z)==0
        XXX2b(1:maxtime,time+saddle2_num)=px;
        YYY2b(1:maxtime,time+saddle2_num)=py;
        ZZZ2b(1:maxtime,time+saddle2_num)=pz;
    end
    
    
    while object(sample_x, sample_y, sample_z)==1
        number=number+1;
        point=[px, py, pz];
        [rx, ry, rz]=skel_computeRPF(point, boundary_pos, num_boundary, RPF_order); 
        AE=sqrt(rx^2+ry^2+rz^2);
        px=px+step*(rx/AE);
        py=py+step*(ry/AE);
        pz=pz+step*(rz/AE);
        
        if px<=size(voxel,1) && py<=size(voxel,2) && pz<=size(voxel,3) && px>=1 && py>=1 && pz>=1
            if object(floor(px+0.5),floor(py+0.5),floor(pz+0.5))==0
                XXX2b(number:maxtime,time+saddle2_num)=px;
                YYY2b(number:maxtime,time+saddle2_num)=py;
                ZZZ2b(number:maxtime,time+saddle2_num)=pz;
                break
            end
        end
        
        if px>size(voxel,1) || py>size(voxel,2) || pz>size(voxel,3) || px<=0 || py<=0 || pz<=0
            XXX2b(number:maxtime,time+saddle2_num)=px;
            YYY2b(number:maxtime,time+saddle2_num)=py;
            ZZZ2b(number:maxtime,time+saddle2_num)=pz;
            break
        end
        
        if number>maxtime
            break
        end
                
        XXX2b(number,time+saddle2_num)=px;
        YYY2b(number,time+saddle2_num)=py;
        ZZZ2b(number,time+saddle2_num)=pz;
        
        sample_x=floor(px+0.5);
        sample_y=floor(py+0.5);
        sample_z=floor(pz+0.5);
        
        if sample_x-1<=0
             sample_x=2;
        end    
        if sample_y-1<=0
             sample_y=2;
        end
        if sample_z-1<=0
             sample_z=2;
        end
        if sample_x+1>size(voxel,1)
             sample_x=size(voxel,1)-1;
        end
        if sample_y+1>size(voxel,2)
             sample_y=size(voxel,2)-1;
        end
        if sample_z+1>size(voxel,3)
             sample_z=size(voxel,3)-1;
        end
        
        full_skeleton(sample_x, sample_y, sample_z)=1;
        
        if  object(sample_x, sample_y, sample_z)==0
            XXX2b(number:maxtime,time+saddle2_num)=px;
            YYY2b(number:maxtime,time+saddle2_num)=py;
            ZZZ2b(number:maxtime,time+saddle2_num)=pz;
        end
    end
    
end