function bbox = createOverlapBox(full_stack_shape_ijk, chunk_shape_ijk, double_pad_depth_plus_one)
    if nargin<1
        full_stack_shape_ijk = [2000 2000 2000];
        chunk_shape_ijk = [200 200 200];
        double_pad_depth_plus_one = 15;
    end
    numDims = length(full_stack_shape_ijk) ;
    if length(double_pad_depth_plus_one)<numDims ,
        double_pad_depth_plus_one = ones(1,numDims) * double_pad_depth_plus_one ;
    end
    pad_depth_ijk = (double_pad_depth_plus_one-1)/2 ;
    [BB{1:numDims}]=deal([]);
    for i=1:numDims
        st = 1:(chunk_shape_ijk(i)-double_pad_depth_plus_one(i)):full_stack_shape_ijk(i)-pad_depth_ijk(i)-1;
        en = min(full_stack_shape_ijk(i),st+(chunk_shape_ijk(i)-1));
        BB{i} = [st(:) en(:)];
    end
    numb = cellfun(@(x) size(x,1),BB);
    if length(numb)==2
        numb = [numb 1];
    end
    [aa,bb,cc] = ndgrid(1:numb(1),1:numb(2),1:numb(3));
    perms = [aa(:),bb(:),cc(:)];
    numblocks = size(perms,1);

    bbox = zeros(numblocks,2*numDims);
    for i=1:numblocks
        bbox(i,:) = [BB{1}(perms(i,1),:) BB{2}(perms(i,2),:) BB{3}(perms(i,3),:)];
    end
end
