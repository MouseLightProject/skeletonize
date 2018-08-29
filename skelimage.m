function skelimage()
%%
% Io = squeeze(h5read(myh5,myh5prob,starts,datasiz));
% function emptyvol = skelimage(Io,BB,brainSize)

% calculated_parameters = '/nrs/mouselight/SAMPLES/2017-10-31/calculated_parameters.jl';
transform_parameters =  '/nrs/mouselight/SAMPLES/2017-09-25/transform.txt';
configfile = 'configfile.cfg';
h5file = '/nrs/mouselight/cluster/classifierOutputs/2017-09-25/20170925_prob0/20170925_prob0_lev-6_chunk-111_111_masked-0.h5';
h5datanama = '/prob0';
%%
params_config = configparser(configfile);
params_config.thr=25;

mkdir(fullfile(params_config.outfolder,'full'))
mkdir(fullfile(params_config.outfolder,'frags'))

params_trans = configparser(transform_parameters);
params_trans.level = params_trans.nl-1;
params_trans.voxres = [params_trans.sx params_trans.sy params_trans.sz]/2^(params_trans.level)/1e3; % in um

% opt.params = params_trans;
% opt.inputh5 = h5file;
% opt.thr = 25;
% opt.sizethreshold = 100;
% opt.outfolder = './skelfulltest'

% %%
% h5datanama = '/prob0';
% try
%     [brainSize,RR,chunk_dims,rank] = h5parser(opt.inputh5,h5datanama);
% catch
%     
% end

%% load BB
%%%%%%%%%%%%%%%%
% read calculated_parameters.jl file
% 0) distribute over lowest level
% 1) read octtree
% 2) read padding
% 3) run 
% grabData(inputfolder,params)

%% skeletonization
W = [100 100 100];
pixloc = um2pix(params_trans,[74167.3, 14663.3, 34743.9])
starts = pixloc - W;
datasiz = 2*W+1;
Io = squeeze(h5read(h5file,h5datanama,starts,datasiz));


[skel,A,subs,edges_] = skeletonimage(Io,params_config);
% figure, imshow3D(Io)
%%
figure, imshow(max(Io,[],3)',[])
hold on
gplot3(A,subs(:,[1 2 3]))

%% reconstruction
G = graph(A);
params_config.params = params_trans;
workflow1(G,subs,params_config)

%%
subs = subs + ones(size(subs,1),1)*BB(1:2:end)-1; % convert to original subs
inds(:,ii) = sub2ind(brainSize,subs(:,1),subs(:,2),subs(:,3)); % convert to original inds

%%
if ~isdeployed& opt.viz
    clear Subs1 Subs2
    [Subs1(:,1),Subs1(:,2),Subs1(:,3)]=ind2sub(brainSize,inds(:,1));
    [Subs2(:,1),Subs2(:,2),Subs2(:,3)]=ind2sub(brainSize,inds(:,2));
    ii = 1:1:size(Subs1,1);
    X = [Subs1(ii,1) Subs2(ii,1) NaN(length(ii),1)]';
    Y = [Subs1(ii,2) Subs2(ii,2) NaN(length(ii),1)]';
    Z = [Subs1(ii,3) Subs2(ii,3) NaN(length(ii),1)]';
    figure,
    plot3(X(:),Y(:),Z(:),'r--')
    axis equal tight
end
% %% connectivity graph
% connG = sparse(edges_(:,1),edges_(:,2),1,max(edges_(:)),max(edges_(:)));

%%
if isdeployed | 1
    %%
    fileID = fopen(outfile,'w');
    if size(inds,2)==2
        fprintf(fileID,'%d %d\n',inds');
    else
        fprintf(fileID,'%d %d %.2f\n',inds');
    end
    fclose(fileID);
end
end
function params = readparams(calculated_parameters)
fid=fopen(calculated_parameters);
clear params
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    strs = strsplit(tline,' ');
    if strcmp(strs{2},'nlevels')
        params.nlevel = eval(strs{4});
    elseif strcmp(strs{2},'shape_leaf_px')
        params.shape_leaf_px = eval(strs{4});
    elseif strcmp(strs{2},'voxelsize_used_um')
        params.voxelsize_used_um = eval(strs{4});
    elseif strcmp(strs{2},'origin_nm')
        params.origin_nm = eval(strs{4});
    else
        
    end
    disp(tline)
end
fclose(fid);
end

function [skel,A,subs_,edges_] = skeletonimage(Io,opt)
% given input image
% Io = squeeze(h5read(myh5,myh5prob,starts+[500 300 100],datasiz));
% figure, imshow(squeeze(max(Io,[],3))',[]),
if ~any(Io(:))
    if isdeployed |1
        %% touch file
        inds = [];
        %         fileID = fopen(outfile,'w');
        %         fclose(fileID);
    end
    return
end
probThr = opt.thr;
fullh = opt.fullh;
% smooth image
Io = smooth3(Io,'gaussian',[3 3 1]);
Io = Io>probThr;
%%
if ~any(Io(:))
    if isdeployed |1
        %% touch file
        inds = [];
        %         fileID = fopen(outfile,'w');
        %         fclose(fileID);
    end
    return
end
%% cleanup image
s  = regionprops(Io, 'centroid','PixelIdxList','Area');
if 0
    % fast
    Iout = zeros(size(Io),'single');
    for ii=1:length(s)
        if s(ii).Area>10
            Iout(s(ii).PixelIdxList)=Io(s(ii).PixelIdxList);
        end
    end
    Io = Iout; clear Iout;
else
    % memory efficient
    for ii=1:length(s)
        if s(ii).Area<opt.sizethreshold
            Io(s(ii).PixelIdxList)=0;
        end
    end
    
end
%%
if ~any(Io(:))
    if isdeployed |1
        %% touch file
        inds = [];
        %         fileID = fopen(outfile,'w');
        %         fclose(fileID);
    end
    return
end
%%
% binarize it before skeletionization
Io = Io>0;
% run skeletonization
% if size(Io) is big limit memory by using less number of nodes
skel = block3d({Io},[200 200 200],fullh,1,@Skeleton3D,[]);
% estimate radius
skel = padarray(skel,ones(1,3),0,'both');
Io = padarray(Io,ones(1,3),0,'both');
bIo=bwdist(~Io);
radskel = double(bIo.*single(skel));

%%
% get the edge pairs
dims = size(skel);
skelinds = find(skel);
if isempty(skelinds)
    % touch file
    if isdeployed |1
        inds = [];
        %         fileID = fopen(outfile,'w');
        %         fclose(fileID);
    end
    return
end
%%
nout = length(dims);
%Compute linear indices
k = [1 cumprod(dims(1:end-1))];
x = [-1:1];
per = zeros(nout^nout,nout);
siz = nout*ones(1,nout);
for i=1:nout
    s = ones(1,nout);
    s(i) = numel(x);
    x = reshape(x,s);
    s = siz;
    s(i) = 1;
    dum = repmat(x,s);
    per(:,i) = dum(:);
end
idxneig = per*k';
idxneig((numel(idxneig)+1)/2)=[];

%% get edge pairs
E = [];
it = 1;
for idx = skelinds(:)'
    inds = idx + idxneig;
    hits = inds(skel(inds));
    rad = radskel(idx);
    % crop back to original size
    E{it} = [[idx*ones(length(hits),1) hits(:) rad*ones(length(hits),1)]]';
    it = it+1;
end
edges = [E{:}]'; clear E
%% map onto original graph
for ii=1:2
    [xx,yy,zz] = ind2sub(dims,edges(:,ii)); % subs on appended crop
    subs = [xx(:),yy(:),zz(:)]-1; % to compansate crop;
    edges(:,ii) = sub2ind(dims-1,subs(:,1),subs(:,2),subs(:,3));
end
% [xx,yy,zz] = ind2sub(dims,edges(:,1)); % subs on appended crop
% subs = [xx(:),yy(:),zz(:)]-1; % to compansate crop;

clear subs_
[keepthese,ia,ic] = unique(edges(:,[1 2]));
[subs_(:,1),subs_(:,2),subs_(:,3)] = ind2sub(dims-1,keepthese);
edges_ = reshape(ic,[],2);
weights_ = edges(ia,3:end);
A = sparse(edges_(:,1),edges_(:,2),1,max(edges_(:,2)),max(edges_(:,2)));
skel = skel(2:end-1,2:end-1,2:end-1);
end
