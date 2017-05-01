function [inds] = octSkel(myh5,myh5prob,BB,outfile,configfile)
%CLUSTER_SKELH5 reads a series of data from octree format and and 
% skeletonize it then writes it in a text file for cluster jobs or 
% returns an array for local tasks
%
% [OUTPUTARGS] = CLUSTER_SKELH5(INPUTARGS)
%
% Inputs:
%     myh5: input hdf5 file path
%     BB: [Xstart Xend Ystart Yend Zstart Zend] Bounding box to be croped
%     (probThr): scalar threshold value
%
% Outputs:
%     (Optional) list of edge pairs
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here
% $Author: base $	$Date: 2016/03/23 11:04:02 $	$Revision: 0.1 $
% Copyright: HHMI 2016
%%
% ox: 70429094
% oy: 12477293
% oz: 31173241
% sx: 16014.695422535211
% sy: 16061.106481481482
% sz: 67494.7875
% nl: 7
if ~isdeployed
    addpath(genpath('./common'))
end
inputfolder = '/nrs/mouselight/cluster/classifierOutputs/2016-10-31/render/161031prob'
skelrange = []; % empty means everything
fileext = '.h5'
seqfile = '/nrs/mouselight/cluster/classifierOutputs/2016-10-31/render/20161031-seq0.txt'
opt = configparser(fullfile(inputfolder,'transform.txt'));
opt.inputfolder = inputfolder;
opt.ext = fileext;
opt.seqtemp = seqfile;
% if isfield(opt,'level')
if isfield(opt,'nl')
    opt.level=opt.nl-1;
else
    error('No level was found')
end
% set block size to a fraction of image size to maximize speed
[~,~,fileext] = fileparts(opt.ext);
if strcmp(fileext,'.h5')
    myh5 = dir(fullfile(opt.inputfolder,'*.h5'));
    info = h5info(fullfile(opt.inputfolder,myh5(1).name));
    imgsiz = info.Datasets.Dataspace.Size;
    opt.imgsiz = imgsiz;
    blocksize = imgsiz/2;
    while any(blocksize>[128 128 128])
        div = blocksize>[128 128 128];
        blocksize = blocksize./(div+1);
    end
    outsiz = opt.imgsiz*2^(opt.level);
else
    mytif = dir(fullfile(opt.inputfolder,'*.tif'));
    info = imfinfo(fullfile(opt.inputfolder,mytif(1).name), 'tif');
    imgsiz = double([info(1).Width info(1).Height length(info)]);
    opt.imgsiz = imgsiz;
    blocksize = imgsiz/2;
    while any(blocksize>[128 128 128])
        div = blocksize>[128 128 128];
        blocksize = blocksize./(div+1);
    end
    outsiz = opt.imgsiz*2^(opt.level);
end
%%
% get sequence
args.level = opt.level;
args.ext = opt.ext;
if exist(opt.seqtemp, 'file') == 2
    % load file directly
else
    mkdir(fileparts(opt.seqtemp))
    args.fid = fopen(opt.seqtemp,'w');
    recdir(opt.inputfolder,args)
end
fid=fopen(opt.seqtemp,'r');
myfiles = textscan(fid,'%s');
myfiles = myfiles{1};
fclose(fid);
%%
RR = zeros(size(myfiles,1),6);
SS = zeros(size(myfiles,1),3);
parfor idx=1:size(myfiles,1)
    %%
    tmps = strsplit(myfiles{idx},filesep);
    seq = [tmps{end-opt.level:end-1}];
    idxtile = str2num(seq);
    lenseq = ceil(log10(idxtile));
    st = [0 0 0];
    bin = opt.imgsiz'*2.^[lenseq-1:-1:0];
    indxyz = zeros(lenseq,3);
    for iseq = 1:lenseq
        is = seq(iseq);
        temp = fliplr(dec2bin(str2num(is)-1,3)); % coordinates wrto xyz
        for ii=1:3 % x y z
            indxyz(iseq,ii) = str2num(temp(ii));
            st(ii) = st(ii) + bin(ii,iseq)*indxyz(iseq,ii);
        end
    end
    RR(idx,:)=[st st+opt.imgsiz];
    SS(idx,:) = [rem(floor(idxtile/10^(opt.level-1)),10) rem(floor(idxtile/10^(opt.level-2)),10) rem(floor(idxtile/10^(opt.level-3)),10)];
end
%%
if nargin<1
    %%
    configfile = 'cmp3_config_skelh5.cfg';
    opt = configparser(configfile);

    myh5 = opt.inputh5;%'/nobackup2/mouselight/cluster/stitching_experiments/renderedvolumes/GN1_tp1_nd4_overlapcut-hdf5_lev-5.h5'
    myh5prob = opt.h5prob;%'/prob1'
    outfile = 'test-1601.txt';
    if opt.brainSize
        brainSize = opt.brainSize;
    else
        inputinfo = h5info(myh5); % opt.inputh5 is redundant for cluster usage
        if length(inputinfo.Groups)>1
            brainSize = inputinfo.Datasets(1).Dataspace.Size;
        else
            brainSize = inputinfo.Datasets.Dataspace.Size;
        end
    end
    opt.brainSize=brainSize;
%%
    probThr = opt.probThr;
    fullh = opt.fullh;

    BB = [1 500 1 500 1 500] + [7489 7489 3193 3193 1 1];
    k=1.0e3;
    bbox = createOverlapBox(brainSize,[k k k],fullh);
    [aa,idx]=min(pdist2(BB,bbox))
    
    BB = bbox(idx,:); % make sure BB is a multiple of chunksize 
    %%
    cluster_skelh5(myh5,myh5prob,BB,outfile,configfile)
else
    opt = configparser(configfile);
    probThr = opt.probThr;
    fullh = opt.fullh;
    % check if brainsize is provided
    if opt.brainSize
        brainSize = opt.brainSize;
    else
        inputinfo = h5info(myh5); % opt.inputh5 is redundant for cluster usage
        brainSize = inputinfo.Datasets.Dataspace.Size;
        opt.brainSize=brainSize;
    end
end

if isdeployed
    inds = [];
    BB = eval(BB);
else
    BB = eval(BB);   
end

%%
starts = BB(1:2:end);
ends = BB(2:2:end);
datasiz = ends-starts+1;

Io = squeeze(h5read(myh5,myh5prob,starts,datasiz));
% Io = squeeze(h5read(myh5,myh5prob,starts+[500 300 100],datasiz));
% figure, imshow(squeeze(max(Io,[],3))',[]),
if ~any(Io(:))
    if isdeployed |1
        %% touch file
        fileID = fopen(outfile,'w');
        fclose(fileID);
    end
    return
end

% smooth image
Io = smooth3(Io,'gaussian',[3 3 1]);
Io = Io>probThr;
%%
if ~any(Io(:))
    if isdeployed |1
        %% touch file
        fileID = fopen(outfile,'w');
        fclose(fileID);
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
        fileID = fopen(outfile,'w');
        fclose(fileID);
    end
    return
end
%%
% binarize it before skeletionization
Io = Io>0;
% run skeletonization
% if size(Io) is big limit memory by using less number of nodes
skel = block3d({Io},[200 200 200],fullh,1,@Skeleton3D,[]);
skel = padarray(skel,ones(1,3),0,'both');
%%
% get the edge pairs
dims = size(skel);
skelinds = find(skel);
if isempty(skelinds)
    % touch file
    if isdeployed |1
        fileID = fopen(outfile,'w');
        fclose(fileID);
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
    % crop back to original size
    E{it} = [[idx*ones(length(hits),1) hits(:)]]';
    it = it+1;
end
edges = [E{:}]'; clear E

%% map onto original graph
inds = zeros(size(edges));
for ii=1:2
    [xx,yy,zz] = ind2sub(dims,edges(:,ii)); % subs on appended crop
    subs = [xx(:),yy(:),zz(:)]-1; % to compansate crop;
    subs = subs + ones(size(subs,1),1)*BB(1:2:end)-1; % convert to original subs
    inds(:,ii) = sub2ind(brainSize,subs(:,1),subs(:,2),subs(:,3)); % convert to original inds
end
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
    fprintf(fileID,'%d %d\n',inds');
    fclose(fileID);
end
end
function runlocal
%%
addpath(genpath('./common'))
clear all
clc
configfile = fullfile(pwd,'./config_files/20150619_oct12_config_skelh5.cfg');
opt = configparser(configfile);

myh5 = opt.inputh5;
inputinfo = h5info(myh5);
numGroups = length(inputinfo.Datasets);
idxGroup = 1;

%
if numGroups>1
    myh5prob = ['/',inputinfo.Datasets(idxGroup).Name];
    cropSize = 10*inputinfo.Datasets(idxGroup).ChunkSize;
    % to get %10 overlap overhead use multiple of 10
    fullh = inputinfo.Datasets(idxGroup).ChunkSize; % add 1 to make it odd (heuristic)
    RR = h5read(myh5,sprintf('%s/ROI',inputinfo.Groups(idxGroup).Name));
    brainSize = inputinfo.Datasets(idxGroup).Dataspace.MaxSize;
else
    myh5prob = ['/',inputinfo.Datasets.Name];
    cropSize = 10*inputinfo.Datasets.ChunkSize;
    % to get %10 overlap overhead use multiple of 10
    fullh = inputinfo.Datasets.ChunkSize; % add 1 to make it odd (heuristic)
    RR = h5read(myh5,sprintf('%s/ROI',inputinfo.Groups.Name));
    brainSize = inputinfo.Datasets.Dataspace.MaxSize;
end
%
[aa,bb,cc]=fileparts(myh5);
% outfolder = '/nobackup2/mouselight/cluster/GN1_autorecon_05/'
outfolder = opt.outfolder;
mkdir(outfolder)
%%
%
% rmdir(outfolder)
s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
compiledfunc = '/groups/mousebrainmicro/home/base/CODE/MATLAB/recontree/compiled/compiledfiles_skelh5/cluster_skelh5'
%find number of random characters to choose from
numRands = length(s);
%specify length of random string to generate
sLength = 10;
%-o /dev/null
% chunck data
%
% fullh = opt.fullh; % 15
%

bbox = createOverlapBox(brainSize,cropSize,fullh);
% bbox = createOverlapBox(brainSize,[cropSize cropSize cropSize],fullh);
%
BBoxes = RR(:,[1 4 2 5 3 6])+1;

X = BBoxes(:,1:2);
Y = BBoxes(:,3:4);
Z = BBoxes(:,5:6);
XYZ = unique([X(:),Y(:),Z(:)],'rows');
in = inhull([bbox(:,1:2:end);bbox(:,2:2:end)],XYZ);
in = any(reshape(in,[],2),2);

finished = zeros(1,size(bbox,1));
if 1 % check any missing file
    myfiles = dir([outfolder,'*.txt']);
    for ii=1:length(myfiles)
        rt=strsplit(myfiles(ii).name,'idx-');
        finished(str2num(rt{2}(1:5))) = 1;
    end
end
%%
for idx = 1:size(bbox,1)
    %%
    %generate random string
    BB = sprintf('[%d %d %d %d %d %d]',bbox(idx,:));
    %% check if BB is outsize of BBoxes
    if ~in(idx) | finished(idx)% skip
        (idx)
        continue
    end
    %%
    outfile = fullfile(outfolder,sprintf('%s_idx-%05d_stxyzendxyz-%d_%d_%d_%d_%d_%d.txt',bb,idx,bbox(idx,1:2:end),bbox(idx,2:2:end)));
    cluster_skelh5(myh5,myh5prob,BB,outfile,configfile)
end
end

function deployment
%qsub -pe batch 4 -l short=true -N tile_test -j y -o ~/logs -b y -cwd -V './compiledfiles_mytest/mytest > output_mytest.log'
%%
% mcc -m -R -nojvm -v cluster_skelh5.m -d ./compiled/compiledfiles_skelh5  -a ./common
%%
addpath(genpath('./common'))
clear all
clc
numcores = 16;
mysh = '20150619_oct12_skel_missing.sh';
configfile = fullfile(pwd,'./config_files/20150619_oct12_config_skelh5.cfg');
opt = configparser(configfile);
%
% myh5 = '/srv/data/probGN1_lvl-5.h5'
% myh5 = '/tier2/mousebrainmicro/mousebrainmicro/cluster/hdf5test/merge_probGN1_lvl-5.h5'
% myh5prob='/renderedVolume'
% myh5 = '/data3/renderedData/2015-07-11/2015-07-11-G3457_lev-3.h5'
myh5 = opt.inputh5;
inputinfo = h5info(myh5);
numGroups = length(inputinfo.Datasets);
idxGroup = 1;

%
if numGroups>1
    myh5prob = ['/',inputinfo.Datasets(idxGroup).Name];
    cropSize = 10*inputinfo.Datasets(idxGroup).ChunkSize;
    % to get %10 overlap overhead use multiple of 10
    fullh = inputinfo.Datasets(idxGroup).ChunkSize; % add 1 to make it odd (heuristic)
    RR = h5read(myh5,sprintf('%s/ROI',inputinfo.Groups(idxGroup).Name));
    brainSize = inputinfo.Datasets(idxGroup).Dataspace.MaxSize;
else
    myh5prob = ['/',inputinfo.Datasets.Name];
    cropSize = 10*inputinfo.Datasets.ChunkSize;
    % to get %10 overlap overhead use multiple of 10
    fullh = inputinfo.Datasets.ChunkSize; % add 1 to make it odd (heuristic)
    RR = h5read(myh5,sprintf('%s/ROI',inputinfo.Groups.Name));
    brainSize = inputinfo.Datasets.Dataspace.MaxSize;
end
%
[aa,bb,cc]=fileparts(myh5);
% outfolder = '/nobackup2/mouselight/cluster/GN1_autorecon_05/'
outfolder = opt.outfolder;
mkdir(outfolder)
%%
%
% rmdir(outfolder)
s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
compiledfunc = '/groups/mousebrainmicro/home/base/CODE/MATLAB/recontree/compiled/compiledfiles_skelh5/cluster_skelh5'
%find number of random characters to choose from
numRands = length(s);
%specify length of random string to generate
sLength = 10;
%-o /dev/null
% chunck data
%
% fullh = opt.fullh; % 15
%

bbox = createOverlapBox(brainSize,cropSize,fullh);
% bbox = createOverlapBox(brainSize,[cropSize cropSize cropSize],fullh);
%
if 1 % gets the BBs used for creating h5
    BBoxes = RR(:,[1 4 2 5 3 6])+1;
else
    load BBoxes
end
X = BBoxes(:,1:2);
Y = BBoxes(:,3:4);
Z = BBoxes(:,5:6);
XYZ = unique([X(:),Y(:),Z(:)],'rows');
in = inhull([bbox(:,1:2:end);bbox(:,2:2:end)],XYZ);
in = any(reshape(in,[],2),2);
%
timelim = 5*60
finished = zeros(1,size(bbox,1));

if 1 % check any missing file
    myfiles = dir([outfolder,'*.txt']);
    for ii=1:length(myfiles)
        rt=strsplit(myfiles(ii).name,'idx-');
        finished(str2num(rt{2}(1:5))) = 1;
    end
end
%%
fid = fopen(mysh,'w');
for idx = 1:size(bbox,1)
    %%
    %generate random string
    BB = bbox(idx,:);
    %% check if BB is outsize of BBoxes
    if ~in(idx) | finished(idx)% skip
        (idx)
        continue
    end
    %%
    randString = s( ceil(rand(1,sLength)*numRands) );
    outfile = fullfile(outfolder,sprintf('%s_idx-%05d_stxyzendxyz-%d_%d_%d_%d_%d_%d.txt',bb,idx,BB(1:2:end),BB(2:2:end)));
    name = sprintf('skel_%05d-%s',idx,randString);
    args = sprintf('''%s %s %s "[%d,%d,%d,%d,%d,%d]" %s %s> output.log''',compiledfunc,myh5,myh5prob,(BB),outfile,configfile);
    mysub = sprintf('qsub -pe batch %d -l h_rt=%d -N %s -j y -o ~/logs -b y -cwd -V %s\n',numcores,timelim,name,args);
    fwrite(fid,mysub);
end
unix(sprintf('chmod +x %s',mysh));
%%
% test
cluster_skelh5('/nobackup2/mouselight/cluster/2016-10-25/20161025_prob0/20161025_prob0_lev-6_chunk-367_888.h5',...
'/prob0','[8961,9660,19366,20505,4626,5120]',...
'/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/2016-10-25/prob0/20161025_prob0_lev-6_chunk-367_888/20161025_prob0_lev-6_chunk-367_888_idx-00036_stxyzendxyz-8961_19366_4626_9660_20505_5120.txt ',...
'/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/2016-10-25/skel_prob0_cfgfiles/367.cfg')
%%
cluster_skelh5('/nrs/mouselight/cluster/classifierOutputs/2015-06-19/150619prob_octants12_prob0_lev-5.h5',...
    '/prob0','[12241,13040,5176,6325,811,1710]',...
    '/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/2015-06-19/oct12/prob0/150619prob_octants12_prob0_lev-5_idx-01158_stxyzendxyz-12241_5176_811_13040_6325_1710.txt',...
    '/groups/mousebrainmicro/home/base/CODE/MATLAB/recontree/config_files/20150619_oct12_config_skelh5.cfg')
%%
% (myh5,myh5prob,BB,outfile,configfile)
% mcc -m -R -nojvm -v cluster_skelh5.m -d ./compiled/compiledfiles_skelh5  -a ./common
qsub -pe batch 4 -l short=true -N skel_00046-T9wrzWzvoB -j y -o ~/logs -b y -cwd -V './compiled/compiledfiles_skelh5/cluster_skelh5 /nobackup2/mouselight/cluster/stitching_experiments/renderedvolumes/GN1_tp1_nd4_minopt_lev-5.h5 /prob1 "[10369,11088,1027,2166,1,1040]" /groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/2015-06-19/GN1/prob1/GN1_tp1_nd4_minopt_lev-5_idx-00046_stxyzendxyz-10369_1027_1_11088_2166_1040.txt ./config_files/cmp3_config_skelh5.cfg> output.log'



end
































