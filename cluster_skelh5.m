function cluster_skelh5(myh5, myh5prob, BB, outfile, sizethreshold, probThr, fullh)
    %CLUSTER_SKELH5 crop a section from h5 and skeletonize it then writes it in
    %a text file for cluster jobs or returns an array for local tasks
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
    if ~isdeployed
        addpath(genpath('./common'))
    end

    %opt = configparser(configfile);
    %probThr = opt.probThr;
    %fullh = opt.fullh;
    do_viz = false ;    

    %     cropSize = 10*chunk_dims;%inputinfo.Datasets.ChunkSize;
    %     % to get %10 overlap overhead use multiple of 10
    %     fullh = chunk_dims; % add 1 to make it odd (heuristic)

    brainSize = h5parser(myh5, myh5prob);

    if isdeployed ,
        BB = eval(BB);
        sizethreshold = eval(sizethreshold) ;
        probThr = eval(probThr) ;
        fullh = eval(fullh) ;
    end

    %%
    starts = BB(1:2:end);
    ends = BB(2:2:end);
    datasiz = ends-starts+1;

    Io = squeeze(h5read(myh5,myh5prob,starts,datasiz));
    % Io = squeeze(h5read(myh5,myh5prob,starts+[500 300 100],datasiz));
    % figure, imshow(squeeze(max(Io,[],3))',[]),
    if ~any(Io(:))
        if isdeployed || true
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
        if isdeployed || true
            %% touch file
            fileID = fopen(outfile,'w');
            fclose(fileID);
        end
        return
    end
    %% cleanup image
    s  = regionprops(Io, 'centroid','PixelIdxList','Area');
%     if 0
%         % fast
%         Iout = zeros(size(Io),'single');
%         for ii=1:length(s)
%             if s(ii).Area>10
%                 Iout(s(ii).PixelIdxList)=Io(s(ii).PixelIdxList);
%             end
%         end
%         Io = Iout; clear Iout;
%     else
    % memory efficient
    for ii=1:length(s)
        if s(ii).Area<sizethreshold
            Io(s(ii).PixelIdxList)=0;
        end
    end
%     end
    %%
    if ~any(Io(:))
        if isdeployed || true
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
    % Heuristic: 0 out boundary pixels to prevent replicating skels in the
    % overlaping region
    s = round((fullh+1)/2);
    skel(1:s,:,:) = 0;
    skel(end-s+1:end,:,:) = 0;
    skel(:,1:s,:) = 0;
    skel(:,end-s+1:end,:) = 0;
    skel(:,:,1:s) = 0;
    skel(:,:,end-s+1:end) = 0;
    %%
    % estimate radius
    Io = padarray(Io,ones(1,3),0,'both');
    bIo=bwdist(~Io);
    radskel = double(bIo.*single(skel));

    %%
    % get the edge pairs
    dims = size(skel);
    skelinds = find(skel);
    if isempty(skelinds)
        % touch file
        if isdeployed || true
            fileID = fopen(outfile,'w');
            fclose(fileID);
        end
        return
    end
    %%
    nout = length(dims);
    %Compute linear indices
    k = [1 cumprod(dims(1:end-1))];
    x = (-1:1) ;
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
        E{it} = [idx*ones(length(hits),1) hits(:) rad*ones(length(hits),1)]';  %#ok<AGROW>
        it = it+1;
    end
    edges = [E{:}]'; clear E

    %% map onto original graph
    inds = edges;%zeros(size(edges));
    for ii=1:2
        [xx,yy,zz] = ind2sub(dims,edges(:,ii)); % subs on appended crop
        subs = [xx(:),yy(:),zz(:)]-1; % to compansate crop;
        subs = subs + ones(size(subs,1),1)*BB(1:2:end)-1; % convert to original subs
        inds(:,ii) = sub2ind(brainSize,subs(:,1),subs(:,2),subs(:,3)); % convert to original inds
    end
    %%
    if ~isdeployed && do_viz
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
    if isdeployed || true
        %%
        fileID = fopen(outfile,'w');
        if size(inds,2)==2
            fprintf(fileID,'%d %d\n',inds');
        else
            fprintf(fileID,'%d %d %.2f\n',inds');
        end
        fclose(fileID);
    end
end  % function

% function runlocal(configfile,myh5prob)
% %%
% % cluster_skelh5('/nrs/mouselight/cluster/classifierOutputs/2015-06-19/150619prob_octants12_prob0_lev-5.h5',...
% %     '/prob0','[12241,13040,5176,6325,811,1710]',...
% %     '/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/2015-06-19/oct12/prob0/150619prob_octants12_prob0_lev-5_idx-01158_stxyzendxyz-12241_5176_811_13040_6325_1710.txt',...
% %     '/groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/skeletonize/config_files/20150619_octant12_prob0_config_skelh5.cfg')
% 
% addpath(genpath('./common'))
% % clear all
% %%
% clc
% opt = configparser(configfile);
% 
% myh5 = opt.inputh5;
% if nargin<3
%     myh5prob = opt.h5prob
% end
% if 1
% %     fid = H5F.open(myh5);
% %     dset_id = H5D.open(fid,myh5prob);
% %     space = H5D.get_space(dset_id);
% %     [~,dims] = H5S.get_simple_extent_dims(space);
% %     H5S.close(space);
% %     
% %     dcpl = H5D.get_create_plist(dset_id);
% %     [rank,chunk_dims] = H5P.get_chunk(dcpl);
% %     H5P.close(dcpl);
% % 
% %     brainSize = dims([3 2 1]);
% %     chunk_dims = chunk_dims([3 2 1]);
% %     RR = h5read(myh5,[myh5prob,'_props/ROI']);
% %     
% %     cropSize = 10*chunk_dims;%inputinfo.Datasets.ChunkSize;
% %     % to get %10 overlap overhead use multiple of 10
% %     fullh = chunk_dims; % add 1 to make it odd (heuristic)
% %     
% %     H5D.close(dset_id);
% %     H5F.close(fid);
%     [brainSize,RR,chunk_dims,rank] = h5parser(myh5,myh5prob);
%     cropSize = round(1000./chunk_dims).*chunk_dims;
%     %cropSize = 10*chunk_dims;%inputinfo.Datasets.ChunkSize;
%     % to get %10 overlap overhead use multiple of 10
%     fullh = chunk_dims; % add 1 to make it odd (heuristic)
% else
%     inputinfo = h5info(myh5)
%     numGroups = length(inputinfo.Datasets);
%     idxGroup = 1;
%     
%     %
%     if numGroups>1
%         myh5prob = ['/',inputinfo.Datasets(idxGroup).Name];
%         cropSize = 10*inputinfo.Datasets(idxGroup).ChunkSize;
%         % to get %10 overlap overhead use multiple of 10
%         fullh = inputinfo.Datasets(idxGroup).ChunkSize; % add 1 to make it odd (heuristic)
%         RR = h5read(myh5,sprintf('%s/ROI',inputinfo.Groups(idxGroup).Name));
%         brainSize = inputinfo.Datasets(idxGroup).Dataspace.MaxSize;
%     else
%         myh5prob = ['/',inputinfo.Datasets.Name];
%         cropSize = 10*inputinfo.Datasets.ChunkSize;
%         % to get %10 overlap overhead use multiple of 10
%         fullh = inputinfo.Datasets.ChunkSize; % add 1 to make it odd (heuristic)
%         RR = h5read(myh5,sprintf('%s/ROI',inputinfo.Groups.Name));
%         brainSize = inputinfo.Datasets.Dataspace.MaxSize;
%     end
% end%
% [aa,bb,cc]=fileparts(myh5);
% % outfolder = '/nobackup2/mouselight/cluster/GN1_autorecon_05/'
% outfolder = opt.outfolder;
% mkdir(outfolder)
% unix(sprintf('umask g+rxw %s',outfolder))
% unix(sprintf('chmod g+rwx %s',outfolder));
% 
% %%
% %
% % rmdir(outfolder)
% s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
% %compiledfunc = '/groups/mousebrainmicro/home/base/CODE/MATLAB/recontree/compiled/compiledfiles_skelh5/cluster_skelh5'
% %find number of random characters to choose from
% numRands = length(s);
% %specify length of random string to generate
% sLength = 10;
% %-o /dev/null
% % chunck data
% %
% % fullh = opt.fullh; % 15
% %
% 
% bbox = createOverlapBox(brainSize,cropSize,fullh);
% % bbox = createOverlapBox(brainSize,[cropSize cropSize cropSize],fullh);
% %
% BBoxes = RR(:,[1 4 2 5 3 6])+1;
% 
% X = BBoxes(:,1:2);
% Y = BBoxes(:,3:4);
% Z = BBoxes(:,5:6);
% XYZ = unique([X(:),Y(:),Z(:)],'rows');
% in = inhull([bbox(:,1:2:end);bbox(:,2:2:end)],XYZ);
% in = any(reshape(in,[],2),2);
% 
% finished = zeros(1,size(bbox,1));
% if 1 % check any missing file
%     myfiles = dir([outfolder,'*.txt']);
%     for ii=1:length(myfiles)
%         rt=strsplit(myfiles(ii).name,'idx-');
%         finished(str2num(rt{2}(1:5))) = 1;
%     end
% end
% sum(finished)
% %%
% for idx = 1:size(bbox,1)
%     %%
%     %generate random string
%     BB = sprintf('[%d %d %d %d %d %d]',bbox(idx,:));
%     %% check if BB is outsize of BBoxes
%     if ~in(idx) | finished(idx)% skip
%         (idx)
%         continue
%     end
%     %%
%     outfile = fullfile(outfolder,sprintf('%s_idx-%05d_stxyzendxyz-%d_%d_%d_%d_%d_%d.txt',bb,idx,bbox(idx,1:2:end),bbox(idx,2:2:end)));
%     cluster_skelh5(myh5,myh5prob,BB,outfile,configfile)
% end
% end  % runlocal()


% function deployment(configfile,mysh,myh5prob)
% %qsub -pe batch 4 -l short=true -N tile_test -j y -o ~/logs -b y -cwd -V './compiledfiles_mytest/mytest > output_mytest.log'
% %%
% % mcc -m -R -nojvm -v cluster_skelh5.m -d ./compiled/compiledfiles_skelh5  -a ./common
% % mcc -m -R -nojvm -v /groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/skeletonize/cluster_skelh5.m -d /groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/skeletonization -a /groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/skeletonize/common
% %%
% addpath(genpath('./common'))
% % clear all
% clc
% numcores = 8;
% % mysh = '20150619_oct12config_skelh5_miss.sh';
% opt = configparser(configfile);
% %
% % myh5 = '/srv/data/probGN1_lvl-5.h5'
% % myh5 = '/tier2/mousebrainmicro/mousebrainmicro/cluster/hdf5test/merge_probGN1_lvl-5.h5'
% % myh5prob='/renderedVolume'
% % myh5 = '/data3/renderedData/2015-07-11/2015-07-11-G3457_lev-3.h5'
% myh5 = opt.inputh5;
% if nargin<3
%     myh5prob = '/prob0'
% end
% if 1
%     % likely breakpoint location
%     [brainSize,RR,chunk_dims,rank] = h5parser(myh5,myh5prob);
%     % get a multiple of chunksize that is around 1000^3
%     cropSize = round(1000./chunk_dims).*chunk_dims;
%     % cropSize = 10*chunk_dims;%inputinfo.Datasets.ChunkSize;
%     % to get %10 overlap overhead use multiple of 10
%     fullh = chunk_dims; % add 1 to make it odd (heuristic)
% end
% %
% [aa,bb,cc]=fileparts(myh5);
% % outfolder = '/nobackup2/mouselight/cluster/GN1_autorecon_05/'
% outfolder = opt.outfolder;
% mkdir(outfolder)
% unix(sprintf('umask g+rxw %s',outfolder))
% unix(sprintf('chmod g+rwx %s',outfolder));
% %%
% %
% % rmdir(outfolder)
% s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
% %old = 0;
% % if old
% %     compiledfunc = '/groups/mousebrainmicro/home/base/CODE/MATLAB/recontree/compiled/compiledfiles_skelh5/cluster_skelh5'
% % else
% %     compiledfunc = '/groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/skeletonization/cluster_skelh5'
% % end
% script_file_path = mfilename('fullpath') ;
% script_folder_path = fileparts(script_file_path) ;
% compiledfunc = fullfile(script_folder_path, 'compiled', 'run_cluster_skelh5_at_janelia.sh') ;
% 
% 
% %find number of random characters to choose from
% numRands = length(s);
% %specify length of random string to generate
% sLength = 10;
% %-o /dev/null
% 
% bbox = createOverlapBox(brainSize,cropSize,fullh);
% % bbox = createOverlapBox(brainSize,[cropSize cropSize cropSize],fullh);
% %
% if 1 % gets the BBs used for creating h5
%     BBoxes = RR(:,[1 4 2 5 3 6])+1;
% else
%     load BBoxes
% end
% X = BBoxes(:,1:2);
% Y = BBoxes(:,3:4);
% Z = BBoxes(:,5:6);
% XYZ = unique([X(:),Y(:),Z(:)],'rows');
% in = inhull([bbox(:,1:2:end);bbox(:,2:2:end)],XYZ);
% in = any(reshape(in,[],2),2);
% %
% timelim = 10*60
% finished = zeros(1,size(bbox,1));
% 
% if 1 % check any missing file
%     myfiles = dir([outfolder,'*.txt']);
%     for ii=1:length(myfiles)
%         rt=strsplit(myfiles(ii).name,'idx-');
%         finished(str2num(rt{2}(1:5))) = 1;
%     end
% end
% sum(~finished)
% %%
% % likely breakpoint location
% iter=0;
% fid = fopen(mysh,'w');
% for idx = 1:size(bbox,1)
%     %%
%     %generate random string
%     BB = bbox(idx,:);
%     %% check if BB is outsize of BBoxes
%     if ~in(idx) | finished(idx)% skip
%         continue
%     end
%     %%
%     randString = s( ceil(rand(1,sLength)*numRands) );
%     outfile = fullfile(outfolder,sprintf('%s_idx-%05d_stxyzendxyz-%d_%d_%d_%d_%d_%d.txt',bb,idx,BB(1:2:end),BB(2:2:end)));
%     name = sprintf('skel_%05d-%s',idx,randString);
%     argsout = sprintf('''%s %s %s "[%d,%d,%d,%d,%d,%d]" %s %s''',compiledfunc,myh5,myh5prob,(BB),outfile,configfile);
%     % mysub = sprintf('qsub -pe batch %d -l d_rt=%d -N %s -j y -o /dev/null -b y -cwd -V %s\n',numcores,timelim,name,args);
%     mysub = sprintf('bsub -P mouselight -n%d -R"affinity[core(1)]" -We %d -J %s -o %s %s\n',numcores,timelim/60,name,'/dev/null',argsout);
%     fwrite(fid,mysub);
%     iter=iter+1;
%     sprintf('iter: %d',iter)
% %     % for debugging
% %     if iter==100 ,
% %         break
% %     end
% end
% fclose(fid) ;
% unix(sprintf('chmod +x %s',mysh));
% % %%
% % % test
% %%
% % cluster_skelh5('/nrs/mouselight/cluster/classifierOutputs/2015-06-19_backup/20150619_oct12_prob0_FC/150619prob_octants12_prob0_FC_lev-5_chunk-111_111_masked-0.h5',...
% % '/prob0','[10801,12000,6193,7052,1585,2464]',...
% % './test.txt',...
% % '/groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/skeletonize/config_files/20150619_octant12_prob0_config_skelh5.cfg')
% % /groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/2015-06-19/octant12/prob0/150619prob_octants12_prob0_FC_lev-5_chunk-111_111_masked-0_idx-01113_stxyzendxyz-10801_6193_1585_12000_7052_2464.txt
% % %%
% % cluster_skelh5('/nrs/mouselight/cluster/classifierOutputs/2015-06-19/150619prob_octants12_prob0_lev-5.h5',...
% %     '/prob0','[12241,13040,5176,6325,811,1710]',...
% %     '/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/2015-06-19/oct12/prob0/150619prob_octants12_prob0_lev-5_idx-01158_stxyzendxyz-12241_5176_811_13040_6325_1710.txt',...
% %     '/groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/skeletonize/config_files/20150619_octant12_prob0_config_skelh5.cfg')
% % %%
% % % (myh5,myh5prob,BB,outfile,configfile)
% % % mcc -m -R -nojvm -v cluster_skelh5.m -d ./compiled/compiledfiles_skelh5  -a ./common
% % qsub -pe batch 4 -l short=true -N skel_00046-T9wrzWzvoB -j y -o ~/logs -b y -cwd -V './compiled/compiledfiles_skelh5/cluster_skelh5 /nobackup2/mouselight/cluster/stitching_experiments/renderedvolumes/GN1_tp1_nd4_minopt_lev-5.h5 /prob1 "[10369,11088,1027,2166,1,1040]" /groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/2015-06-19/GN1/prob1/GN1_tp1_nd4_minopt_lev-5_idx-00046_stxyzendxyz-10369_1027_1_11088_2166_1040.txt ./config_files/cmp3_config_skelh5.cfg> output.log'
% %
% 
% 
% end
































