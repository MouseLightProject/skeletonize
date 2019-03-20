function cluster_skelh5_2018_12_01_take_2()
    sample='20181201' ;
    myh5prob = '/prob0' ;
    exp = sprintf('%s_%s',sample,myh5prob(2:end));
    configfile = fullfile(pwd,sprintf('./config_files/%s_config_skelh5_take_2.cfg',exp));
    mysh = sprintf('%s_%s.sh',datestr(now,'ddmmyyHHMM'),exp) ;
    generate_shell_script_that_bsubs_jobs(configfile,mysh,myh5prob)
end



function generate_shell_script_that_bsubs_jobs(configfile,mysh,myh5prob)
    %qsub -pe batch 4 -l short=true -N tile_test -j y -o ~/logs -b y -cwd -V './compiledfiles_mytest/mytest > output_mytest.log'
    %%
    % mcc -m -R -nojvm -v cluster_skelh5.m -d ./compiled/compiledfiles_skelh5  -a ./common
    % mcc -m -R -nojvm -v /groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/skeletonize/cluster_skelh5.m -d /groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/skeletonization -a /groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/skeletonize/common
    %%
    %addpath(genpath('./common'))
    % clear all
    %clc
    numcores = 8;
    % mysh = '20150619_oct12config_skelh5_miss.sh';
    opt = configparser(configfile);
    %
    % myh5 = '/srv/data/probGN1_lvl-5.h5'
    % myh5 = '/tier2/mousebrainmicro/mousebrainmicro/cluster/hdf5test/merge_probGN1_lvl-5.h5'
    % myh5prob='/renderedVolume'
    % myh5 = '/data3/renderedData/2015-07-11/2015-07-11-G3457_lev-3.h5'
    myh5 = opt.inputh5;
%     if nargin<3 ,
%         myh5prob = '/prob0' ;
%     end
    % likely breakpoint location
    [brainSize,RR,chunk_dims] = h5parser(myh5,myh5prob);
    % get a multiple of chunksize that is around 1000^3
    cropSize = round(1000./chunk_dims).*chunk_dims;
    % cropSize = 10*chunk_dims;%inputinfo.Datasets.ChunkSize;
    % to get %10 overlap overhead use multiple of 10
    fullh = chunk_dims; % add 1 to make it odd (heuristic)
    %
    [~,bb,~]=fileparts(myh5);
    % outfolder = '/nobackup2/mouselight/cluster/GN1_autorecon_05/'
    outfolder = opt.outfolder;
    if ~exist(outfolder, 'dir') ,
        mkdir(outfolder)
    end
    unix(sprintf('umask g+rxw %s',outfolder))
    unix(sprintf('chmod g+rwx %s',outfolder));
    %%
    %
    % rmdir(outfolder)
    s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
    %old = 0;
    % if old
    %     compiledfunc = '/groups/mousebrainmicro/home/base/CODE/MATLAB/recontree/compiled/compiledfiles_skelh5/cluster_skelh5'
    % else
    %     compiledfunc = '/groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/skeletonization/cluster_skelh5'
    % end
    script_file_path = mfilename('fullpath') ;
    script_folder_path = fileparts(script_file_path) ;
    compiledfunc = fullfile(script_folder_path, 'compiled', 'run_cluster_skelh5_at_janelia.sh') ;


    %find number of random characters to choose from
    numRands = length(s);
    %specify length of random string to generate
    sLength = 10;
    %-o /dev/null

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
    total_number_of_boxes = sum(in) ;
    fprintf('Total number of boxes: %d\n', total_number_of_boxes) ;
    
    %
    timelim = 10*60 ;
    finished = zeros(1,size(bbox,1));

    if 1 % check any missing file
        myfiles = dir([outfolder,'*.txt']);
        for ii=1:length(myfiles)
            rt=strsplit(myfiles(ii).name,'idx-');
            finished(str2num(rt{2}(1:5))) = 1;   %#ok<ST2NM>
        end
    end
    sum(~finished)
    %%
    % likely breakpoint location
    iter=0;
    fid = fopen(mysh,'w');
    for idx = 1:size(bbox,1)
        %%
        %generate random string
        BB = bbox(idx,:);
        %% check if BB is outsize of BBoxes
        if ~in(idx) || finished(idx) ,  % skip
            continue
        end
        %%
        randString = s( ceil(rand(1,sLength)*numRands) );
        outfile = fullfile(outfolder,sprintf('%s_idx-%05d_stxyzendxyz-%d_%d_%d_%d_%d_%d.txt',bb,idx,BB(1:2:end),BB(2:2:end)));
        name = sprintf('skel_%05d-%s',idx,randString);
        argsout = sprintf('''%s %s %s "[%d,%d,%d,%d,%d,%d]" %s %s''',compiledfunc,myh5,myh5prob,(BB),outfile,configfile);
        % mysub = sprintf('qsub -pe batch %d -l d_rt=%d -N %s -j y -o /dev/null -b y -cwd -V %s\n',numcores,timelim,name,args);
        mysub = sprintf('bsub -P mouselight -n%d -R"affinity[core(1)]" -We %d -J %s -o %s %s\n',numcores,timelim/60,name,'/dev/null',argsout);
        fwrite(fid,mysub);
        iter=iter+1;
        sprintf('iter: %d',iter)
    end
    fclose(fid) ;
    unix(sprintf('chmod +x %s',mysh));
end
