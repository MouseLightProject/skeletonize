function submit_bsub_jobs_for_skeletonization(sample_date, sizethreshold, probThr, fullh)
    sample_date
    sizethreshold
    probThr
    fullh
    
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
    %opt = configparser(configuration_file_name);
    %output_folder_path = opt.outfolder;
    output_folder_path = sprintf('/nrs/mouselight/cluster/classifierOutputs/%s/skeletonization', sample_date) ;
    %p_map_dataset_path = opt.h5prob ;
    p_map_dataset_path = '/prob0' ;    
    %whole_brain_p_map_h5_file_path = opt.inputh5 ;
    whole_brain_p_map_h5_file_path = ...
        sprintf('/nrs/mouselight/cluster/classifierOutputs/%s/whole-brain-p-map-as-h5/%s-whole-brain-p-map.h5', sample_date, sample_date) ;

    %bsub_script_file_path = determine_bsub_script_file_path(output_folder_path) ;
    
    %
    % myh5 = '/srv/data/probGN1_lvl-5.h5'
    % myh5 = '/tier2/mousebrainmicro/mousebrainmicro/cluster/hdf5test/merge_probGN1_lvl-5.h5'
    % myh5prob='/renderedVolume'
    % myh5 = '/data3/renderedData/2015-07-11/2015-07-11-G3457_lev-3.h5'
%     if nargin<3 ,
%         myh5prob = '/prob0' ;
%     end
    % likely breakpoint location
    [brainSize,RR,chunk_dims] = h5parser(whole_brain_p_map_h5_file_path, p_map_dataset_path);
    % get a multiple of chunksize that is around 1000^3
    cropSize = round(1000./chunk_dims).*chunk_dims;
    % cropSize = 10*chunk_dims;%inputinfo.Datasets.ChunkSize;
    % to get %10 overlap overhead use multiple of 10
    fullh_for_create_ovelap_box = chunk_dims; % add 1 to make it odd (heuristic)
    %
    [~,whole_brain_p_map_h5_base_name,~] = fileparts(whole_brain_p_map_h5_file_path) ;
    % outfolder = '/nobackup2/mouselight/cluster/GN1_autorecon_05/'
    if ~exist(output_folder_path, 'dir') ,
        mkdir(output_folder_path) ;
    end
    unix(sprintf('umask g+rxw %s',output_folder_path)) ;
    unix(sprintf('chmod g+rwx %s',output_folder_path)) ;
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
    %script_file_path = mfilename('fullpath') ;
    %script_folder_path = fileparts(script_file_path) ;
    %compiled_function_front_end_script_path = fullfile(script_folder_path, 'compiled', 'run_cluster_skelh5_at_janelia.sh') ;


    %find number of random characters to choose from
    numRands = length(s);
    %specify length of random string to generate
    sLength = 10;
    %-o /dev/null

    bbox = createOverlapBox(brainSize,cropSize,fullh_for_create_ovelap_box);
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
    time_limit_in_seconds = 10*60 ;
    is_finished = false(1,size(bbox,1));

    % check any missing file
    output_text_file_template_path = fullfile(output_folder_path, '*.txt') ;
    extant_output_text_file_names = simple_dir(output_text_file_template_path) ;
    for ii = 1:length(extant_output_text_file_names) ,
        extant_output_text_file_name = extant_output_text_file_names{ii} ;
        bbox_index = bounding_box_index_from_file_name(extant_output_text_file_name) ;
        is_finished(bbox_index) = true ;
    end
    %sum(~is_finished)
    %%
    % likely breakpoint location
    jobs_submitted_count=0;
    %fid = fopen(bsub_script_file_path,'w');
    for bounding_box_index = 1:size(bbox,1)
        BB = bbox(bounding_box_index, :) ;
        % check if BB is outsize of BBoxes
        if ~in(bounding_box_index) || is_finished(bounding_box_index) ,  % skip
            continue
        end
        random_string_for_job_name = s( ceil(rand(1,sLength)*numRands) ) ;
        output_file_name = sprintf('%s_idx-%05d_stxyzendxyz-%d_%d_%d_%d_%d_%d.txt', ...
                                   whole_brain_p_map_h5_base_name, ...
                                   bounding_box_index, ...
                                   BB(1:2:end), ...
                                   BB(2:2:end)) ;
        output_file_path = fullfile(output_folder_path, output_file_name) ;
        job_name = sprintf('skel-%05d-%s',bounding_box_index,random_string_for_job_name);
        
        matlab_command_template = ...
            'try; modpath; cluster_skelh5(''%s'', ''%s'', %s, ''%s'', ''%s'', %.18g, %.18g, %.18g); catch err; fprintf(2, ''%%s\\n'', err.getReport()); quit(1); end; quit(0);' ;
        matlab_command = sprintf(matlab_command_template, ...
                                 whole_brain_p_map_h5_file_path, ...
                                 p_map_dataset_path, ...
                                 mat2str(BB), ...
                                 output_file_path, ...
                                 sizethreshold, ...
                                 probThr, ...
                                 fullh) ;
        matlab_command_line_as_tokens = { '/misc/local/matlab-2018b/bin/matlab',  '-nodisplay',  '-r', matlab_command } ;  
        %matlab_command_line_as_string = bash_command_line_from_token_list(matlab_command_line_as_tokens) ;
        % mysub = sprintf('qsub -pe batch %d -l d_rt=%d -N %s -j y -o /dev/null -b y -cwd -V %s\n',numcores,timelim,name,args);
        %stdout_file_name = sprintf('%s.stdout', job_name) ;
        %stderr_file_name = sprintf('%s.stderr', job_name) ;        
        stdout_file_name = sprintf('/dev/null') ;
        stderr_file_name = sprintf('/dev/null') ;
%         bsub_command_line = sprintf('bsub -P mouselight -n%d -We %d -J %s -o %s -e %s "%s"', ...
%                                     numcores, ...
%                                     time_limit_in_seconds/60, ...
%                                     job_name, ...
%                                     stdout_file_name, ...
%                                     stderr_file_name, ...
%                                     raw_command_line);
        bsub_command_line_as_tokens = { 'bsub', ...
                                        '-P', 'mouselight', ...
                                        '-n', num2str(numcores), ...
                                        '-We', num2str(time_limit_in_seconds/60), ...
                                        '-J', job_name, ...
                                        '-o', stdout_file_name, ...
                                        '-e', stderr_file_name, ...
                                        matlab_command_line_as_tokens } ;
        bsub_command_line_as_string = bash_command_line_from_token_list(bsub_command_line_as_tokens) ;                            
        fprintf('%s\n', bsub_command_line_as_string) ;                   
        return_code = system(bsub_command_line_as_string) ;  % submit the job!
        %return_code = 0 ;  % for testing
        if return_code ~= 0 ,
            fprintf('There was a problem with the submission of the job to create file %s\n', output_file_path) ;
        end
        %fwrite(fid,mysub);
        jobs_submitted_count = jobs_submitted_count + 1 ;
        %fprintf('jobs submitted: %d\n', jobs_submitted_count) ;
        maybe_fprint_dot(false, jobs_submitted_count) ;
    end
    maybe_fprint_dot(true, jobs_submitted_count) ;  % fprintf() final newline
    %fclose(fid) ;
    %unix(sprintf('chmod +x %s',bsub_script_file_path));
    fprintf('%d jobs submitted\n', jobs_submitted_count) ;
end

function maybe_fprint_dot(is_done, count)
    if is_done ,
        % print a newline, but only if with didn't just print a newline for
        % the final iteration
        if ~mod(count, 100) ,
            fprintf('\n') ;
        end
    else
        fprintf('.') ;
        if mod(count, 100) ,
            fprintf('\n') ;
        end
    end
end

% function result = determine_bsub_script_file_path(output_folder_path)
%     did_determine_result = false ;
%     max_round_count = 20 ;
%     for round = 1:max_round_count ,
%         putative_bsub_script_file_name = sprintf('submit_skeletonization_jobs_round_%02d.sh', round) ;
%         putative_bsub_script_file_path = fullfile(output_folder_path, putative_bsub_script_file_name) ;
%         if ~exist(putative_bsub_script_file_path, 'file') ,
%             did_determine_result = true ;
%             result = putative_bsub_script_file_path ;
%             break
%         end
%     end
%     if ~did_determine_result ,
%         error('You''ve already done %d rounds of skeletonization.  Maybe something is wrong?', max_round_count) ;
%     end
% end
