function submit_bsub_jobs_for_skeletonization_from_paths(output_folder_path, whole_brain_p_map_h5_file_path, sizethreshold, probThr, fullh)
    output_folder_path
    whole_brain_p_map_h5_file_path
    sizethreshold
    probThr
    fullh
    
    core_count_for_each_bjob = 8;
    p_map_dataset_path = '/prob0' ;    
    
    [full_stack_shape_ijk, has_ladenness_info, laden_octree_chunk_bounding_boxes_ijk0, h5_chunk_shape_ijk] = ...
        h5parser_new(whole_brain_p_map_h5_file_path, p_map_dataset_path) ;
    
    % get a multiple of h5_chunk_shape_ijk that is around 1000^3
    analysis_chunk_shape_ijk = round(1000./h5_chunk_shape_ijk).*h5_chunk_shape_ijk;
    
    % to get %10 overlap overhead use multiple of 10
    fullh_for_create_ovelap_box = h5_chunk_shape_ijk ; % add 1 to make it odd (heuristic)
    
    [~,whole_brain_p_map_h5_base_name,~] = fileparts(whole_brain_p_map_h5_file_path) ;
    
    if ~exist(output_folder_path, 'dir') ,
        mkdir(output_folder_path) ;
    end
    
    unix(sprintf('umask g+rxw %s',output_folder_path)) ;
    unix(sprintf('chmod g+rwx %s',output_folder_path)) ;
    
    s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';

    %find number of random characters to choose from
    numRands = length(s);
    sLength = 10;

    % Construct a set of overlapped chunks for analysis
    padded_analysis_chunk_bounding_boxes_ijk1 = createOverlapBox(full_stack_shape_ijk, analysis_chunk_shape_ijk, fullh_for_create_ovelap_box) ;
    analysis_chunk_count = size(padded_analysis_chunk_bounding_boxes_ijk1, 1) ;
    fprintf('Count of analysis chunks: %d\n', analysis_chunk_count) ;
    
    % Figure out which of the analysis chunks overlap with the laden octree
    % chunks, erring on the side of caution.
    if has_ladenness_info , 
        BBoxes = laden_octree_chunk_bounding_boxes_ijk0(:,[1 4 2 5 3 6])+1 ;
        X = BBoxes(:,1:2);
        Y = BBoxes(:,3:4);
        Z = BBoxes(:,5:6);
        laden_octree_chunk_bounding_box_corners_ijk1 = unique([X(:),Y(:),Z(:)],'rows') ;
        padded_analysis_chunk_bounding_box_corners_ijk1 = ...
            [ padded_analysis_chunk_bounding_boxes_ijk1(:,1:2:end) ; ...
              padded_analysis_chunk_bounding_boxes_ijk1(:,2:2:end) ] ;
        is_padded_analysis_chunk_bounding_box_corner_in_laden_hull = ...
            inhull(padded_analysis_chunk_bounding_box_corners_ijk1, ...
                   laden_octree_chunk_bounding_box_corners_ijk1) ;
        is_analysis_chunk_laden = any(reshape(is_padded_analysis_chunk_bounding_box_corner_in_laden_hull, [analysis_chunk_count 2]),2) ;
    else
        is_analysis_chunk_laden = true(analysis_chunk_count, 1) ;
    end
    laden_analysis_chunk_count = sum(is_analysis_chunk_laden) ;
    fprintf('Count of laden analysis chunks: %d\n', laden_analysis_chunk_count) ;
    
    % Figure out which ones have already be analysed
    is_analysis_chunk_already_skeletonized = false(analysis_chunk_count, 1) ;
    output_text_file_template_path = fullfile(output_folder_path, '*.txt') ;
    extant_output_text_file_names = simple_dir(output_text_file_template_path) ;
    for ii = 1:length(extant_output_text_file_names) ,
        extant_output_text_file_name = extant_output_text_file_names{ii} ;
        bbox_index = bounding_box_index_from_file_name(extant_output_text_file_name) ;
        is_analysis_chunk_already_skeletonized(bbox_index) = true ;
    end
    fprintf('Count of already-skeletonized analysis chunks: %d\n', sum(is_analysis_chunk_already_skeletonized)) ;

    % Determine which chunks need to be skeletonized
    does_analysis_chunk_need_to_be_run = is_analysis_chunk_laden & ~is_analysis_chunk_already_skeletonized ;
    fprintf('Count of analysis chunks to be run: %d\n', sum(does_analysis_chunk_need_to_be_run)) ;
    
    % Submit the jobs
    time_limit_in_seconds = 10*60 ;
    jobs_submitted_count = 0 ;
    for analysis_chunk_index = 1:analysis_chunk_count ,
        if ~does_analysis_chunk_need_to_be_run(analysis_chunk_index) ,  % skip
            continue
        end
        padded_analysis_chunk_bounding_box_ijk1 = padded_analysis_chunk_bounding_boxes_ijk1(analysis_chunk_index, :) ;
        random_string_for_job_name = s( ceil(rand(1,sLength)*numRands) ) ;
        output_file_name = sprintf('%s_idx-%05d_stxyzendxyz-%d_%d_%d_%d_%d_%d.txt', ...
                                   whole_brain_p_map_h5_base_name, ...
                                   analysis_chunk_index, ...
                                   padded_analysis_chunk_bounding_box_ijk1(1:2:end), ...
                                   padded_analysis_chunk_bounding_box_ijk1(2:2:end)) ;
        output_file_path = fullfile(output_folder_path, output_file_name) ;
        job_name = sprintf('skel-%05d-%s', analysis_chunk_index, random_string_for_job_name) ;
        
        matlab_command_template = ...
            'try; modpath; cluster_skelh5(''%s'', ''%s'', %s, ''%s'', %.18g, %.18g, %.18g); catch err; fprintf(2, ''%%s\\n'', err.getReport()); quit(1); end; quit(0);' ;
        matlab_command = sprintf(matlab_command_template, ...
                                 whole_brain_p_map_h5_file_path, ...
                                 p_map_dataset_path, ...
                                 mat2str(padded_analysis_chunk_bounding_box_ijk1), ...
                                 output_file_path, ...
                                 sizethreshold, ...
                                 probThr, ...
                                 fullh) ;
        matlab_command_line_as_tokens = { '/misc/local/matlab-2018b/bin/matlab',  '-nodisplay',  '-r', matlab_command } ;  
        %matlab_command_line_as_string = bash_command_line_from_token_list(matlab_command_line_as_tokens) ;
        % mysub = sprintf('qsub -pe batch %d -l d_rt=%d -N %s -j y -o /dev/null -b y -cwd -V %s\n',numcores,timelim,name,args);
        stdout_file_name = sprintf('%s.out.txt', job_name) ;
        stderr_file_name = sprintf('%s.out.txt', job_name) ;        
        %stdout_file_name = sprintf('/dev/null') ;
        %stderr_file_name = sprintf('/dev/null') ;
%         bsub_command_line = sprintf('bsub -P mouselight -n%d -We %d -J %s -o %s -e %s "%s"', ...
%                                     numcores, ...
%                                     time_limit_in_seconds/60, ...
%                                     job_name, ...
%                                     stdout_file_name, ...
%                                     stderr_file_name, ...
%                                     raw_command_line);
        bsub_command_line_as_tokens = { 'bsub', ...
                                        '-P', 'mouselight', ...
                                        '-n', num2str(core_count_for_each_bjob), ...
                                        '-We', num2str(time_limit_in_seconds/60), ...
                                        '-J', job_name, ...
                                        '-oo', stdout_file_name, ...
                                        '-eo', stderr_file_name, ...
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
