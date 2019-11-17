function collect_skeleton_txts_given_sample_date(sample_date)
    whole_brain_h5_p_map_file_path = ...
        sprintf('/nrs/mouselight/cluster/classifierOutputs/%s/20180801_prob0/20180801_prob0_lev-6_chunk-111_111_masked-0.h5', ...
                sample_date) ;
    reconstructions_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/skeletonization', ...
                                          sample_date) ;
    skeletonization_folder_path = fullfile(reconstructions_folder_path, 'skeletonization') ;
    graph_as_mat_file_path = fullfile(reconstructions_folder_path, 'skeleton-graph.mat') ;
    
    brainSize = h5parser(whole_brain_h5_p_map_file_path, '/prob0') ;


    if exist(graph_as_mat_file_path, 'file') && ~options.do_force_computations ,
        % do nothing
        %load(graph_file_path, 'skeleton_graph', 'skeleton_ijks') ;        
    else        
        [skeleton_ijks,~,A,~] = skel2graph(skeletonization_folder_path, brainSize) ;
        skeleton_graph = graph(max(A,A')) ;
        if ~exist(reconstructions_folder_path, 'file') ,
            mkdir(reconstructions_folder_path) ;
        end
        save(graph_as_mat_file_path, 'skeleton_graph', 'skeleton_ijks', '-v7.3') ;
    end
end
