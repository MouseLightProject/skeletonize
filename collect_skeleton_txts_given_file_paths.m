function collect_skeleton_txts_given_file_paths(graph_as_mat_file_path, ...
                                                graph_as_single_text_file_path, ...
                                                whole_brain_h5_p_map_file_path, ...
                                                skeletonization_folder_path)
    
    if ~exist(graph_as_mat_file_path, 'file') || ~exist(graph_as_single_text_file_path, 'file') ,
        whole_brain_stack_shape = h5parser_new(whole_brain_h5_p_map_file_path, '/prob0') ;
        [skeleton_ijks,~,A,~] = skel2graph(skeletonization_folder_path, whole_brain_stack_shape) ;
    end        
    
    if ~exist(graph_as_mat_file_path, 'file') ,
        output_folder_path = fileparts(graph_as_mat_file_path) ;  % the folder that will hold the output files
        if ~exist(output_folder_path, 'file') ,
            mkdir(output_folder_path) ;
        end
        skeleton_graph = graph(A) ;
        save(graph_as_mat_file_path, 'skeleton_graph', 'skeleton_ijks', '-v7.3') ;
    end

    if ~exist(graph_as_single_text_file_path, 'file') ,
        output_folder_path = fileparts(graph_as_single_text_file_path) ;  % the folder that will hold the output files
        if ~exist(output_folder_path, 'file') ,
            mkdir(output_folder_path) ;
        end
        neighbors_from_node_index = neighbors_from_node_index_from_adjacency(A) ;
        save_skeleton_graph_as_text(graph_as_single_text_file_path, skeleton_ijks, neighbors_from_node_index) ;
    end

%     if ~exist(graph_as_single_binary_file_path, 'file') ,
%         output_folder_path = fileparts(graph_as_single_binary_file_path) ;  % the folder that will hold the output files
%         if ~exist(output_folder_path, 'file') ,
%             mkdir(output_folder_path) ;
%         end
%         save_skeleton_graph_as_binary(graph_as_single_binary_file_path, skeleton_ijks, neighbors_from_node_index) ;
%     end    
end
