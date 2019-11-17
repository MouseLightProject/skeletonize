sample_date = '2018-12-01' ;
old_reconstructions_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/old-stuff-for-0p25-um-render', ...
                                          sample_date) ;
old_whole_brain_h5_p_map_file_path = ...
    fullfile(old_reconstructions_folder_path, 'whole-brain-p-map.h5') ;
%old_skeletonization_folder_path = fullfile(old_reconstructions_folder_path, 'skeletonization') ;
old_graph_as_mat_file_path = fullfile(old_reconstructions_folder_path, 'skeleton-graph.mat') ;

old_render_path = sprintf('/nrs/mouselight/SAMPLES/%s-at-0p25-um', sample_date) ;
new_render_path = sprintf('/nrs/mouselight/SAMPLES/%s', sample_date) ;

new_reconstructions_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s', ...
                                          sample_date) ;
new_graph_as_single_text_file_path = fullfile(new_reconstructions_folder_path, 'skeleton-graph-converted-from-0p25-um.txt') ;

[new_shape_xyz, new_origin, new_spacing, new_jaws_origin] = load_sample_shape_origin_and_spacing(new_render_path)
%[new_shape_xyz, new_origin, new_spacing, new_jaws_origin] = load_sample_shape_origin_and_spacing_old(new_render_path)
[old_shape_xyz, old_origin, old_spacing, old_jaws_origin] = load_sample_shape_origin_and_spacing(old_render_path)

% whole_brain_stack_shape = h5parser(old_whole_brain_h5_p_map_file_path, '/prob0') ;
% [old_skeleton_ijks,~,A,~] = skel2graph(old_skeletonization_folder_path, whole_brain_stack_shape) ;
load(old_graph_as_mat_file_path, 'skeleton_graph', 'skeleton_ijks') ;
old_skeleton_ijks = skeleton_ijks ;
clear('skeleton_ijks') ;
A = skeleton_graph.adjacency() ;
clear skeleton_graph ;

skeleton_xyzs = old_origin + old_spacing .* old_skeleton_ijks ;

new_skeleton_ijks_unrounded = (skeleton_xyzs - new_origin) ./ new_spacing ;

neighbors_from_node_index = neighbors_from_node_index_from_adjacency(A) ;
save_unrounded_skeleton_graph_as_text(new_graph_as_single_text_file_path, new_skeleton_ijks_unrounded, neighbors_from_node_index) ;




