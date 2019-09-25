output_folder_path = '/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/2018-08-01/prob0' ;
whole_brain_p_map_h5_file_path = '/nrs/mouselight/cluster/classifierOutputs/2018-08-01/20180801_prob0/20180801_prob0_lev-6_chunk-111_111_masked-0.h5' ;
sizethreshold = 300 ;
probThr = 100 ;
fullh = 15 ;
submit_bsub_jobs_for_skeletonization_from_paths(output_folder_path, whole_brain_p_map_h5_file_path, sizethreshold, probThr, fullh) ;

