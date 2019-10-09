sample_date = '2019-09-06' ;
output_folder_path = sprintf('/nrs/mouselight-v/cluster/classifierOutputs/%s/skeletonization', sample_date) ;
whole_brain_p_map_h5_file_path = ...
    sprintf('/nrs/mouselight-v/cluster/classifierOutputs/%s/whole-brain-p-map-as-h5/whole-brain-p-map.h5', sample_date) ;

sizethreshold = 300 ;
probThr = 100 ;
fullh = 15 ;

submit_bsub_jobs_for_skeletonization_from_paths(output_folder_path, whole_brain_p_map_h5_file_path, sizethreshold, probThr, fullh) ;
