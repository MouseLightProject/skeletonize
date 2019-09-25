function submit_bsub_jobs_for_skeletonization(sample_date, sizethreshold, probThr, fullh)
    output_folder_path = sprintf('/nrs/mouselight/cluster/classifierOutputs/%s/skeletonization', sample_date) ;
    whole_brain_p_map_h5_file_path = ...
        sprintf('/nrs/mouselight/cluster/classifierOutputs/%s/whole-brain-p-map-as-h5/%s-whole-brain-p-map.h5', sample_date, sample_date) ;

    submit_bsub_jobs_for_skeletonization_from_paths(output_folder_path, whole_brain_p_map_h5_file_path, sizethreshold, probThr, fullh)
end
