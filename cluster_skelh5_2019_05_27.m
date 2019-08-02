function cluster_skelh5_2019_05_27()
    this_file_path = mfilename('fullpath') ;
    this_folder_path = fileparts(this_file_path) ;
    configuration_file_path = fullfile(this_folder_path, 'config_files', '20190527_prob0_config_skelh5.cfg') ;
    submit_bsub_jobs_for_skeletonization(configuration_file_path) ;
end
