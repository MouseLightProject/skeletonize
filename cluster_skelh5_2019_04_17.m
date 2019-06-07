function cluster_skelh5_2019_04_17()
    this_file_path = mfilename('fullpath') ;
    this_folder_path = fileparts(this_file_path) ;
    configuration_file_path = fullfile(this_folder_path, 'config_files', '20190417_prob0_config_skelh5.cfg') ;
    generate_shell_script_that_bsubs_jobs(configuration_file_path) ;
end
