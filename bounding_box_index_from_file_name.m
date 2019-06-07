function result = bounding_box_index_from_file_name(file_name)
    file_name_parts_before_and_after_idx = strsplit(file_name,'idx-') ;
    file_name_part_after_idx = file_name_parts_before_and_after_idx{2} ;
    result_as_string = file_name_part_after_idx(1:5) ;  % the five characters after "idx-" are the bounding box index
    result = str2double(result_as_string) ;
end
