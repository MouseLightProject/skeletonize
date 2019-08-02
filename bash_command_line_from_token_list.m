function result  = bash_command_line_from_token_list(token_list)
    % token_list is a cell array, each element either a string or another
    % token_list.
    %
    % Returns a string that, when used in bash, will be parsed into the
    % given token list.
    result = reshape('', [1 0]) ;
    n = length(token_list) ;
    for i = 1 : n ,
        element = token_list{i} ;
        if ischar(element) ,
            element_as_string = escape_string(element) ;
            subresult = element_as_string ;
        else
            % must be another token list
            element_as_string = bash_command_line_from_token_list(element) ;
            subresult = escape_string(element_as_string) ;
        end
        if i==1 ,
            result = subresult ;            
        else
            result = horzcat(result, ' ', subresult) ;  %#ok<AGROW>
        end
    end
end


function result = escape_string(str)
    % Escape a string so that it remains a single token when parsed by
    % bash.
    % Basically, we replace each single quote in the string with '\'' (four
    % characters), and then add a single quote at the begining and end.    
    result = horzcat('''', replace(str, '''', '''\'''''), '''') ;
end
