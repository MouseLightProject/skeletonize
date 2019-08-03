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
    % We first check if it's a "safe" string, in which case is can be
    % returned as-is.  (Note that this is only a subset of the strings that
    % could be returned as-is.  But it will cover the most common cases.)
    % Otherwise, we replace each single quote in the string with '\'' (four
    % characters), and then add a single quote at the begining and end.    
    % Note that we could skip the safety-checking part, and single-quote
    % everything, but it's easier to read and scan for issues if we return
    % a lot of strings unchanged.
    if is_string_safe(str) ,
        result = str ;
    else
        result = horzcat('''', replace(str, '''', '''\'''''), '''') ;
    end
end


function result = is_string_safe(str)
    result = (length(str) == length(regexp(str, '[a-zA-Z0-9_\-/]'))) ;
end
