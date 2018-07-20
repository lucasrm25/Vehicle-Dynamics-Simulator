function [ name_last_version ] = get_last_version( path, type )

%   The first argument is the folder path where the version files are
%   The second argument is 'end' or 'beginning'

    files = dir(path);

    names   = {};
    for i=1:numel(files)
        names{i} = files(i).name;

         temp = regexp(names{i},['\d+'],'match');
         if ~isempty(temp)
            if strcmp(type,'beginning')
                numbers(i) = str2num(temp{1});
            else
                numbers(i) = str2num(temp{end});
            end
         else
             numers(i) = 0;
         end
    end
    
    [M,I] = max(numbers);
    
    name_last_version = names{I};

end

