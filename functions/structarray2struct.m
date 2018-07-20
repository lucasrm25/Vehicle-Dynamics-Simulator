function [ structOut ] = structarray2struct( structarray )

    fieldNames = fieldnames(structarray);
    matrix = cellfun(@double,struct2cell(structarray));
    structOut = struct();
    for i=1:numel(fieldNames)
        structOut.(fieldNames{i}) = squeeze(matrix(i,:,:));
    end

end

